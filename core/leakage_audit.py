"""Modulo A -- Audit del leakage per paziente.

Confronta una valutazione onesta (split raggruppato per paziente) con un controllo
negativo (split casuale sulle singole cellule, che ignora il paziente e quindi lascia
"gemelle" della stessa cellula in train e test). Il divario fra le due misura quanto
della performance riportata da una pipeline che non raggruppa per paziente sarebbe
un artefatto.

Logica di preprocessing e correzione statistica portate da pdac-ml/src/scpipe.py e
pdac-ml/src/11_model_zoo.py (generalizzate: nomi di colonna configurabili, non piu'
legate al dataset PDAC originale).
"""

from __future__ import annotations

from dataclasses import dataclass, field

import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse as sp
from sklearn.base import BaseEstimator, TransformerMixin, clone
from sklearn.ensemble import HistGradientBoostingClassifier, RandomForestClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import f1_score
from sklearn.model_selection import LeaveOneGroupOut, StratifiedGroupKFold, StratifiedKFold
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler

from core.stats import NadeauBengioResult, nadeau_bengio_test, paired_wilcoxon, wilcoxon_min_pvalue

ALPHA = 0.05


class LogCPM(BaseEstimator, TransformerMixin):
    """CPM + log1p per cellula. Stateless (non guarda le altre cellule): non puo'
    causare leakage, ma resta in Pipeline per garantire l'ordine delle operazioni."""

    def __init__(self, target_sum: float = 1e4):
        self.target_sum = target_sum

    def fit(self, X, y=None):  # noqa: ARG002
        return self

    def transform(self, X):
        counts = np.asarray(X.sum(axis=1)).ravel()
        counts[counts == 0] = 1.0
        X = sp.csr_matrix(X).multiply(self.target_sum / counts[:, None])
        X = sp.csr_matrix(X)
        X.data = np.log1p(X.data)
        return X


class ToDense(BaseEstimator, TransformerMixin):
    """HistGradientBoostingClassifier non accetta input sparso."""

    def fit(self, X, y=None):  # noqa: ARG002
        return self

    def transform(self, X):
        return X.toarray() if sp.issparse(X) else np.asarray(X)


def _pipeline(estimator: BaseEstimator, needs_dense: bool = False) -> Pipeline:
    steps = [("logcpm", LogCPM()), ("scale", StandardScaler(with_mean=False))]
    if needs_dense:
        steps.append(("densify", ToDense()))
    steps.append(("clf", estimator))
    return Pipeline(steps)


def _reference_model(seed: int) -> LogisticRegression:
    return LogisticRegression(C=1.0, class_weight="balanced", max_iter=1000, random_state=seed)


def _model_registry(seed: int) -> dict[str, tuple[BaseEstimator, bool]]:
    """Modello -> (istanza, richiede input denso)."""
    return {
        "logreg": (_reference_model(seed), False),
        "random_forest": (RandomForestClassifier(
            n_estimators=300, class_weight="balanced_subsample", n_jobs=-1, random_state=seed), False),
        "hist_gb": (HistGradientBoostingClassifier(class_weight="balanced", random_state=seed), True),
    }


def _macro_f1(y_true, y_pred, classes: np.ndarray) -> float:
    return float(f1_score(y_true, y_pred, average="macro", labels=classes, zero_division=0))


@dataclass(frozen=True)
class SchemeSummary:
    scheme: str
    fold_scores: list[float]
    n_train: list[int]
    n_test: list[int]

    @property
    def mean(self) -> float:
        return float(np.mean(self.fold_scores))

    @property
    def std(self) -> float:
        return float(np.std(self.fold_scores, ddof=1)) if len(self.fold_scores) > 1 else 0.0


def _run_scheme(X, y: np.ndarray, groups: np.ndarray | None, scheme: str,
                 n_folds: int, seed: int) -> SchemeSummary:
    classes = np.array(sorted(np.unique(y)))
    if scheme == "grouped":
        n_splits = min(n_folds, len(np.unique(groups)))
        cv = StratifiedGroupKFold(n_splits=n_splits, shuffle=True, random_state=seed)
        split_iter = cv.split(X, y, groups=groups)
    elif scheme == "random":
        n_splits = min(n_folds, int(pd.Series(y).value_counts().min()))
        cv = StratifiedKFold(n_splits=n_splits, shuffle=True, random_state=seed)
        split_iter = cv.split(X, y)
    else:
        raise ValueError(scheme)

    scores, n_train, n_test = [], [], []
    for tr, te in split_iter:
        pipe = _pipeline(_reference_model(seed))
        pipe.fit(X[tr], y[tr])
        pred = pipe.predict(X[te])
        scores.append(_macro_f1(y[te], pred, classes))
        n_train.append(len(tr))
        n_test.append(len(te))
    return SchemeSummary(scheme=scheme, fold_scores=scores, n_train=n_train, n_test=n_test)


@dataclass(frozen=True)
class PairwiseComparison:
    model: str
    mean_macro_f1: float
    delta_vs_best: float
    p_nadeau_bengio: float
    p_wilcoxon: float
    wilcoxon_floor: float
    divergent: bool
    note: str


@dataclass(frozen=True)
class ModelComparisonResult:
    scheme: str
    best_model: str
    scores: dict[str, SchemeSummary]
    comparisons: list[PairwiseComparison]


def _compare_models(X, y: np.ndarray, groups: np.ndarray, seed: int) -> ModelComparisonResult:
    classes = np.array(sorted(np.unique(y)))
    logo = LeaveOneGroupOut()
    registry = _model_registry(seed)

    scores: dict[str, SchemeSummary] = {}
    for name, (estimator, needs_dense) in registry.items():
        fold_scores, n_train, n_test = [], [], []
        for tr, te in logo.split(X, y, groups=groups):
            pipe = _pipeline(clone(estimator), needs_dense)
            pipe.fit(X[tr], y[tr])
            pred = pipe.predict(X[te])
            fold_scores.append(_macro_f1(y[te], pred, classes))
            n_train.append(len(tr))
            n_test.append(len(te))
        scores[name] = SchemeSummary(scheme="fold_lopo", fold_scores=fold_scores,
                                      n_train=n_train, n_test=n_test)

    best_model = max(scores, key=lambda m: scores[m].mean)
    best = scores[best_model]

    comparisons = []
    for name, summary in scores.items():
        if name == best_model:
            continue
        nb: NadeauBengioResult = nadeau_bengio_test(
            best.fold_scores, summary.fold_scores, best.n_train, best.n_test)
        p_w = paired_wilcoxon(best.fold_scores, summary.fold_scores)
        d = np.array(best.fold_scores) - np.array(summary.fold_scores)
        k_nonzero = int(np.sum(~np.isclose(d, 0)))
        floor = wilcoxon_min_pvalue(k_nonzero)

        sig_nb = nb.p_value < ALPHA
        sig_w = (not np.isnan(p_w)) and p_w < ALPHA
        divergent = sig_nb != sig_w
        note = ""
        if divergent:
            if sig_nb and not sig_w and p_w >= floor - 1e-12 and floor > ALPHA:
                note = (f"Nadeau-Bengio trova una differenza significativa, Wilcoxon no -- ma con "
                        f"{k_nonzero} coppie non nulle il p-value minimo raggiungibile da Wilcoxon "
                        f"e' {floor:.3f} > 0.05: strutturalmente incapace di essere significativo "
                        f"con cosi' pochi fold, indipendentemente dalla dimensione dell'effetto.")
            else:
                note = "Le due conclusioni divergono a soglia 0.05: da riportare, non da ignorare."
        comparisons.append(PairwiseComparison(
            model=name, mean_macro_f1=summary.mean, delta_vs_best=summary.mean - best.mean,
            p_nadeau_bengio=nb.p_value, p_wilcoxon=p_w, wilcoxon_floor=floor,
            divergent=divergent, note=note,
        ))

    return ModelComparisonResult(scheme="fold_lopo", best_model=best_model,
                                  scores=scores, comparisons=comparisons)


@dataclass(frozen=True)
class LeakageAuditResult:
    task: str
    n_cells: int
    n_patients: int
    n_classes: int
    n_folds_used: int
    grouped: SchemeSummary
    random: SchemeSummary
    gap: float
    std_ratio: float
    narrative: str
    model_comparison: ModelComparisonResult | None = field(default=None)


def _narrative(grouped: SchemeSummary, random_: SchemeSummary, gap: float, std_ratio: float) -> str:
    direzione = "sovrastima" if gap > 0 else ("sottostima" if gap < 0 else "non altera")
    parts = [
        f"Split per paziente (valutazione onesta): macro-F1 = {grouped.mean:.3f} ± {grouped.std:.3f} "
        f"su {len(grouped.fold_scores)} fold.",
        f"Split casuale sulle cellule (controllo negativo, NON una valutazione valida perche' "
        f"mette cellule dello stesso paziente sia in train sia in test): macro-F1 = "
        f"{random_.mean:.3f} ± {random_.std:.3f}.",
        f"Divario sulla media: {gap:+.3f} -- lo split casuale {direzione} l'accuratezza reale "
        f"del modello di circa {abs(gap):.3f} punti di macro-F1.",
    ]
    if std_ratio > 0 and np.isfinite(std_ratio):
        if std_ratio > 1.05:
            parts.append(
                f"La deviazione standard fra fold e' {std_ratio:.1f}x piu' ampia nello split onesto: "
                f"lo split casuale sottostima l'incertezza reale sulla performance di quel fattore."
            )
        elif std_ratio < 0.95:
            parts.append(
                f"La deviazione standard fra fold e' piu' ampia nello split casuale "
                f"({1/std_ratio:.1f}x): in questo caso specifico non mostra il pattern atteso, "
                f"da verificare con piu' fold o piu' pazienti."
            )
        else:
            parts.append("Le due deviazioni standard sono simili in questo dataset.")
    return " ".join(parts)


def run_leakage_audit(
    adata: ad.AnnData,
    target_col: str,
    patient_col: str,
    n_folds: int = 5,
    min_patients_for_model_comparison: int = 8,
    seed: int = 0,
) -> LeakageAuditResult:
    """Esegue il Modulo A su un task di classificazione a scelta dell'utente.

    ``target_col`` e ``patient_col`` sono nomi di colonna in ``adata.obs``, configurabili
    dall'utente (non fissi come nello script originale). Le cellule con target mancante
    sono escluse.
    """
    obs = adata.obs
    if target_col not in obs.columns:
        raise ValueError(f"colonna target '{target_col}' non trovata in adata.obs")
    if patient_col not in obs.columns:
        raise ValueError(f"colonna paziente '{patient_col}' non trovata in adata.obs")

    mask = obs[target_col].notna().values
    if mask.sum() < len(obs):
        adata = adata[mask]
        obs = adata.obs

    y = obs[target_col].astype(str).values
    groups = obs[patient_col].astype(str).values
    X = adata.X
    n_patients = len(np.unique(groups))
    n_classes = len(np.unique(y))

    if n_classes < 2:
        raise ValueError(f"servono almeno 2 classi in '{target_col}', trovate {n_classes}")
    if n_patients < 2:
        raise ValueError(f"servono almeno 2 pazienti in '{patient_col}', trovati {n_patients}")

    grouped = _run_scheme(X, y, groups, "grouped", n_folds, seed)
    random_ = _run_scheme(X, y, None, "random", n_folds, seed)

    gap = random_.mean - grouped.mean
    std_ratio = (grouped.std / random_.std) if random_.std > 0 else float("inf")
    narrative = _narrative(grouped, random_, gap, std_ratio)

    model_comparison = None
    if n_patients >= min_patients_for_model_comparison:
        model_comparison = _compare_models(X, y, groups, seed)

    return LeakageAuditResult(
        task=target_col, n_cells=adata.n_obs,
        n_patients=n_patients, n_classes=n_classes,
        n_folds_used=len(grouped.fold_scores), grouped=grouped, random=random_,
        gap=gap, std_ratio=std_ratio, narrative=narrative, model_comparison=model_comparison,
    )
