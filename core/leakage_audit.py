"""Modulo A -- Audit del leakage per paziente.

Confronta una valutazione onesta (split raggruppato per paziente) con un controllo
negativo (split casuale sulle singole cellule, che ignora il paziente e quindi lascia
"gemelle" della stessa cellula in train e test). Il divario fra le due misura quanto
della performance riportata da una pipeline che non raggruppa per paziente sarebbe
un artefatto.

Motore pensato per dati reali (decine di migliaia di cellule, ~20-30k geni):
- la matrice resta SPARSA (CSR) dall'inizio alla fine: nessuna conversione globale in denso;
- dentro ogni fold, e stimati SOLO sul training: log-CPM (per cellula), selezione dei 2000
  geni a varianza più alta (HVG), standardizzazione senza centratura e, per i modelli
  non lineari, TruncatedSVD a 50 componenti (come pdac-ml/src/11_model_zoo.py);
- training limitato a ``max_train_cells`` cellule per fold con sottocampionamento
  stratificato per paziente (seme fisso, dichiarato nell'output);
- macro-F1 calcolata sulle sole classi presenti nel fold di test; le classi assenti sono
  registrate e dichiarate (prima venivano contate come F1 = 0: un fold perfetto valeva 0.667).

Logica di preprocessing e correzione statistica portate da pdac-ml/src/scpipe.py e
pdac-ml/src/11_model_zoo.py.
"""

from __future__ import annotations

import time
import warnings
from collections.abc import Callable
from dataclasses import dataclass, field

import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse as sp
from sklearn.base import BaseEstimator, TransformerMixin, clone
from sklearn.decomposition import TruncatedSVD
from sklearn.ensemble import HistGradientBoostingClassifier, RandomForestClassifier
from sklearn.exceptions import ConvergenceWarning
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import f1_score
from sklearn.model_selection import LeaveOneGroupOut, StratifiedGroupKFold, StratifiedKFold
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler

from core.stats import NadeauBengioResult, nadeau_bengio_test, paired_wilcoxon, wilcoxon_min_pvalue

ALPHA = 0.05
N_HVG_DEFAULT = 2000
N_SVD_DEFAULT = 50
MAX_TRAIN_CELLS_DEFAULT = 20_000
N_TOP_GENES_XAI = 50
MIN_PATIENTS_FOR_MODEL_COMPARISON_DEFAULT = 8
# soglia del divario oltre la quale il report segnala "leakage rilevabile"
GAP_ALERT = 0.05

Progress = Callable[[str], None]


# --------------------------------------------------------------------------- #
# Trasformatori (tutti sparsi)
# --------------------------------------------------------------------------- #
class LogCPM(BaseEstimator, TransformerMixin):
    """CPM + log1p per cellula. Stateless (non guarda le altre cellule): non può
    causare leakage, ma resta in Pipeline per garantire l'ordine delle operazioni."""

    def __init__(self, target_sum: float = 1e4):
        self.target_sum = target_sum

    def fit(self, X, y=None):  # noqa: ARG002
        return self

    def transform(self, X):
        X = sp.csr_matrix(X, dtype=np.float64)
        counts = np.asarray(X.sum(axis=1)).ravel()
        counts[counts == 0] = 1.0
        X = sp.csr_matrix(sp.diags(self.target_sum / counts) @ X)
        X.data = np.log1p(X.data)
        return X


class HVGSelector(BaseEstimator, TransformerMixin):
    """Seleziona gli ``n_top`` geni con varianza più alta (sui valori log-CPM), stimata SOLO
    sulle cellule passate a ``fit`` (il training del fold): nessuna informazione del test
    entra nella scelta. In caso di parità vince l'indice di gene più basso (deterministico).
    Lavora sulla matrice sparsa, senza densificarla."""

    def __init__(self, n_top: int = N_HVG_DEFAULT):
        self.n_top = n_top

    def fit(self, X, y=None):  # noqa: ARG002
        X = sp.csr_matrix(X)
        mean = np.asarray(X.mean(axis=0)).ravel()
        sq = np.asarray(X.multiply(X).mean(axis=0)).ravel()
        var = sq - mean ** 2
        k = min(self.n_top, X.shape[1])
        order = np.argsort(-var, kind="stable")
        self.genes_ = np.sort(order[:k])
        return self

    def transform(self, X):
        return sp.csr_matrix(X)[:, self.genes_]


def _pipeline(estimator: BaseEstimator, reduce: bool, n_hvg: int, n_svd: int, seed: int) -> Pipeline:
    steps = [("logcpm", LogCPM()), ("hvg", HVGSelector(n_top=n_hvg)),
             ("scale", StandardScaler(with_mean=False))]
    if reduce:
        steps.append(("svd", TruncatedSVD(n_components=n_svd, random_state=seed)))
    steps.append(("clf", estimator))
    return Pipeline(steps)


def _reference_model(seed: int) -> LogisticRegression:
    return LogisticRegression(C=1.0, class_weight="balanced", max_iter=1000, random_state=seed)


def _model_registry(seed: int) -> dict[str, tuple[BaseEstimator, bool]]:
    """Modello -> (istanza, usa la riduzione SVD). La regressione logistica lavora sui geni
    HVG (sparsi, interpretabili); random forest e gradient boosting sulle componenti SVD."""
    return {
        "logreg": (_reference_model(seed), False),
        "random_forest": (RandomForestClassifier(
            n_estimators=300, class_weight="balanced_subsample", n_jobs=-1, random_state=seed), True),
        "hist_gb": (HistGradientBoostingClassifier(class_weight="balanced", random_state=seed), True),
    }


# --------------------------------------------------------------------------- #
# Metrica e sottocampionamento
# --------------------------------------------------------------------------- #
def fold_macro_f1(y_true: np.ndarray, y_pred: np.ndarray,
                  all_classes: np.ndarray) -> tuple[float, list[str]]:
    """Macro-F1 sulle sole classi presenti nel fold di test (``y_true``). Una classe del
    dataset assente dal fold non entra nella media (non vale F1 = 0) ed è restituita nella
    lista delle classi assenti, da dichiarare. Le predizioni di una classe assente restano
    errori: abbassano il recall delle classi presenti."""
    present = np.unique(y_true)
    score = float(f1_score(y_true, y_pred, average="macro", labels=present, zero_division=0))
    absent = sorted(set(map(str, all_classes)) - set(map(str, present)))
    return score, absent


def cap_training(tr: np.ndarray, groups: np.ndarray, max_cells: int, seed: int) -> np.ndarray:
    """Se il training supera ``max_cells`` cellule, ne estrae ``max_cells`` con allocazione
    proporzionale per paziente (almeno 1 cellula per paziente), con seme fisso."""
    if max_cells is None or len(tr) <= max_cells:
        return tr
    rng = np.random.default_rng(seed)
    g = groups[tr]
    pats, counts = np.unique(g, return_counts=True)
    quota = np.maximum(1, np.floor(counts / counts.sum() * max_cells).astype(int))
    keep = []
    for p, q in zip(pats, quota):
        idx = tr[g == p]
        keep.append(rng.choice(idx, size=min(q, len(idx)), replace=False))
    return np.sort(np.concatenate(keep))


# --------------------------------------------------------------------------- #
# Schemi di valutazione
# --------------------------------------------------------------------------- #
@dataclass(frozen=True)
class SchemeSummary:
    scheme: str
    fold_scores: list[float]
    n_train: list[int]
    n_test: list[int]
    absent_classes: list[list[str]] = field(default_factory=list)
    n_train_before_cap: list[int] = field(default_factory=list)
    top_genes: list[list[str]] = field(default_factory=list)
    not_converged: list[bool] = field(default_factory=list)

    @property
    def mean(self) -> float:
        return float(np.mean(self.fold_scores))

    @property
    def std(self) -> float:
        return float(np.std(self.fold_scores, ddof=1)) if len(self.fold_scores) > 1 else 0.0

    @property
    def n_folds_with_absent_classes(self) -> int:
        return sum(bool(a) for a in self.absent_classes)

    @property
    def n_not_converged(self) -> int:
        return sum(self.not_converged)

    @property
    def capped(self) -> bool:
        return any(b > a for a, b in zip(self.n_train, self.n_train_before_cap))


class _Clock:
    """Avanzamento: numero di addestramenti fatti/previsti, tempo trascorso e stima."""

    def __init__(self, total: int, progress: Progress | None):
        self.total, self.done, self.t0, self.progress = total, 0, time.time(), progress

    def step(self, what: str) -> None:
        self.done += 1
        if self.progress is None:
            return
        el = time.time() - self.t0
        eta = el / self.done * (self.total - self.done)
        self.progress(f"[{self.done}/{self.total}] {what} -- trascorsi {el:.0f} s, "
                      f"stima residua {eta:.0f} s")


def _top_genes(pipe: Pipeline, gene_names: np.ndarray, k: int) -> list[str]:
    clf = pipe.named_steps["clf"]
    coef = np.abs(np.atleast_2d(clf.coef_)).max(axis=0)
    hvg = pipe.named_steps["hvg"].genes_
    order = np.argsort(-coef, kind="stable")[:k]
    return [str(gene_names[hvg[i]]) for i in order]


def _fit_fold(X, y, groups, tr, te, i, classes, make_pipe, max_train, seed, gene_names):
    """Un fold: sottocampionamento del training, fit, macro-F1, geni principali. Funzione pura
    (stessi argomenti -> stesso risultato), quindi eseguibile in parallelo senza cambiare i numeri."""
    n_before = len(tr)
    tr = cap_training(tr, groups, max_train, seed + i)
    present = np.unique(y[tr])
    if len(present) < 2:
        raise ValueError(
            f"nel fold {i + 1} il training contiene una sola classe ('{present[0]}'): lasciando fuori "
            f"i pazienti del test non resta nessuna cellula delle altre classi. Succede quando una "
            f"classe compare in un solo paziente: servono almeno 2 pazienti per classe.")
    pipe = make_pipe()
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always", ConvergenceWarning)
        pipe.fit(X[tr], y[tr])
    not_converged = any(issubclass(w.category, ConvergenceWarning) for w in caught)
    s, a = fold_macro_f1(y[te], pipe.predict(X[te]), classes)
    top = _top_genes(pipe, gene_names, N_TOP_GENES_XAI) if gene_names is not None else None
    return s, a, len(tr), len(te), n_before, top, not_converged


def _run_folds(X, y, groups, splits, classes, make_pipe, label: str, clock: _Clock,
               max_train: int, seed: int, gene_names: np.ndarray | None,
               n_jobs: int = 1) -> SchemeSummary:
    """Esegue i fold, in parallelo se ``n_jobs`` > 1 (stessi risultati: ogni fold è
    deterministico e indipendente dagli altri). L'avanzamento segue l'ordine dei fold."""
    args = [(X, y, groups, tr, te, i, classes, make_pipe, max_train, seed, gene_names)
            for i, (tr, te) in enumerate(splits)]
    if n_jobs > 1 and len(splits) > 1:
        from joblib import Parallel, delayed
        results = Parallel(n_jobs=min(n_jobs, len(splits)), return_as="generator")(
            delayed(_fit_fold)(*a) for a in args)
    else:
        results = (_fit_fold(*a) for a in args)
    scores, n_train, n_test, absent, before, tops, nc = [], [], [], [], [], [], []
    for i, (s, a, ntr, nte, nb, top, not_conv) in enumerate(results):
        nc.append(bool(not_conv))
        scores.append(s)
        absent.append(a)
        n_train.append(ntr)
        n_test.append(nte)
        before.append(nb)
        if top is not None:
            tops.append(top)
        clock.step(f"{label}, fold {i + 1}/{len(splits)}")
    return SchemeSummary(scheme=label, fold_scores=scores, n_train=n_train, n_test=n_test,
                         absent_classes=absent, n_train_before_cap=before, top_genes=tops,
                         not_converged=nc)


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


def _compare(scores: dict[str, SchemeSummary]) -> ModelComparisonResult:
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
                        f"è {floor:.3f} > 0.05: strutturalmente incapace di essere significativo "
                        f"con così pochi fold, indipendentemente dalla dimensione dell'effetto.")
            else:
                note = "Le due conclusioni divergono a soglia 0.05: da riportare, non da ignorare."
        comparisons.append(PairwiseComparison(
            model=name, mean_macro_f1=summary.mean, delta_vs_best=summary.mean - best.mean,
            p_nadeau_bengio=nb.p_value, p_wilcoxon=p_w, wilcoxon_floor=floor,
            divergent=divergent, note=note))
    return ModelComparisonResult(scheme="fold_lopo", best_model=best_model, scores=scores,
                                 comparisons=comparisons)


# --------------------------------------------------------------------------- #
# Stabilità delle spiegazioni (descrittiva)
# --------------------------------------------------------------------------- #
@dataclass(frozen=True)
class XaiStability:
    """Jaccard fra i ``k`` geni con coefficiente più grande della regressione logistica nei
    diversi fold. Misura DESCRITTIVA, senza inferenza: se rimuovere pazienti cambia i geni
    scelti più di quanto lo cambi rimuovere cellule a caso, l'informazione è organizzata per
    paziente."""
    k: int
    grouped_matrix: np.ndarray
    random_matrix: np.ndarray
    grouped_mean: float
    random_mean: float


def _jaccard_matrix(sets: list[list[str]]) -> np.ndarray:
    n = len(sets)
    m = np.ones((n, n))
    for i in range(n):
        for j in range(i + 1, n):
            a, b = set(sets[i]), set(sets[j])
            m[i, j] = m[j, i] = len(a & b) / len(a | b) if a | b else float("nan")
    return m


def _off_diag_mean(m: np.ndarray) -> float:
    n = m.shape[0]
    return float(m[~np.eye(n, dtype=bool)].mean()) if n > 1 else float("nan")


# --------------------------------------------------------------------------- #
# Risultato e narrativa
# --------------------------------------------------------------------------- #
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
    xai: XaiStability | None = field(default=None)
    settings: dict = field(default_factory=dict)
    elapsed_seconds: float = 0.0


def _narrative(grouped: SchemeSummary, random_: SchemeSummary, gap: float, std_ratio: float) -> str:
    direzione = "sovrastima" if gap > 0 else ("sottostima" if gap < 0 else "non altera")
    parts = [
        f"Split per paziente (valutazione onesta): macro-F1 = {grouped.mean:.3f} ± {grouped.std:.3f} "
        f"su {len(grouped.fold_scores)} fold.",
        f"Split casuale sulle cellule (controllo negativo, NON una valutazione valida perché "
        f"mette cellule dello stesso paziente sia in train sia in test): macro-F1 = "
        f"{random_.mean:.3f} ± {random_.std:.3f}.",
        f"Divario sulla media: {gap:+.3f} -- lo split casuale {direzione} l'accuratezza reale "
        f"del modello di circa {abs(gap):.3f} punti di macro-F1.",
    ]
    if std_ratio > 0 and np.isfinite(std_ratio):
        if std_ratio > 1.05:
            parts.append(
                f"La deviazione standard fra fold è {std_ratio:.1f}x più ampia nello split onesto: "
                f"lo split casuale sottostima l'incertezza reale sulla performance di quel fattore."
            )
        elif std_ratio < 0.95:
            parts.append(
                f"La deviazione standard fra fold è più ampia nello split casuale "
                f"({1/std_ratio:.1f}x): in questo caso specifico non mostra il pattern atteso, "
                f"da verificare con più fold o più pazienti."
            )
        else:
            parts.append("Le due deviazioni standard sono simili in questo dataset.")
    if grouped.n_folds_with_absent_classes:
        parts.append(
            f"In {grouped.n_folds_with_absent_classes} fold su {len(grouped.fold_scores)} dello "
            f"split per paziente alcune classi non compaiono nel test: la macro-F1 di quei fold "
            f"è calcolata sulle sole classi presenti (le classi assenti sono elencate nel report).")
    nc = grouped.n_not_converged + random_.n_not_converged
    if nc:
        parts.append(
            f"In {nc} fold su {len(grouped.fold_scores) + len(random_.fold_scores)} il classificatore "
            f"non ha raggiunto la convergenza entro il numero massimo di iterazioni: i punteggi di "
            f"quei fold vanno letti con cautela (succede tipicamente quando il segnale è debole o "
            f"assente).")
    if grouped.capped or random_.capped:
        parts.append(
            f"Il training è stato limitato a {max(grouped.n_train + random_.n_train):,} cellule "
            f"per fold con un sottocampionamento stratificato per paziente (seme fisso).")
    return " ".join(parts)


def run_leakage_audit(
    adata: ad.AnnData,
    target_col: str,
    patient_col: str,
    n_folds: int = 5,
    min_patients_for_model_comparison: int = MIN_PATIENTS_FOR_MODEL_COMPARISON_DEFAULT,
    seed: int = 0,
    benchmark: bool = True,
    max_train_cells: int | None = MAX_TRAIN_CELLS_DEFAULT,
    n_hvg: int = N_HVG_DEFAULT,
    n_svd: int = N_SVD_DEFAULT,
    progress: Progress | None = None,
    n_jobs: int = -1,
) -> LeakageAuditResult:
    """Esegue il Modulo A su un task di classificazione a scelta dell'utente.

    ``benchmark=False`` (flag ``--rapido`` della CLI) salta il confronto fra modelli e usa solo
    il classificatore di riferimento. ``n_jobs``: processi per i fold della regressione
    logistica (-1 = tutti i core meno uno); non cambia i risultati. Le cellule con target o
    paziente mancante sono escluse.
    """
    t0 = time.time()
    if n_jobs is None or n_jobs < 1:
        import os
        n_jobs = max(1, (os.cpu_count() or 1) - 1)
    obs = adata.obs
    if target_col not in obs.columns:
        raise ValueError(f"colonna target '{target_col}' non trovata in adata.obs")
    if patient_col not in obs.columns:
        raise ValueError(f"colonna paziente '{patient_col}' non trovata in adata.obs")

    mask = obs[target_col].notna().values & obs[patient_col].notna().values
    y = obs[target_col].astype(str).values[mask]
    groups = obs[patient_col].astype(str).values[mask]
    X = sp.csr_matrix(adata.X[mask] if mask.sum() < len(obs) else adata.X)
    if X.nnz and (X.data < 0).any():
        raise ValueError("la matrice contiene valori negativi: il Modulo A richiede conteggi grezzi")
    sample = X.data[: min(X.nnz, 1_000_000)]
    if sample.size and not np.allclose(sample, np.round(sample)):
        raise ValueError(
            "la matrice non contiene conteggi interi (sembra già normalizzata o log-trasformata): "
            "il Modulo A applica la propria normalizzazione e richiede conteggi grezzi. Usa il layer "
            "dei conteggi (es. adata.layers['counts']) come adata.X.")
    gene_names = np.asarray(adata.var_names.astype(str))
    n_patients = len(np.unique(groups))
    classes = np.array(sorted(np.unique(y)))
    if len(classes) < 2:
        raise ValueError(f"servono almeno 2 classi in '{target_col}', trovate {len(classes)}")
    if n_patients < 2:
        raise ValueError(f"servono almeno 2 pazienti in '{patient_col}', trovati {n_patients}")
    min_class = int(pd.Series(y).value_counts().min())
    if min_class < 2:
        raise ValueError("almeno una classe ha una sola cellula: lo split casuale stratificato non è definito")

    # la SVD non può avere più componenti dei geni selezionati
    n_svd = max(1, min(n_svd, min(n_hvg, X.shape[1]) - 1))
    n_grouped = min(n_folds, n_patients)
    n_random = min(n_folds, min_class)
    do_benchmark = benchmark and n_patients >= min_patients_for_model_comparison
    n_models = len(_model_registry(seed))
    clock = _Clock(n_grouped + n_random + (n_models * n_patients if do_benchmark else 0), progress)
    if progress:
        progress(f"Modulo A: {len(y):,} cellule, {X.shape[1]:,} geni, {n_patients} pazienti, "
                 f"{len(classes)} classi. Addestramenti previsti: {clock.total}"
                 + ("" if do_benchmark else " (confronto fra modelli escluso)") + ".")

    ref = lambda: _pipeline(_reference_model(seed), False, n_hvg, n_svd, seed)  # noqa: E731
    sgkf = StratifiedGroupKFold(n_splits=n_grouped, shuffle=True, random_state=seed)
    grouped = _run_folds(X, y, groups, list(sgkf.split(X, y, groups=groups)), classes, ref,
                         "split per paziente", clock, max_train_cells, seed, gene_names, n_jobs)
    skf = StratifiedKFold(n_splits=n_random, shuffle=True, random_state=seed)
    random_ = _run_folds(X, y, groups, list(skf.split(X, y)), classes, ref,
                         "split casuale", clock, max_train_cells, seed, gene_names, n_jobs)

    gap = random_.mean - grouped.mean
    std_ratio = (grouped.std / random_.std) if random_.std > 0 else float("inf")

    gm, rm = _jaccard_matrix(grouped.top_genes), _jaccard_matrix(random_.top_genes)
    xai = XaiStability(k=N_TOP_GENES_XAI, grouped_matrix=gm, random_matrix=rm,
                       grouped_mean=_off_diag_mean(gm), random_mean=_off_diag_mean(rm))

    model_comparison = None
    if do_benchmark:
        logo = list(LeaveOneGroupOut().split(X, y, groups=groups))
        scores = {}
        for name, (est, reduce) in _model_registry(seed).items():
            make = lambda est=est, reduce=reduce: _pipeline(clone(est), reduce, n_hvg, n_svd, seed)  # noqa: E731
            # random forest e gradient boosting usano già tutti i core al loro interno
            scores[name] = _run_folds(X, y, groups, logo, classes, make, f"confronto modelli: {name}",
                                      clock, max_train_cells, seed, None,
                                      n_jobs if name == "logreg" else 1)
        model_comparison = _compare(scores)

    return LeakageAuditResult(
        task=target_col, n_cells=int(len(y)), n_patients=n_patients, n_classes=len(classes),
        n_folds_used=len(grouped.fold_scores), grouped=grouped, random=random_,
        gap=gap, std_ratio=std_ratio, narrative=_narrative(grouped, random_, gap, std_ratio),
        model_comparison=model_comparison, xai=xai,
        settings={"n_hvg": n_hvg, "n_svd": n_svd, "max_train_cells": max_train_cells,
                  "benchmark": bool(do_benchmark), "seed": seed, "n_jobs": n_jobs,
                  "min_patients_for_model_comparison": min_patients_for_model_comparison},
        elapsed_seconds=time.time() - t0,
    )
