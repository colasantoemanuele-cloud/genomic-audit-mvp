"""Stimatori statistici validati: correzione di Nadeau-Bengio e cluster bootstrap.

Logica portata da pdac-ml/src/scpipe.py (nadeau_bengio_ttest), pdac-ml/src/07_robustness.py
(cluster_bootstrap sui pazienti) e pdac-ml/src/10_loco.py (cluster_boot). Entrambi gli
stimatori esistevano gia' validati nel progetto originale: qui sono generalizzati (nomi di
gruppo/colonna configurabili) ma la matematica non e' stata reimplementata da zero.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.stats import t as student_t
from scipy.stats import wilcoxon as _scipy_wilcoxon


@dataclass(frozen=True)
class NadeauBengioResult:
    mean_diff: float
    t_stat: float
    p_value: float
    df: int
    k_folds: int


def nadeau_bengio_test(
    scores_a: np.ndarray,
    scores_b: np.ndarray,
    n_train,
    n_test,
) -> NadeauBengioResult:
    """t-test appaiato corretto per la correlazione fra fold (Nadeau & Bengio, 2003).

    Il t-test standard sui fold di una cross-validation e' anticonservativo: i training
    set dei fold si sovrappongono, quindi le differenze di score non sono osservazioni
    indipendenti e la varianza del t-test ordinario e' sottostimata. La correzione
    moltiplica la varianza campionaria dei delta per (1/k + n_test/n_train), dove il
    termine 1/k e' quello del t-test appaiato standard e n_test/n_train e' il termine
    correttivo per la sovrapposizione dei training set.

    ``n_train`` e ``n_test`` possono essere scalari (gia' mediati) o array-like con un
    valore per fold: in quel caso si usa la media, perche' con StratifiedGroupKFold le
    dimensioni dei fold non sono garantite identiche.
    """
    a = np.asarray(scores_a, dtype=float)
    b = np.asarray(scores_b, dtype=float)
    if a.shape != b.shape:
        raise ValueError(f"scores_a e scores_b devono avere la stessa forma: {a.shape} vs {b.shape}")

    d = a - b
    k = len(d)
    n_train_mean = float(np.mean(n_train))
    n_test_mean = float(np.mean(n_test))

    if k < 2 or np.allclose(d, 0.0):
        return NadeauBengioResult(mean_diff=float(d.mean()) if k else 0.0, t_stat=0.0,
                                   p_value=1.0, df=max(k - 1, 0), k_folds=k)

    var_corretta = d.var(ddof=1) * (1.0 / k + n_test_mean / max(n_train_mean, 1.0))
    denom = np.sqrt(var_corretta)
    if denom == 0:
        return NadeauBengioResult(mean_diff=float(d.mean()), t_stat=0.0, p_value=1.0,
                                   df=k - 1, k_folds=k)

    t_stat = float(d.mean() / denom)
    p_value = float(2 * (1 - student_t.cdf(abs(t_stat), df=k - 1)))
    return NadeauBengioResult(mean_diff=float(d.mean()), t_stat=t_stat, p_value=p_value,
                               df=k - 1, k_folds=k)


def paired_wilcoxon(scores_a: np.ndarray, scores_b: np.ndarray) -> float:
    """Wilcoxon signed-rank appaiato su due vettori di score fold-per-fold. NaN se tutte
    le differenze sono nulle (test non definito)."""
    a, b = np.asarray(scores_a, dtype=float), np.asarray(scores_b, dtype=float)
    try:
        return float(_scipy_wilcoxon(a, b).pvalue)
    except ValueError:
        return float("nan")


def wilcoxon_min_pvalue(k_nonzero: int) -> float:
    """P-value minimo raggiungibile dal test di Wilcoxon esatto con k coppie non nulle,
    quando tutte le differenze hanno lo stesso segno: 2 / 2^k. Con pochi fold (k<=5-6)
    questo minimo puo' restare sopra 0.05, rendendo il test strutturalmente incapace di
    dichiarare significativita' anche di fronte a un effetto reale e consistente."""
    if k_nonzero < 1:
        return 1.0
    return 2.0 / (2.0 ** k_nonzero)


@dataclass(frozen=True)
class BootstrapResult:
    mean: float
    ci_low: float
    ci_high: float
    n_groups: int
    n_obs: int
    sufficient: bool


def cluster_bootstrap(
    values: np.ndarray,
    groups: np.ndarray,
    n_boot: int = 2000,
    seed: int = 0,
    ci_percentiles: tuple[float, float] = (2.5, 97.5),
    min_groups: int = 2,
) -> BootstrapResult:
    """Cluster bootstrap: ricampiona con reinserimento i GRUPPI (es. pazienti), non le
    singole osservazioni. Necessario quando le osservazioni sono annidate in gruppi
    (es. piu' cellule/cloni per paziente) e non sono quindi indipendenti fra loro: un
    bootstrap "ingenuo" sulle singole osservazioni produce intervalli di confidenza
    artificialmente stretti.

    Se il numero di gruppi e' inferiore a ``min_groups``, l'intervallo non e'
    affidabile: ``sufficient`` sara' False e ci_low/ci_high saranno NaN, cosi' il
    chiamante puo' mostrare "numerosita' insufficiente" invece di un numero fuorviante.
    """
    v = np.asarray(values, dtype=float)
    g = np.asarray(groups)
    if len(v) != len(g):
        raise ValueError("values e groups devono avere la stessa lunghezza")

    unique_groups = np.unique(g)
    n_groups = len(unique_groups)
    mean = float(v.mean()) if len(v) else float("nan")

    if n_groups < min_groups:
        return BootstrapResult(mean=mean, ci_low=float("nan"), ci_high=float("nan"),
                                n_groups=n_groups, n_obs=len(v), sufficient=False)

    by_group = {grp: v[g == grp] for grp in unique_groups}
    rng = np.random.default_rng(seed)
    boot_means = np.empty(n_boot)
    for i in range(n_boot):
        draw = rng.choice(unique_groups, size=n_groups, replace=True)
        boot_means[i] = np.concatenate([by_group[grp] for grp in draw]).mean()

    lo, hi = np.percentile(boot_means, ci_percentiles)
    return BootstrapResult(mean=mean, ci_low=float(lo), ci_high=float(hi),
                            n_groups=n_groups, n_obs=len(v), sufficient=True)
