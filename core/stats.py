"""Stimatori statistici: correzione di Nadeau-Bengio, intervallo cluster-robusto sui pazienti,
cluster bootstrap percentile (solo per riprodurre numeri precedenti).

- ``nadeau_bengio_test``: portato da pdac-ml/src/scpipe.py (nadeau_bengio_ttest).
- ``cluster_t_interval``: intervallo usato per OGNI intervallo esposto dal Modulo B. Intervallo t
  sulla media pesata per osservazioni (cellule, cloni, coppie) con varianza cluster-robusta CR2 e
  gradi di libertà di Bell-McCaffrey; per i tassi (valori in [0, 1]) è calcolato su scala logit.
  Sostituisce il cluster bootstrap percentile (revisione, A1), che con pochi gruppi è
  anticonservativo: 14,5% di intervalli che mancano il valore vero con 5 gruppi nel test di
  calibrazione (tests/test_stats_calibration.py, revisione/log/a1_calibrazione.txt), contro il 5%
  nominale; nelle simulazioni esplorative 18% con 5 gruppi e 9% con 10 (revisione/REGISTRO.md).
- ``cluster_bootstrap``: bootstrap percentile sui gruppi, portato da pdac-ml/src/07_robustness.py e
  10_loco.py. Resta solo per riprodurre gli intervalli della tesi (validation/) e per la frazione di
  CD8 sperimentale (core/cd8_propagation.py, che ne replica lo schema di ricampionamento). Non è
  usato per gli intervalli esposti dal Modulo B.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.stats import t as student_t
from scipy.stats import wilcoxon as _scipy_wilcoxon


@dataclass(frozen=True)
class NadeauBengioResult:
    """Esito del t-test corretto di Nadeau-Bengio: differenza media, statistica, p-value, gradi di libertà e fold."""
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

    ``n_train`` e ``n_test`` possono essere scalari (già mediati) o array-like con un
    valore per fold: il termine correttivo è il rapporto fra le MEDIE, mean(n_test) /
    mean(n_train), non la media dei rapporti per fold. Lo strumento lo usa solo sul confronto
    fra modelli con LeaveOneGroupOut (un fold per paziente, k >= 8 fold di dimensioni diverse):
    è il regime calibrato in tests/test_stats_calibration.py.
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


def wilcoxon_test(scores_a: np.ndarray, scores_b: np.ndarray) -> tuple[float, str]:
    """Wilcoxon signed-rank appaiato fold per fold. Ritorna (p-value, metodo).

    Le differenze nulle sono escluse (metodo di Wilcoxon). Il metodo è scelto qui, in modo
    esplicito e uguale per ogni versione di scipy: distribuzione esatta se non ci sono pareggi
    fra i valori assoluti delle differenze e al massimo 50 coppie non nulle; con pareggi,
    permutazione esatta dei segni fino a 13 coppie; oltre, approssimazione normale. Solo nel
    primo caso il p-value non può scendere sotto ``wilcoxon_min_pvalue``; con l'approssimazione
    normale sì. (nan, "non definito") se tutte le differenze sono nulle."""
    d = np.asarray(scores_a, dtype=float) - np.asarray(scores_b, dtype=float)
    d = d[~np.isclose(d, 0.0)]
    n = len(d)
    if n == 0:
        return float("nan"), "non definito"
    ties = len(np.unique(np.round(np.abs(d), 12))) < n
    if not ties and n <= 50:
        return float(_scipy_wilcoxon(d, method="exact").pvalue), "esatto"
    if n <= 13:
        from scipy.stats import PermutationMethod
        return float(_scipy_wilcoxon(d, method=PermutationMethod()).pvalue), "permutazione esatta"
    return float(_scipy_wilcoxon(d, method="asymptotic").pvalue), "approssimazione normale"


def paired_wilcoxon(scores_a: np.ndarray, scores_b: np.ndarray) -> float:
    """Solo il p-value di ``wilcoxon_test``. NaN se tutte le differenze sono nulle."""
    return wilcoxon_test(scores_a, scores_b)[0]


def holm_adjust(pvalues: list[float]) -> list[float]:
    """P-value corretti con il metodo di Holm (1979) per confronti multipli, nello stesso
    ordine dell'input. I NaN restano NaN e non contano nel numero di confronti."""
    p = np.asarray(pvalues, dtype=float)
    out = np.full(len(p), np.nan)
    ok = np.flatnonzero(~np.isnan(p))
    m = len(ok)
    if m == 0:
        return out.tolist()
    order = ok[np.argsort(p[ok], kind="stable")]
    running = 0.0
    for rank, idx in enumerate(order):
        running = max(running, min(1.0, (m - rank) * p[idx]))
        out[idx] = running
    return out.tolist()


def wilcoxon_min_pvalue(k_nonzero: int) -> float:
    """P-value minimo raggiungibile dal test di Wilcoxon esatto (o dalla permutazione esatta dei
    segni) con k coppie non nulle, quando tutte le differenze hanno lo stesso segno: 2 / 2^k.
    Con k <= 5 il minimo resta sopra 0.05: il test non può dichiarare significatività, qualunque
    sia l'effetto. Non vale quando ``wilcoxon_test`` usa l'approssimazione normale (pareggi e più
    di 13 coppie)."""
    if k_nonzero < 1:
        return 1.0
    return 2.0 / (2.0 ** k_nonzero)


@dataclass(frozen=True)
class ClusterInterval:
    """Stima puntuale (media pesata per osservazioni) e intervallo al 95% sui gruppi.

    ``sufficient`` è False quando l'intervallo non è prodotto (troppi pochi gruppi, variabilità
    nulla fra gruppi, tasso pari a 0 o 1 su scala logit): ci_low e ci_high sono allora NaN e
    ``reason`` dice perché. ``method`` e ``df`` descrivono come è stato calcolato."""
    mean: float
    ci_low: float
    ci_high: float
    n_groups: int
    n_obs: int
    sufficient: bool
    method: str = ""
    df: float | None = None
    reason: str = ""


# Nome storico, usato dal codice e dai test precedenti alla revisione.
BootstrapResult = ClusterInterval


def bell_mccaffrey_df(sizes: np.ndarray) -> float:
    """Gradi di libertà di Bell e McCaffrey (2002) per la varianza CR2 della media pesata
    (regressione sulla sola intercetta), con modello di lavoro a errori indipendenti e
    omoschedastici. Con gruppi di uguale dimensione vale G - 1; con dimensioni diverse è minore
    (Imbens e Kolesár 2016). Formula: (tr M)^2 / tr(M^2), con M = A^(1/2) K A^(1/2),
    K = diag(n) - n n'/N e A = diag(1 / ((1 - n_g/N) N^2))."""
    n = np.asarray(sizes, dtype=float)
    N = n.sum()
    a = 1.0 / ((1.0 - n / N) * N ** 2)
    K = np.diag(n) - np.outer(n, n) / N
    r = np.sqrt(a)
    M = r[:, None] * K * r[None, :]
    return float(np.trace(M) ** 2 / np.sum(M * M))


def cluster_t_interval(
    values: np.ndarray,
    groups: np.ndarray,
    min_groups: int,
    scale: str = "identity",
    level: float = 0.95,
) -> ClusterInterval:
    """Intervallo per la media pesata per osservazioni quando le osservazioni sono annidate in
    gruppi (cellule, cloni o coppie annidati nei pazienti) e il numero di gruppi è piccolo.

    Stima puntuale: media di tutte le osservazioni (pesata per numero di osservazioni del
    gruppo, come prima della revisione). Varianza cluster-robusta CR2 (Bell e McCaffrey 2002)
    della regressione sulla sola intercetta: V = sum_g S_g^2 / ((1 - n_g/N) N^2), con S_g la
    somma dei residui del gruppo g; quantile della t di Student con i gradi di libertà di
    Bell-McCaffrey (``bell_mccaffrey_df``). ``scale="logit"`` (per i tassi, valori in [0, 1]):
    l'intervallo è calcolato per logit(media) con il metodo delta e riportato sulla scala
    originale, così resta in (0, 1) e tiene conto dell'asimmetria dei tassi vicini a zero.

    Calibrazione (tests/test_stats_calibration.py, bande dichiarate prima dell'esecuzione):
    5, 8, 10 e 15 gruppi di dimensioni diverse, con valori normali, con eterogeneità fra
    gruppi e con tassi binari modellati su GSE278694. Sotto ``min_groups`` gruppi nessun
    intervallo (``sufficient=False``)."""
    if scale not in ("identity", "logit"):
        raise ValueError(f"scale deve essere 'identity' o 'logit', non '{scale}'")
    v = np.asarray(values, dtype=float)
    g = np.asarray(groups)
    if len(v) != len(g):
        raise ValueError("values e groups devono avere la stessa lunghezza")
    if np.isnan(v).any():
        raise ValueError("values contiene valori mancanti (NaN)")
    method = "t cluster-robusto CR2, gdl Bell-McCaffrey" + (", scala logit" if scale == "logit" else "")
    if len(v) == 0:
        return ClusterInterval(float("nan"), float("nan"), float("nan"), 0, 0, False, method,
                               None, "nessuna osservazione")
    if scale == "logit" and ((v < 0).any() or (v > 1).any()):
        raise ValueError("scale='logit' richiede valori in [0, 1] (tassi)")
    uniques, codes = np.unique(g, return_inverse=True)
    G = len(uniques)
    mean = float(v.mean())

    def refused(reason: str, df: float | None = None) -> ClusterInterval:
        return ClusterInterval(mean, float("nan"), float("nan"), G, len(v), False, method, df, reason)

    if G < max(min_groups, 2):
        return refused(f"{G} pazienti, ne servono almeno {max(min_groups, 2)}")
    sizes = np.bincount(codes).astype(float)
    sums = np.bincount(codes, weights=v)
    N = sizes.sum()
    S = sums - mean * sizes
    var = float(np.sum(S ** 2 / (1.0 - sizes / N)) / N ** 2)
    df = bell_mccaffrey_df(sizes)
    if not var > 0:
        return refused("nessuna variabilità fra pazienti: intervallo non definito", df)
    q = float(student_t.ppf(0.5 + level / 2, df))
    se = float(np.sqrt(var))
    if scale == "logit":
        if not 0.0 < mean < 1.0:
            return refused(f"tasso pari a {mean:g}: intervallo su scala logit non definito", df)
        lg = np.log(mean / (1.0 - mean))
        sl = se / (mean * (1.0 - mean))
        lo, hi = 1.0 / (1.0 + np.exp(-(lg - q * sl))), 1.0 / (1.0 + np.exp(-(lg + q * sl)))
    else:
        lo, hi = mean - q * se, mean + q * se
    return ClusterInterval(mean, float(lo), float(hi), G, len(v), True, method, df)


def cluster_bootstrap(
    values: np.ndarray,
    groups: np.ndarray,
    n_boot: int = 2000,
    seed: int = 0,
    ci_percentiles: tuple[float, float] = (2.5, 97.5),
    min_groups: int = 2,
) -> BootstrapResult:
    """Cluster bootstrap percentile: ricampiona con reinserimento i GRUPPI (es. pazienti), non
    le singole osservazioni (un bootstrap sulle singole osservazioni darebbe intervalli troppo
    stretti quando le osservazioni sono annidate nei gruppi).

    NON è usato per gli intervalli esposti: con pochi gruppi l'intervallo percentile è troppo
    stretto (14,5% di intervalli che mancano il valore vero con 5 gruppi nel test di calibrazione;
    18% con 5 e 9% con 10 nelle simulazioni esplorative di revisione/REGISTRO.md). Resta per riprodurre gli intervalli della tesi e per
    la frazione di CD8 sperimentale. Con meno di ``min_groups`` gruppi ``sufficient`` è False e
    ci_low/ci_high sono NaN.
    """
    v = np.asarray(values, dtype=float)
    g = np.asarray(groups)
    if len(v) != len(g):
        raise ValueError("values e groups devono avere la stessa lunghezza")

    unique_groups = np.unique(g)
    n_groups = len(unique_groups)
    mean = float(v.mean()) if len(v) else float("nan")

    if n_groups < min_groups:
        return ClusterInterval(mean=mean, ci_low=float("nan"), ci_high=float("nan"),
                               n_groups=n_groups, n_obs=len(v), sufficient=False,
                               method="bootstrap percentile sui gruppi",
                               reason=f"{n_groups} gruppi, ne servono almeno {min_groups}")

    by_group = {grp: v[g == grp] for grp in unique_groups}
    rng = np.random.default_rng(seed)
    boot_means = np.empty(n_boot)
    for i in range(n_boot):
        draw = rng.choice(unique_groups, size=n_groups, replace=True)
        boot_means[i] = np.concatenate([by_group[grp] for grp in draw]).mean()

    lo, hi = np.percentile(boot_means, ci_percentiles)
    return ClusterInterval(mean=mean, ci_low=float(lo), ci_high=float(hi), n_groups=n_groups,
                           n_obs=len(v), sufficient=True, method="bootstrap percentile sui gruppi")
