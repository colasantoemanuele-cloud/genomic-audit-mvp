"""Verifica di calibrazione su dati simulati per gli stimatori di core/stats.py.

OBBLIGATORIO: nessuno stimatore statistico va esposto nell'interfaccia prima che questi
test passino (vedi prompt di progetto, sezione "Requisito trasversale").

Per ciascuno stimatore (intervallo t cluster-robusto, esposto; correzione di Nadeau-Bengio, Wilcoxon e
t non corretto in LeaveOneGroupOut, non esposti perché fuori banda o non usati; cluster bootstrap
percentile, solo per riprodurre i numeri della tesi):
  - sotto ipotesi nulla (nessun effetto reale simulato): su >=200 repliche, la frazione
    di risultati dichiarati significativi a soglia 0.05 deve cadere approssimativamente
    fra 2% e 8% (calibrazione: né zero, né esplosa);
  - sotto un effetto noto e sostanziale (iniettato deliberatamente): lo stimatore deve
    rilevarlo nella maggioranza delle repliche (potenza).

I semi casuali sono fissi ovunque: i test sono deterministici, non flaky.

Regime d'uso (revisione, A1-A2). Lo strumento usa Nadeau-Bengio SOLO nel confronto fra modelli,
con LeaveOneGroupOut (un fold per paziente, k >= 8 fold di dimensioni diverse) e il rapporto
mean(n_test)/mean(n_train), ma è fuori banda e il confronto è esposto solo in forma descrittiva;
l'intervallo cluster-robusto sui pazienti con 5-15 pazienti (sotto 5 nessun intervallo). Le
calibrazioni in quei regimi sono in fondo al file. I test con DummyClassifier a k = 3 e k = 5
restano come controlli fuori dal regime d'uso: due DummyClassifier predicono a caso, senza usare
il training, quindi le loro differenze fra fold NON sono correlate dalla sovrapposizione dei
training set, cioè manca proprio il fenomeno che la correzione compensa.
"""

from __future__ import annotations

import numpy as np
import pytest
from sklearn.dummy import DummyClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import f1_score
from sklearn.model_selection import StratifiedKFold

from core.stats import (
    bell_mccaffrey_df,
    cluster_bootstrap,
    cluster_t_interval,
    nadeau_bengio_test,
    wilcoxon_min_pvalue,
)

ALPHA = 0.05
CALIBRATION_LOW, CALIBRATION_HIGH = 0.02, 0.08


# --------------------------------------------------------------------------- #
# Nadeau-Bengio: ipotesi nulla
# --------------------------------------------------------------------------- #
def _nb_null_replicate(rep: int, seed0: int, n: int, k: int) -> float:
    """Un confronto fra due DummyClassifier(strategy='stratified') su una CV a k fold: nessuna
    vera differenza di performance (entrambi ignorano le feature). Poiché i DummyClassifier non
    usano il training, le differenze fra fold non sono correlate dalla sovrapposizione dei
    training set: questo scenario NON contiene il fenomeno che Nadeau-Bengio corregge."""
    rng = np.random.default_rng(seed0 + rep)
    y = rng.integers(0, 2, n)
    X = np.zeros((n, 1))
    skf = StratifiedKFold(n_splits=k, shuffle=True, random_state=rep)
    fa, fb, n_train, n_test = [], [], [], []
    for tr, te in skf.split(X, y):
        clf_a = DummyClassifier(strategy="stratified", random_state=2 * rep)
        clf_b = DummyClassifier(strategy="stratified", random_state=2 * rep + 1)
        clf_a.fit(X[tr], y[tr])
        clf_b.fit(X[tr], y[tr])
        fa.append(f1_score(y[te], clf_a.predict(X[te]), average="macro", zero_division=0))
        fb.append(f1_score(y[te], clf_b.predict(X[te]), average="macro", zero_division=0))
        n_train.append(len(tr))
        n_test.append(len(te))
    p = nadeau_bengio_test(fa, fb, n_train, n_test).p_value
    return p


@pytest.mark.slow
def test_nadeau_bengio_calibration_null():
    """k=3 con DummyClassifier: controllo fuori dal regime d'uso (vedi il docstring del file)."""
    R, seed0, n, k = 1000, 55_003, 90, 3
    n_sig = sum(_nb_null_replicate(rep, seed0, n, k) < ALPHA for rep in range(R))
    frac = n_sig / R
    assert CALIBRATION_LOW <= frac <= CALIBRATION_HIGH, (
        f"Nadeau-Bengio non calibrato sotto H0: {n_sig}/{R} = {frac:.3f} significativi, "
        f"atteso in [{CALIBRATION_LOW}, {CALIBRATION_HIGH}]"
    )


@pytest.mark.slow
@pytest.mark.xfail(
    reason=(
        "Fuori dal regime d'uso (lo strumento usa Nadeau-Bengio solo con LeaveOneGroupOut, k >= 8). "
        "Con due DummyClassifier a k=5 la correzione è conservativa sotto H0: 18/1000=1.8% (R=1000) e "
        "38/2000=1.9% (R=2000), sotto la banda [2%, 8%], che NON è stata allargata. Causa: i "
        "DummyClassifier predicono a caso senza usare il training, quindi le differenze fra fold non "
        "sono correlate e il termine n_test/n_train gonfia una varianza che non andava corretta. Non "
        "dipende dai gradi di libertà. strict=True: uno xpass fa fallire la suite."
    ),
    strict=True,
)
def test_nadeau_bengio_calibration_null_k5_default():
    """k=5 con DummyClassifier: controllo fuori dal regime d'uso (il nome del test è storico;
    k=5 NON è il regime in cui lo strumento usa Nadeau-Bengio). Banda [2%, 8%]."""
    R, seed0, n, k = 1000, 55_005, 90, 5
    n_sig = sum(_nb_null_replicate(rep, seed0, n, k) < ALPHA for rep in range(R))
    frac = n_sig / R
    assert CALIBRATION_LOW <= frac <= CALIBRATION_HIGH, (
        f"Nadeau-Bengio (k=5, default) non calibrato sotto H0: {n_sig}/{R} = {frac:.3f}, "
        f"atteso in [{CALIBRATION_LOW}, {CALIBRATION_HIGH}]"
    )


# --------------------------------------------------------------------------- #
# Nadeau-Bengio: potenza
# --------------------------------------------------------------------------- #
@pytest.mark.slow
def test_nadeau_bengio_power():
    R, seed0, n, p, k, wscale = 300, 51_005, 120, 8, 5, 8.0
    n_sig = 0
    for rep in range(R):
        rng = np.random.default_rng(seed0 + rep)
        X = rng.normal(size=(n, p))
        w = rng.normal(size=p) * wscale
        prob = 1 / (1 + np.exp(-(X @ w)))
        y = (rng.random(n) < prob).astype(int)
        skf = StratifiedKFold(n_splits=k, shuffle=True, random_state=rep)
        fa, fb, n_train, n_test = [], [], [], []
        for tr, te in skf.split(X, y):
            clf_a = LogisticRegression(max_iter=1000).fit(X[tr], y[tr])
            clf_b = DummyClassifier(strategy="stratified", random_state=rep).fit(X[tr], y[tr])
            fa.append(f1_score(y[te], clf_a.predict(X[te]), average="macro", zero_division=0))
            fb.append(f1_score(y[te], clf_b.predict(X[te]), average="macro", zero_division=0))
            n_train.append(len(tr))
            n_test.append(len(te))
        p_val = nadeau_bengio_test(fa, fb, n_train, n_test).p_value
        if p_val < ALPHA:
            n_sig += 1
    frac = n_sig / R
    assert frac > 0.5, f"potenza insufficiente: {n_sig}/{R} = {frac:.3f} (atteso maggioranza)"


def test_wilcoxon_min_pvalue_formula():
    assert wilcoxon_min_pvalue(5) == pytest.approx(2 / 32)
    assert wilcoxon_min_pvalue(1) == pytest.approx(1.0)
    assert wilcoxon_min_pvalue(0) == 1.0


def test_nadeau_bengio_zero_variance_returns_p_one():
    """Se le due sequenze di score sono identiche fold per fold (nessuna differenza
    osservabile), il test deve restituire p=1.0 e non dividere per zero."""
    scores = [0.7, 0.71, 0.69, 0.72, 0.70]
    res = nadeau_bengio_test(scores, scores, n_train=[80] * 5, n_test=[20] * 5)
    assert res.p_value == 1.0


# --------------------------------------------------------------------------- #
# Cluster bootstrap: ipotesi nulla
# --------------------------------------------------------------------------- #
def _simulate_grouped(rng: np.random.Generator, n_groups: int, effect: float) -> tuple[np.ndarray, np.ndarray]:
    """Popolazione annidata in gruppi (es. pazienti): dimensione di gruppo variabile,
    come nel caso reale in cui ogni paziente contribuisce un numero diverso di cloni."""
    groups, values = [], []
    for gi in range(n_groups):
        n_obs = rng.integers(3, 15)
        groups += [gi] * n_obs
        values += list(rng.normal(effect, 1.0, n_obs))
    return np.array(values), np.array(groups)


@pytest.mark.slow
def test_cluster_bootstrap_calibration_null():
    R, seed0, n_groups, n_boot = 1000, 90_000, 25, 500
    n_sig = 0
    for rep in range(R):
        rng = np.random.default_rng(seed0 + rep)
        values, groups = _simulate_grouped(rng, n_groups, effect=0.0)
        res = cluster_bootstrap(values, groups, n_boot=n_boot, seed=rep)
        assert res.sufficient
        if res.ci_low > 0 or res.ci_high < 0:
            n_sig += 1
    frac = n_sig / R
    assert CALIBRATION_LOW <= frac <= CALIBRATION_HIGH, (
        f"cluster bootstrap non calibrato sotto H0: {n_sig}/{R} = {frac:.3f} significativi, "
        f"atteso in [{CALIBRATION_LOW}, {CALIBRATION_HIGH}]"
    )


@pytest.mark.slow
def test_cluster_bootstrap_power():
    R, seed0, n_groups, n_boot, effect = 300, 95_050, 25, 500, 0.5
    n_sig = 0
    for rep in range(R):
        rng = np.random.default_rng(seed0 + rep)
        values, groups = _simulate_grouped(rng, n_groups, effect=effect)
        res = cluster_bootstrap(values, groups, n_boot=n_boot, seed=rep)
        if res.sufficient and (res.ci_low > 0 or res.ci_high < 0):
            n_sig += 1
    frac = n_sig / R
    assert frac > 0.5, f"potenza insufficiente: {n_sig}/{R} = {frac:.3f} (atteso maggioranza)"


def test_cluster_bootstrap_insufficient_groups_flagged():
    """Con meno di min_groups gruppi, l'IC non va riportato: sufficient=False, non un
    numero fuorviante calcolato su troppo pochi pazienti."""
    rng = np.random.default_rng(0)
    values = rng.normal(0, 1, 12)
    groups = np.array([0, 0, 0, 0, 1, 1, 1, 1, 2, 2, 2, 2])  # solo 3 gruppi
    res = cluster_bootstrap(values, groups, n_boot=200, min_groups=5)
    assert not res.sufficient
    assert np.isnan(res.ci_low) and np.isnan(res.ci_high)
    assert res.n_groups == 3


def test_cluster_bootstrap_naive_vs_cluster_ci_width():
    """Controllo di contrasto (non di calibrazione): un bootstrap 'ingenuo' che
    ricampiona le osservazioni singole invece dei gruppi produce un IC più stretto di
    quello a cluster, quando le osservazioni sono fortemente correlate entro gruppo --
    coerente con 07_robustness.py, dove l'IC a cluster risultava più largo."""
    rng = np.random.default_rng(1)
    n_groups = 15
    groups, values = [], []
    for gi in range(n_groups):
        group_mean = rng.normal(0, 1.5)  # forte eterogeneità fra pazienti
        n_obs = rng.integers(10, 20)
        groups += [gi] * n_obs
        values += list(rng.normal(group_mean, 0.2, n_obs))  # poca variabilità entro paziente
    values, groups = np.array(values), np.array(groups)

    cluster_res = cluster_bootstrap(values, groups, n_boot=1000, seed=0)

    rng_naive = np.random.default_rng(0)
    naive_boot = np.array([
        rng_naive.choice(values, len(values), replace=True).mean() for _ in range(1000)
    ])
    naive_lo, naive_hi = np.percentile(naive_boot, [2.5, 97.5])

    cluster_width = cluster_res.ci_high - cluster_res.ci_low
    naive_width = naive_hi - naive_lo
    assert cluster_width > naive_width


# --------------------------------------------------------------------------- #
# Intervallo cluster-robusto (revisione, A1): calibrazione nel regime d'uso
# --------------------------------------------------------------------------- #
# Bande dichiarate PRIMA di eseguire questi test (non si allargano):
#   - sotto H0, frazione di intervalli al 95% che NON contengono il valore vero in [2%, 8%],
#     per ogni scenario e per 5, 8, 10 e 15 gruppi (il regime d'uso del Modulo B: 4-11
#     pazienti su GSE278694), 1000 repliche per caso. Lo stimatore è in forma chiusa: 1000
#     repliche costano pochi secondi e riducono l'errore Monte Carlo (d.s. 0,7 punti al 5%);
#   - potenza: con un effetto noto, l'intervallo esclude lo zero nella maggioranza delle
#     repliche (10 gruppi);
#   - il percentile (cluster_bootstrap) nello scenario normale a 5 gruppi resta fuori banda
#     (> 8%): documenta il difetto che ha motivato la sostituzione.
# Scenari, con dimensioni dei gruppi diverse:
#   normale     3-14 osservazioni per gruppo, N(0, 1) (lo scenario della revisione);
#   eterogeneo  3-14 osservazioni, effetto di gruppo N(0, 0,7^2) più N(0, 1);
#   cellula     tasso d'errore per cellula come su GSE278694 (scala logit): 40-1000 cellule per
#               paziente, tasso per paziente Beta con media 0,09 e d.s. 0,06, cellule raggruppate
#               in cloni di 1-10 cellule con tasso di clone Beta(4p, 4(1-p));
#   clone       tasso per clone (scala logit): 5-60 unità per paziente, Beta media 0,18 d.s. 0,125;
#   eccesso     eccesso di discordanza: 13-285 coppie per paziente, valori in [-1, 1] (metà zeri,
#               il resto N(0, 0,28^2)), media di paziente N(0, 0,05^2), troncati a [-1, 1].
CI_REPS = 1000
CI_GROUPS = (5, 8, 10, 15)
CI_SCENARIOS = ("normale", "eterogeneo", "cellula", "clone", "eccesso")
CI_SCALE = {"cellula": "logit", "clone": "logit"}


def _beta_params(m: float, s: float) -> tuple[float, float]:
    k = m * (1 - m) / s ** 2 - 1
    return m * k, (1 - m) * k


def _ci_scenario(rng: np.random.Generator, G: int, kind: str, effect: float = 0.0):
    """(valori, gruppi) per uno scenario; la media vera è ``_ci_truth(kind) + effect``."""
    if kind in ("normale", "eterogeneo"):
        sizes = rng.integers(3, 15, G)
        shift = rng.normal(0, 0.7, G) if kind == "eterogeneo" else np.zeros(G)
        groups = np.repeat(np.arange(G), sizes)
        return rng.normal(effect, 1.0, sizes.sum()) + shift[groups], groups
    if kind == "cellula":
        sizes = rng.integers(40, 1001, G)
        p = rng.beta(*_beta_params(0.09, 0.06), G)
        vals, groups = [], []
        for g, (n, pg) in enumerate(zip(sizes, p)):
            k = rng.integers(1, 11, n)  # più cloni di quanti servono, poi si tronca a n cellule
            q = rng.beta(4 * pg + 1e-3, 4 * (1 - pg) + 1e-3, n)
            cell_q = np.repeat(q, k)[:n]
            vals.append((rng.random(n) < cell_q).astype(float))
            groups.append(np.full(n, g))
        return np.concatenate(vals), np.concatenate(groups)
    if kind == "clone":
        sizes = rng.integers(5, 61, G)
        p = rng.beta(*_beta_params(0.18, 0.125), G)
        groups = np.repeat(np.arange(G), sizes)
        return (rng.random(sizes.sum()) < p[groups]).astype(float), groups
    if kind == "eccesso":
        sizes = rng.integers(13, 286, G)
        groups = np.repeat(np.arange(G), sizes)
        n = sizes.sum()
        x = np.where(rng.random(n) < 0.5, 0.0, rng.normal(0, 0.28, n)) + rng.normal(0, 0.05, G)[groups]
        return np.clip(x + effect, -1, 1), groups
    raise ValueError(kind)


_TRUTH: dict[str, float] = {}


def _ci_truth(kind: str) -> float:
    """Media vera dello scenario sotto H0 (per 'eccesso' la troncatura a [-1, 1] la sposta di
    pochissimo da zero: calcolata una volta con 4 milioni di estrazioni)."""
    if kind not in _TRUTH:
        if kind == "cellula":
            _TRUTH[kind] = 0.09
        elif kind == "clone":
            _TRUTH[kind] = 0.18
        elif kind == "eccesso":
            rng = np.random.default_rng(123)
            m = 4_000_000
            x = np.where(rng.random(m) < 0.5, 0.0, rng.normal(0, 0.28, m)) + rng.normal(0, 0.05, m)
            _TRUTH[kind] = float(np.clip(x, -1, 1).mean())
        else:
            _TRUTH[kind] = 0.0
    return _TRUTH[kind]


@pytest.mark.slow
@pytest.mark.parametrize("n_groups", CI_GROUPS)
@pytest.mark.parametrize("kind", CI_SCENARIOS)
def test_cluster_t_interval_calibration_null(kind, n_groups):
    truth = _ci_truth(kind)
    seed0 = 3_000_000 + 100_000 * CI_SCENARIOS.index(kind) + 1_000 * n_groups
    miss = 0
    for rep in range(CI_REPS):
        values, groups = _ci_scenario(np.random.default_rng(seed0 + rep), n_groups, kind)
        r = cluster_t_interval(values, groups, min_groups=5, scale=CI_SCALE.get(kind, "identity"))
        assert r.sufficient, r.reason
        miss += not (r.ci_low <= truth <= r.ci_high)
    frac = miss / CI_REPS
    print(f"\n[A1 {kind}, {n_groups} gruppi] fuori dall'IC95%: {miss}/{CI_REPS} = {frac:.3f} "
          f"(banda [{CALIBRATION_LOW}, {CALIBRATION_HIGH}])")
    assert CALIBRATION_LOW <= frac <= CALIBRATION_HIGH


@pytest.mark.slow
def test_cluster_t_interval_power():
    R, n_sig = 300, 0
    for rep in range(R):
        values, groups = _ci_scenario(np.random.default_rng(4_000_000 + rep), 10, "normale", effect=0.5)
        r = cluster_t_interval(values, groups, min_groups=5)
        n_sig += r.ci_low > 0 or r.ci_high < 0
    print(f"\n[A1 potenza] effetto 0,5, 10 gruppi: {n_sig}/{R} = {n_sig / R:.3f}")
    assert n_sig / R > 0.5


@pytest.mark.slow
def test_percentile_bootstrap_is_anticonservative_with_five_groups():
    """Il difetto che ha motivato A1, misurato: il percentile a 5 gruppi è fuori banda."""
    R, miss = 400, 0
    for rep in range(R):
        values, groups = _ci_scenario(np.random.default_rng(5_000_000 + rep), 5, "normale")
        r = cluster_bootstrap(values, groups, n_boot=500, seed=rep)
        miss += not (r.ci_low <= 0.0 <= r.ci_high)
    print(f"\n[A1 percentile, 5 gruppi] fuori dall'IC95%: {miss}/{R} = {miss / R:.3f}")
    assert miss / R > CALIBRATION_HIGH


def test_cluster_t_interval_refuses_below_min_groups_and_degenerate_cases():
    g = np.repeat(np.arange(4), 5)
    r = cluster_t_interval(np.arange(20.0), g, min_groups=5)
    assert not r.sufficient and np.isnan(r.ci_low) and "4 pazienti" in r.reason
    zeros = cluster_t_interval(np.zeros(30), np.repeat(np.arange(6), 5), min_groups=5, scale="logit")
    assert not zeros.sufficient and zeros.mean == 0.0
    with pytest.raises(ValueError, match="logit"):
        cluster_t_interval(np.array([-0.5, 0.2] * 10), np.repeat(np.arange(10), 2), 5, scale="logit")


def test_bell_mccaffrey_df_equal_sizes_is_g_minus_one():
    assert bell_mccaffrey_df(np.array([7, 7, 7, 7, 7])) == pytest.approx(4.0)
    assert bell_mccaffrey_df(np.array([3, 300, 5, 8, 40])) < 4.0


def test_cluster_t_interval_point_estimate_is_pooled_mean():
    values = np.array([1, 0, 0, 1, 1, 1, 0, 0, 0, 0], dtype=float)
    groups = np.array([0, 0, 1, 1, 2, 2, 3, 3, 4, 4])
    r = cluster_t_interval(values, groups, min_groups=5, scale="logit")
    assert r.mean == pytest.approx(0.4) and r.ci_low < 0.4 < r.ci_high


# --------------------------------------------------------------------------- #
# Nadeau-Bengio nel regime d'uso (revisione, A2): LeaveOneGroupOut, 8-20 gruppi
# --------------------------------------------------------------------------- #
# Bande dichiarate PRIMA di eseguire questi test (non si allargano): sotto H0, frazione di
# p < 0.05 in [2%, 8%] su 400 repliche, per k = 8, 10, 15 e 20 gruppi, separatamente per
# Nadeau-Bengio (come lo usa leakage_audit: rapporto fra le medie di n_test e n_train), per il
# Wilcoxon appaiato (core.stats.wilcoxon_test) e per il t appaiato non corretto (solo misura di
# contrasto: lo strumento non lo usa, perché in teoria è anticonservativo quando i training si
# sovrappongono). Uno stimatore che non entra nella banda non viene esposto.
# Scenario: G pazienti con 15-60 cellule ciascuno (dimensioni diverse), prevalenza della classe
# positiva per paziente U(0.3, 0.7), 10 feature con effetto paziente N(0, 0.8^2) su ciascuna; le
# feature 0-4 e 5-9 portano lo stesso segnale (+-0.6 secondo la classe), quindi due regressioni
# logistiche, una sulle prime cinque e una sulle altre cinque, sono ugualmente buone in
# aspettativa (H0 esatta per simmetria) e sono addestrate sugli stessi training sovrapposti di
# LeaveOneGroupOut. Metrica: macro-F1 sulle classi presenti nel paziente di test, come lo strumento.
NB_LOGO_REPS = 400
NB_LOGO_GROUPS = (8, 10, 15, 20)


def _logo_replicate(rep: int, G: int) -> tuple[float, float, float]:
    from scipy.stats import t as student_t

    from core.stats import wilcoxon_test

    rng = np.random.default_rng(6_000_000 + 10_000 * G + rep)
    sizes = rng.integers(15, 61, G)
    groups = np.repeat(np.arange(G), sizes)
    n = sizes.sum()
    y = (rng.random(n) < rng.uniform(0.3, 0.7, G)[groups]).astype(int)
    X = rng.normal(0, 1, (n, 10)) + rng.normal(0, 0.8, (G, 10))[groups]
    X += 0.6 * (2 * y[:, None] - 1)
    fa, fb, n_train, n_test = [], [], [], []
    for g in range(G):
        te = groups == g
        tr = ~te
        labels = np.unique(y[te])
        for cols, out in ((slice(0, 5), fa), (slice(5, 10), fb)):
            clf = LogisticRegression(max_iter=500).fit(X[tr][:, cols], y[tr])
            out.append(f1_score(y[te], clf.predict(X[te][:, cols]), average="macro", labels=labels,
                                zero_division=0))
        n_train.append(int(tr.sum()))
        n_test.append(int(te.sum()))
    p_nb = nadeau_bengio_test(fa, fb, n_train, n_test).p_value
    p_w = wilcoxon_test(fa, fb)[0]
    d = np.array(fa) - np.array(fb)
    t_stat = d.mean() / np.sqrt(d.var(ddof=1) / G)
    p_t = float(2 * (1 - student_t.cdf(abs(t_stat), G - 1)))
    return p_nb, p_w, p_t


_LOGO_CACHE: dict[int, np.ndarray] = {}


def _logo_pvalues(G: int) -> np.ndarray:
    if G not in _LOGO_CACHE:
        _LOGO_CACHE[G] = np.array([_logo_replicate(rep, G) for rep in range(NB_LOGO_REPS)])
    return _LOGO_CACHE[G]


def _logo_rate(G: int, col: int, name: str) -> float:
    p = _logo_pvalues(G)[:, col]
    frac = float(np.mean(np.nan_to_num(p, nan=1.0) < ALPHA))
    print(f"\n[A2 LOGO, {G} gruppi] {name}: p < 0.05 in {frac:.3f} delle {NB_LOGO_REPS} repliche "
          f"(banda [{CALIBRATION_LOW}, {CALIBRATION_HIGH}])")
    return frac


# Esito misurato (revisione/log/a2_calibrazione.txt): Nadeau-Bengio 0.018, 0.007, 0.007, 0.005 per
# 8, 10, 15, 20 gruppi, tutti sotto la banda (conservativo: in LeaveOneGroupOut le differenze fra
# fold sono quasi scorrelate, mentre la correzione assume una correlazione n_test/n_train e quasi
# raddoppia la varianza); Wilcoxon 0.060, 0.083, 0.072, 0.045, fuori banda a 10 gruppi. Le bande
# NON sono state allargate. Conseguenza: nessuno dei due p-value è mostrato da CLI, app o report;
# il confronto fra modelli è descrittivo.
_NB_XFAIL = pytest.mark.xfail(strict=True, reason=(
    "Nadeau-Bengio in LeaveOneGroupOut è conservativo sotto H0 (0.5-1.8% di p < 0.05, banda [2%, 8%]): "
    "la correzione assume fra le differenze dei fold una correlazione n_test/n_train che in LOGO quasi "
    "non c'è. Non esposto."))
_W_XFAIL = pytest.mark.xfail(strict=True, reason=(
    "Wilcoxon in LeaveOneGroupOut a 10 gruppi: 0.083 di p < 0.05 sotto H0, sopra la banda [2%, 8%]. "
    "Non esposto."))


@pytest.mark.slow
@pytest.mark.parametrize("G", [pytest.param(g, marks=_NB_XFAIL) for g in NB_LOGO_GROUPS])
def test_nadeau_bengio_calibration_logo(G):
    assert CALIBRATION_LOW <= _logo_rate(G, 0, "Nadeau-Bengio") <= CALIBRATION_HIGH


@pytest.mark.slow
@pytest.mark.parametrize("G", [pytest.param(g, marks=_W_XFAIL) if g == 10 else g for g in NB_LOGO_GROUPS])
def test_wilcoxon_calibration_logo(G):
    assert CALIBRATION_LOW <= _logo_rate(G, 1, "Wilcoxon") <= CALIBRATION_HIGH


@pytest.mark.slow
@pytest.mark.parametrize("G", NB_LOGO_GROUPS)
def test_uncorrected_t_calibration_logo(G):
    assert CALIBRATION_LOW <= _logo_rate(G, 2, "t non corretto") <= CALIBRATION_HIGH
