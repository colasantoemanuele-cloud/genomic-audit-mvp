"""Calibrazione dell'audit del disegno (Intervento 1) su metadati simulati.

OBBLIGATORIO prima di esporre core/design_audit.py in CLI, app o report.

Bande e tolleranze dichiarate PRIMA di eseguire i test (non vanno allargate per farli
passare: un caso fuori banda si corregge nel metodo o si marca xfail(strict=True)):

  (a) fattori indipendenti: frazione di repliche con almeno un allarme su Cramér V
      (V corretto di Bergsma >= 0.5, solo su tabelle valutabili) in [0%, 5%], su 200
      repliche. Nessun fatto strutturale (annidamento, uno-a-uno, esito determinato)
      deve comparire: il generatore lo rende di fatto impossibile.
  (b) annidamento iniettato: rilevato nel 100% delle 200 repliche (deterministico).
  (c) confondimento completo paziente-esito: rilevato nel 100% delle 200 repliche.
  Potenza del V: detection rate a V iniettato 0.3 / 0.5 / 0.7 / 0.9 riportata come
      curva (stampata), senza banda; si richiede solo che sia monotona non decrescente.
  Scomposizione della varianza (quote 50% paziente / 30% tessuto / 20% residuo,
      12 pazienti x 2 tessuti, 300 geni): in >= 95% di 200 repliche ogni quota stimata
      entro +-0.10 dal vero, e distorsione media (media delle stime - vero) entro +-0.03.
  Intervallo bootstrap delle quote (sui pazienti): copertura dell'IC 95% del valore vero
      in [90%, 99%] su 200 repliche, per ciascuna quota.
  Fattori confusi: la scomposizione rifiuta di produrre quote ("non identificabile").
  Caso di riferimento GSE278694: tumore vs adiacente stimabile con bassa potenza su
      5 pazienti; scRNA vs snRNA non stimabile.

I semi sono fissi: i test sono deterministici.
"""

from __future__ import annotations

import functools

import numpy as np
import pytest

from core.design_audit import (
    CRAMER_V_THRESHOLD,
    run_design_audit,
    variance_decomposition,
    variance_decomposition_from_units,
    variance_shares,
)
from core.synthetic import (
    make_design_sheet_association,
    make_design_sheet_independent,
    make_design_sheet_nested,
    make_design_sheet_outcome_confounded,
    make_gse278694_like_sheet,
    make_pseudobulk_adata,
    make_variance_units,
)

R = 200
V_FALSE_ALARM_BAND = (0.0, 0.05)
SHARE_TOL, SHARE_FRACTION_WITHIN, SHARE_BIAS_TOL = 0.10, 0.95, 0.03
COVERAGE_BAND = (0.90, 0.99)
TRUE_SHARES = {"patient": 0.5, "tissue": 0.3}


# --------------------------------------------------------------------------- #
# (a) fattori indipendenti
# --------------------------------------------------------------------------- #
def test_independent_factors_false_alarm_rate():
    n_v_alarm, n_structural = 0, 0
    for rep in range(R):
        sheet = make_design_sheet_independent(seed=10_000 + rep)
        res = run_design_audit(
            sheet, patient_col="patient",
            technical_cols={"batch": "batch", "chemistry": "chemistry", "protocol": "protocol",
                            "library": "library"},
            outcome_cols=["tissue"],
        )
        if any(p.v_alarm for p in res.pairs):
            n_v_alarm += 1
        if any(f.kind in ("annidamento", "uno-a-uno", "esito-determinato") for f in res.findings):
            n_structural += 1
    frac = n_v_alarm / R
    print(f"\n[a] falsi allarmi Cramér V: {n_v_alarm}/{R} = {frac:.3f} "
          f"(banda {V_FALSE_ALARM_BAND}); repliche con fatti strutturali: {n_structural}/{R}")
    assert V_FALSE_ALARM_BAND[0] <= frac <= V_FALSE_ALARM_BAND[1]
    assert n_structural == 0


# --------------------------------------------------------------------------- #
# (b) annidamento iniettato
# --------------------------------------------------------------------------- #
def test_injected_nesting_always_detected():
    n_detected = 0
    for rep in range(R):
        sheet = make_design_sheet_nested(seed=20_000 + rep)
        res = run_design_audit(sheet, patient_col="patient",
                               technical_cols={"batch": "batch", "library": "library"},
                               outcome_cols=["tissue"])
        if any(f.kind == "annidamento" and f.factors == ("patient", "batch") for f in res.findings):
            n_detected += 1
    print(f"\n[b] annidamento paziente-in-batch rilevato: {n_detected}/{R}")
    assert n_detected == R


# --------------------------------------------------------------------------- #
# (c) confondimento completo paziente-esito
# --------------------------------------------------------------------------- #
def test_complete_patient_outcome_confounding_always_detected():
    n_detected, n_not_estimable = 0, 0
    for rep in range(R):
        sheet = make_design_sheet_outcome_confounded(seed=30_000 + rep)
        res = run_design_audit(sheet, patient_col="patient", tissue_col="tissue",
                               technical_cols={"batch": "batch", "library": "library"},
                               comparisons=[("tissue", "Tumor", "Normal")])
        if any(f.kind == "esito-determinato" and f.factors == ("patient", "tissue")
               for f in res.findings):
            n_detected += 1
        if res.comparisons[0].classification == "non stimabile":
            n_not_estimable += 1
    print(f"\n[c] esito determinato dal paziente rilevato: {n_detected}/{R}; "
          f"confronto tessuto classificato 'non stimabile': {n_not_estimable}/{R}")
    assert n_detected == R
    assert n_not_estimable == R


# --------------------------------------------------------------------------- #
# Potenza del Cramér V (curva, senza banda)
# --------------------------------------------------------------------------- #
def test_cramer_v_power_curve():
    rates = {}
    for v in (0.3, 0.5, 0.7, 0.9):
        n_det = 0
        for rep in range(R):
            sheet = make_design_sheet_association(v, seed=40_000 + rep)
            res = run_design_audit(sheet, patient_col="library",
                                   outcome_cols=["factor_a", "factor_b"])
            pair = next(p for p in res.pairs if {p.factor_a, p.factor_b} == {"factor_a", "factor_b"})
            structural = pair.a_nested_in_b or pair.b_nested_in_a
            if pair.v_alarm or structural:
                n_det += 1
        rates[v] = n_det / R
    print(f"\n[potenza] detection rate (soglia V >= {CRAMER_V_THRESHOLD}, n=60, 3x3): "
          + ", ".join(f"V={v}: {r:.3f}" for v, r in rates.items()))
    vals = list(rates.values())
    assert all(b >= a for a, b in zip(vals, vals[1:]))


# --------------------------------------------------------------------------- #
# Scomposizione della varianza
# --------------------------------------------------------------------------- #
def test_variance_shares_recover_known_components():
    est = {f: [] for f in TRUE_SHARES}
    for rep in range(R):
        Y, units = make_variance_units(seed=50_000 + rep)
        shares = variance_shares(Y, units[["patient", "tissue"]], fixed_cols=["tissue"])
        assert shares is not None
        for f in TRUE_SHARES:
            est[f].append(shares[f])
    for f, truth in TRUE_SHARES.items():
        e = np.array(est[f])
        within = float(np.mean(np.abs(e - truth) <= SHARE_TOL))
        bias = float(e.mean() - truth)
        print(f"\n[quote] {f}: vero {truth}, media stime {e.mean():.3f} (bias {bias:+.3f}), "
              f"entro +-{SHARE_TOL}: {within:.3f}")
        assert within >= SHARE_FRACTION_WITHIN
        assert abs(bias) <= SHARE_BIAS_TOL


@functools.lru_cache(maxsize=1)
def _coverage_counts() -> tuple[dict[str, int], int]:
    covered = {f: 0 for f in TRUE_SHARES}
    n_valid = 0
    for rep in range(R):
        Y, units = make_variance_units(seed=60_000 + rep)
        dec = variance_decomposition_from_units(
            Y, units, patient_col="patient", factor_cols=["patient", "tissue"],
            fixed_cols=["tissue"], n_boot=200, seed=rep)
        assert dec.identifiable
        n_valid += 1
        for s in dec.shares:
            if s.ci_low <= TRUE_SHARES[s.factor] <= s.ci_high:
                covered[s.factor] += 1
    return covered, n_valid


def _assert_coverage(factor: str) -> None:
    covered, n_valid = _coverage_counts()
    cov = covered[factor] / n_valid
    print(f"\n[copertura IC95%] {factor}: {covered[factor]}/{n_valid} = {cov:.3f} (banda {COVERAGE_BAND})")
    assert COVERAGE_BAND[0] <= cov <= COVERAGE_BAND[1]


def test_variance_share_bootstrap_coverage_patient():
    _assert_coverage("patient")


@pytest.mark.xfail(
    strict=True,
    reason=(
        "Copertura dell'IC95% della quota del TESSUTO misurata 200/200 = 1.000, sopra la banda "
        "dichiarata [0.90, 0.99]: intervallo conservativo (troppo largo). Metodo già corretto "
        "due volte (1: divisore per effetti casuali + intervallo 'basic', copertura 0/200 -> "
        "stimatore puntuale non distorto ma intervallo ancora spostato; 2: intervallo normale "
        "stima +- 1.96 SE bootstrap -> paziente 0.965 in banda, tessuto 1.000 fuori). La banda "
        "NON è stata allargata. Conseguenza: la scomposizione della varianza non è esposta in "
        "CLI, app o report. Probabile causa: il bootstrap sui pazienti propaga all'intervallo "
        "del tessuto (effetto fisso) anche la variabilità del denominatore dovuta al paziente."
    ),
)
def test_variance_share_bootstrap_coverage_tissue():
    _assert_coverage("tissue")


def test_pseudobulk_end_to_end_smoke():
    """Non è una calibrazione: verifica che la pseudobulk da conteggi per cellula
    produca quote nell'ordine atteso (paziente > tessuto > 0)."""
    adata = make_pseudobulk_adata(seed=3)
    dec = variance_decomposition(adata, patient_col="patient_id", tissue_col="tissue",
                                 n_hvg=300, n_boot=100, seed=0)
    assert dec.identifiable
    sh = {s.factor: s.share for s in dec.shares}
    print(f"\n[pseudobulk] quote stimate: {sh}")
    assert sh["patient_id"] > sh["tissue"] > 0.05


def test_variance_decomposition_refuses_confounded_factors():
    Y, units = make_variance_units(seed=1)
    # 'library' coincide con la coppia paziente-tessuto: nessun grado di libertà residuo
    units = units.assign(library=units.patient + "-" + units.tissue)
    dec = variance_decomposition_from_units(
        Y, units, patient_col="patient", factor_cols=["patient", "tissue", "library"],
        n_boot=50, seed=0)
    assert not dec.identifiable
    assert dec.shares == []
    assert "non identificabile" in dec.sentence
    # batch che coincide con il paziente: confuso
    units2 = units.assign(batch=units.patient.map(lambda p: f"B{p}"))
    dec2 = variance_decomposition_from_units(
        Y, units2, patient_col="patient", factor_cols=["patient", "tissue", "batch"],
        n_boot=50, seed=0)
    assert not dec2.identifiable
    assert dec2.shares == []


def test_variance_decomposition_refuses_too_few_patients():
    Y, units = make_variance_units(n_patients=4, seed=2)
    dec = variance_decomposition_from_units(
        Y, units, patient_col="patient", factor_cols=["patient", "tissue"], n_boot=50, seed=0)
    assert not dec.identifiable
    assert dec.shares == []


# --------------------------------------------------------------------------- #
# Caso di riferimento GSE278694
# --------------------------------------------------------------------------- #
def test_gse278694_reference_design():
    sheet = make_gse278694_like_sheet()
    res = run_design_audit(
        sheet, patient_col="patient", tissue_col="tissue",
        technical_cols={"protocol": "protocol", "library": "library"},
        comparisons=[("tissue", "Tumor", "Adjacent_normal"), ("protocol", "scRNA", "snRNA")],
    )
    tum_adj, sc_sn = res.comparisons
    print(f"\n[GSE278694] {tum_adj.sentence}\n[GSE278694] {sc_sn.sentence}")
    assert tum_adj.n_units == 5
    assert tum_adj.classification == "stimabile con bassa potenza"
    assert sc_sn.n_units == 0
    assert sc_sn.classification == "non stimabile"
    # libreria = coppia paziente-tessuto; paziente annidato nel protocollo (coorti disgiunte)
    kinds = {(f.kind, f.factors) for f in res.findings}
    assert ("unità-tecnica", ("library",)) in kinds
    assert ("annidamento", ("patient", "protocol")) in kinds
