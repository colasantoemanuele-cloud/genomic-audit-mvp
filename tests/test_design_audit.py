"""Test funzionali dell'audit del disegno (casi limite, non calibrazione: quella e' in
test_design_audit_calibration.py)."""

from __future__ import annotations

import numpy as np
import pandas as pd

from core.design_audit import (
    assess_comparison,
    cramers_v_bergsma,
    run_design_audit,
    sample_sheet_from_obs,
)


def test_cramers_v_perfect_and_null():
    assert cramers_v_bergsma(np.array([[50, 0], [0, 50]])) > 0.95
    assert cramers_v_bergsma(np.array([[25, 25], [25, 25]])) == 0.0
    assert cramers_v_bergsma(np.array([[10, 0]])) is None


def test_small_table_is_not_evaluable():
    sheet = pd.DataFrame({"patient": ["P1", "P2", "P3", "P4"], "batch": ["a", "b", "a", "b"],
                          "tissue": ["T", "T", "N", "N"]})
    res = run_design_audit(sheet, patient_col="patient", technical_cols={"batch": "batch"},
                           outcome_cols=["tissue"])
    pair = next(p for p in res.pairs if {p.factor_a, p.factor_b} == {"batch", "tissue"})
    assert pair.cramer_v is None
    assert "non valutabile" in pair.sentence


def test_batch_perfectly_separating_tissue_makes_comparison_not_estimable():
    rows = []
    for p in range(8):
        rows.append({"patient": f"P{p}", "tissue": "Tumor", "batch": "run1", "library": f"L{p}T"})
        rows.append({"patient": f"P{p}", "tissue": "Normal", "batch": "run2", "library": f"L{p}N"})
    sheet = pd.DataFrame(rows)
    res = run_design_audit(sheet, patient_col="patient", tissue_col="tissue",
                           technical_cols={"batch": "batch", "library": "library"},
                           comparisons=[("tissue", "Tumor", "Normal")])
    c = res.comparisons[0]
    assert c.n_units == 8
    assert c.classification == "non stimabile"
    assert "batch" in c.sentence
    assert any(f.kind == "uno-a-uno" and set(f.factors) == {"tissue", "batch"} for f in res.findings)


def test_paired_comparison_with_enough_patients_is_estimable():
    rows = []
    rng = np.random.default_rng(0)
    for p in range(8):
        for t in ("Tumor", "Normal"):
            rows.append({"patient": f"P{p}", "tissue": t, "batch": f"run{rng.integers(0, 3)}"})
    res = run_design_audit(pd.DataFrame(rows), patient_col="patient", tissue_col="tissue",
                           technical_cols={"batch": "batch"},
                           comparisons=[("tissue", "Tumor", "Normal")])
    c = res.comparisons[0]
    assert c.design.startswith("appaiato")
    assert c.n_units == 8
    assert c.classification == "stimabile"


def test_patient_level_outcome_is_between_patient_comparison():
    rows = []
    for p in range(12):
        resp = "R" if p < 6 else "NR"
        for t in ("Tumor", "PBMC"):
            rows.append({"patient": f"P{p}", "tissue": t, "response": resp})
    res = run_design_audit(pd.DataFrame(rows), patient_col="patient", tissue_col="tissue",
                           outcome_cols=["response"], comparisons=[("response", "R", "NR")])
    c = res.comparisons[0]
    assert c.design == "fra pazienti"
    assert c.n_patients_a == 6 and c.n_patients_b == 6
    assert c.classification == "stimabile"
    assert any(f.kind == "esito-determinato" and f.factors == ("patient", "response")
               for f in res.findings)


def test_missing_level_is_not_estimable():
    sheet = pd.DataFrame({"patient": ["P1", "P2"], "tissue": ["Tumor", "Tumor"]})
    c = assess_comparison(sheet, "tissue", "Tumor", "Normal", "patient",
                          roles={"patient": "patient", "tissue": "tissue"})
    assert c.classification == "non stimabile"
    assert "non compare" in c.sentence


def test_sample_sheet_from_cell_level_obs():
    obs = pd.DataFrame({"patient": ["P1"] * 5 + ["P2"] * 3, "tissue": ["T"] * 8,
                        "celltype": list("abcdeabc")})
    sheet = sample_sheet_from_obs(obs, ["patient", "tissue"])
    assert len(sheet) == 2


def test_protocol_note_present_only_when_protocol_declared():
    sheet = pd.DataFrame({"patient": ["P1", "P2", "P3"], "protocol": ["scRNA", "snRNA", "scRNA"]})
    with_p = run_design_audit(sheet, patient_col="patient", technical_cols={"protocol": "protocol"})
    without = run_design_audit(sheet, patient_col="patient")
    assert any("snRNA" in n for n in with_p.notes)
    assert not any("snRNA" in n for n in without.notes)


def test_every_comparison_reports_units_and_min_pvalue_with_fixed_note():
    """B2: la classe non cambia; il testo riporta unita', p-value minimo e la frase fissa."""
    from core.design_audit import MIN_PVALUE_NOTE
    rows = []
    for p in range(8):
        for t in ("Tumor", "Normal"):
            rows.append({"patient": f"P{p}", "tissue": t, "resp": "R" if p < 4 else "NR"})
    sheet = pd.DataFrame(rows)
    res = run_design_audit(sheet, patient_col="patient", tissue_col="tissue", outcome_cols=["resp"],
                           comparisons=[("tissue", "Tumor", "Normal"), ("resp", "R", "NR")])
    paired, between = res.comparisons
    assert paired.classification == "stimabile"  # 8 coppie: classe invariata
    assert paired.min_pvalue == 2 / 2 ** 8 and paired.min_pvalue_test == "Wilcoxon appaiato"
    assert between.classification == "stimabile con bassa potenza"  # 4 contro 4 < 5
    assert abs(between.min_pvalue - 2 / 70) < 1e-12 and between.min_pvalue_test == "Mann-Whitney"
    for c in (paired, between):
        assert MIN_PVALUE_NOTE in c.sentence
        assert f"Unita' indipendenti: {c.n_units}" in c.sentence
        assert f"{c.min_pvalue:.3f}" in c.sentence


def test_not_estimable_without_units_says_pvalue_undefined():
    sheet = pd.DataFrame({"patient": ["P1", "P2", "P3"], "protocol": ["sc", "sn", "sc"]})
    c = assess_comparison(sheet, "protocol", "sc", "sn", "patient",
                          roles={"patient": "patient", "protocol": "protocol"})
    assert c.classification == "non stimabile" and c.min_pvalue is None
    assert "p-value minimo non definito" in c.sentence


def test_missing_values_are_an_explicit_level_not_silently_dropped():
    """Regressione (validazione esterna, GSE132465): lo stadio manca per i campioni normali.
    Prima i mancanti venivano scartati da crosstab/groupby ma contati in n: l'audit riportava
    'paziente determina lo stadio' (falso: i pazienti appaiati hanno stadio e mancante) e non
    produceva il Cramér V per una tabella valutabile."""
    import io
    rows = []
    for p in range(12):
        rows.append({"patient": f"P{p}", "tissue": "Tumor", "stage": str(1 + p % 3)})
        if p < 6:
            rows.append({"patient": f"P{p}", "tissue": "Normal", "stage": None})
    csv = pd.DataFrame(rows).to_csv(index=False)
    sheet = pd.read_csv(io.StringIO(csv), dtype=str)  # come la CLI
    assert sheet["stage"].isna().sum() == 6
    res = run_design_audit(sheet, patient_col="patient", tissue_col="tissue", outcome_cols=["stage"])
    facts = {(f.kind, f.factors) for f in res.findings}
    assert ("esito-determinato", ("patient", "stage")) not in facts
    pair = next(p for p in res.pairs if {p.factor_a, p.factor_b} == {"tissue", "stage"})
    assert pair.n_levels_b == 4 or pair.n_levels_a == 4  # 3 stadi + "NA"
    assert pair.cramer_v is not None
    assert any("\"NA\"" in n and "'stage' (6 righe)" in n for n in res.notes)


def test_strong_v_with_structural_relation_is_not_described_as_below_threshold():
    """Regressione (validazione esterna, GSE132465): V = 0.95 con un annidamento veniva
    descritto come 'sotto la soglia di 0.5'."""
    rows = []
    for p in range(12):
        rows.append({"patient": f"P{p}", "tissue": "Tumor", "stage": str(1 + p % 3)})
        if p < 6:
            rows.append({"patient": f"P{p}", "tissue": "Normal", "stage": None})
    res = run_design_audit(pd.DataFrame(rows), patient_col="patient", tissue_col="tissue",
                           outcome_cols=["stage"])
    pair = next(p for p in res.pairs if {p.factor_a, p.factor_b} == {"tissue", "stage"})
    assert pair.cramer_v is not None and pair.cramer_v >= 0.5
    assert pair.a_nested_in_b or pair.b_nested_in_a
    assert "sotto la soglia" not in pair.sentence
    assert "strutturale" in pair.sentence
