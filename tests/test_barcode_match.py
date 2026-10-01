"""Pulizia 2 -- controllo del match dei barcode VDJ <-> metadati."""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from core.synthetic import make_tcr_validation_dataset
from core.tcr_validation import run_tcr_validation

MARKERS = {"CD4T": ["CD4"], "CD8T": ["CD8A", "CD8B"]}


def _run(adata, contigs):
    return run_tcr_validation(adata, contigs, patient_col="patient_id", compartment_col="tissue",
                              celltype_col="celltype", barcode_col="barcode", n_boot=200, seed=0,
                              marker_map=MARKERS, reference_compartment="PBMC")


def _numbers(res):
    t = res.marker_error.by_compartment["Tumor"]
    return (res.n_cells_with_tcr, res.discordance.mean_excess, res.discordance.ci_low, t.mean, t.ci_high)


def test_same_suffix_on_both_sides_no_normalization():
    adata, contigs = make_tcr_validation_dataset(n_patients=6, n_clones_per_patient=8, seed=1)
    ref = _run(adata, contigs)
    adata.obs["barcode"] = adata.obs["barcode"].astype(str) + "-1"
    contigs = contigs.assign(barcode=contigs["barcode"] + "-1")
    res = _run(adata, contigs)
    print("\n" + res.barcode_match.sentence)
    assert res.barcode_match.fraction == 1.0
    assert res.barcode_match.normalized_suffixes == ()
    assert "Normalizzato" not in res.barcode_match.sentence
    assert _numbers(res) == _numbers(ref)


def test_suffix_missing_on_one_side_is_normalized_and_declared():
    adata, contigs = make_tcr_validation_dataset(n_patients=6, n_clones_per_patient=8, seed=1)
    ref = _run(adata, contigs)
    contigs = contigs.assign(barcode=contigs["barcode"] + "-1")  # solo il lato VDJ ha il suffisso
    res = _run(adata, contigs)
    print("\n" + res.barcode_match.sentence)
    assert res.barcode_match.fraction == 1.0
    assert res.barcode_match.normalized_suffixes == ("-1",)
    assert "Normalizzato il suffisso -1" in res.barcode_match.sentence
    assert _numbers(res) == _numbers(ref)
    # i flag restano sulle righe giuste
    assert res.cell_flags.equals(ref.cell_flags)


def test_normalization_refused_when_it_would_merge_distinct_barcodes():
    adata, contigs = make_tcr_validation_dataset(n_patients=6, n_clones_per_patient=8, seed=1)
    # tutti i barcode VDJ con -1 e due cellule duplicate con -2: togliere il suffisso
    # fonderebbe barcode diversi, quindi la normalizzazione non va applicata e il
    # match resta allo 0%
    dup = contigs.iloc[:2].assign(barcode=lambda d: d["barcode"] + "-2")
    contigs = contigs.assign(barcode=contigs["barcode"] + "-1")
    contigs = pd.concat([contigs, dup], ignore_index=True)
    with pytest.raises(ValueError, match="ritrovato nei metadati"):
        _run(adata, contigs)


def test_no_match_raises_clear_error_with_examples():
    adata, contigs = make_tcr_validation_dataset(n_patients=6, n_clones_per_patient=8, seed=1)
    rng = np.random.default_rng(0)
    contigs = contigs.assign(barcode=[f"ZZ{rng.integers(1e9)}" for _ in range(len(contigs))])
    with pytest.raises(ValueError) as e:
        _run(adata, contigs)
    msg = str(e.value)
    print("\n" + msg)
    assert "0.0%" in msg and "Esempi VDJ" in msg and "esempi metadati" in msg
