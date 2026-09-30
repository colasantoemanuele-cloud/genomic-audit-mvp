"""Flag per cellula dal Modulo B (Intervento 2).

Requisiti verificati:
  - Coerenza: la media di audit_label_vs_reference sulle cellule valutabili di un
    compartimento coincide ESATTAMENTE con la stima puntuale di marker_error_rate.
  - Regressione: gli output di run_tcr_validation sui fixture sintetici sono identici
    (tolleranza 1e-12) a quelli prodotti dal codice PRIMA della modifica
    (tests/fixtures/tcr_regression_baseline.json, generato sul commit di baseline).
  - Le colonne originali dell'AnnData restano invariate; i flag stanno sulle righe giuste
    anche con l'ordine delle righe rimescolato.
  - NA dove la cellula non era verificabile (nessun TCR, clone senza riferimento,
    etichetta fuori dalla marker_map, compartimento di riferimento).
"""

from __future__ import annotations

import json

import anndata as ad
import numpy as np
import pandas as pd
import pytest

from core.synthetic import make_tcr_validation_dataset
from core.tcr_validation import export_audited, run_tcr_validation
from tests.fixtures.tcr_regression import BASELINE_PATH, CONFIGS, run_config, summarize

MARKERS = {"CD4T": ["CD4"], "CD8T": ["CD8A", "CD8B"]}
TOL = 1e-12


def _run(adata, contigs, **kw):
    return run_tcr_validation(adata, contigs, patient_col="patient_id", compartment_col="tissue",
                              celltype_col="celltype", barcode_col="barcode", n_boot=200, seed=0,
                              marker_map=MARKERS, reference_compartment="PBMC", **kw)


def _assert_close(a, b, path="") -> None:
    if isinstance(a, dict):
        assert set(a) == set(b), f"{path}: chiavi diverse {set(a) ^ set(b)}"
        for k in a:
            _assert_close(a[k], b[k], f"{path}.{k}")
    elif isinstance(a, list):
        assert len(a) == len(b), f"{path}: lunghezze diverse"
        for i, (x, y) in enumerate(zip(a, b)):
            _assert_close(x, y, f"{path}[{i}]")
    elif isinstance(a, float) and isinstance(b, (int, float)):
        assert abs(a - b) <= TOL, f"{path}: {a!r} != {b!r}"
    else:
        assert a == b, f"{path}: {a!r} != {b!r}"


@pytest.mark.parametrize("name", list(CONFIGS))
def test_regression_identical_to_baseline(name):
    baseline = json.loads(BASELINE_PATH.read_text())[name]
    current = json.loads(json.dumps(summarize(run_config(name)[2]), ensure_ascii=False, default=int))
    _assert_close(baseline, current, name)


def test_flags_mean_equals_marker_error_rate_exactly():
    for seed in range(5):
        adata, contigs = make_tcr_validation_dataset(
            n_patients=10, n_clones_per_patient=12, injected_excess=0.3,
            compartments=("PBMC", "Adjacent", "Tumor"), seed=seed)
        res = _run(adata, contigs)
        flags = res.cell_flags
        comp = adata.obs["tissue"].astype(str).values
        for c, br in res.marker_error.by_compartment.items():
            v = flags.loc[comp == c, "audit_label_vs_reference"].dropna().astype(bool)
            assert len(v) == br.n_obs
            assert v.mean() == br.mean, f"seed {seed}, {c}: {v.mean()!r} != {br.mean!r}"


def test_reference_compartment_and_unverifiable_cells_are_na():
    adata, contigs = make_tcr_validation_dataset(n_patients=8, n_clones_per_patient=10, seed=3)
    # meta' delle cellule senza TCR: si tolgono i loro contig
    drop = set(adata.obs["barcode"].iloc[::2])
    contigs = contigs[~contigs["barcode"].isin(drop)]
    res = _run(adata, contigs)
    f = res.cell_flags
    obs = adata.obs
    assert f["audit_label_vs_reference"][obs["tissue"].values == "PBMC"].isna().all()
    no_tcr = obs["barcode"].isin(drop).values
    assert f.loc[no_tcr].isna().all().all()
    # etichetta fuori dalla marker_map (il distrattore 'NK'): mai True/False
    nk = obs["celltype"].values == "NK"
    assert f["audit_label_vs_reference"][nk].isna().all()
    assert 0 < res.flag_coverage < 1


def test_flags_follow_rows_when_order_is_shuffled():
    adata, contigs = make_tcr_validation_dataset(n_patients=8, n_clones_per_patient=10, seed=5)
    adata.obs_names = [f"cell{i}" for i in range(adata.n_obs)]
    ref = _run(adata, contigs).cell_flags
    perm = np.random.default_rng(0).permutation(adata.n_obs)
    shuffled = adata[perm].copy()
    got = _run(shuffled, contigs).cell_flags
    assert list(got.index) == list(shuffled.obs_names)
    pd.testing.assert_frame_equal(got.loc[ref.index], ref)


def test_export_writes_copy_and_leaves_original_untouched(tmp_path):
    adata, contigs = make_tcr_validation_dataset(n_patients=6, n_clones_per_patient=8, seed=6)
    adata.obs_names = [f"cell{i}" for i in range(adata.n_obs)]
    obs_before = adata.obs.copy()
    x_before = adata.X.copy()
    res = _run(adata, contigs)
    h5ad_path, csv_path = export_audited(adata, res, tmp_path / "dati")
    assert h5ad_path.name == "dati_audited.h5ad"
    # l'AnnData in memoria non e' stato toccato
    pd.testing.assert_frame_equal(adata.obs, obs_before)
    assert (adata.X != x_before).nnz == 0
    # la copia ha le colonne originali identiche + le due colonne di flag
    out = ad.read_h5ad(h5ad_path)
    for col in obs_before.columns:
        assert (out.obs[col].astype(str).values == obs_before[col].astype(str).values).all()
    assert list(out.obs.columns) == list(obs_before.columns) + [
        "audit_reference_label", "audit_label_vs_reference"]
    assert "genomic_audit_flags" in out.uns
    csv = pd.read_csv(csv_path, index_col=0)
    assert list(csv.index.astype(str)) == list(adata.obs_names)
    assert csv["audit_label_vs_reference"].isna().sum() == res.cell_flags[
        "audit_label_vs_reference"].isna().sum()


def test_export_refuses_existing_flag_columns(tmp_path):
    adata, contigs = make_tcr_validation_dataset(n_patients=6, n_clones_per_patient=8, seed=7)
    res = _run(adata, contigs)
    adata.obs["audit_reference_label"] = "x"
    with pytest.raises(ValueError):
        export_audited(adata, res, tmp_path / "dati")


def test_no_flags_without_marker_map():
    adata, contigs = make_tcr_validation_dataset(n_patients=6, n_clones_per_patient=8, seed=8)
    res = run_tcr_validation(adata, contigs, patient_col="patient_id", compartment_col="tissue",
                             celltype_col="celltype", barcode_col="barcode", n_boot=100, seed=0)
    assert res.cell_flags is None
