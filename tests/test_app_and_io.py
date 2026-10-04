"""Test di lettura formati, sintesi dei controlli, report Markdown, tabella TCR per etichetta e web app."""

from __future__ import annotations

import gzip
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
import scipy.io
import scipy.sparse as sp

from core.design_audit import run_design_audit
from core.io import load_matrix_market, read_table
from core.report import render_markdown_report
from core.synthetic import make_gse278694_like_sheet, make_leakage_dataset, make_tcr_validation_dataset
from core.tcr_validation import run_tcr_validation, tcr_by_celltype
from core.verdict import GIALLO, GRIGIO, ROSSO, VERDE, design_summary, leakage_summary, tcr_summary

ROOT = Path(__file__).resolve().parents[1]


def _write_10x(folder: Path, gz: bool) -> pd.DataFrame:
    rng = np.random.default_rng(0)
    X = sp.random(30, 12, density=0.3, format="coo", random_state=0, data_rvs=lambda k: rng.integers(1, 9, k))
    op = (lambda p: gzip.open(str(p) + ".gz", "wb")) if gz else (lambda p: open(p, "wb"))
    with op(folder / "matrix.mtx") as fh:
        scipy.io.mmwrite(fh, X)  # geni x cellule (30 geni, 12 cellule)
    bcs = [f"BC{i:02d}-1" for i in range(12)]
    with op(folder / "barcodes.tsv") as fh:
        fh.write(("\n".join(bcs) + "\n").encode())
    with op(folder / "features.tsv") as fh:
        fh.write(("\n".join(f"ENSG{i}\tGENE{i % 25}\tGene Expression" for i in range(30)) + "\n").encode())
    return pd.DataFrame({"cell": bcs[::-1], "patient": ["P1", "P2"] * 6})


@pytest.mark.parametrize("gz", [False, True])
def test_matrix_market_loader_aligns_metadata_and_dedups_genes(tmp_path, gz):
    meta = _write_10x(tmp_path, gz)
    a = load_matrix_market(tmp_path, meta, "cell")
    assert a.shape == (12, 30) and sp.issparse(a.X)
    assert list(a.obs_names) == [f"BC{i:02d}-1" for i in range(12)]
    assert (a.obs["cell"].values == a.obs_names).all()
    assert a.var_names.is_unique


def test_matrix_market_loader_refuses_missing_metadata(tmp_path):
    meta = _write_10x(tmp_path, False).iloc[:5]
    with pytest.raises(ValueError, match="non hanno metadati"):
        load_matrix_market(tmp_path, meta, "cell")


def test_read_table_detects_tsv_and_csv(tmp_path):
    (tmp_path / "a.tsv").write_text("x\ty\n1\t2\n")
    (tmp_path / "a.csv").write_text("x,y\n1,2\n")
    for f in ("a.tsv", "a.csv"):
        t = read_table(tmp_path / f)
        assert list(t.columns) == ["x", "y"] and t.iloc[0, 0] == "1"
    assert list(read_table(b"x\ty\n1\t2\n").columns) == ["x", "y"]


def _state(summary, name_part: str) -> str:
    return next(c.state for c in summary.checks if name_part in c.name)


def test_summary_rules_on_real_module_outputs():
    sheet = make_gse278694_like_sheet()
    d = run_design_audit(sheet, "patient", "tissue", {"protocol": "protocol"},
                         comparisons=[("protocol", "scRNA", "snRNA")])
    assert _state(design_summary(d), "protocol") == ROSSO  # confronto non stimabile
    d2 = run_design_audit(sheet, "patient", "tissue", comparisons=[("tissue", "Tumor", "Adjacent_normal")])
    assert _state(design_summary(d2), "tissue") == GIALLO  # 5 pazienti: bassa potenza
    from core.leakage_audit import run_leakage_audit
    leak = run_leakage_audit(make_leakage_dataset(n_patients=8, cells_per_patient=15, n_genes=150, seed=6),
                             "label", "patient_id", benchmark=False)
    assert _state(leakage_summary(leak), "Confronto fra modelli") == GRIGIO  # non eseguito: mai verde
    adata, contigs = make_tcr_validation_dataset(n_patients=12, injected_excess=0.35, seed=0)
    t = run_tcr_validation(adata, contigs, "patient_id", "tissue", "celltype", "barcode", n_boot=200,
                           marker_map={"CD4T": ["CD4"], "CD8T": ["CD8A", "CD8B"]}, reference_compartment="PBMC")
    assert _state(tcr_summary(t), "discordanza") == GIALLO  # intervallo sopra lo zero: stima con limiti
    assert _state(tcr_summary(t), "[cell] in Tumor") == VERDE  # tasso con intervallo prodotto
    adata0, contigs0 = make_tcr_validation_dataset(n_patients=3, seed=0)
    t0 = run_tcr_validation(adata0, contigs0, "patient_id", "tissue", "celltype", "barcode", n_boot=100)
    assert _state(tcr_summary(t0), "discordanza") == ROSSO  # 3 pazienti: stima non possibile
    assert _state(tcr_summary(t0), "Tasso d'errore") == GRIGIO  # non calcolato senza mappa dei marcatori


def test_markdown_report_contains_verdicts_definitions_and_limits():
    d = run_design_audit(make_gse278694_like_sheet(), "patient", "tissue",
                         comparisons=[("tissue", "Tumor", "Adjacent_normal")])
    md = render_markdown_report(design_result=d, dataset_name="prova")
    assert "## Sintesi dei controlli" in md and "giallo: stima con limiti" in md and "Regola Disegno" in md
    assert "Verdetto" not in md
    assert "## Limiti dichiarati" in md and "un solo dataset reale" in md


def test_tcr_by_celltype_is_descriptive_table():
    adata, contigs = make_tcr_validation_dataset(n_patients=4, seed=1)
    tab = tcr_by_celltype(adata, contigs, "patient_id", "tissue", "celltype", "barcode")
    assert set(tab.index) == {"CD4T", "CD8T", "NK"}
    assert (tab["con_TRB"] <= tab["cellule"]).all()
    assert bool(tab.loc["NK", "etichetta T"]) is False and bool(tab.loc["CD8T", "etichetta T"]) is True


def test_demo_data_present_and_raw_counts():
    import anndata as ad
    a = ad.read_h5ad(ROOT / "data" / "demo" / "GSE125449_demo.h5ad")
    assert a.obs["patient"].nunique() == 10 and a.n_obs == 1861
    assert np.allclose(a.X.data, np.round(a.X.data))
    for f in ("GSE125449_samples.csv", "GSE132465_samples.csv", "GSE131907_samples.csv"):
        assert (ROOT / "data" / "demo" / f).exists()


def test_streamlit_config_keeps_everything_local():
    cfg = (ROOT / ".streamlit" / "config.toml").read_text()
    assert "gatherUsageStats = false" in cfg and 'address = "localhost"' in cfg


def test_web_app_demo_runs_end_to_end():
    from streamlit.testing.v1 import AppTest
    at = AppTest.from_file(str(ROOT / "app.py"), default_timeout=600).run()
    at.checkbox[0].check().run()  # Modulo A rapido, per contenere i tempi del test
    next(b for b in at.button if b.label == "Carica la demo").click().run()
    assert not at.exception, at.exception
    text = " ".join(m.value for m in at.markdown)
    assert "Confronti: che cosa il disegno permette di stimare" in " ".join(s.value for s in at.subheader)
    assert "badge b-" in text
    assert any("Sintesi dei controlli" in t.label for t in at.tabs)
    assert not any("Verdetto" in t.label for t in at.tabs)
