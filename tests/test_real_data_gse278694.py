"""Test sui dati reali di GSE278694 (cartella pdac-ml, sola lettura).

Saltati con motivo se i dati non sono presenti. Verificano che la convenzione "per clone"
del Modulo B riproduca esattamente i numeri della tesi (10_loco.py, schema blood):
Tumor 0.1951 [0.1230; 0.2737], Adjacent_normal 0.0593, differenza Tumor - Adjacent
+0.0795 su 88 cloni. Il test carica ~6 GB in memoria (circa 1 minuto).
"""

from __future__ import annotations

import re
from pathlib import Path

import anndata as ad
import pandas as pd
import pytest

from core.tcr_validation import parse_vdj_contigs, run_tcr_validation

PDAC = Path("/home/nemo/Uni/PROGETTO_ML/pdac-ml")
H5AD = PDAC / "data/interim/sc_raw.h5ad"
TCR = PDAC / "data/raw/tcr"
TISSUE = {"Tumor": "Tumor", "Normal": "Adjacent_normal", "PBMC": "PBMC"}

pytestmark = pytest.mark.skipif(not (H5AD.exists() and TCR.exists()),
                                reason="dati reali di GSE278694 (pdac-ml) non presenti su questa macchina")


def _vdj_files(tmp: Path, strip_suffix: bool) -> list[tuple[Path, str, str]]:
    files = []
    for p in sorted(TCR.glob("*_filtered_contig_annotations.csv.gz")):
        m = re.match(r"GSM\d+_PA_(\d+)-(\w+)-VDJ_filtered_contig_annotations\.csv\.gz", p.name)
        patient, comp = f"PA{int(m.group(1)):02d}", TISSUE[m.group(2)]
        if strip_suffix:
            df = pd.read_csv(p, low_memory=False)
            df["barcode"] = df["barcode"].astype(str).str.replace(r"-\d+$", "", regex=True)
            q = tmp / p.name.replace(".csv.gz", ".csv")
            df.to_csv(q, index=False)
            files.append((q, patient, comp))
        else:
            files.append((p, patient, comp))
    return files


@pytest.fixture(scope="module")
def real(tmp_path_factory):
    adata = ad.read_h5ad(H5AD)
    contigs = parse_vdj_contigs(_vdj_files(tmp_path_factory.mktemp("vdj"), strip_suffix=True))
    res = run_tcr_validation(
        adata, contigs, patient_col="patients", compartment_col="tissue",
        celltype_col="all_celltype", barcode_col="barcode", n_boot=2000, seed=0,
        marker_map={"CD4T": ["CD4"], "CD8T": ["CD8A", "CD8B"]}, reference_compartment="PBMC",
        convention="both", clone_error_labels=["NK"], clone_marker_priority=["CD8T", "CD4T"])
    return adata, res


def test_clone_convention_reproduces_thesis_exactly(real):
    _, res = real
    clone = res.conventions["clone"]
    t, a = clone.by_compartment["Tumor"], clone.by_compartment["Adjacent_normal"]
    d, n_cl = clone.paired_differences[("Tumor", "Adjacent_normal")]
    print(f"\n[reale, per clone] Tumor {t.mean:.4f} [{t.ci_low:.4f}; {t.ci_high:.4f}] {t.n_groups} pz | "
          f"Adjacent {a.mean:.4f} ({a.n_groups} pz) | Tumor-Adj {d.mean:+.4f} su {n_cl} cloni")
    assert abs(t.mean - 0.1951219512195122) < 1e-12
    assert abs(t.ci_low - 0.12300170068027211) < 1e-12 and abs(t.ci_high - 0.27365049812586306) < 1e-12
    assert abs(a.mean - 0.059322033898305086) < 1e-12
    assert abs(d.mean - 0.07954545454545454) < 1e-12 and n_cl == 88


def test_cell_convention_unchanged_on_real_data(real):
    _, res = real
    cell = res.conventions["cell"].by_compartment["Tumor"]
    print(f"\n[reale, per cellula] Tumor {cell.mean:.4f} [{cell.ci_low:.4f}; {cell.ci_high:.4f}]")
    assert abs(cell.mean - 0.0925967851482907) < 1e-12
    assert cell is res.marker_error.by_compartment["Tumor"]
