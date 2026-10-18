"""Test sui dati reali di GSE278694 (cartella pdac-ml, sola lettura).

Saltati con motivo se i dati non sono presenti. Sono test di regressione sui numeri reali del
Modulo B con il codice attuale, con le opzioni della tesi (--clone-error-labels NK
--clone-marker-priority CD8T,CD4T) e con i file VDJ originali (barcode con "-1").

Rispetto ai numeri della tesi (Tumor 0.1951 [0.1230; 0.2737], Adjacent 0.0593, differenza +0.0795
su 88 cloni) la stima per clone è diversa per due regole rese esplicite nella revisione: i pareggi
dell'etichetta di maggioranza (11 unità su 288 nel tumore; prima decisi dall'ordine delle righe) e la scelta
fra contig TRB a pari UMI (346 cellule; prima decisa dall'ordine dell'ordinamento). L'intervallo è
quello cluster-robusto (revisione A1), non più il percentile. La riproduzione della tesi con le sue
regole di allora è in validation/gse278694_moduleB.py. Il test carica ~6 GB (circa 1 minuto).
"""

from __future__ import annotations

import os
from pathlib import Path

import anndata as ad
import pytest

from core.io import read_vdj_manifest
from core.tcr_validation import parse_vdj_contigs, run_tcr_validation

# Dati reali di GSE278694: non inclusi nella repo. Cartella indicata da GSE278694_DIR, altrimenti
# una cartella "pdac-ml" accanto a questa repo. Se i dati mancano, i test vengono saltati.
PDAC = Path(os.environ.get("GSE278694_DIR", Path(__file__).resolve().parents[2] / "pdac-ml"))
H5AD = PDAC / "data/interim/sc_raw.h5ad"
TCR = PDAC / "data/raw/tcr"

pytestmark = pytest.mark.skipif(not (H5AD.exists() and TCR.exists()),
                                reason="dati reali di GSE278694 (pdac-ml) non presenti su questa macchina")


@pytest.fixture(scope="module")
def real(tmp_path_factory):
    from validation.make_vdj_manifest import build_manifest
    man = tmp_path_factory.mktemp("vdj") / "manifest.csv"
    assert build_manifest(man, PDAC) == 30
    adata = ad.read_h5ad(H5AD)
    contigs = parse_vdj_contigs(read_vdj_manifest(man))
    res = run_tcr_validation(
        adata, contigs, patient_col="patients", compartment_col="tissue",
        celltype_col="all_celltype", barcode_col="barcode",
        marker_map={"CD4T": ["CD4"], "CD8T": ["CD8A", "CD8B"]}, reference_compartment="PBMC",
        convention="both", clone_error_labels=["NK"], clone_marker_priority=["CD8T", "CD4T"],
        compartment_order=["Tumor", "Adjacent_normal"])
    return adata, res


def test_clone_convention_on_real_data(real):
    _, res = real
    clone = res.conventions["clone"]
    t, a = clone.by_compartment["Tumor"], clone.by_compartment["Adjacent_normal"]
    d, n_cl = clone.paired_differences[("Tumor", "Adjacent_normal")]
    print(f"\n[reale, per clone] Tumor {t.mean:.4f} [{t.ci_low:.4f}; {t.ci_high:.4f}] {t.n_groups} pz | "
          f"Adjacent {a.mean:.4f} ({a.n_groups} pz) | Tumor-Adj {d.mean:+.4f} su {n_cl} cloni")
    assert abs(t.mean - 0.18576388888888892) < 1e-12 and t.n_obs == 288 and t.n_groups == 10
    assert abs(t.ci_low - 0.11295427591468271) < 1e-9 and abs(t.ci_high - 0.2901543690548171) < 1e-9
    assert abs(a.mean - 0.0635593220338983) < 1e-12 and not a.sufficient and a.n_groups == 4
    assert abs(d.mean - 0.06439393939393939) < 1e-12 and n_cl == 88 and not d.sufficient


def test_cell_convention_and_discordance_on_real_data(real):
    _, res = real
    cell = res.conventions["cell"].by_compartment["Tumor"]
    print(f"\n[reale, per cellula] Tumor {cell.mean:.4f} [{cell.ci_low:.4f}; {cell.ci_high:.4f}]")
    assert abs(cell.mean - 0.0930285196921684) < 1e-12 and cell.n_obs == 4418
    assert abs(cell.ci_low - 0.043398641533973964) < 1e-9 and abs(cell.ci_high - 0.18824547127832675) < 1e-9
    assert cell is res.marker_error.by_compartment["Tumor"]
    x = res.discordance
    assert abs(x.mean_excess - 0.058874659699858695) < 1e-12 and x.n_pairs == 1124 and x.n_patients == 11
    assert abs(x.ci_low - 0.01984978980133284) < 1e-9 and abs(x.ci_high - 0.09789952959838455) < 1e-9
    assert abs(x.equal_weights.mean - 0.04260169851300013) < 1e-12


def test_barcode_suffix_normalized_and_declared_on_real_data(real):
    _, res = real
    m = res.barcode_match
    print("\n[reale] " + m.sentence)
    assert m.normalized_suffixes == ("-1",)
    assert "Normalizzato il suffisso -1" in m.sentence
    assert m.n_matched == 105_809 and m.fraction > 0.5
