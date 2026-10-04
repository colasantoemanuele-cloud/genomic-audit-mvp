"""Pulizia 4 -- diagnostica su P(chiamata CD8 | vera CD4) = 0.45 nel tumore (GSE278694 reale).

Verifica se quel valore dipende da come il sangue assegna l'identità del clone. Ricalcola
la matrice di confusione del tumore con soglie di riferimento più o meno stringenti
(margine 0.20 / 0.40, minimo di cellule nel sangue 3 / 5), SENZA cambiare i default
dell'MVP: usa le stesse funzioni (assign_reference_identity, cd8_fraction_intervals) con
parametri diversi. Per ogni combinazione riporta anche la frazione dei cloni "CD4" di
riferimento con almeno una cellula CD8A+ nel sangue.

Uso:
    python -m validation.gse278694_cd8_diagnostic --vdj-manifest <manifest.csv> --out <dir>
"""

from __future__ import annotations

import os
import argparse
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd

from core.cd8_propagation import cd8_fraction_intervals
from core.tcr_validation import (
    assign_reference_identity,
    build_clonotypes,
    compute_marker_positivity,
    match_barcodes,
    parse_vdj_contigs,
)

# Dati reali di GSE278694: non inclusi nella repo. Cartella indicata da GSE278694_DIR, altrimenti
# una cartella "pdac-ml" accanto a questa repo.
H5AD = Path(os.environ.get("GSE278694_DIR", Path(__file__).resolve().parents[2] / "pdac-ml")) / "data/interim/sc_raw.h5ad"
MARKERS = {"CD4T": ["CD4"], "CD8T": ["CD8A", "CD8B"]}


def load_cells(manifest: Path) -> pd.DataFrame:
    """Stessa preparazione di run_tcr_validation, più la positività al solo CD8A."""
    adata = ad.read_h5ad(H5AD)
    obs = adata.obs[["patients", "tissue", "all_celltype", "barcode"]].astype(str).reset_index(drop=True)
    obs.columns = ["patient", "compartment", "celltype", "barcode"]
    pos = compute_marker_positivity(adata, {**MARKERS, "CD8A_only": ["CD8A"]}).reset_index(drop=True)
    obs = pd.concat([obs, pos], axis=1)
    del adata
    m = pd.read_csv(manifest)
    clones = build_clonotypes(parse_vdj_contigs([(Path(r.path), str(r.patient), str(r.compartment))
                                                 for r in m.itertuples()]))
    obs, clones, bm = match_barcodes(obs, clones)
    print(bm.sentence)
    return obs.merge(clones, on=["patient", "compartment", "barcode"], how="inner")


def diagnose(cells: pd.DataFrame, min_cells: int, min_margin: float) -> dict:
    ref = assign_reference_identity(cells, list(MARKERS), reference_compartment="PBMC",
                                    min_cells=min_cells, min_margin=min_margin)
    j = cells.merge(ref[["patient", "clone_id", "reference_identity"]], on=["patient", "clone_id"], how="left")
    res = cd8_fraction_intervals(j, "patient", "compartment", "celltype", "reference_identity", "Tumor",
                                 n_boot=500, seed=0)
    # cloni CD4 di riferimento con almeno una cellula CD8A+ nel sangue
    blood = cells[cells.compartment == "PBMC"].merge(ref, on=["patient", "clone_id"])
    cd4 = blood[blood.reference_identity == "CD4T"].groupby(["patient", "clone_id"]).pos_CD8A_only.any()
    # cloni di riferimento effettivamente presenti fra le cellule usate nel tumore
    tum = j[(j.compartment == "Tumor") & j.reference_identity.notna()]
    m = res.matrix
    return {
        "margine": min_margin, "min cellule sangue": min_cells,
        "cloni di riferimento (tutti)": len(ref),
        "cloni CD4 / CD8 di riferimento": f"{int((ref.reference_identity == 'CD4T').sum())} / "
                                          f"{int((ref.reference_identity == 'CD8T').sum())}",
        "cloni usati nel tumore": int(tum.groupby(["patient", "clone_id"]).ngroups),
        "cellule tumore con riferimento": res.n_reference_cells,
        "pazienti": res.n_reference_patients,
        "P(chiamata CD8 | vera CD4)": None if m is None else round(float(m.iloc[0, 1]), 3),
        "P(chiamata CD4 | vera CD8)": None if m is None else round(float(m.iloc[1, 0]), 3),
        "J": None if res.youden_j is None else round(res.youden_j, 3),
        "cloni CD4 rif. con >=1 cellula CD8A+ nel sangue": f"{cd4.mean():.1%} ({int(cd4.sum())}/{len(cd4)})",
    }


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--vdj-manifest", type=Path, required=True)
    ap.add_argument("--out", type=Path, required=True)
    args = ap.parse_args()
    cells = load_cells(args.vdj_manifest)
    rows = [diagnose(cells, mc, mm) for mm in (0.20, 0.40) for mc in (3, 5)]
    tab = pd.DataFrame(rows)
    tab.to_csv(args.out / "pulizia4_cd8_diagnostic.csv", index=False)
    with pd.option_context("display.width", 250, "display.max_columns", 20):
        print(tab.to_string(index=False))


if __name__ == "__main__":
    np.seterr(all="ignore")
    main()
