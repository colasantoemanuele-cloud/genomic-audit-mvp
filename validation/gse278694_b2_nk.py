"""Revisione B2 -- di quanto scende la differenza fra la convenzione per cellula e quella per clone
(tesi) togliendo le etichette NK dalle etichette contate come errore, su GSE278694. È una
descrizione, non un'attribuzione di causa: la scomposizione non è additiva e la quota dipende
dall'ordine (aggiungendo le NK alla convenzione per cellula, passo S2 -> S3 di
validation/gse278694_moduleB.py, la quota è diversa).

Quota = (per clone con NK come errore - per clone senza etichette fuori mappa) /
        (per clone con NK come errore - per cellula), nel tumore, con la regola di positività
esclusiva della tesi (CD8T > CD4T). Riporta anche la convenzione per clone di default (tutte le
etichette fuori mappa contate come errore, anche mieloidi e B).

Uso: python -m validation.gse278694_b2_nk --vdj-manifest validation/results/vdj_manifest_gse278694.csv \
         --out validation/results
"""

from __future__ import annotations

import argparse
import json
import os
from pathlib import Path

import anndata as ad

from core.io import read_vdj_manifest
from core.tcr_validation import parse_vdj_contigs, run_tcr_validation

PDAC = Path(os.environ.get("GSE278694_DIR", Path(__file__).resolve().parents[2] / "pdac-ml"))
MARKERS = {"CD4T": ["CD4"], "CD8T": ["CD8A", "CD8B"]}


def main() -> None:
    """Calcola la quota NK della differenza fra convenzioni e scrive b2_quota_nk.json."""
    ap = argparse.ArgumentParser()
    ap.add_argument("--vdj-manifest", type=Path, required=True)
    ap.add_argument("--out", type=Path, required=True)
    args = ap.parse_args()
    adata = ad.read_h5ad(PDAC / "data/interim/sc_raw.h5ad")
    contigs = parse_vdj_contigs(read_vdj_manifest(args.vdj_manifest))
    out = {}
    for name, labels in (("per_clone_NK_errore", ["NK"]), ("per_clone_nessuna_fuori_mappa", []),
                         ("per_clone_tutte_fuori_mappa", None)):
        r = run_tcr_validation(adata, contigs, "patients", "tissue", "all_celltype", "barcode",
                               marker_map=MARKERS, reference_compartment="PBMC", clone_error_labels=labels,
                               clone_marker_priority=["CD8T", "CD4T"])
        out[name] = r.conventions["clone"].by_compartment["Tumor"].mean
        out["per_cellula"] = r.conventions["cell"].by_compartment["Tumor"].mean
    out["quota_NK_della_differenza"] = ((out["per_clone_NK_errore"] - out["per_clone_nessuna_fuori_mappa"])
                                        / (out["per_clone_NK_errore"] - out["per_cellula"]))
    (args.out / "b2_quota_nk.json").write_text(json.dumps(out, indent=1))
    print(json.dumps(out, indent=1))


if __name__ == "__main__":
    main()
