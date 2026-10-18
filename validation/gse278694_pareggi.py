"""Revisione A8 -- quante scelte di GSE278694 dipendevano dall'ordine delle righe.

1. Contig TRB a pari UMI massime nella stessa cellula: prima della revisione parse_vdj_contigs
   sceglieva secondo l'ordine dell'ordinamento (non stabile); ora vince il contig con più reads,
   poi la CDR3 minore in ordine alfabetico.
2. Pareggi dell'etichetta di maggioranza nelle unità della convenzione per clone (opzioni della
   tesi): prima della revisione decideva il primo valore di value_counts (ordine delle righe); ora
   l'unità è esclusa, come nella versione finale della tesi, e come sensibilità contata come errore
   o come corretta (``--tie-rule``).

Uso: python -m validation.gse278694_pareggi --vdj-manifest validation/results/vdj_manifest_gse278694.csv \
         --out validation/results
"""

from __future__ import annotations

import argparse
import json
import os
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd

from core.io import read_vdj_manifest
from core.tcr_validation import (
    _KEY,
    assign_reference_identity,
    build_clonotypes,
    compute_marker_positivity,
    match_barcodes,
    parse_vdj_contigs,
)

PDAC = Path(os.environ.get("GSE278694_DIR", Path(__file__).resolve().parents[2] / "pdac-ml"))
MARKERS = {"CD4T": ["CD4"], "CD8T": ["CD8A", "CD8B"]}


def contig_ties(manifest: Path) -> dict:
    """Cellule con due o più CDR3 TRB a pari UMI massime, e quante restano pari anche nelle reads."""
    frames = []
    for path, patient, comp, *_ in read_vdj_manifest(manifest):
        df = pd.read_csv(path, low_memory=False)
        for col in ("is_cell", "high_confidence", "productive", "full_length"):
            if col in df.columns:
                df = df[df[col].astype(str).str.lower().isin(["true", "1"])]
        df = df[(df.chain == "TRB") & df.cdr3_nt.notna() & (df.cdr3_nt.astype(str) != "None")]
        frames.append(df.assign(patient=patient, compartment=comp))
    c = pd.concat(frames)
    key = ["patient", "compartment", "barcode"]
    top = c[c.umis == c.groupby(key).umis.transform("max")]
    tied = top.groupby(key).cdr3_nt.nunique()
    tied = tied[tied > 1].index
    t = top.set_index(key).loc[tied]
    still = int(t.groupby(level=[0, 1, 2]).apply(lambda x: x[x.reads == x.reads.max()].cdr3_nt.nunique() > 1).sum())
    return {"cellule_con_TRB": int(c.groupby(key).ngroups), "pari_UMI": int(len(tied)),
            "pari_anche_nelle_reads": still}


def majority_ties(manifest: Path) -> dict:
    """Unità per clone (opzioni della tesi: NK come errore, positività CD8T > CD4T, >= 3 cellule) in
    cui l'etichetta di maggioranza è a pari merito, per compartimento giudicato."""
    adata = ad.read_h5ad(PDAC / "data/interim/sc_raw.h5ad")
    obs = pd.DataFrame({"patient": adata.obs["patients"].astype(str).values,
                        "compartment": adata.obs["tissue"].astype(str).values,
                        "celltype": adata.obs["all_celltype"].astype(str).values,
                        "barcode": adata.obs["barcode"].astype(str).values, "library": ""})
    obs = pd.concat([obs, compute_marker_positivity(adata, MARKERS).reset_index(drop=True)], axis=1)
    del adata
    obs, clones, _ = match_barcodes(obs, build_clonotypes(parse_vdj_contigs(read_vdj_manifest(manifest))))
    x = obs.merge(clones, on=_KEY)
    x = x[x.celltype.isin(["CD4T", "CD8T", "NK"])].copy()
    taken = np.zeros(len(x), dtype=bool)
    for lab in ("CD8T", "CD4T"):
        x[f"pos_{lab}"] = x[f"pos_{lab}"].values & ~taken
        taken |= x[f"pos_{lab}"].values
    ref = assign_reference_identity(x, list(MARKERS), reference_compartment="PBMC")
    t = x[x.compartment != "PBMC"].merge(ref[["patient", "clone_id", "reference_identity"]], on=["patient", "clone_id"])
    out = {}
    for (comp,), g in t.groupby(["compartment"]):
        n_units = n_ties = 0
        for _, u in g.groupby(["patient", "clone_id"]):
            if len(u) < 3:
                continue
            n_units += 1
            vc = u.celltype.value_counts()
            n_ties += int(len(vc) > 1 and vc.iloc[0] == vc.iloc[1])
        out[comp] = {"unita": n_units, "pareggi": n_ties}
    out["totale"] = {"unita": sum(v["unita"] for v in out.values()), "pareggi": sum(v["pareggi"] for v in out.values())}
    return out


def main() -> None:
    """Conta le scelte che dipendevano dall'ordine e scrive a8_pareggi.json."""
    ap = argparse.ArgumentParser()
    ap.add_argument("--vdj-manifest", type=Path, required=True)
    ap.add_argument("--out", type=Path, required=True)
    args = ap.parse_args()
    out = {"contig_TRB": contig_ties(args.vdj_manifest), "maggioranza_per_clone": majority_ties(args.vdj_manifest)}
    (args.out / "a8_pareggi.json").write_text(json.dumps(out, indent=1))
    print(json.dumps(out, indent=1))


if __name__ == "__main__":
    main()
