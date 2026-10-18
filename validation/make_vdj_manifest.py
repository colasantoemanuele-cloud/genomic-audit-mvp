"""Manifest VDJ di GSE278694 per `audit-sc tcr`, dai file originali di Cell Ranger.

I file GEO si chiamano GSM<n>_PA_<paziente>-<Tumor|Normal|PBMC>-VDJ_filtered_contig_annotations.csv.gz;
il manifest usa gli stessi identificativi dei metadati di sc_raw.h5ad (PA01..PA15, Tumor,
Adjacent_normal, PBMC). I percorsi sono scritti relativi alla cartella del manifest, come li legge
la CLI. I barcode VDJ hanno il suffisso "-1" e quelli dell'AnnData no: lo strumento normalizza il
suffisso da solo e lo dichiara.

Uso: python -m validation.make_vdj_manifest --out validation/results/vdj_manifest_gse278694.csv
"""

from __future__ import annotations

import argparse
import os
import re
from pathlib import Path

import pandas as pd

PDAC = Path(os.environ.get("GSE278694_DIR", Path(__file__).resolve().parents[2] / "pdac-ml"))
TISSUE = {"Tumor": "Tumor", "Normal": "Adjacent_normal", "PBMC": "PBMC"}


def build_manifest(out: Path, pdac: Path = PDAC) -> int:
    """Scrive il manifest in ``out`` (percorsi relativi alla sua cartella); ritorna il numero di file."""
    rows = []
    for p in sorted((pdac / "data/raw/tcr").glob("*_filtered_contig_annotations.csv.gz")):
        m = re.match(r"GSM\d+_PA_(\d+)-(\w+)-VDJ_filtered_contig_annotations\.csv\.gz", p.name)
        rel = os.path.relpath(p.resolve(), out.resolve().parent)
        rows.append({"path": rel, "patient": f"PA{int(m.group(1)):02d}", "compartment": TISSUE[m.group(2)]})
    out.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows).to_csv(out, index=False)
    return len(rows)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", type=Path, required=True)
    args = ap.parse_args()
    print(f"{build_manifest(args.out)} file VDJ -> {args.out}")


if __name__ == "__main__":
    main()
