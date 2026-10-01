"""Dati della demo della web app (data/demo/), ricostruibili dai file GEO scaricati.

- GSE125449_demo.h5ad: lo stesso sottoinsieme di build_small.py (10 pazienti su 19, al massimo
  200 cellule per paziente, seme 0), scritto con compressione gzip;
- *_samples.csv: tabelle dei campioni prodotte da parse_geo.py.

Uso: python -m validation_esterna.parse_geo && python -m validation_esterna.build_small &&
     python -m validation_esterna.build_demo
"""

from __future__ import annotations

import shutil
from pathlib import Path

import anndata as ad

HERE = Path(__file__).parent
DEMO = HERE.parent / "data" / "demo"


def main() -> None:
    DEMO.mkdir(parents=True, exist_ok=True)
    a = ad.read_h5ad(HERE / "data" / "GSE125449_small.h5ad")
    a.write_h5ad(DEMO / "GSE125449_demo.h5ad", compression="gzip")
    for g in ("GSE125449", "GSE132465", "GSE131907"):
        shutil.copy(HERE / "results" / f"{g}_samples.csv", DEMO / f"{g}_samples.csv")
    for p in sorted(DEMO.iterdir()):
        print(f"{p.stat().st_size:>12,}  {p.name}")


if __name__ == "__main__":
    main()
