"""Coorte ridotta del Modulo A (CRITERI.md, sezione 5): 10 pazienti su 19 di GSE125449
(default_rng(0), senza reinserimento) e al massimo 200 cellule per paziente
(default_rng(0)). Nessun filtro su geni o classi.

Uso: python -m validation_esterna.build_small
"""

from __future__ import annotations

from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd

DATA = Path(__file__).parent / "data"


def main() -> None:
    a = ad.read_h5ad(DATA / "GSE125449.h5ad")
    rng = np.random.default_rng(0)
    patients = sorted(a.obs["patient"].astype(str).unique())
    chosen = sorted(rng.choice(patients, size=10, replace=False).tolist())
    idx = []
    for p in chosen:
        rows = np.flatnonzero(a.obs["patient"].astype(str).values == p)
        idx += sorted(rng.choice(rows, size=min(200, len(rows)), replace=False).tolist())
    small = a[np.array(idx)].copy()
    small.write_h5ad(DATA / "GSE125449_small.h5ad")
    print(f"pazienti scelti: {chosen}")
    print(f"{small.n_obs} cellule, {small.n_vars} geni, {small.obs.patient.nunique()} pazienti")
    print(pd.crosstab(small.obs.patient, small.obs.Type).to_string())


if __name__ == "__main__":
    main()
