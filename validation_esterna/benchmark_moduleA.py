"""Benchmark dei tempi del Modulo A (obiettivo: esecuzione completa, confronto fra modelli
incluso, entro 10 minuti per ~25.000 cellule e 20 pazienti su CPU locale).

- ``sintetico``: 25.000 cellule, 20 pazienti, 20.000 geni, matrice sparsa (~5% non zeri)
  con impronta genica per paziente e segnale di classe (4 classi);
- ``reale``: GSE125449 completo (9.946 cellule, 18.372 geni, 19 pazienti).

Uso: python -m validation_esterna.benchmark_moduleA {sintetico|reale}
"""

from __future__ import annotations

import json
import sys
import time
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse as sp

from core.leakage_audit import run_leakage_audit

RES = Path(__file__).parent / "results"


def synthetic(n_cells=25_000, n_patients=20, n_genes=20_000, n_classes=4, density=0.05, seed=0) -> ad.AnnData:
    rng = np.random.default_rng(seed)
    patient = rng.integers(0, n_patients, n_cells)
    label = rng.integers(0, n_classes, n_cells)
    X = sp.random(n_cells, n_genes, density=density, format="csr", random_state=seed,
                  data_rvs=lambda k: rng.poisson(2.0, k) + 1).astype(np.float32)
    X = X.tolil()
    for c in range(n_classes):  # segnale di classe: 30 geni per classe
        genes = np.arange(c * 30, c * 30 + 30)
        rows = np.flatnonzero(label == c)
        for g in genes:
            X[rows, g] = rng.poisson(6.0, len(rows)) + 1
    for p in range(n_patients):  # impronta per paziente: 50 geni
        genes = 1000 + p * 50 + np.arange(50)
        rows = np.flatnonzero(patient == p)
        for g in genes:
            X[rows, g] = rng.poisson(8.0, len(rows)) + 1
    obs = pd.DataFrame({"patient": [f"P{p:02d}" for p in patient], "label": [f"C{c}" for c in label]})
    return ad.AnnData(X=X.tocsr(), obs=obs, var=pd.DataFrame(index=[f"g{i}" for i in range(n_genes)]))


def main() -> None:
    which = sys.argv[1]
    if which == "sintetico":
        t = time.time()
        adata = synthetic()
        print(f"generazione: {time.time() - t:.0f} s")
        target, patient = "label", "patient"
    else:
        adata = ad.read_h5ad(Path(__file__).parent / "data" / "GSE125449.h5ad")
        target, patient = "Type", "patient"
    t = time.time()
    r = run_leakage_audit(adata, target_col=target, patient_col=patient, seed=0,
                          progress=lambda m: print(m, flush=True))
    el = time.time() - t
    out = {"dataset": which, "cellule": r.n_cells, "geni": int(adata.n_vars), "pazienti": r.n_patients,
           "secondi": round(el, 1), "entro_10_minuti": el <= 600,
           "raggruppato": r.grouped.mean, "casuale": r.random.mean, "divario": r.gap,
           "logo": {m: s.mean for m, s in r.model_comparison.scores.items()} if r.model_comparison else None,
           "fold_con_classi_assenti": r.grouped.n_folds_with_absent_classes,
           "xai_jaccard": [r.xai.grouped_mean, r.xai.random_mean]}
    (RES / f"benchmark_moduleA_{which}.json").write_text(json.dumps(out, indent=1))
    print(json.dumps(out, indent=1))


if __name__ == "__main__":
    main()
