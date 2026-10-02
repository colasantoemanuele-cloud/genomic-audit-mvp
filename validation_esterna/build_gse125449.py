"""Costruisce validation_esterna/data/GSE125449.h5ad dai file GEO (Set1 + Set2).

- Conteggi grezzi interi dai file .mtx (geni x cellule, trasposti in cellule x geni).
- Geni: intersezione degli ID Ensembl dei due set (dichiarata nell'output), nello stesso ordine.
- obs: `Sample` e `Type` da *_samples.txt (etichetta degli autori, non ricalcolata), `patient`
  = token "LCP<id>" del nome del campione (stessa regex di parse_geo.py), `set` = Set1/Set2.
- Nessun filtro sulle cellule e nessuna normalizzazione.

Uso: python -m validation_esterna.build_gse125449
"""

from __future__ import annotations

from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scipy.io
import scipy.sparse as sp

DATA = Path(__file__).parent / "data"


def load_set(name: str) -> ad.AnnData:
    X = scipy.io.mmread(DATA / f"GSE125449_{name}_matrix.mtx.gz").T.tocsr()
    genes = pd.read_csv(DATA / f"GSE125449_{name}_genes.tsv.gz", sep="\t", header=None, names=["id", "symbol"])
    bc = pd.read_csv(DATA / f"GSE125449_{name}_barcodes.tsv.gz", sep="\t", header=None)[0].astype(str)
    samples = pd.read_csv(DATA / f"GSE125449_{name}_samples.txt.gz", sep="\t")
    assert bc.is_unique, f"{name}: barcode duplicati"
    assert samples["Cell Barcode"].is_unique and set(samples["Cell Barcode"]) == set(bc), \
        f"{name}: i barcode della matrice e di samples.txt non coincidono"
    obs = samples.set_index("Cell Barcode").loc[bc.values].reset_index()
    obs["patient"] = obs["Sample"].str.extract(r"_(LCP\d+)$")[0]
    assert obs["patient"].notna().all()
    obs["set"] = name
    obs.index = [f"{name}_{b}" for b in obs["Cell Barcode"]]
    assert np.allclose(X.data, np.round(X.data)), f"{name}: conteggi non interi"
    return ad.AnnData(X=X.astype(np.float32), obs=obs, var=genes.set_index("id"))


def main() -> None:
    a1, a2 = load_set("Set1"), load_set("Set2")
    common = a1.var_names.intersection(a2.var_names)
    print(f"geni Set1 {a1.n_vars}, Set2 {a2.n_vars}, in comune {len(common)}")
    adata = ad.concat([a1[:, common], a2[:, common]], join="inner")
    adata.var = a1.var.loc[common]
    assert sp.issparse(adata.X)
    adata.write_h5ad(DATA / "GSE125449.h5ad")
    o = adata.obs
    print(f"{adata.n_obs} cellule, {adata.n_vars} geni, {o.patient.nunique()} pazienti, "
          f"{o.Sample.nunique()} campioni")
    print(o.Type.value_counts().to_string())
    print(pd.crosstab(o.patient, o.Type).to_string())


if __name__ == "__main__":
    main()
