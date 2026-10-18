import anndata as ad, numpy as np, pandas as pd, scipy.sparse as sp

parts, genes = [], []
for chunk in pd.read_csv("../GSE144735_processed_KUL3_CRC_10X_raw_UMI_count_matrix.txt.gz", sep="\t", index_col=0, chunksize=2000):  # 2000 geni per volta
    parts.append(sp.csr_matrix(chunk.to_numpy(dtype=np.float32)))
    genes.extend(chunk.index.astype(str))
cells = chunk.columns.astype(str)
X = sp.vstack(parts).T.tocsr()                                  # cellule x geni
obs = pd.read_csv("../GSE144735_processed_KUL3_CRC_10X_annotation.txt.gz", sep="\t", index_col=0, dtype=str)
obs = obs.loc[cells]                                            # stesse cellule, stesso ordine
a = ad.AnnData(X=X, obs=obs, var=pd.DataFrame(index=genes))
a.var_names_make_unique()
a.write_h5ad("studio.h5ad")
