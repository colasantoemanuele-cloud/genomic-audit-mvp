"""Lettura dei formati reali usati dalla web app e dalla CLI.

- metadati per campione o per cellula: CSV o TSV (separatore riconosciuto dall'estensione o
  dal contenuto), letti come stringhe;
- matrice di conteggi: AnnData (.h5ad) oppure cartella 10x Matrix Market
  (matrix.mtx[.gz], barcodes.tsv[.gz], features.tsv[.gz] o genes.tsv[.gz]) con una tabella di
  metadati per cellula che contiene la colonna dei barcode.
Ogni incoerenza (file mancante, barcode non trovati, dimensioni diverse) solleva un errore
esplicito: nessuna riga viene scartata in silenzio.
"""

from __future__ import annotations

import io
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scipy.io
import scipy.sparse as sp


def read_table(source, name: str | None = None) -> pd.DataFrame:
    """CSV o TSV come stringhe. ``source``: percorso o bytes/file caricato."""
    if isinstance(source, (str, Path)):
        path = Path(source)
        sep = "\t" if path.suffix.lower() in (".tsv", ".txt") or path.name.lower().endswith((".tsv.gz", ".txt.gz")) else None
        return pd.read_csv(path, sep=sep, dtype=str, engine="python" if sep is None else "c")
    data = source if isinstance(source, bytes) else source.getvalue()
    head = data[:4096].decode("utf-8", errors="ignore")
    sep = "\t" if head.count("\t") > head.count(",") else ","
    return pd.read_csv(io.BytesIO(data), sep=sep, dtype=str)


def _first(folder: Path, names: list[str]) -> Path:
    for n in names:
        for cand in (folder / n, folder / f"{n}.gz"):
            if cand.exists():
                return cand
    raise FileNotFoundError(f"nella cartella {folder} manca uno di: {', '.join(names)} (anche .gz)")


def load_matrix_market(folder: str | Path, cell_metadata: pd.DataFrame, barcode_col: str) -> ad.AnnData:
    """Cartella 10x -> AnnData (cellule x geni, conteggi sparsi). Le righe di ``cell_metadata``
    vengono allineate ai barcode della matrice; i barcode della matrice senza metadati sono un
    errore, non vengono scartati."""
    folder = Path(folder)
    X = scipy.io.mmread(_first(folder, ["matrix.mtx"]))
    X = sp.csr_matrix(X).T.tocsr()  # 10x: geni x cellule
    bc = pd.read_csv(_first(folder, ["barcodes.tsv"]), sep="\t", header=None, dtype=str)[0]
    feat = pd.read_csv(_first(folder, ["features.tsv", "genes.tsv"]), sep="\t", header=None, dtype=str)
    if X.shape != (len(bc), len(feat)):
        raise ValueError(f"dimensioni incoerenti: matrice {X.shape}, barcode {len(bc)}, geni {len(feat)}")
    if barcode_col not in cell_metadata.columns:
        raise ValueError(f"colonna dei barcode '{barcode_col}' assente dai metadati per cellula")
    meta = cell_metadata.drop_duplicates(barcode_col).set_index(barcode_col)
    missing = ~bc.isin(meta.index)
    if missing.any():
        raise ValueError(f"{int(missing.sum())} barcode della matrice su {len(bc)} non hanno metadati "
                         f"(es. {bc[missing].head(3).tolist()}): controlla la colonna '{barcode_col}'")
    obs = meta.loc[bc.values].reset_index()
    obs.index = bc.values
    symbols = feat[1] if feat.shape[1] > 1 else feat[0]
    var = pd.DataFrame({"gene_id": feat[0].values}, index=pd.Index(symbols.values).astype(str))
    var.index = pd.Index(_dedup(var.index.tolist()))
    return ad.AnnData(X=X.astype(np.float32), obs=obs, var=var)


def _dedup(names: list[str]) -> list[str]:
    seen: dict[str, int] = {}
    out = []
    for n in names:
        k = seen.get(n, 0)
        out.append(n if k == 0 else f"{n}-{k}")
        seen[n] = k + 1
    return out
