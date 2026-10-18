"""Lettura dei formati reali usati dalla web app e dalla CLI.

- metadati per campione o per cellula: CSV o TSV (separatore riconosciuto dall'estensione o
  dal contenuto), letti come stringhe;
- matrice di conteggi: AnnData (.h5ad) oppure cartella 10x Matrix Market
  (matrix.mtx[.gz], barcodes.tsv[.gz], features.tsv[.gz] o genes.tsv[.gz]) con una tabella di
  metadati per cellula che contiene la colonna dei barcode.
- manifest dei file VDJ di Cell Ranger (``read_vdj_manifest``).
Ogni incoerenza (file mancante, barcode non trovati o ripetuti, dimensioni diverse, colonne con lo
stesso nome) solleva un errore esplicito: nessuna riga viene scartata in silenzio.
"""

from __future__ import annotations

import io
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scipy.io
import scipy.sparse as sp


def read_table(source) -> pd.DataFrame:
    """CSV o TSV come stringhe. ``source``: percorso o bytes/file caricato. Errore chiaro se il
    file non esiste o non è una tabella leggibile, e se ha colonne con lo stesso nome."""
    try:
        if isinstance(source, (str, Path)):
            path = Path(source)
            if not path.exists():
                raise FileNotFoundError(f"file non trovato: {path}")
            low = path.name.lower()
            sep = "\t" if path.suffix.lower() in (".tsv", ".txt") or low.endswith((".tsv.gz", ".txt.gz")) else None
            df = pd.read_csv(path, sep=sep, dtype=str, engine="python" if sep is None else "c")
        else:
            data = source if isinstance(source, bytes) else source.getvalue()
            head = data[:4096].decode("utf-8", errors="ignore")
            sep = "\t" if head.count("\t") > head.count(",") else ","
            df = pd.read_csv(io.BytesIO(data), sep=sep, dtype=str)
    except (pd.errors.ParserError, pd.errors.EmptyDataError, UnicodeDecodeError) as e:
        raise ValueError(f"tabella illeggibile ({e.__class__.__name__}): {e}") from e
    dup = [c for c in df.columns if "." in str(c) and str(c).rsplit(".", 1)[0] in df.columns
           and str(c).rsplit(".", 1)[1].isdigit()]
    if dup:
        raise ValueError(f"la tabella ha colonne con lo stesso nome (rinominate da pandas in {dup}): "
                         "ogni colonna deve avere un nome unico")
    return df


def _first(folder: Path, names: list[str]) -> Path:
    """Primo file esistente fra ``names`` (anche con suffisso .gz) nella cartella."""
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
    dup = cell_metadata[barcode_col].duplicated(keep=False)
    if dup.any():
        raise ValueError(f"{int(dup.sum())} righe dei metadati per cellula hanno un barcode ripetuto "
                         f"nella colonna '{barcode_col}' (es. {cell_metadata.loc[dup, barcode_col].head(3).tolist()}): "
                         "ogni barcode deve comparire una volta sola")
    if bc.duplicated().any():
        raise ValueError("barcodes.tsv contiene barcode ripetuti: la matrice non identifica le cellule")
    meta = cell_metadata.set_index(barcode_col)
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
    """Nomi resi unici aggiungendo -1, -2, ... alle ripetizioni (simboli genici duplicati)."""
    seen: dict[str, int] = {}
    out = []
    for n in names:
        k = seen.get(n, 0)
        out.append(n if k == 0 else f"{n}-{k}")
        seen[n] = k + 1
    return out


def read_vdj_manifest(path: Path) -> list[tuple]:
    """Manifest VDJ -> lista per ``parse_vdj_contigs``. Colonne path, patient, compartment e,
    facoltativa, library. I percorsi relativi sono rispetto alla cartella del manifest."""
    m = read_table(path)
    for col in ("path", "patient", "compartment"):
        if col not in m.columns:
            raise ValueError(f"il manifest VDJ deve avere le colonne path,patient,compartment (manca '{col}')")
    if m[["path", "patient", "compartment"]].isna().any().any():
        raise ValueError("il manifest VDJ ha celle vuote nelle colonne path, patient o compartment")
    base = Path(path).resolve().parent
    files = []
    for r in m.itertuples(index=False):
        p = Path(r.path)
        p = p if p.is_absolute() else base / p
        if not p.exists():
            raise ValueError(f"file VDJ del manifest non trovato: {p}")
        lib = getattr(r, "library", None) if "library" in m.columns else None
        lib = None if lib is None or pd.isna(lib) or str(lib).strip() == "" else str(lib)
        files.append((p, str(r.patient), str(r.compartment)) + ((lib,) if lib is not None else ()))
    return files
