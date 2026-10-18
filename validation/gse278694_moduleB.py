"""Verifica di fedeltà del porting del Modulo B su GSE278694 (non una validazione del metodo).

Legge (sola lettura) i dati di pdac-ml e confronta le definizioni del Modulo B con quelle della
tesi (07_robustness.py, 10_loco.py, versione da cui il codice è stato portato). Non tara nessun
parametro: parte dalle definizioni dello strumento e le sostituisce UNA ALLA VOLTA con quelle della
tesi, fino a riprodurre il numero della tesi. Per riprodurre ESATTAMENTE quei numeri usa le regole
della tesi di allora: il bootstrap percentile (core.stats.cluster_bootstrap, seme 0, 2000 repliche)
e, per i pareggi dell'etichetta di maggioranza, il primo valore di value_counts (dipende dall'ordine
delle righe). Lo strumento oggi usa invece l'intervallo t cluster-robusto e una regola dei pareggi
che non dipende dall'ordine: questi numeri servono a mostrare che il porting segue le regole della
tesi, non che le regole siano giuste. Gli intervalli con meno di 5 pazienti qui riportati
(min_groups=3) sono quelli della tesi: lo strumento li rifiuta, e con 4 pazienti esistono solo 35
ricampionamenti distinti.

Uso:
    python -m validation.gse278694_moduleB --vdj-manifest <manifest.csv> --out <dir>

``manifest.csv``: path,patient,compartment dei CSV VDJ originali (validation.make_vdj_manifest);
il suffisso "-1" dei barcode è normalizzato come nello strumento (match_barcodes).
"""

from __future__ import annotations

import argparse
import json
import os
from math import comb
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd

from core.io import read_vdj_manifest
from core.stats import cluster_bootstrap
from core.tcr_validation import build_clonotypes, match_barcodes, pairwise_excess_discordance, parse_vdj_contigs

# Dati reali di GSE278694: non inclusi nella repo. Cartella indicata da GSE278694_DIR, altrimenti
# una cartella "pdac-ml" accanto a questa repo.
H5AD = Path(os.environ.get("GSE278694_DIR", Path(__file__).resolve().parents[2] / "pdac-ml")) / "data/interim/sc_raw.h5ad"
N_BOOT, SEED = 2000, 0
LYMPHOID = ("CD4T", "CD8T", "NK")


def load_cells(manifest: Path) -> pd.DataFrame:
    """Una riga per cellula con TCR (TRB): stessa ingestione e stessa merge dell'MVP, più le
    conte grezze di CD4, CD8A, CD8B."""
    adata = ad.read_h5ad(H5AD)
    obs = adata.obs[["patients", "tissue", "all_celltype", "barcode"]].astype(str).reset_index(drop=True)
    obs.columns = ["patient", "compartment", "celltype", "barcode"]
    for g in ("CD4", "CD8A", "CD8B"):
        i = adata.var_names.get_loc(g)
        obs[g] = np.asarray(adata.X[:, i].todense()).ravel()
    del adata
    contigs = parse_vdj_contigs(read_vdj_manifest(manifest))
    obs, clones, _ = match_barcodes(obs, build_clonotypes(contigs))
    return obs.merge(clones, on=["patient", "compartment", "library", "barcode"], how="inner")


def discordance(cells: pd.DataFrame, lymphoid_only: bool) -> dict:
    """Eccesso di discordanza con il bootstrap percentile della tesi."""
    x = cells[cells.celltype.isin(LYMPHOID)] if lymphoid_only else cells
    pw = pairwise_excess_discordance(x)
    b = cluster_bootstrap(pw.excess.values, pw.patient.values, n_boot=N_BOOT, seed=SEED, min_groups=5)
    return {"cellule": len(x), "coppie": len(pw), "pazienti": int(pw.patient.nunique()),
            "d_within": float(pw.d_within.mean()), "eccesso": b.mean, "ic": [b.ci_low, b.ci_high]}


def error_units(cells: pd.DataFrame, *, lymphoid_only: bool, exclusive_cd4: bool, nk_as_error: bool,
                min_target: int, unit: str) -> pd.DataFrame:
    """Unità di errore (cellule o coppie clone-compartimento) con l'identità dal solo sangue.
    Con tutte le opzioni a False/0/'cell' coincide con marker_error_rate dell'MVP."""
    x = cells[cells.celltype.isin(LYMPHOID)].copy() if lymphoid_only else cells.copy()
    x["cd8_pos"] = (x.CD8A > 0) | (x.CD8B > 0)
    x["cd4_pos"] = (x.CD4 > 0) & ~x.cd8_pos if exclusive_cd4 else (x.CD4 > 0)
    ref = x[x.compartment == "PBMC"].groupby(["patient", "clone_id"], observed=True).agg(
        n=("cd8_pos", "size"), f8=("cd8_pos", "mean"), f4=("cd4_pos", "mean")).reset_index()
    ref = ref[(ref.n >= 3) & ((ref.f8 - ref.f4).abs() >= 0.20)]
    ref["identity"] = np.where(ref.f8 > ref.f4, "CD8T", "CD4T")
    t = x[x.compartment != "PBMC"].merge(ref[["patient", "clone_id", "identity"]], on=["patient", "clone_id"])
    size = t.groupby(["patient", "clone_id", "compartment"], observed=True).celltype.transform("size")
    t = t[size >= min_target]
    allowed = set(LYMPHOID) if nk_as_error else {"CD4T", "CD8T"}
    if unit == "cell":
        t = t[t.celltype.isin(allowed)]
        return pd.DataFrame({"patient": t.patient, "clone": t.clone_id, "tissue": t.compartment,
                             "error": (t.celltype != t.identity).astype(float)})
    rows = []
    for (p, c, tis), g in t.groupby(["patient", "clone_id", "compartment"], observed=True):
        lab = g.celltype.value_counts().index[0]
        if lab in allowed:
            rows.append({"patient": p, "clone": c, "tissue": tis, "error": float(lab != g.identity.iloc[0])})
    return pd.DataFrame(rows)


def _note(n_patients: int) -> dict:
    """Annotazione degli intervalli calcolati come nella tesi (min_groups=3): sotto 5 pazienti lo
    strumento non li produce; con n pazienti il bootstrap ha C(2n-1, n) ricampionamenti distinti."""
    if n_patients >= 5:
        return {}
    return {"nota": (f"intervallo della tesi con {n_patients} pazienti: lo strumento lo rifiuta (servono "
                     f"almeno 5); il bootstrap sui pazienti ha solo {comb(2 * n_patients - 1, n_patients)} "
                     f"ricampionamenti distinti")}


def summarize(u: pd.DataFrame) -> dict:
    """Tasso d'errore per compartimento e differenza appaiata, con il bootstrap percentile della tesi."""
    out = {}
    for tis in ("Tumor", "Adjacent_normal"):
        g = u[u.tissue == tis]
        b3 = cluster_bootstrap(g.error.values, g.patient.values, n_boot=N_BOOT, seed=SEED, min_groups=3)
        out[tis] = {"errore": b3.mean, "ic": [b3.ci_low, b3.ci_high], "unita": len(g),
                    "pazienti": b3.n_groups, "ic_prodotto_dallo_strumento": b3.n_groups >= 5,
                    **_note(b3.n_groups)}
    per = u.groupby(["patient", "clone", "tissue"]).error.mean().unstack("tissue")
    d = per.dropna(subset=["Tumor", "Adjacent_normal"]).reset_index()
    diff = (d.Tumor - d.Adjacent_normal).values
    b = cluster_bootstrap(diff, d.patient.values, n_boot=N_BOOT, seed=SEED, min_groups=3)
    out["Tumor-Adjacent"] = {"diff": b.mean, "ic": [b.ci_low, b.ci_high], "cloni": len(d),
                             "pazienti": b.n_groups, **_note(b.n_groups)}
    return out


LADDER = [
    ("S0 MVP (per cellula, tutte le cellule con TCR, CD4>0 non esclusivo, etichette CD4T/CD8T)",
     dict(lymphoid_only=False, exclusive_cd4=False, nk_as_error=False, min_target=1, unit="cell")),
    ("S1 + solo cellule etichettate CD4T/CD8T/NK (anche nel riferimento)",
     dict(lymphoid_only=True, exclusive_cd4=False, nk_as_error=False, min_target=1, unit="cell")),
    ("S2 + CD4 positivo solo se CD8A e CD8B negativi (regola esclusiva della tesi)",
     dict(lymphoid_only=True, exclusive_cd4=True, nk_as_error=False, min_target=1, unit="cell")),
    ("S3 + etichetta NK contata come errore",
     dict(lymphoid_only=True, exclusive_cd4=True, nk_as_error=True, min_target=1, unit="cell")),
    ("S4 + almeno 3 cellule del clone anche nel compartimento giudicato",
     dict(lymphoid_only=True, exclusive_cd4=True, nk_as_error=True, min_target=3, unit="cell")),
    ("S5 + unità = coppia clone-compartimento con l'etichetta maggioritaria (= tesi)",
     dict(lymphoid_only=True, exclusive_cd4=True, nk_as_error=True, min_target=3, unit="clone")),
]


def main() -> None:
    """Esegue la scala delle definizioni e scrive a1_moduleB.json."""
    ap = argparse.ArgumentParser()
    ap.add_argument("--vdj-manifest", type=Path, required=True)
    ap.add_argument("--out", type=Path, required=True)
    args = ap.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)
    cells = load_cells(args.vdj_manifest)
    res = {"cellule_con_TRB": len(cells),
           "discordanza": {"MVP (tutte le cellule con TCR)": discordance(cells, False),
                           "tesi (solo CD4T/CD8T/NK)": discordance(cells, True)},
           "errore_solo_sangue": {}}
    for name, kw in LADDER:
        res["errore_solo_sangue"][name] = summarize(error_units(cells, **kw))
    (args.out / "a1_moduleB.json").write_text(json.dumps(res, indent=1, ensure_ascii=False))
    print(json.dumps(res, indent=1, ensure_ascii=False))


if __name__ == "__main__":
    main()
