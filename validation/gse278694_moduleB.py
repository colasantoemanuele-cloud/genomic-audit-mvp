"""Verifica di fedeltà del porting del Modulo B su GSE278694 (non una validazione del metodo).

Legge (sola lettura) i dati di pdac-ml e confronta le definizioni del Modulo B con quelle della
tesi. Riferimento: versione finale di progetto-ml (ramo ``revisione``, commit 05947d0),
``src/10_loco.py`` schema ``blood`` e le tabelle ``results/tables/error_rate_blood.csv``,
``error_rate_tie_sensitivity.csv`` ed ``error_rate_clone_ab.csv``, lette da qui se presenti.
Non tara nessun parametro: parte dalle definizioni dello strumento e le sostituisce UNA ALLA VOLTA
con quelle della tesi (S0-S5), fino alla definizione della tesi.

Regola dei pareggi della tesi, la stessa dello strumento: se due o più etichette hanno lo stesso
conteggio massimo, l'unità clone-compartimento è esclusa; come sensibilità è contata come errore o
come corretta. La stima puntuale e il numero di cloni devono coincidere con la tesi; l'intervallo
no, perché il metodo è diverso (qui: ``ic`` = bootstrap percentile sui pazienti, seme 0, 2000
repliche, come la tesi prima della revisione; ``ic_strumento`` = intervallo t cluster-robusto
dello strumento; la tesi oggi usa un bootstrap «expanded percentile»). Il vecchio valore della
tesi, 0,1951 [0,1230; 0,2737] su 288 unità, che risolveva i pareggi con l'ordine delle righe, è
SUPERATO.

Le unità che differiscono dalla tesi (``arbitration_blood.csv``) sono elencate con le cellule che
hanno due contig TRB a pari UMI: lo strumento sceglie il contig con più reads (poi la CDR3 minore
in ordine alfabetico), la tesi quello che un ordinamento non stabile lascia per primo
(``sort_values("umis")`` + ``drop_duplicates``). ``thesis_contig_choice`` ricostruisce quella
scelta e la sostituisce a quella dello strumento, lasciando invariato il resto del codice: se i
numeri tornano uguali a quelli della tesi, la differenza viene solo dalla scelta dei contig.

Gli intervalli con meno di 5 pazienti qui riportati (min_groups=3) sono quelli della tesi: lo
strumento li rifiuta, e con 4 pazienti esistono solo 35 ricampionamenti distinti.

Uso:
    python -m validation.gse278694_moduleB --vdj-manifest <manifest.csv> --out <dir>

``manifest.csv``: path,patient,compartment dei CSV VDJ originali (validation.make_vdj_manifest);
il suffisso "-1" dei barcode è normalizzato come nello strumento (match_barcodes).
"""

from __future__ import annotations

import argparse
import json
import os
import re
from math import comb
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd

from core.io import read_vdj_manifest
from core.stats import cluster_bootstrap, cluster_t_interval
from core.tcr_validation import (
    _KEY,
    _strip,
    build_clonotypes,
    match_barcodes,
    pairwise_excess_discordance,
    parse_vdj_contigs,
)

# Dati reali di GSE278694: non inclusi nella repo. Cartella indicata da GSE278694_DIR, altrimenti
# una cartella "pdac-ml" accanto a questa repo.
PDAC = Path(os.environ.get("GSE278694_DIR", Path(__file__).resolve().parents[2] / "pdac-ml"))
H5AD = PDAC / "data/interim/sc_raw.h5ad"
TESI = PDAC / "results/tables"
N_BOOT, SEED = 2000, 0
LYMPHOID = ("CD4T", "CD8T", "NK")
TISSUE = {"Tumor": "Tumor", "Normal": "Adjacent_normal", "PBMC": "PBMC"}


def load_cells(manifest: Path) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Una riga per cellula con TCR (TRB): stessa ingestione e stessa merge dello strumento, più
    le conte grezze di CD4, CD8A, CD8B e la CDR3 della TRA scelta (per il clone TRA+TRB). Restituisce
    anche le cellule con due o più contig TRB a pari UMI massime (chiave con barcode normalizzato)."""
    adata = ad.read_h5ad(H5AD)
    obs = adata.obs[["patients", "tissue", "all_celltype", "barcode"]].astype(str).reset_index(drop=True)
    obs.columns = ["patient", "compartment", "celltype", "barcode"]
    for g in ("CD4", "CD8A", "CD8B"):
        i = adata.var_names.get_loc(g)
        obs[g] = np.asarray(adata.X[:, i].todense()).ravel()
    del adata
    files = read_vdj_manifest(manifest)
    contigs = parse_vdj_contigs(files)
    obs, clones, _ = match_barcodes(obs, build_clonotypes(contigs))
    cells = obs.merge(clones, on=_KEY, how="inner")
    tra = contigs[contigs.chain == "TRA"].assign(barcode=lambda d: _strip(d.barcode))
    cells = cells.merge(tra[_KEY + ["cdr3_nt"]].rename(columns={"cdr3_nt": "cdr3_TRA"}), on=_KEY, how="left")
    return cells, _trb_umi_ties(files)


def _trb_umi_ties(files: list) -> pd.DataFrame:
    """Cellule (paziente, compartimento, barcode senza suffisso) con due o più CDR3 TRB a pari UMI."""
    frames = []
    for path, patient, comp, *_ in files:
        df = pd.read_csv(path, low_memory=False)
        for col in ("is_cell", "high_confidence", "productive", "full_length"):
            if col in df.columns:
                df = df[df[col].astype(str).str.lower().isin(["true", "1"])]
        df = df[(df.chain == "TRB") & df.cdr3_nt.notna() & (df.cdr3_nt.astype(str) != "None")]
        frames.append(df.assign(patient=patient, compartment=comp, barcode=_strip(df.barcode.astype(str))))
    c = pd.concat(frames)
    key = ["patient", "compartment", "barcode"]
    top = c[c.umis == c.groupby(key).umis.transform("max")]
    n = top.groupby(key).cdr3_nt.nunique()
    return top.set_index(key).loc[n[n > 1].index].reset_index()[key + ["cdr3_nt", "umis", "reads"]]


def discordance(cells: pd.DataFrame, lymphoid_only: bool) -> dict:
    """Eccesso di discordanza con il bootstrap percentile della tesi."""
    x = cells[cells.celltype.isin(LYMPHOID)] if lymphoid_only else cells
    pw = pairwise_excess_discordance(x)
    b = cluster_bootstrap(pw.excess.values, pw.patient.values, n_boot=N_BOOT, seed=SEED, min_groups=5)
    return {"cellule": len(x), "coppie": len(pw), "pazienti": int(pw.patient.nunique()),
            "d_within": float(pw.d_within.mean()), "eccesso": b.mean, "ic": [b.ci_low, b.ci_high]}


def error_units(cells: pd.DataFrame, *, lymphoid_only: bool, exclusive_cd4: bool, nk_as_error: bool,
                min_target: int, unit: str, tie_rule: str = "exclude") -> pd.DataFrame:
    """Unità di errore (cellule o coppie clone-compartimento) con l'identità dal solo sangue.
    Con tutte le opzioni a False/0/'cell' coincide con marker_error_rate dello strumento. Per le
    unità clone-compartimento, pareggi dell'etichetta di maggioranza come in 10_loco.py: esclusi
    (``tie_rule="exclude"``) oppure contati come errore ("error") o come corretti ("correct");
    un'unità in pareggio è tenuta se almeno una delle etichette a pari merito è ammessa."""
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
        vc = g.celltype.value_counts()
        top = list(vc.index[vc.values == vc.values.max()])
        if not set(top) & allowed:
            continue
        if len(top) == 1:
            err = float(top[0] != g.identity.iloc[0])
        elif tie_rule == "exclude":
            err = float("nan")
        else:
            err = 1.0 if tie_rule == "error" else 0.0
        rows.append({"patient": p, "clone": c, "tissue": tis, "error": err, "tie": len(top) > 1})
    return pd.DataFrame(rows, columns=["patient", "clone", "tissue", "error", "tie"])


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
    n_ties = u.groupby("tissue").tie.sum() if "tie" in u else pd.Series(dtype=int)
    u = u.dropna(subset=["error"])
    for tis in ("Tumor", "Adjacent_normal"):
        g = u[u.tissue == tis]
        b3 = cluster_bootstrap(g.error.values, g.patient.values, n_boot=N_BOOT, seed=SEED, min_groups=3)
        t = cluster_t_interval(g.error.values, g.patient.values, 5, scale="logit")
        out[tis] = {"errore": b3.mean, "ic": [b3.ci_low, b3.ci_high], "unita": len(g),
                    "pareggi": int(n_ties.get(tis, 0)), "pazienti": b3.n_groups,
                    "ic_prodotto_dallo_strumento": b3.n_groups >= 5,
                    "ic_strumento": [t.ci_low, t.ci_high] if t.sufficient else None, **_note(b3.n_groups)}
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
    ("S5 + unità = coppia clone-compartimento con l'etichetta maggioritaria, pareggi esclusi (= tesi)",
     dict(lymphoid_only=True, exclusive_cd4=True, nk_as_error=True, min_target=3, unit="clone")),
]
S5 = LADDER[-1][1]
SENSIBILITA = {
    "S5, pareggi contati come errore": dict(S5, tie_rule="error"),
    "S5, pareggi contati come corretti": dict(S5, tie_rule="correct"),
}
# righe delle tabelle della tesi (results/tables) corrispondenti a S5 e alle sensibilità
TESI_RIGHE = {
    "S5": ("error_rate_tie_sensitivity.csv", "ties_excluded"),
    "S5, pareggi contati come errore": ("error_rate_tie_sensitivity.csv", "ties_as_error"),
    "S5, pareggi contati come corretti": ("error_rate_tie_sensitivity.csv", "ties_as_correct"),
}


def thesis_value(name: str) -> dict | None:
    """Tumor e Tumor-Adjacent della tesi (schema blood) per S5 e le sensibilità; None se le tabelle
    di progetto-ml non sono presenti."""
    f, version = TESI_RIGHE[name]
    if not (TESI / f).exists():
        return None
    t = pd.read_csv(TESI / f)
    t = t[(t.scheme == "blood") & (t.version == version)].set_index("quantity")
    r, d = t.loc["error_Tumor"], t.loc["diff_Tumor_minus_Adjacent_normal"]
    return {"Tumor": {"errore": float(r.value), "unita": int(r.n_units), "pareggi": int(r.n_ties),
                      "pazienti": int(r.n_patients), "ic_tesi": [float(r.lo), float(r.hi)]},
            "Tumor-Adjacent": {"diff": float(d.value), "cloni": int(d.n_units), "pazienti": int(d.n_patients)},
            "fonte": f"progetto-ml results/tables/{f}, scheme=blood, version={version}"}


def compare(ours: dict, tesi: dict | None) -> dict:
    """Coincidenza della stima puntuale (a 1e-12) e del numero di unità/cloni con la tesi."""
    if tesi is None:
        return {"tesi": None, "nota": "tabelle di progetto-ml non presenti: confronto non eseguito"}
    t, d = ours["Tumor"], ours["Tumor-Adjacent"]
    return {"tesi": tesi,
            "coincide_Tumor": abs(t["errore"] - tesi["Tumor"]["errore"]) < 1e-12
            and t["unita"] == tesi["Tumor"]["unita"] and t["pazienti"] == tesi["Tumor"]["pazienti"],
            "coincide_Tumor-Adjacent": abs(d["diff"] - tesi["Tumor-Adjacent"]["diff"]) < 1e-12
            and d["cloni"] == tesi["Tumor-Adjacent"]["cloni"]}


def unit_differences(cells: pd.DataFrame, ties: pd.DataFrame) -> list[dict]:
    """Unità di S5 presenti solo nello strumento o solo nella tesi (arbitration_blood.csv), con le
    cellule dell'unità che hanno due contig TRB a pari UMI."""
    f = TESI / "arbitration_blood.csv"
    if not f.exists():
        return []
    ours = error_units(cells, **S5)[["patient", "clone", "tissue", "tie"]]
    th = pd.read_csv(f)[["patient", "clone", "tissue", "tie"]]
    m = ours.merge(th, on=["patient", "clone", "tissue"], how="outer", indicator=True, suffixes=("", "_tesi"))
    out = []
    for _, r in m[m._merge != "both"].iterrows():
        g = cells[(cells.patient == r.patient) & (cells.compartment == r.tissue)]
        same = g[g.clone_id == r.clone]
        tied = ties.merge(g[["patient", "compartment", "barcode"]].drop_duplicates(),
                          on=["patient", "compartment", "barcode"])
        tied = tied[tied.barcode.isin(same.barcode) | tied.cdr3_nt.eq(r.clone)]
        out.append({"paziente": r.patient, "clone": r.clone, "compartimento": r.tissue,
                    "presente_in": "strumento" if r._merge == "left_only" else "tesi",
                    "pareggio_maggioranza": bool(r.tie if r._merge == "left_only" else r.tie_tesi),
                    "etichette_nello_strumento": same.celltype.value_counts().to_dict(),
                    "cellule_con_contig_TRB_a_pari_UMI": tied.to_dict("records")})
    return out


def thesis_contig_choice(cells: pd.DataFrame) -> pd.DataFrame:
    """Le stesse cellule con la CDR3 TRB e TRA scelte come in 10_loco.load_tcr della tesi:
    ordinamento per UMI decrescenti (non stabile) e prima riga per (paziente, tessuto, barcode,
    catena), file letti in ordine alfabetico."""
    frames = []
    for p in sorted((PDAC / "data/raw/tcr").glob("*_filtered_contig_annotations.csv.gz")):
        m = re.match(r"GSM\d+_PA_(\d+)-(\w+)-VDJ_filtered_contig_annotations\.csv\.gz", p.name)
        if not m:
            continue
        df = pd.read_csv(p, compression="gzip", low_memory=False)
        for col in ("is_cell", "high_confidence", "productive", "full_length"):
            if col in df.columns:
                df = df[df[col].astype(str).str.lower().isin(["true", "1"])]
        df = df[df.chain.isin(["TRA", "TRB"])]
        df = df[df.cdr3_nt.notna() & (df.cdr3_nt.astype(str) != "None")]
        df["patient"] = f"PA{int(m.group(1)):02d}"
        df["compartment"] = TISSUE[m.group(2)]
        df["barcode"] = df.barcode.astype(str).str.extract(r"([ACGT]{16})", expand=False)
        frames.append(df[["patient", "compartment", "barcode", "chain", "cdr3_nt", "umis"]])
    c = pd.concat(frames, ignore_index=True).dropna(subset=["barcode"])
    c = c.sort_values("umis", ascending=False).drop_duplicates(["patient", "compartment", "barcode", "chain"])
    th = c.pivot_table(index=["patient", "compartment", "barcode"], columns="chain", values="cdr3_nt",
                       aggfunc="first").reset_index()
    j = cells.merge(th, on=["patient", "compartment", "barcode"], how="left")
    return j.assign(clone_id=j.TRB, cdr3_TRA=j.TRA).drop(columns=["TRA", "TRB"])


def clone_ab(cells: pd.DataFrame) -> pd.DataFrame:
    """Clone definito da TRA+TRB (solo cellule con entrambe le catene), come clone_ab della tesi."""
    x = cells[cells.cdr3_TRA.notna()]
    return x.assign(clone_id=x.cdr3_TRA + "|" + x.clone_id)


def main() -> None:
    """Esegue la scala delle definizioni e scrive a1_moduleB.json."""
    ap = argparse.ArgumentParser()
    ap.add_argument("--vdj-manifest", type=Path, required=True)
    ap.add_argument("--out", type=Path, required=True)
    args = ap.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)
    cells, ties = load_cells(args.vdj_manifest)
    res = {"cellule_con_TRB": len(cells),
           "discordanza": {"MVP (tutte le cellule con TCR)": discordance(cells, False),
                           "tesi (solo CD4T/CD8T/NK)": discordance(cells, True)},
           "errore_solo_sangue": {}}
    for name, kw in LADDER:
        res["errore_solo_sangue"][name] = summarize(error_units(cells, **kw))
    s5 = res["errore_solo_sangue"][LADDER[-1][0]]
    res["fedelta_S5"] = {"S5": compare(s5, thesis_value("S5"))}
    for name, kw in SENSIBILITA.items():
        r = summarize(error_units(cells, **kw))
        res["fedelta_S5"][name] = {"strumento": r, **compare(r, thesis_value(name))}
    res["fedelta_S5"]["unita_diverse_dalla_tesi"] = unit_differences(cells, ties)
    ab = summarize(error_units(clone_ab(cells), **S5))
    tesi_ab = None
    if (TESI / "error_rate_clone_ab.csv").exists():
        t = pd.read_csv(TESI / "error_rate_clone_ab.csv")
        t = t[(t.clone_def == "clone_ab") & (t.scheme == "blood") & (t.tissue == "Tumor")].iloc[0]
        tesi_ab = {"errore": float(t.error), "unita": int(t.n), "pareggi": int(t.n_ties), "pazienti": int(t.n_pat),
                   "ic_tesi": [float(t.lo), float(t.hi)],
                   "fonte": "progetto-ml results/tables/error_rate_clone_ab.csv, clone_ab, blood, Tumor"}
    res["fedelta_S5"]["clone TRA+TRB"] = {
        "strumento": ab, "tesi": tesi_ab,
        "coincide_Tumor": None if tesi_ab is None else bool(
            abs(ab["Tumor"]["errore"] - tesi_ab["errore"]) < 1e-12 and ab["Tumor"]["unita"] == tesi_ab["unita"])}
    # stesso codice, contig scelti come nella tesi: isola l'effetto della scelta dei contig a pari UMI
    th = thesis_contig_choice(cells)
    n_trb = int((th.clone_id != cells.clone_id).sum())
    n_tra = int((th.cdr3_TRA.fillna("") != cells.cdr3_TRA.fillna("")).sum())
    out = {}
    for name, x, rule in [("S5", th, "exclude"), ("S5, pareggi contati come errore", th, "error"),
                          ("S5, pareggi contati come corretti", th, "correct"), ("clone TRA+TRB", clone_ab(th), "exclude")]:
        t = summarize(error_units(x.dropna(subset=["clone_id"]), **dict(S5, tie_rule=rule)))["Tumor"]
        ref = tesi_ab if name == "clone TRA+TRB" else (thesis_value(name) or {}).get("Tumor")
        out[name] = {"errore": t["errore"], "unita": t["unita"], "pareggi": t["pareggi"], "pazienti": t["pazienti"],
                     "coincide_con_la_tesi": None if ref is None else bool(
                         abs(t["errore"] - ref["errore"]) < 1e-12 and t["unita"] == ref["unita"])}
    res["fedelta_S5"]["con la scelta dei contig della tesi"] = {
        "cellule_con_CDR3_TRB_diversa": n_trb, "cellule_con_CDR3_TRA_diversa": n_tra, "Tumor": out}
    (args.out / "a1_moduleB.json").write_text(json.dumps(res, indent=1, ensure_ascii=False))
    print(json.dumps(res, indent=1, ensure_ascii=False))


if __name__ == "__main__":
    main()
