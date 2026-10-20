"""Verifiche di eleggibilità su GSE200996 che richiedono solo i metadati (prima dell'archivio RAW).

Legge da ``meta_dir`` i metadati per cellula degli autori (CD4/CD8/CD45 × PBMC/tumor), le due
chiavi degli hashtag (``hashtag_info_sc_GEX_PBMC``: corsie non selezionate; ``Hashtag_info_sorted_
sc_GEX_PBMC``: corsie selezionate CD38+ HLA-DR+) e le series matrix. Verifica:
- la parte di E1 che dipende dagli hashtag: pazienti con cellule di sangue assegnate a un solo
  paziente e con cellule di tumore dopo la terapia, entrambe con etichetta CD4 o CD8;
- E2: etichetta CD4/CD8 per cellula in entrambi i compartimenti (appartenenza al file CD4 o CD8);
- E5: il sangue usato contiene sia CD4 sia CD8 per paziente, e viene solo da corsie non selezionate.

Regola di identificazione delle cellule (dichiarata qui e in build_gse200996.py): l'ID di cellula
degli autori è ``<barcode 10x>_<paziente>_<prelievo>``, con prelievo B1/B2/B3 per il sangue e
Pre-Tx/Post-Tx per il tumore. Nel sangue la coppia paziente-prelievo identifica un hashtag in una
sola corsia non selezionata (verificato qui sotto), quindi la corsia di ogni cellula è determinata.

Scrive ``results/metadati.json``. Uso:
python -I validation_esterna_B/GSE200996/verifica_metadati.py ../dati_esterni/GSE200996/meta
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import pandas as pd

HERE = Path(__file__).resolve().parent
MIN_PATIENTS = 5


def read_meta(meta_dir: Path, subset: str, tissue: str) -> pd.DataFrame:
    """Metadati per cellula degli autori, con barcode, paziente e prelievo ricavati dall'ID."""
    m = pd.read_csv(meta_dir / f"GSE200996_{subset}.{tissue}.single.cell.meta.data.txt.gz", sep="\t",
                    index_col=0, low_memory=False)
    parts = m.index.to_series().str.extract(r"^([ACGT]{16})_(P\d+)_(.+)$")
    if parts.isna().any().any():
        raise ValueError(f"{subset}.{tissue}: ID di cellula fuori formato, es. {m.index[parts[0].isna()][:3].tolist()}")
    m = m.assign(barcode=parts[0].values, pid=parts[1].values, prelievo=parts[2].values)
    if not ((m["pid"] == m["Patient_ID"].astype(str)) & (m["prelievo"] == m["Stage"].astype(str))).all():
        raise ValueError(f"{subset}.{tissue}: paziente o prelievo dell'ID diversi dalle colonne Patient_ID/Stage")
    return m


def hashtag_key(path: Path) -> pd.DataFrame:
    """Chiave degli hashtag: una riga per (corsia, hashtag, campione)."""
    k = pd.read_csv(path)
    k = k[k["id"].notna()]
    lanes = [c for c in k.columns if c not in ("id", "name", "read", "pattern", "sequence", "feature_type")]
    rows = [{"corsia": lane, "hashtag": r["id"], "campione": r[lane]}
            for lane in lanes for _, r in k.iterrows() if isinstance(r[lane], str) and r[lane]]
    return pd.DataFrame(rows)


def main(meta_dir: Path) -> None:
    out: dict = {"accessione": "GSE200996"}
    key = hashtag_key(meta_dir / "GSE200996_hashtag_info_sc_GEX_PBMC.csv.gz")
    key_sorted = hashtag_key(meta_dir / "GSE200996_Hashtag_info_sorted_sc_GEX_PBMC.csv.gz")
    dup = key[key.duplicated("campione", keep=False)]
    if not dup.empty:
        raise ValueError(f"campioni con più hashtag o corsie: {dup.to_dict(orient='records')}")
    key = key.assign(paziente=key["campione"].str.split("-").str[0], prelievo=key["campione"].str.split("-").str[1])
    lane_of = key.groupby("paziente")["corsia"].nunique()
    out["hashtag"] = {
        "corsie_non_selezionate": sorted(key["corsia"].unique()),
        "campioni_per_corsia": {c: sorted(g["campione"]) for c, g in key.groupby("corsia")},
        "pazienti_in_piu_corsie": sorted(lane_of[lane_of > 1].index),
        "corsie_selezionate": sorted(key_sorted["corsia"].unique()),
        "pazienti_nelle_corsie_selezionate": sorted(key_sorted["campione"].unique()),
    }

    meta = {(s, t): read_meta(meta_dir, s, t) for s in ("CD4", "CD8", "CD45") for t in ("PBMC", "tumor")}
    for t in ("PBMC", "tumor"):
        cd4, cd8, cd45 = (meta[(s, t)] for s in ("CD4", "CD8", "CD45"))
        both = cd4.index.intersection(cd8.index)
        t_cells = cd45[cd45["CellType_ID"] == "T cell"]
        lab = cd4.index.union(cd8.index)
        out[f"etichette_{t}"] = {
            "cellule_CD45": len(cd45), "tipi_CD45": {str(k): int(v) for k, v in cd45["CellType_ID"].value_counts().items()},
            "CD4": len(cd4), "CD8": len(cd8), "in_entrambi": len(both),
            "etichettate_non_in_CD45": int((~lab.isin(cd45.index)).sum()),
            "etichettate_per_tipo_CD45": {str(k): int(v) for k, v in cd45.loc[cd45.index.isin(lab), "CellType_ID"].value_counts().items()},
            "T_CD45_senza_etichetta": int((~t_cells.index.isin(lab)).sum()),
        }
    # sangue: ogni coppia paziente-prelievo dei metadati deve esistere nella chiave non selezionata
    pb = pd.concat([meta[("CD4", "PBMC")], meta[("CD8", "PBMC")], meta[("CD45", "PBMC")]])
    pairs = set(zip(pb["pid"], pb["prelievo"]))
    in_key = set(zip(key["paziente"], key["prelievo"]))
    out["coppie_paziente_prelievo_sangue_fuori_chiave"] = sorted(f"{p}-{b}" for p, b in pairs - in_key)

    per = {}
    for s in ("CD4", "CD8"):
        for t in ("PBMC", "tumor"):
            m = meta[(s, t)]
            for (p, st), n in m.groupby(["pid", "prelievo"]).size().items():
                per.setdefault(p, {}).setdefault(f"{s}_{t}_{st}", int(n))
    pazienti = {}
    for p in sorted(per):
        v = per[p]
        blood = {s: sum(n for k, n in v.items() if k.startswith(f"{s}_PBMC")) for s in ("CD4", "CD8")}
        post = {s: v.get(f"{s}_tumor_Post-Tx", 0) for s in ("CD4", "CD8")}
        pre = {s: v.get(f"{s}_tumor_Pre-Tx", 0) for s in ("CD4", "CD8")}
        pazienti[p] = {"conteggi": v, "sangue": blood, "tumore_post": post, "tumore_pre": pre,
                       "corsia_sangue": key.loc[key["paziente"] == p, "corsia"].iloc[0] if p in set(key["paziente"]) else None,
                       "e1_metadati": sum(blood.values()) > 0 and sum(post.values()) > 0,
                       "e5_cd4_e_cd8_nel_sangue": blood["CD4"] > 0 and blood["CD8"] > 0}
    eleg = [p for p, v in pazienti.items() if v["e1_metadati"]]
    out["pazienti"] = pazienti
    out["E1_hashtag"] = {"pazienti": eleg, "n": len(eleg), "richiesti": MIN_PATIENTS, "soddisfatta": len(eleg) >= MIN_PATIENTS}
    out["E2"] = {"soddisfatta": all(out[f"etichette_{t}"]["CD4"] > 0 and out[f"etichette_{t}"]["CD8"] > 0 for t in ("PBMC", "tumor"))}
    e5_pat = [p for p in eleg if pazienti[p]["e5_cd4_e_cd8_nel_sangue"]]
    out["E5"] = {"pazienti_eleggibili_con_CD4_e_CD8_nel_sangue": e5_pat,
                 "soddisfatta": len(e5_pat) == len(eleg) and not out["coppie_paziente_prelievo_sangue_fuori_chiave"]}
    out["pazienti_con_tumore_pre"] = [p for p in eleg if sum(pazienti[p]["tumore_pre"].values()) > 0]

    (HERE / "results").mkdir(exist_ok=True)
    (HERE / "results" / "metadati.json").write_text(json.dumps(out, indent=2, ensure_ascii=False) + "\n")
    for t in ("PBMC", "tumor"):
        print(t, out[f"etichette_{t}"])
    print("hashtag:", {k: v for k, v in out["hashtag"].items() if k != "campioni_per_corsia"})
    print("coppie paziente-prelievo del sangue fuori chiave:", out["coppie_paziente_prelievo_sangue_fuori_chiave"])
    for p in eleg:
        v = pazienti[p]
        print(f"  {p} corsia {v['corsia_sangue']}: sangue {v['sangue']}, tumore post {v['tumore_post']}, pre {v['tumore_pre']}")
    print(f"E1 (parte hashtag): {len(eleg)} pazienti; E2: {out['E2']['soddisfatta']}; E5: {out['E5']['soddisfatta']}; "
          f"tumore pre-terapia per {out['pazienti_con_tumore_pre']}")


if __name__ == "__main__":
    main(Path(sys.argv[1]))
