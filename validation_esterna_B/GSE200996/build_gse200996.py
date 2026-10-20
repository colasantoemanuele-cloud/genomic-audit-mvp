"""Costruisce gli input del Modulo B per GSE200996 (validation_esterna_B/GSE200996/CRITERI.md).

Uscite in ``out_dir`` (fuori dal repository):
- ``GSE200996_principale.h5ad``: PBMC (tutti i prelievi B1-B3 delle corsie non selezionate) + Tumor
  dopo la terapia (Post-Tx), solo cellule etichettate CD4T/CD8T secondo MAPPA_ETICHETTE.md;
- ``GSE200996_descrittivo.h5ad``: PBMC + Tumor prima della terapia (Pre-Tx), per i pazienti che
  hanno sia la matrice sia il TCR single-cell del tumore prima della terapia;
- ``vdj/`` e ``manifest_principale.csv``, ``manifest_descrittivo.csv``: manifest VDJ nel formato
  di docs/GUIDA_UTENTE.md §3.3 (path,patient,compartment).
Scrive anche ``results/build.json`` con i conteggi per campione (condizioni E1, E3, E4) e li stampa.

Regole di abbinamento dei barcode (condizione E3):
- ID di cellula degli autori = ``<barcode 16 basi>_<paziente>_<prelievo>``; nella matrice 10x e nel
  file VDJ dello stesso campione la cellula ha barcode ``<barcode 16 basi>-1``.
- Tumore: un campione (matrice .h5 e file VDJ) per paziente e prelievo; cellula = (campione, barcode).
- Sangue: un campione per corsia con più pazienti e prelievi (hashtag). La corsia di una cellula è
  quella della chiave ``GSE200996_hashtag_info_sc_GEX_PBMC.csv.gz`` che contiene la sua coppia
  paziente-prelievo (ogni paziente sta in una sola corsia). Il file VDJ della corsia è diviso per
  paziente: un contig va al paziente a cui gli autori hanno assegnato il suo barcode in quella corsia
  (file CD45.PBMC, cioè tutte le cellule demultiplate, non solo le T etichettate); i contig dei
  barcode non assegnati (doppietti, negativi, cellule tolte dal controllo di qualità) sono scartati e
  contati. Le righe tenute sono copiate senza modifiche (stesse colonne): nessuna informazione persa
  per le cellule assegnate (E4).
- Nella matrice si tengono solo le feature «Gene Expression»; le altre (Antibody Capture) sono
  elencate in build.json (B-E5).

Uso: python -I validation_esterna_B/GSE200996/build_gse200996.py <meta_dir> <raw_dir> <out_dir>
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import anndata as ad
import h5py
import numpy as np
import pandas as pd
import scipy.sparse as sp

HERE = Path(__file__).resolve().parent
LABELS = {"CD4": "CD4T", "CD8": "CD8T"}


def read_meta(meta_dir: Path, subset: str, tissue: str) -> pd.DataFrame:
    """Metadati degli autori con barcode (con suffisso -1), paziente e prelievo dall'ID."""
    m = pd.read_csv(meta_dir / f"GSE200996_{subset}.{tissue}.single.cell.meta.data.txt.gz", sep="\t",
                    index_col=0, low_memory=False)
    parts = m.index.to_series().str.extract(r"^([ACGT]{16})_(P\d+)_(.+)$")
    assert not parts.isna().any().any(), f"ID fuori formato in {subset}.{tissue}"
    return m.assign(barcode=parts[0].values + "-1", patient=parts[1].values, timepoint=parts[2].values)


def read_10x_h5(path: Path, keep: pd.Index) -> tuple[sp.csr_matrix, list[str], dict]:
    """Righe (barcode) ``keep`` della matrice 10x, solo feature Gene Expression. Controlla che i
    conteggi siano interi non negativi. Restituisce matrice cellule × geni, simboli e informazioni."""
    with h5py.File(path, "r") as f:
        g = f["matrix"]
        data, indices, indptr = g["data"][:], g["indices"][:], g["indptr"][:]
        shape = tuple(g["shape"][:])
        barcodes = pd.Index([b.decode() for b in g["barcodes"][:]])
        names = [x.decode() for x in g["features"]["name"][:]]
        ftype = np.array([x.decode() for x in g["features"]["feature_type"][:]])
    m = sp.csc_matrix((data, indices, indptr), shape=shape)  # geni × barcode
    pos = barcodes.get_indexer(keep)
    missing = int((pos < 0).sum())
    sub = m[:, pos[pos >= 0]].T.tocsr()
    gex = np.where(ftype == "Gene Expression")[0]
    sub = sub[:, gex]
    info = {"dtype": str(data.dtype), "interi_non_negativi": bool(np.all(data >= 0) and np.all(np.mod(data, 1) == 0)),
            "barcode_nella_matrice": int(len(barcodes)), "cellule_cercate": int(len(keep)), "cellule_non_trovate": missing,
            "feature_per_tipo": {t: int((ftype == t).sum()) for t in np.unique(ftype)},
            "feature_non_GEX": [n for n, t in zip(names, ftype) if t != "Gene Expression"]}
    return sub.astype(np.float32), [names[i] for i in gex], info


def productive_cells(vdj: pd.DataFrame, chain: str | None) -> set[str]:
    """Barcode con almeno un contig is_cell, high_confidence, full_length, productive (della catena
    indicata, o TRA/TRB se None) e CDR3 nucleotidica presente."""
    ok = vdj
    for c in ("is_cell", "high_confidence", "full_length", "productive"):
        if c in ok.columns:
            ok = ok[ok[c].astype(str).str.lower().isin(["true", "1"])]
    ok = ok[ok["chain"].isin([chain] if chain else ["TRA", "TRB"])]
    ok = ok[ok["cdr3_nt"].notna() & (ok["cdr3_nt"].astype(str) != "None")]
    return set(ok["barcode"].astype(str))


def main(meta_dir: Path, raw_dir: Path, out_dir: Path) -> None:
    # percorsi assoluti: nel manifest i percorsi relativi sarebbero letti rispetto alla sua cartella
    meta_dir, raw_dir, out_dir = meta_dir.resolve(), raw_dir.resolve(), out_dir.resolve()
    out_dir.mkdir(parents=True, exist_ok=True)
    (out_dir / "vdj").mkdir(exist_ok=True)
    elig = json.loads((HERE / "results" / "metadati.json").read_text(encoding="utf-8"))
    patients = elig["E1_hashtag"]["pazienti"]
    files = {p.name: p for p in raw_dir.iterdir()}

    def find(*tokens: str) -> Path | None:
        hits = [p for n, p in files.items() if all(t in n for t in tokens)]
        assert len(hits) <= 1, (tokens, hits)
        return hits[0] if hits else None

    lab = {(s, t): read_meta(meta_dir, s, t) for s in ("CD4", "CD8") for t in ("PBMC", "tumor")}
    cd45 = {t: read_meta(meta_dir, "CD45", t) for t in ("PBMC", "tumor")}
    key = pd.read_csv(meta_dir / "GSE200996_hashtag_info_sc_GEX_PBMC.csv.gz")
    key = key[key["id"].notna()]
    lanes = [c for c in key.columns if c not in ("id", "name", "read", "pattern", "sequence", "feature_type")]
    lane_of = {str(v).split("-")[0]: lane for lane in lanes for v in key[lane].dropna()}

    labelled = []
    for (s, t), m in lab.items():
        labelled.append(m.assign(label=LABELS[s], tissue=t))
    labelled = pd.concat(labelled)
    assert not labelled.index.duplicated().any()

    samples, blocks, genes_ref = [], [], None
    vdj_rows: list[dict] = []

    def add_block(m: pd.DataFrame, compartment: str, h5: Path, sample: str, lane: str) -> None:
        nonlocal genes_ref
        X, genes, info = read_10x_h5(h5, pd.Index(m["barcode"]))
        assert info["cellule_non_trovate"] == 0, (sample, info["cellule_non_trovate"])
        if genes_ref is None:
            genes_ref = genes
        assert genes == genes_ref, f"geni diversi in {h5.name}"
        obs = pd.DataFrame({"patient": m["patient"].values, "compartment": compartment,
                            "timepoint": m["timepoint"].values, "label": m["label"].values,
                            "barcode": m["barcode"].values, "sample": sample, "lane": lane}, index=m.index)
        blocks.append((X, obs))
        samples.append({"campione": sample, "file": h5.name, "compartimento": compartment, **info})

    # ---- tumore: un campione per paziente e prelievo
    for stage, tag in (("Post-Tx", "post-Tx"), ("Pre-Tx", "pre-Tx")):
        for p in patients:
            m = labelled[(labelled.tissue == "tumor") & (labelled.patient == p) & (labelled.timepoint == stage)]
            h5 = find(f"_{p}_{tag}_GEX_sc_tumor.h5")
            vdj = find(f"_{p}_{tag}_TCR_sc_tumor.csv.gz")
            if m.empty or h5 is None:
                continue
            if stage == "Pre-Tx" and vdj is None:
                vdj_rows.append({"campione": f"{p}_{stage}", "nota": "matrice senza TCR single-cell: escluso"})
                continue
            add_block(m, "Tumor", h5, f"{p}_{stage}", "")
            v = pd.read_csv(vdj, low_memory=False)
            allc = cd45["tumor"][(cd45["tumor"].patient == p) & (cd45["tumor"].timepoint == stage)]
            vdj_rows.append({"campione": f"{p}_{stage}", "patient": p, "compartment": "Tumor", "timepoint": stage,
                             "path": str(vdj), "righe_vdj": len(v), "barcode_vdj": int(v["barcode"].nunique()),
                             "cellule_CD45": len(allc), "T_CD45_non_etichettate": int(
                                 ((allc["CellType_ID"] == "T cell") & ~allc.index.isin(labelled.index)).sum()),
                             "etichettate": len(m),
                             "etichettate_con_TRB_produttiva": int(m["barcode"].isin(productive_cells(v, "TRB")).sum()),
                             "etichettate_con_TCR_produttivo": int(m["barcode"].isin(productive_cells(v, None)).sum())})

    # ---- sangue: un campione per corsia, diviso per paziente
    for lane in sorted({lane_of[p] for p in patients}):
        h5 = find(f"_{lane}_GEX_sc_PBMC.h5")
        vdjf = find(f"_{lane}_TCR_sc_PBMC.csv.gz")
        assert h5 is not None and vdjf is not None and "sorted" not in h5.name, lane
        lane_pat = [p for p, ln in lane_of.items() if ln == lane]
        m = labelled[(labelled.tissue == "PBMC") & labelled.patient.isin(lane_pat) & labelled.patient.isin(patients)]
        add_block(m, "PBMC", h5, f"PBMC_{lane}", lane)
        v = pd.read_csv(vdjf, low_memory=False)
        allc = cd45["PBMC"][cd45["PBMC"].patient.isin(lane_pat)]
        assert not allc["barcode"].duplicated().any(), f"barcode ripetuti nella corsia {lane}"
        owner = dict(zip(allc["barcode"], allc["patient"]))
        v = v.assign(_p=v["barcode"].map(owner))
        cells_vdj = v.drop_duplicates("barcode")
        isc = cells_vdj["is_cell"].astype(str).str.lower().isin(["true", "1"]) if "is_cell" in v.columns else True
        samples[-1]["vdj_corsia"] = {"barcode": int(len(cells_vdj)), "barcode_is_cell": int(np.sum(isc)),
                                    "barcode_assegnati": int(cells_vdj["_p"].notna().sum()),
                                    "barcode_is_cell_non_assegnati": int((cells_vdj["_p"].isna() & isc).sum())}
        for p in lane_pat:
            if p not in patients:
                continue
            vp = v[v["_p"] == p].drop(columns="_p")
            path = out_dir / "vdj" / f"{p}_PBMC_filtered_contig_annotations.csv"
            vp.to_csv(path, index=False)
            mp = m[m.patient == p]
            ap = allc[allc.patient == p]
            vdj_rows.append({"campione": f"{p}_PBMC", "patient": p, "compartment": "PBMC", "timepoint": "B1-B3",
                             "lane": lane, "path": str(path), "righe_vdj": len(vp), "barcode_vdj": int(vp["barcode"].nunique()),
                             "cellule_CD45": len(ap), "cellule_per_prelievo": {k: int(n) for k, n in ap["timepoint"].value_counts().items()},
                             "T_CD45_non_etichettate": int(((ap["CellType_ID"] == "T cell") & ~ap.index.isin(labelled.index)).sum()),
                             "etichettate": len(mp),
                             "etichettate_con_TRB_produttiva": int(mp["barcode"].isin(productive_cells(vp, "TRB")).sum()),
                             "etichettate_con_TCR_produttivo": int(mp["barcode"].isin(productive_cells(vp, None)).sum())})

    X = sp.vstack([b[0] for b in blocks]).tocsr()
    obs = pd.concat([b[1] for b in blocks])
    full = ad.AnnData(X=X, obs=obs, var=pd.DataFrame(index=pd.Index(genes_ref, name="gene")))
    full.var_names_make_unique()
    vr = pd.DataFrame([r for r in vdj_rows if "path" in r])
    pre_pat = sorted(set(vr.loc[(vr.compartment == "Tumor") & (vr.timepoint == "Pre-Tx"), "patient"]))
    outputs = {}
    for name, stage, pats in (("principale", "Post-Tx", patients), ("descrittivo", "Pre-Tx", pre_pat)):
        keep = full.obs["patient"].isin(pats) & ((full.obs["compartment"] == "PBMC") | (full.obs["timepoint"] == stage))
        a = full[keep.values].copy()
        for c in ("patient", "compartment", "timepoint", "label", "sample", "lane"):
            a.obs[c] = a.obs[c].astype(str)
        a.write_h5ad(out_dir / f"GSE200996_{name}.h5ad", compression="gzip")
        man = vr[vr.patient.isin(pats) & ((vr.compartment == "PBMC") | (vr.timepoint == stage))]
        man[["path", "patient", "compartment"]].to_csv(out_dir / f"manifest_{name}.csv", index=False)
        outputs[name] = {"cellule": int(a.n_obs), "geni": int(a.n_vars), "pazienti": sorted(set(a.obs.patient)),
                         "per_compartimento_ed_etichetta": {f"{c}|{lab_}": int(n) for (c, lab_), n in
                                                             a.obs.groupby(["compartment", "label"]).size().items()},
                         "file_vdj": len(man)}

    def pct(r: dict, k: str) -> float:
        return r[k] / r["etichettate"] if r.get("etichettate") else float("nan")

    per_comp = {}
    for comp in ("PBMC", "Tumor"):
        sub = vr[(vr.compartment == comp) & ((vr.timepoint != "Pre-Tx"))]
        per_comp[comp] = {"etichettate": int(sub["etichettate"].sum()),
                          "con_TRB_produttiva": int(sub["etichettate_con_TRB_produttiva"].sum()),
                          "con_TCR_produttivo": int(sub["etichettate_con_TCR_produttivo"].sum())}
    res = {"campioni_matrici": samples, "campioni_vdj": vdj_rows, "uscite": outputs, "E3_principale": per_comp,
           "pazienti_pre_tx_descrittivo": pre_pat}
    (HERE / "results" / "build.json").write_text(json.dumps(res, indent=2, ensure_ascii=False, default=str) + "\n")

    print("campione | compartimento | cellule CD45 assegnate | T etichettate | T CD45 non etichettate (escluse) | "
          "con TRB produttiva | con TCR produttivo")
    for r in vdj_rows:
        if "path" not in r:
            print(f"{r['campione']}: {r['nota']}")
            continue
        print(f"{r['campione']:12s} {r['compartment']:6s} {r['cellule_CD45']:6d} {r['etichettate']:6d} "
              f"{r['T_CD45_non_etichettate']:6d} {r['etichettate_con_TRB_produttiva']:6d} ({pct(r, 'etichettate_con_TRB_produttiva'):.1%}) "
              f"{r['etichettate_con_TCR_produttivo']:6d} ({pct(r, 'etichettate_con_TCR_produttivo'):.1%})")
    for s in samples:
        print(f"[matrice] {s['file']}: dtype {s['dtype']}, interi {s['interi_non_negativi']}, feature {s['feature_per_tipo']}, "
              f"non GEX {s['feature_non_GEX'][:8]}" + (f", VDJ corsia {s['vdj_corsia']}" if "vdj_corsia" in s else ""))
    for comp, v in per_comp.items():
        print(f"E3 {comp} (analisi principale): {v['con_TRB_produttiva']}/{v['etichettate']} "
              f"({v['con_TRB_produttiva'] / v['etichettate']:.1%}) con TRB produttiva; "
              f"{v['con_TCR_produttivo'] / v['etichettate']:.1%} con TCR produttivo")
    print(json.dumps(outputs, indent=1))


if __name__ == "__main__":
    main(Path(sys.argv[1]), Path(sys.argv[2]), Path(sys.argv[3]))
