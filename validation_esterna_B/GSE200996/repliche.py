"""B-E3, B-E4 e B-E6 su GSE200996 (CRITERI.md) in un solo processo, più i numeri dello strumento per B-E2.

Carica una volta l'h5ad principale e i file VDJ del manifest e costruisce una volta la tabella
cellula-clone-etichetta con le funzioni pubbliche del Modulo B usate da ``run_tcr_validation``
(chiamata da cli.py): ``parse_vdj_contigs``, ``build_clonotypes``, ``compute_marker_positivity``,
``match_barcodes``. Ogni replica modifica solo quella tabella (compartimento, etichetta o clone) e
richiama ``pairwise_excess_discordance``, ``clone_level_units`` e ``cluster_t_interval`` con i default
della 0.9.0. Controllo: sui dati originali i numeri coincidono esattamente con ``run_tcr_validation``
(stesso percorso della riga di comando) e, arrotondati, con l'output della CLI di B-E1.

Semi: B-E3 permutazioni con semi 0-19 (``numpy.random.default_rng(seme)``); B-E4 un solo ordinamento
casuale delle unità candidate con seme 0, il 10% sono le prime round(0,10 N) e il 20% le prime
round(0,20 N) (dosi annidate).

Variante TRA+TRB di B-E6: clone = coppia (CDR3 nucleotidica TRB, CDR3 nucleotidica TRA), con il contig
per catena scelto da ``parse_vdj_contigs``; le cellule senza TRA escono da questa variante.

Uso: python -I validation_esterna_B/GSE200996/repliche.py <h5ad> <manifest> <marcatori.json>
"""

from __future__ import annotations

import json
import re
import sys
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[1]))  # radice del repository (con -I non è nel percorso)

import anndata as ad  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

from core.io import read_vdj_manifest  # noqa: E402
from core.stats import cluster_t_interval  # noqa: E402
from core.tcr_validation import (  # noqa: E402
    _KEY,
    MIN_PATIENTS_FOR_CI_DEFAULT,
    assign_reference_identity,
    build_clonotypes,
    clone_level_units,
    compute_marker_positivity,
    match_barcodes,
    pairwise_excess_discordance,
    parse_vdj_contigs,
    run_tcr_validation,
)

REF, JUDGED = "PBMC", "Tumor"
MINP = MIN_PATIENTS_FOR_CI_DEFAULT


def ci(res) -> dict:
    """ClusterInterval -> dizionario JSON."""
    f = lambda x: None if x is None or (isinstance(x, float) and np.isnan(x)) else float(x)  # noqa: E731
    return {"mean": f(res.mean), "ci_low": f(res.ci_low), "ci_high": f(res.ci_high), "df": f(res.df),
            "n_groups": int(res.n_groups), "n_obs": int(res.n_obs), "sufficient": bool(res.sufficient),
            "reason": res.reason}


def excess(cells: pd.DataFrame) -> dict:
    """Eccesso pesato e a pesi uguali per la coppia PBMC-Tumor con intervalli."""
    pw = pairwise_excess_discordance(cells)
    pw = pw[(pw.comp_a == REF) & (pw.comp_b == JUDGED)]
    return {"n_pairs": len(pw), "n_patients": int(pw.patient.nunique()),
            "weighted": ci(cluster_t_interval(pw["excess"].values, pw["patient"].values, MINP)),
            "equal": ci(cluster_t_interval(pw["excess_equal"].values, pw["patient"].values, MINP))}


def clone_rate(cells: pd.DataFrame, labels: list[str], tie_rule: str = "exclude") -> tuple[dict, pd.DataFrame]:
    """Tasso d'errore per clone nel tumore (convenzione della tesi) con le unità."""
    u = clone_level_units(cells, labels, REF, tie_rule=tie_rule)
    u = u[u.compartment == JUDGED]
    used = u.dropna(subset=["error"])
    res = cluster_t_interval(used["error"].values.astype(float), used["patient"].values, MINP, scale="logit")
    return ({"n_units": len(used), "n_ties": int(u["tie"].sum()), "n_errors": int(used["error"].sum()),
             "n_patients": int(used.patient.nunique()), "rate": ci(res)}, u)


def main(h5ad: Path, manifest: Path, marker_path: Path) -> None:
    t0 = time.time()
    marker_map = json.loads(marker_path.read_text())
    labels = list(marker_map)
    adata = ad.read_h5ad(h5ad)
    contigs = parse_vdj_contigs(read_vdj_manifest(manifest))

    # --- percorso della riga di comando (controllo)
    tool = run_tcr_validation(adata, contigs, patient_col="patient", compartment_col="compartment",
                              celltype_col="label", barcode_col="barcode", marker_map=marker_map,
                              reference_compartment=REF)

    # --- tabella cellula-clone-etichetta, costruita una volta
    contigs = contigs.assign(library="")
    obs = pd.DataFrame({"patient": adata.obs["patient"].astype(str).values,
                        "compartment": adata.obs["compartment"].astype(str).values,
                        "celltype": adata.obs["label"].astype(str).values,
                        "barcode": adata.obs["barcode"].astype(str).values, "library": ""})
    obs = pd.concat([obs, compute_marker_positivity(adata, marker_map).reset_index(drop=True)], axis=1)
    obs_m, clones, _ = match_barcodes(obs, build_clonotypes(contigs))
    cells = obs_m.merge(clones, on=_KEY, how="inner")
    tra = contigs[contigs.chain == "TRA"][_KEY + ["cdr3_nt"]].rename(columns={"cdr3_nt": "tra"})
    cells_ab = cells.merge(tra, on=_KEY, how="inner")
    cells_ab = cells_ab.assign(clone_id=cells_ab["clone_id"] + "|" + cells_ab["tra"]).drop(columns="tra")

    base_ex = excess(cells)
    base_rate, base_units = clone_rate(cells, labels)

    # --- controllo contro run_tcr_validation (uguaglianza esatta) e contro la CLI (arrotondata)
    d, cl = tool.discordance, tool.conventions["clone"].by_compartment[JUDGED]
    eq = {
        "eccesso_pesato": [d.mean_excess, d.ci_low, d.ci_high] == [base_ex["weighted"][k] for k in ("mean", "ci_low", "ci_high")],
        "eccesso_pesi_uguali": [d.equal_weights.mean, d.equal_weights.ci_low, d.equal_weights.ci_high]
        == [base_ex["equal"][k] for k in ("mean", "ci_low", "ci_high")],
        "coppie": d.n_pairs == base_ex["n_pairs"],
        "tasso_clone": [cl.mean, cl.ci_low, cl.ci_high, cl.n_obs, cl.n_groups]
        == [base_rate["rate"][k] for k in ("mean", "ci_low", "ci_high", "n_obs", "n_groups")],
        "pareggi": tool.conventions["clone"].n_ties.get(JUDGED, 0) == base_rate["n_ties"],
    }
    cli_out = (HERE / "results" / "BE1_cli_stdout.txt").read_text(encoding="utf-8")
    m = re.search(r"'Tumor': ([\d.]+) \(IC95% \[([\d.]+), ([\d.]+)\], (\d+) pazienti, (\d+) unità; unità in pareggio: (\d+)\)", cli_out)
    eq["cli_tasso_clone_arrotondato"] = m is not None and m.groups() == (
        f"{base_rate['rate']['mean']:.3f}", f"{base_rate['rate']['ci_low']:.3f}", f"{base_rate['rate']['ci_high']:.3f}",
        str(base_rate["rate"]["n_groups"]), str(base_rate["n_units"]), str(base_rate["n_ties"]))
    eq["cli_eccesso_arrotondato"] = (f"eccesso di discordanza di {base_ex['weighted']['mean']:+.3f} (IC95% "
                                     f"[{base_ex['weighted']['ci_low']:+.3f}, {base_ex['weighted']['ci_high']:+.3f}], su "
                                     f"{base_ex['n_pairs']:,} coppie") in cli_out

    # --- numeri dello strumento per B-E2
    strumento = {"clone_tumor": base_rate, "excess_PBMC_Tumor": base_ex,
                 "n_cells_with_trb": int(tool.n_cells_with_tcr), "n_reference_clones": int(tool.marker_error.n_resolved_clones)}

    # --- B-E3: permutazione del compartimento dentro ogni clone di ogni paziente
    be3 = []
    for seed in range(20):
        rng = np.random.default_rng(seed)
        perm = cells.copy()
        perm["compartment"] = perm.groupby(["patient", "clone_id"], observed=True)["compartment"].transform(
            lambda s: s.values[rng.permutation(len(s))])
        ex = excess(perm)
        w = ex["weighted"]
        be3.append({"seme": seed, **ex, "contiene_zero": bool(w["sufficient"] and w["ci_low"] <= 0 <= w["ci_high"])})
    means = [r["weighted"]["mean"] for r in be3]
    be3_sum = {"media": float(np.mean(means)), "min": float(np.min(means)), "max": float(np.max(means)),
               "contiene_zero": int(sum(r["contiene_zero"] for r in be3)),
               "superato": bool(-0.01 <= np.mean(means) <= 0.01 and sum(r["contiene_zero"] for r in be3) >= 18)}

    # --- B-E4: errori iniettati
    ident = assign_reference_identity(cells, labels, reference_compartment=REF).rename(columns={"clone_id": "clone"})
    cand = base_units.merge(ident, on=["patient", "clone"])
    cand = cand[(cand.reference_identity == "CD4T") & (cand.error == 0) & (~cand.tie)][["patient", "clone"]]
    cand = cand.sort_values(["patient", "clone"]).reset_index(drop=True)
    order = np.random.default_rng(0).permutation(len(cand))
    be4 = {"candidate": len(cand), "valutabile": len(cand) >= 20, "dosi": []}
    if be4["valutabile"]:
        for dose in (0.10, 0.20):
            k = int(round(dose * len(cand)))
            inj = cand.iloc[order[:k]]
            key = set(zip(inj.patient, inj.clone))
            mod = cells.copy()
            hit = (mod.compartment == JUDGED) & pd.Series(list(zip(mod.patient, mod.clone_id))).isin(key).values
            mod.loc[hit, "celltype"] = "CD8T"
            r, _ = clone_rate(mod, labels)
            ex = excess(mod)
            be4["dosi"].append({
                "dose": dose, "unita_iniettate": k, "cellule_rietichettate": int(hit.sum()),
                "errori_prima": base_rate["n_errors"], "errori_dopo": r["n_errors"],
                "aumento_errori": r["n_errors"] - base_rate["n_errors"], "pareggi_prima": base_rate["n_ties"],
                "pareggi_dopo": r["n_ties"], "unita_prima": base_rate["n_units"], "unita_dopo": r["n_units"],
                "tasso_dopo": r["rate"], "eccesso_pesi_uguali_prima": base_ex["equal"]["mean"],
                "eccesso_pesi_uguali_dopo": ex["equal"]["mean"], "eccesso_pesi_uguali_dopo_ci": ex["equal"],
                "eccesso_pesato_prima": base_ex["weighted"]["mean"], "eccesso_pesato_dopo": ex["weighted"]["mean"],
                "superato": bool(r["n_errors"] - base_rate["n_errors"] == k and ex["equal"]["mean"] > base_ex["equal"]["mean"])})
        be4["superato"] = all(d["superato"] for d in be4["dosi"])

    # --- B-E6: sensibilità
    variants = {"TRB, pareggi esclusi (default)": base_rate}
    variants["TRA+TRB, pareggi esclusi"] = clone_rate(cells_ab, labels)[0]
    variants["TRB, pareggi = errore"] = clone_rate(cells, labels, "error")[0]
    variants["TRB, pareggi = corretti"] = clone_rate(cells, labels, "correct")[0]
    b = base_rate["rate"]["mean"]
    be6 = {"varianti": variants, "differenze": {k: float(v["rate"]["mean"] - b) for k, v in variants.items()},
           "ci_sotto_5_pazienti": [k for k, v in variants.items() if v["rate"]["sufficient"] and v["rate"]["n_groups"] < 5],
           "cellule_TRA_TRB": len(cells_ab), "cellule_TRB": len(cells)}
    be6["superato"] = bool(max(abs(x) for x in be6["differenze"].values()) <= 0.05 and not be6["ci_sotto_5_pazienti"])

    # --- descrittivi sui dati originali: errori per direzione
    dirs = base_units[base_units.error == 1].merge(ident, on=["patient", "clone"])
    tumor_cells = cells[cells.compartment == JUDGED]
    maj = (tumor_cells.groupby(["patient", "clone_id"])["celltype"].agg(lambda s: s.value_counts().idxmax())
           .rename("maggioranza").reset_index().rename(columns={"clone_id": "clone"}))
    dirs = dirs.merge(maj, on=["patient", "clone"])
    descr = {"unita_per_identita": {str(k): int(v) for k, v in base_units.merge(ident, on=["patient", "clone"])
                                    .groupby("reference_identity").size().items()},
             "errori_per_direzione": {f"{a}->{m_}": int(n) for (a, m_), n in dirs.groupby(["reference_identity", "maggioranza"]).size().items()},
             "errori_per_paziente": {str(k): int(v) for k, v in dirs.groupby("patient").size().items()},
             "cloni_con_riferimento": int(len(ident)), "cellule_con_trb": int(len(cells))}

    out = {"controllo_su_dati_originali": eq, "strumento_BE2": strumento, "BE3": {"repliche": be3, "sintesi": be3_sum},
           "BE4": be4, "BE6": be6, "descrittivi": descr, "secondi": round(time.time() - t0, 1)}
    (HERE / "results" / "repliche.json").write_text(json.dumps(out, indent=2, ensure_ascii=False) + "\n")
    (HERE / "results" / "BE2_strumento.json").write_text(json.dumps(strumento, indent=2) + "\n")
    print("controllo:", eq)
    print("B-E3:", be3_sum)
    print("B-E4:", {k: v for k, v in be4.items() if k != "dosi"}, [(d["unita_iniettate"], d["aumento_errori"], d["eccesso_pesi_uguali_prima"], d["eccesso_pesi_uguali_dopo"]) for d in be4["dosi"]])
    print("B-E6:", be6["differenze"], "superato", be6["superato"])
    print("descrittivi:", descr, f"{out['secondi']} s")


if __name__ == "__main__":
    main(Path(sys.argv[1]), Path(sys.argv[2]), Path(sys.argv[3]))
