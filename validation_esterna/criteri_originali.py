"""Valuta il codice ATTUALE con le tre versioni dei criteri (REPORT.md, «Criteri originali e
modifiche»). Non modifica lo strumento: cambia solo la regola di valutazione.

Versioni dei criteri di CRITERI.md:
- O  = sezioni 1-4 (commit 0b5a847), scritte prima di toccare i dati esterni;
- M1 = O + sezione 5 (commit b08b2d2): Modulo A valutato sulla coorte ridotta;
- M2 = O + sezione 6 (commit 0662a40): Modulo A sul dataset completo, con la pipeline
  dichiarata del nuovo motore in A2.

File di ingresso, tutti in results/ (comandi dalla radice della repo):
  python -m validation_esterna.indipendenti.leakage_check --pipeline originale \
      --h5ad validation_esterna/data/GSE125449.h5ad --target-col Type --patient-col patient \
      --out validation_esterna/results/A2_originale_completo.json
  python -m validation_esterna.indipendenti.leakage_check --pipeline originale \
      --h5ad data/demo/GSE125449_demo.h5ad --target-col Type --patient-col patient \
      --out validation_esterna/results/A2_originale_ridotta.json
  python -m validation_esterna.indipendenti.leakage_check --pipeline sezione6 \
      --h5ad data/demo/GSE125449_demo.h5ad --target-col Type --patient-col patient \
      --out validation_esterna/results/A2_sezione6_ridotta.json
  python -m validation_esterna.run_moduleA base    data/demo/GSE125449_demo.h5ad _ridotta
  python -m validation_esterna.run_moduleA perm    data/demo/GSE125449_demo.h5ad _ridotta
  python -m validation_esterna.run_moduleA formati data/demo/GSE125449_demo.h5ad _ridotta
  python cli.py leakage --h5ad data/demo/GSE125449_demo.h5ad --target-col Type \
      --patient-col patient --out validation_esterna/results/leakage_ridotta.html \
      > validation_esterna/results/A1_ridotta_cli.txt; echo $? > validation_esterna/results/A1_ridotta_exit.txt

Uso: python -m validation_esterna.criteri_originali
"""

from __future__ import annotations

import json
from pathlib import Path

R = Path(__file__).parent / "results"
TOL = 0.02  # tolleranza di A2, uguale in tutte le versioni dei criteri


def _j(name: str) -> dict:
    return json.loads((R / name).read_text())


def _a2(tool: dict, ind: dict, metric: str) -> dict:
    """Differenze assolute fra le medie dello strumento e quelle dello script indipendente."""
    rows = {k: {"strumento": tool[k], "indipendente": ind[k][f"media_{metric}"],
                "diff": abs(tool[k] - ind[k][f"media_{metric}"]),
                # informativo: la stessa differenza con le altre due definizioni di macro-F1
                "diff_altre_metriche": {m: abs(tool[k] - ind[k][f"media_{m}"])
                                        for m in ("classi_del_test", "standard", "tutte_le_classi") if m != metric}}
            for k in tool}
    return {"metrica_indipendente": metric, "confronti": rows,
            "esito": all(r["diff"] <= TOL for r in rows.values())}


def _a4c(tool_logo: float, ind: dict) -> dict:
    """La media dello strumento deve coincidere con la macro-F1 sulle classi presenti e non con
    quella che conta come zero le classi assenti."""
    lg = ind["logo_logreg"]
    n_abs = sum(bool(x) for x in lg["classi_assenti_e_non_predette"])
    return {"fold_con_classi_assenti": n_abs, "strumento": tool_logo,
            "classi_presenti": lg["media_classi_del_test"], "assenti_come_zero": lg["media_tutte_le_classi"],
            "esito": abs(tool_logo - lg["media_classi_del_test"]) <= TOL
            and (n_abs == 0 or abs(tool_logo - lg["media_tutte_le_classi"]) > TOL)}


def main() -> None:
    bm, fm, pm = _j("benchmark_moduleA_reale.json"), _j("moduleA_formati.json"), _j("moduleA_perm.json")
    ind6, ind_o = _j("A2_indipendente.json"), _j("A2_originale_completo.json")
    full_tool = {"raggruppato": sum(fm["fold_raggruppati"]["riferimento"]) / len(fm["fold_raggruppati"]["riferimento"]),
                 "logo_logreg": bm["logo"]["logreg"]}

    bs, fs, ps = _j("moduleA_base_ridotta.json"), _j("moduleA_formati_ridotta.json"), _j("moduleA_perm_ridotta.json")
    s_o, s_6 = _j("A2_originale_ridotta.json"), _j("A2_sezione6_ridotta.json")
    small_tool = {"raggruppato": bs["grouped_mean"], "logo_logreg": bs["logo_mean"]["logreg"]}
    a1_small = int((R / "A1_ridotta_exit.txt").read_text().strip()) == 0

    def a3(p):
        return {"media": p["media_raggruppata"], "atteso": p["atteso_1_su_K"],
                "tag_leakage": p["permutazioni_con_tag_leakage"], "esito": bool(p["A3_i"] and p["A3_ii"])}

    def a4ab(f):
        return bool(f["A4a_densa_identica"] and f["A4b_stringhe_identica"] and f["A4b_categoriche_identica"])

    out = {
        "O": {"criteri": "sezioni 1-4 (0b5a847)", "dataset": "GSE125449 completo (19 pazienti)",
              "A1": True, "A2": _a2(full_tool, ind_o, "standard"), "A3": a3(pm), "A4ab": a4ab(fm),
              "A4c": _a4c(bm["logo"]["logreg"], ind6)},
        "M1": {"criteri": "sezioni 1-5 (b08b2d2)", "dataset": "GSE125449 ridotto (10 pazienti, 1.861 cellule)",
               "A1": a1_small, "A2": _a2(small_tool, s_o, "standard"), "A3": a3(ps), "A4ab": a4ab(fs),
               "A4c": _a4c(bs["logo_mean"]["logreg"], s_6)},
        "M2": {"criteri": "sezioni 1-4 e 6 (0662a40)", "dataset": "GSE125449 completo (19 pazienti)",
               "A1": True,
               "A2": _a2({**full_tool, "casuale": bm["casuale"]}, ind6, "classi_del_test"),
               "A3": a3(pm), "A4ab": a4ab(fm), "A4c": _a4c(bm["logo"]["logreg"], ind6)},
    }
    # informativo: la stessa coorte ridotta valutata con la pipeline della sezione 6
    out["M1"]["A2_pipeline_sezione6"] = _a2({**small_tool, "casuale": bs["random_mean"]}, s_6, "classi_del_test")
    for v in out.values():
        v["tutti_i_criteri_A"] = bool(v["A1"] and v["A2"]["esito"] and v["A3"]["esito"] and v["A4ab"]
                                      and v["A4c"]["esito"])
    (R / "criteri_originali.json").write_text(json.dumps(out, indent=1, ensure_ascii=False))
    for k, v in out.items():
        a2 = "; ".join(f"{n}: {r['strumento']:.4f} contro {r['indipendente']:.4f} (diff {r['diff']:.4f})"
                       for n, r in v["A2"]["confronti"].items())
        print(f"{k} [{v['criteri']}] {v['dataset']}: A1={v['A1']} A2={v['A2']['esito']} ({a2}) "
              f"A3={v['A3']['esito']} A4ab={v['A4ab']} A4c={v['A4c']['esito']} -> tutti={v['tutti_i_criteri_A']}")


if __name__ == "__main__":
    main()
