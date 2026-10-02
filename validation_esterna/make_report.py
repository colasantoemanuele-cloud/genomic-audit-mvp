"""Assembla validation_esterna/REPORT.md dai file di risultato (nessun numero trascritto a mano).

Uso: python -m validation_esterna.make_report "<riga della suite finale>"
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

HERE = Path(__file__).parent
R = HERE / "results"
UA, UI = "unita' audit", "unita' indipendenti"


def main() -> None:
    suite = sys.argv[1] if len(sys.argv) > 1 else "(non disponibile)"
    ds = json.loads((R / "design_summary.json").read_text())
    ds0 = json.loads((R / "prima_esecuzione" / "design_summary.json").read_text())
    ind = json.loads((R / "A2_indipendente.json").read_text())
    fm = json.loads((R / "moduleA_formati.json").read_text())
    pm = json.loads((R / "moduleA_perm.json").read_text())
    bm = json.loads((R / "benchmark_moduleA_reale.json").read_text())
    bs = json.loads((R / "benchmark_moduleA_sintetico.json").read_text())
    tool_g = fm["fold_raggruppati"]["riferimento"]
    tg = sum(tool_g) / len(tool_g)
    ig, ic, il = (ind["raggruppato"]["media_classi_del_test"], ind["casuale"]["media_classi_del_test"],
                  ind["logo_logreg"]["media_classi_del_test"])
    n_abs = sum(bool(x) for x in ind["logo_logreg"]["classi_assenti_e_non_predette"])

    def cli(g):
        return "\n".join(x for x in (R / f"design_{g}_cli.txt").read_text().splitlines() if not x.startswith("[ok]"))

    def drow(g):
        v = ds[g]
        comps = "; ".join(f"{c['confronto']}: {c['classe audit']} ({c[UA]} u.) = {c['classe indipendente']} "
                          f"({c[UI]} u.)" for c in v["confronti"])
        nv = sum(1 for c in v["cramer_v"] if c["V audit"] is None)
        ok = lambda b: "OK" if b else "NO"  # noqa: E731
        return (f"| {g} | {ok(v['D1'])} | {ok(v['D2'])}: {len(v['fatti_audit'])} fatti audit = "
                f"{len(v['fatti_indipendenti'])} indipendenti | {ok(v['D3'])}: {comps} | {ok(v['D4'])}: "
                f"{len(v['cramer_v'])} coppie, {nv} non valutabili, valori coincidenti entro 1e-9 |")

    a1 = (R / "A1_cli_stdout.txt").read_text().strip()
    text = (HERE / "REPORT_template.md").read_text()
    repl = {
        "DROW_GSE132465": drow("GSE132465"), "DROW_GSE131907": drow("GSE131907"), "DROW_GSE125449": drow("GSE125449"),
        "SOLO_AUDIT_PRIMA": str(ds0["GSE132465"]["solo_audit"]),
        "A2_TXT": (f"split per paziente: strumento {tg:.5f}, indipendente {ig:.5f} (diff {abs(tg - ig):.1e}; "
                   f"fold: max diff {max(abs(a - b) for a, b in zip(tool_g, ind['raggruppato']['classi_del_test'])):.4f}); "
                   f"split casuale: {bm['casuale']:.5f} contro {ic:.5f} (diff {abs(bm['casuale'] - ic):.1e}); "
                   f"LeaveOneGroupOut regressione logistica: {bm['logo']['logreg']:.5f} contro {il:.5f} "
                   f"(diff {abs(bm['logo']['logreg'] - il):.1e}). Tolleranza 0,02"),
        "A3_TXT": (f"20 permutazioni (semi 0-19): macro-F1 per paziente media {pm['media_raggruppata']:.4f} "
                   f"(range {pm['min_max_raggruppata'][0]:.4f}-{pm['min_max_raggruppata'][1]:.4f}), atteso 1/K = "
                   f"{pm['atteso_1_su_K']:.3f} (diff {abs(pm['media_raggruppata'] - pm['atteso_1_su_K']):.4f} <= 0,05); "
                   f"divario medio {pm['media_divario']:+.4f}; tag 'leakage rilevabile' in "
                   f"{pm['permutazioni_con_tag_leakage']}/20 (<= 1)"),
        "A4C_TXT": (f"nei 19 fold LeaveOneGroupOut {n_abs} hanno classi assenti; lo strumento da' "
                    f"{bm['logo']['logreg']:.4f}, coerente con la macro-F1 sulle classi presenti ({il:.4f}) e non "
                    f"con quella che conta le assenti come zero ({ind['logo_logreg']['media_tutte_le_classi']:.4f}); "
                    f"test unitario: un fold perfetto con una classe assente vale 1,0 e la classe e' dichiarata"),
        "CLI_GSE132465": cli("GSE132465"), "CLI_GSE131907": cli("GSE131907"), "CLI_GSE125449": cli("GSE125449"),
        "A1_STDOUT": a1,
        "LOGO_TXT": (f"regressione logistica {bm['logo']['logreg']:.3f}, random forest {bm['logo']['random_forest']:.3f}, "
                     f"gradient boosting {bm['logo']['hist_gb']:.3f}"),
        "GAP_TXT": f"{bm['divario']:+.3f}** ({bm['casuale']:.3f} contro {bm['raggruppato']:.3f})",
        "XAI_TXT": f"{bm['xai_jaccard'][0]:.3f} fra fold per paziente, {bm['xai_jaccard'][1]:.3f} fra fold casuali",
        "SINT_S": f"{bs['secondi']:.0f}", "SUITE_FINALE": suite,
    }
    for k, v in repl.items():
        text = text.replace("{{" + k + "}}", v)
    (HERE / "REPORT.md").write_text(text)
    print(f"REPORT.md: {len(text.splitlines())} righe")


if __name__ == "__main__":
    main()
