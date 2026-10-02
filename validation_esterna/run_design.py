"""Esecuzione dell'audit del disegno sui dataset esterni e confronto con la verifica
indipendente (criteri D1-D4 di CRITERI.md).

D1: la CLI viene lanciata come processo separato con i soli flag documentati.
D2-D4: i fatti strutturati dell'audit (core.design_audit, letto come fa la CLI: CSV come
stringhe) vengono confrontati con quelli di validation_esterna/indipendenti/design_check.py,
che non importa nulla da core/.

Uso: python -m validation_esterna.run_design
"""

from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path

import pandas as pd

from core.design_audit import run_design_audit
from validation_esterna.indipendenti import design_check as ind

ROOT = Path(__file__).resolve().parents[1]
RES = Path(__file__).parent / "results"

CONFIGS = {
    "GSE132465": dict(
        patient="patient_id", tissue="tissue type",
        technical={"library": "gsm", "batch": "platform"}, outcomes=["tumor stage", "region"],
        comparisons=[("tissue type", "Colorectal cancer", "Normal mucosa"), ("tumor stage", "2", "3")]),
    "GSE131907": dict(
        patient="patient id", tissue="tissue origin abbrevation",
        technical={"library": "gsm", "batch": "platform"}, outcomes=["tumor stage"],
        comparisons=[("tissue origin abbrevation", "tLung", "nLung"),
                     ("tissue origin abbrevation", "mLN", "nLN"),
                     ("tissue origin abbrevation", "tLung", "mBrain"), ("tumor stage", "I", "IV")]),
    "GSE125449": dict(
        patient="patient_from_title", tissue=None,
        technical={"library": "gsm", "batch": "platform"}, outcomes=["cancer type"],
        comparisons=[("cancer type", "Hepatocellular carcinoma", "Intrahepatic cholangiocarcinoma"),
                     ("platform", "GPL18573", "GPL20301")]),
}


def cli_args(gse: str, cfg: dict) -> list[str]:
    a = [sys.executable, "cli.py", "design", "--meta", str(RES / f"{gse}_samples.csv"),
         "--patient-col", cfg["patient"]]
    if cfg["tissue"]:
        a += ["--tissue-col", cfg["tissue"]]
    for role, col in cfg["technical"].items():
        a += ["--technical", f"{role}={col}"]
    for o in cfg["outcomes"]:
        a += ["--outcome-col", o]
    for f, x, y in cfg["comparisons"]:
        a += ["--compare", f"{f}:{x}:{y}"]
    return a + ["--out", str(RES / f"design_{gse}.html")]


def tool_facts(res) -> set[tuple]:
    out = set()
    for f in res.findings:
        if f.kind == "uno-a-uno":
            out.add(("uno-a-uno",) + tuple(f.factors))
        elif f.kind in ("annidamento", "esito-determinato", "unita'-tecnica"):
            out.add((f.kind,) + tuple(f.factors))
    return out


def main() -> None:
    summary = {}
    for gse, cfg in CONFIGS.items():
        # D1
        proc = subprocess.run(cli_args(gse, cfg), cwd=ROOT, capture_output=True, text=True)
        (RES / f"design_{gse}_cli.txt").write_text(proc.stdout + ("\n[stderr]\n" + proc.stderr if proc.returncode else ""))
        d1 = proc.returncode == 0 and "Traceback" not in proc.stderr

        roles = {cfg["patient"]: "patient"}
        if cfg["tissue"]:
            roles[cfg["tissue"]] = "tissue"
        for role, col in cfg["technical"].items():
            roles.setdefault(col, role)
        for o in cfg["outcomes"]:
            roles.setdefault(o, "outcome")
        sheet = pd.read_csv(RES / f"{gse}_samples.csv", dtype=str)
        res = run_design_audit(sheet, patient_col=cfg["patient"], tissue_col=cfg["tissue"],
                               technical_cols=cfg["technical"], outcome_cols=cfg["outcomes"],
                               comparisons=cfg["comparisons"])
        # verifica indipendente
        df = ind.read_sheet(RES / f"{gse}_samples.csv", list(roles))
        outcome_like = ([cfg["tissue"]] if cfg["tissue"] else []) + cfg["outcomes"]
        f_ind = ind.structural_facts(df, roles, outcome_like, cfg["patient"], cfg["tissue"])
        f_tool = tool_facts(res)
        # D3
        comp_rows = []
        for c, (f, x, y) in zip(res.comparisons, cfg["comparisons"]):
            cls_i, units_i = ind.classify(df, roles, cfg["patient"], f, x, y)
            comp_rows.append({"confronto": f"{f}: {x} vs {y}", "classe audit": c.classification,
                              "classe indipendente": cls_i, "unita' audit": c.n_units,
                              "unita' indipendenti": units_i,
                              "coincide": c.classification == cls_i and c.n_units == units_i})
        # D4
        v_ind = ind.cramer_v_expected(df, roles)
        v_rows, d4 = [], True
        for p in res.pairs:
            key = frozenset((p.factor_a, p.factor_b))
            vi = v_ind.get(key, "assente")
            ok = (vi is None and p.cramer_v is None) or (
                isinstance(vi, float) and p.cramer_v is not None and abs(vi - p.cramer_v) < 1e-9)
            d4 &= ok
            v_rows.append({"coppia": f"{p.factor_a} x {p.factor_b}", "V audit": p.cramer_v,
                           "V indipendente": vi, "coincide": ok})
        d4 &= len(v_ind) == len(res.pairs)
        summary[gse] = {
            "D1": d1, "exit_code": proc.returncode,
            "D2": f_ind == f_tool,
            "fatti_audit": sorted(map(list, f_tool)), "fatti_indipendenti": sorted(map(list, f_ind)),
            "solo_audit": sorted(map(list, f_tool - f_ind)), "solo_indipendenti": sorted(map(list, f_ind - f_tool)),
            "D3": all(r["coincide"] for r in comp_rows), "confronti": comp_rows,
            "D4": d4, "cramer_v": v_rows,
        }
        print(f"== {gse}: D1={d1} D2={f_ind == f_tool} D3={summary[gse]['D3']} D4={d4} "
              f"| fatti audit {len(f_tool)}, indipendenti {len(f_ind)}")
    (RES / "design_summary.json").write_text(json.dumps(summary, indent=1, ensure_ascii=False, default=str))


if __name__ == "__main__":
    main()
