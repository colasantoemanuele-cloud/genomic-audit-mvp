"""Parti A3 e A4 -- frazione di CD8 e copertura dei flag su GSE278694 reale.

Rilegge gli ``obs`` reali (sola lettura) e il CSV dei flag scritto da ``cli.py tcr
--export-flags``, ricalcola la propagazione con gli stessi parametri della CLI (2000
repliche, seme 0: stessi numeri dell'output CLI) e scrive le tabelle.

Uso:
    python -m validation.gse278694_cd8_flags --flags <sc_raw_audit_flags.csv> --out <dir>
"""

from __future__ import annotations

import os
import argparse
import json
from pathlib import Path

import anndata as ad
import pandas as pd

from core.cd8_propagation import cd8_fraction_intervals

# Dati reali di GSE278694: non inclusi nella repo. Cartella indicata da GSE278694_DIR, altrimenti
# una cartella "pdac-ml" accanto a questa repo.
H5AD = Path(os.environ.get("GSE278694_DIR", Path(__file__).resolve().parents[2] / "pdac-ml")) / "data/interim/sc_raw.h5ad"


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--flags", type=Path, required=True)
    ap.add_argument("--out", type=Path, required=True)
    args = ap.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)

    obs = ad.read_h5ad(H5AD, backed="r").obs[["patients", "tissue", "all_celltype"]].astype(str)
    flags = pd.read_csv(args.flags, index_col=0)
    assert list(flags.index) == list(obs.index), "flag non allineati alle righe"
    df = obs.join(flags)

    # A4: copertura dei flag, totale e per compartimento
    cov = df.groupby("tissue").agg(
        cellule=("all_celltype", "size"),
        con_riferimento=("audit_reference_label", lambda s: int(s.notna().sum())),
        valutabili=("audit_label_vs_reference", lambda s: int(s.notna().sum())),
        discordanti=("audit_label_vs_reference", lambda s: int((s == True).sum())),  # noqa: E712
    )
    cov.loc["Totale"] = cov.sum()
    cov["% valutabili"] = (100 * cov.valutabili / cov.cellule).round(2)
    t = df[df.all_celltype.isin(["CD4T", "CD8T"])]
    cov_t = t.groupby("tissue").audit_label_vs_reference.apply(lambda s: round(100 * s.notna().mean(), 2))
    cov["% valutabili fra CD4T+CD8T"] = cov_t
    cov.loc["Totale", "% valutabili fra CD4T+CD8T"] = round(100 * t.audit_label_vs_reference.notna().mean(), 2)
    cov.to_csv(args.out / "a4_flag_coverage.csv")
    print(cov.to_string())

    # A3: frazione di CD8 nel tumore
    res = cd8_fraction_intervals(df, "patients", "tissue", "all_celltype", "audit_reference_label",
                                 "Tumor", n_boot=2000, seed=0)
    rows = []
    for p in res.patients:
        row = {"paziente": p.patient, "n_CD4T+CD8T": p.n_cd4_called + p.n_cd8_called,
               "riportata": p.reported, "IC_conteggi": f"{p.naive_low:.2f}-{p.naive_high:.2f}"}
        for k, iv in sorted(p.scenarios.items()):
            row[f"{k:g}x"] = (f"{iv.low:.2f}-{iv.high:.2f}" if iv.low is not None
                              else f"non prodotto: {iv.refused_reason}")
        rows.append(row)
    tab = pd.DataFrame(rows)
    tab.to_csv(args.out / "a3_cd8_fraction.csv", index=False)
    summary = {
        "J": res.youden_j, "cellule_di_riferimento": res.n_reference_cells,
        "pazienti_di_riferimento": res.n_reference_patients,
        "matrice": res.matrix.round(4).to_dict() if res.matrix is not None else None,
        "intervalli_prodotti_per_scenario": {
            f"{k:g}x": sum(p.scenarios[k].low is not None for p in res.patients) for k in (0.5, 1.0, 2.0)},
        "pazienti": len(res.patients),
        "rifiuti": sorted({p.scenarios[k].refused_reason for p in res.patients for k in (0.5, 1.0, 2.0)
                           if p.scenarios[k].low is None}),
    }
    (args.out / "a3_summary.json").write_text(json.dumps(summary, indent=1, ensure_ascii=False))
    print(res.matrix.round(3).to_string() if res.matrix is not None else "matrice non stimata")
    print(json.dumps(summary, indent=1, ensure_ascii=False))
    print(tab.to_string(index=False))


if __name__ == "__main__":
    main()
