#!/usr/bin/env python3
"""Esecuzione da riga di comando del motore analitico (core/), senza Streamlit.

Esempi:
    python cli.py demo --out results/demo_report.html
    python cli.py leakage --h5ad dati.h5ad --target-col tissue --patient-col patient_id \
        --out results/leakage_report.html
    python cli.py tcr --h5ad dati.h5ad --vdj-manifest vdj_manifest.csv \
        --patient-col patient_id --compartment-col tissue --celltype-col celltype \
        --barcode-col barcode --out results/tcr_report.html

``vdj_manifest.csv`` ha tre colonne: path,patient,compartment -- una riga per ogni CSV
Cell Ranger da includere (i CSV Cell Ranger non contengono queste informazioni al loro
interno, e non c'e' una convenzione di nome file universale per dedurle).
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import anndata as ad
import pandas as pd

from core.leakage_audit import run_leakage_audit
from core.report import save_report
from core.synthetic import make_leakage_dataset, make_tcr_validation_dataset, write_vdj_csvs
from core.tcr_validation import parse_vdj_contigs, run_tcr_validation


def _cmd_leakage(args: argparse.Namespace) -> None:
    adata = ad.read_h5ad(args.h5ad)
    result = run_leakage_audit(
        adata, target_col=args.target_col, patient_col=args.patient_col,
        n_folds=args.n_folds, min_patients_for_model_comparison=args.min_patients_for_model_comparison,
        seed=args.seed,
    )
    print(result.narrative)
    out = save_report(args.out, leakage_result=result, dataset_name=Path(args.h5ad).stem)
    print(f"[ok] report scritto in {out}")


def _cmd_tcr(args: argparse.Namespace) -> None:
    adata = ad.read_h5ad(args.h5ad)
    manifest = pd.read_csv(args.vdj_manifest)
    for col in ("path", "patient", "compartment"):
        if col not in manifest.columns:
            sys.exit(f"vdj_manifest.csv deve avere le colonne path,patient,compartment (manca '{col}')")
    files = [(Path(r.path), str(r.patient), str(r.compartment)) for r in manifest.itertuples()]
    contigs = parse_vdj_contigs(files)

    marker_map = json.loads(Path(args.marker_map).read_text()) if args.marker_map else None

    result = run_tcr_validation(
        adata, contigs, patient_col=args.patient_col, compartment_col=args.compartment_col,
        celltype_col=args.celltype_col, barcode_col=args.barcode_col,
        n_boot=args.n_boot, seed=args.seed, marker_map=marker_map,
        reference_compartment=args.reference_compartment,
    )
    print(result.narrative)
    out = save_report(args.out, tcr_result=result, dataset_name=Path(args.h5ad).stem)
    print(f"[ok] report scritto in {out}")


def _cmd_demo(args: argparse.Namespace) -> None:
    out_dir = Path(args.data_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    leakage_adata = make_leakage_dataset(seed=args.seed)
    leakage_result = run_leakage_audit(
        leakage_adata, target_col="label", patient_col="patient_id", seed=args.seed,
    )
    print("=== Modulo A (dataset sintetico) ===")
    print(leakage_result.narrative)

    tcr_adata, contigs = make_tcr_validation_dataset(seed=args.seed)
    write_vdj_csvs(contigs, out_dir / "vdj")
    marker_map = {"CD4T": ["CD4"], "CD8T": ["CD8A", "CD8B"]}
    tcr_result = run_tcr_validation(
        tcr_adata, contigs, patient_col="patient_id", compartment_col="tissue",
        celltype_col="celltype", barcode_col="barcode", seed=args.seed,
        marker_map=marker_map, reference_compartment="PBMC",
    )
    print("\n=== Modulo B (dataset sintetico) ===")
    print(tcr_result.narrative)

    out = save_report(args.out, leakage_result=leakage_result, tcr_result=tcr_result,
                       dataset_name="demo sintetico")
    print(f"\n[ok] report scritto in {out}")
    print(f"[ok] file VDJ di esempio scritti in {out_dir / 'vdj'}")


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="command", required=True)

    p_demo = sub.add_parser("demo", help="genera dati sintetici ed esegue entrambi i moduli")
    p_demo.add_argument("--out", default="results/demo_report.html", type=Path)
    p_demo.add_argument("--data-dir", default="data/synthetic", type=Path)
    p_demo.add_argument("--seed", type=int, default=0)
    p_demo.set_defaults(func=_cmd_demo)

    p_leak = sub.add_parser("leakage", help="Modulo A: audit del leakage per paziente")
    p_leak.add_argument("--h5ad", required=True, type=Path)
    p_leak.add_argument("--target-col", required=True)
    p_leak.add_argument("--patient-col", required=True)
    p_leak.add_argument("--n-folds", type=int, default=5)
    p_leak.add_argument("--min-patients-for-model-comparison", type=int, default=8)
    p_leak.add_argument("--seed", type=int, default=0)
    p_leak.add_argument("--out", default="results/leakage_report.html", type=Path)
    p_leak.set_defaults(func=_cmd_leakage)

    p_tcr = sub.add_parser("tcr", help="Modulo B: validazione dell'annotazione via TCR")
    p_tcr.add_argument("--h5ad", required=True, type=Path)
    p_tcr.add_argument("--vdj-manifest", required=True, type=Path,
                        help="CSV con colonne path,patient,compartment")
    p_tcr.add_argument("--patient-col", required=True)
    p_tcr.add_argument("--compartment-col", required=True)
    p_tcr.add_argument("--celltype-col", required=True)
    p_tcr.add_argument("--barcode-col", required=True)
    p_tcr.add_argument("--marker-map", type=Path, default=None,
                        help="JSON {\"etichetta\": [\"gene1\", \"gene2\"]}, opzionale")
    p_tcr.add_argument("--reference-compartment", default=None,
                        help="richiesto insieme a --marker-map per il tasso d'errore per compartimento")
    p_tcr.add_argument("--n-boot", type=int, default=2000)
    p_tcr.add_argument("--seed", type=int, default=0)
    p_tcr.add_argument("--out", default="results/tcr_report.html", type=Path)
    p_tcr.set_defaults(func=_cmd_tcr)

    args = ap.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
