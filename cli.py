#!/usr/bin/env python3
"""Esecuzione da riga di comando del motore analitico (core/), senza Streamlit.

Esempi:
    python cli.py demo --out results/demo_report.html
    python cli.py design --meta campioni.csv --patient-col patient --tissue-col tissue \
        --technical batch=run --technical protocol=protocol \
        --compare tissue:Tumor:Adjacent_normal --out results/design_report.html
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

from core.cd8_propagation import cd8_fraction_intervals, format_cd8_text
from core.design_audit import TECHNICAL_ROLES, run_design_audit, sample_sheet_from_obs
from core.leakage_audit import run_leakage_audit
from core.report import save_report
from core.synthetic import (
    make_gse278694_like_sheet,
    make_leakage_dataset,
    make_tcr_validation_dataset,
    write_vdj_csvs,
)
from core.tcr_validation import export_audited, parse_vdj_contigs, run_tcr_validation


def _parse_technical(items: list[str]) -> dict[str, str]:
    out = {}
    for item in items:
        if "=" not in item:
            sys.exit(f"--technical vuole ruolo=colonna, ricevuto '{item}'")
        role, col = item.split("=", 1)
        if role not in TECHNICAL_ROLES:
            sys.exit(f"ruolo tecnico '{role}' non riconosciuto; ammessi: {', '.join(TECHNICAL_ROLES)}")
        out[role] = col
    return out


def _parse_comparisons(items: list[str]) -> list[tuple[str, str, str]]:
    out = []
    for item in items:
        parts = item.split(":")
        if len(parts) != 3:
            sys.exit(f"--compare vuole colonna:livello_a:livello_b, ricevuto '{item}'")
        out.append((parts[0], parts[1], parts[2]))
    return out


def _print_design(result) -> None:
    print(result.narrative)
    for pr in result.pairs:
        print(" -", pr.sentence)
    for note in result.notes:
        print(" *", note)


def _cmd_design(args: argparse.Namespace) -> None:
    technical = _parse_technical(args.technical)
    comparisons = _parse_comparisons(args.compare)
    cols = [args.patient_col] + ([args.tissue_col] if args.tissue_col else []) \
        + list(technical.values()) + list(args.outcome_col)
    if args.meta:
        sheet = pd.read_csv(args.meta, dtype=str)
        name = Path(args.meta).stem
    else:
        adata = ad.read_h5ad(args.h5ad, backed="r")
        sheet = sample_sheet_from_obs(adata.obs, cols)
        name = Path(args.h5ad).stem
    result = run_design_audit(
        sheet, patient_col=args.patient_col, tissue_col=args.tissue_col,
        technical_cols=technical, outcome_cols=list(args.outcome_col), comparisons=comparisons,
        min_units=args.min_units,
    )
    _print_design(result)
    out = save_report(args.out, design_result=result, dataset_name=name)
    print(f"[ok] report scritto in {out}")


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
    if result.flag_coverage is not None:
        print(f"Flag per cellula valutabili (audit_label_vs_reference non NA): "
              f"{result.flag_coverage:.1%} delle cellule.")
    if args.export_flags:
        if result.cell_flags is None:
            sys.exit("--export-flags richiede --marker-map e --reference-compartment")
        prefix = Path(args.out).parent / Path(args.h5ad).stem
        h5ad_path, csv_path = export_audited(adata, result, prefix)
        print(f"[ok] copia con i flag scritta in {h5ad_path} (il file originale non e' toccato)")
        print(f"[ok] flag per cellula scritti in {csv_path}")
    cd8 = None
    if args.cd8_compartment:
        if result.cell_flags is None:
            sys.exit("--cd8-compartment richiede --marker-map e --reference-compartment")
        obs = adata.obs.join(result.cell_flags)
        cd8 = cd8_fraction_intervals(
            obs, patient_col=args.patient_col, compartment_col=args.compartment_col,
            celltype_col=args.celltype_col, reference_col="audit_reference_label",
            target_compartment=args.cd8_compartment, cd4_label=args.cd4_label,
            cd8_label=args.cd8_label, n_boot=args.n_boot, seed=args.seed)
        print(format_cd8_text(cd8))
    out = save_report(args.out, tcr_result=result, dataset_name=Path(args.h5ad).stem,
                      cd8_result=cd8)
    print(f"[ok] report scritto in {out}")


def _cmd_demo(args: argparse.Namespace) -> None:
    out_dir = Path(args.data_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    design_sheet = make_gse278694_like_sheet()
    design_result = run_design_audit(
        design_sheet, patient_col="patient", tissue_col="tissue",
        technical_cols={"protocol": "protocol", "library": "library"},
        comparisons=[("tissue", "Tumor", "Adjacent_normal"), ("protocol", "scRNA", "snRNA")],
    )
    print("=== Audit del disegno (struttura sintetica modellata su GSE278694) ===")
    _print_design(design_result)
    print()

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
    cd8_result = cd8_fraction_intervals(
        tcr_adata.obs.join(tcr_result.cell_flags), patient_col="patient_id",
        compartment_col="tissue", celltype_col="celltype", reference_col="audit_reference_label",
        target_compartment="Tumor", seed=args.seed)
    print("\n=== Frazione di CD8 nel tumore (dataset sintetico) ===")
    print(format_cd8_text(cd8_result))

    out = save_report(args.out, leakage_result=leakage_result, tcr_result=tcr_result,
                       dataset_name="demo sintetico", design_result=design_result,
                       cd8_result=cd8_result)
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

    p_des = sub.add_parser("design", help="audit del disegno e del confondimento (solo metadati)")
    src = p_des.add_mutually_exclusive_group(required=True)
    src.add_argument("--meta", type=Path, help="CSV dei metadati, una riga per campione/libreria")
    src.add_argument("--h5ad", type=Path, help="AnnData: si usa solo adata.obs, ridotto ai campioni")
    p_des.add_argument("--patient-col", required=True)
    p_des.add_argument("--tissue-col", default=None)
    p_des.add_argument("--technical", action="append", default=[],
                       help=f"ruolo=colonna, ruoli: {', '.join(TECHNICAL_ROLES)} (ripetibile)")
    p_des.add_argument("--outcome-col", action="append", default=[], help="colonna di esito (ripetibile)")
    p_des.add_argument("--compare", action="append", default=[],
                       help="colonna:livello_a:livello_b (ripetibile)")
    p_des.add_argument("--min-units", type=int, default=5)
    p_des.add_argument("--out", default="results/design_report.html", type=Path)
    p_des.set_defaults(func=_cmd_design)

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
    p_tcr.add_argument("--export-flags", action="store_true",
                        help="scrive <nome>_audited.h5ad (copia con i flag) e <nome>_audit_flags.csv "
                             "nella cartella di --out; richiede --marker-map")
    p_tcr.add_argument("--cd8-compartment", default=None,
                        help="compartimento (es. Tumor) su cui stimare l'intervallo della frazione di "
                             "CD8 per paziente; richiede --marker-map")
    p_tcr.add_argument("--cd4-label", default="CD4T")
    p_tcr.add_argument("--cd8-label", default="CD8T")
    p_tcr.add_argument("--n-boot", type=int, default=2000)
    p_tcr.add_argument("--seed", type=int, default=0)
    p_tcr.add_argument("--out", default="results/tcr_report.html", type=Path)
    p_tcr.set_defaults(func=_cmd_tcr)

    args = ap.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
