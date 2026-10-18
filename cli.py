#!/usr/bin/env python3
"""Esecuzione da riga di comando del motore analitico (core/). `serve` avvia la web app locale.

Strumento per la ricerca. Non è un dispositivo medico e non va usato per decisioni cliniche
su singoli pazienti.

Esempi:
    audit-sc demo --out results/demo_report.html
    audit-sc design --meta campioni.csv --patient-col patient --tissue-col tissue \
        --technical batch=run --technical protocol=protocol \
        --compare tissue:Tumor:Adjacent_normal --out results/design_report.html
    audit-sc leakage --h5ad dati.h5ad --target-col tissue --patient-col patient_id \
        --out results/leakage_report.html
    audit-sc tcr --h5ad dati.h5ad --vdj-manifest vdj_manifest.csv \
        --patient-col patient_id --compartment-col tissue --celltype-col celltype \
        --barcode-col barcode --out results/tcr_report.html

``vdj_manifest.csv`` ha le colonne path,patient,compartment (e, se ci sono più librerie per la
stessa coppia paziente-compartimento, library): una riga per ogni CSV Cell Ranger da includere.
I percorsi relativi sono letti rispetto alla cartella del manifest.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import anndata as ad

from core.cd8_propagation import cd8_fraction_intervals, format_cd8_text
from core.design_audit import (
    LIBRARY_COLLAPSE_NOTE,
    TECHNICAL_ROLES,
    run_design_audit,
    sample_sheet_from_obs,
)
from core.io import read_table, read_vdj_manifest
from core.leakage_audit import run_leakage_audit
from core.report import save_report
from core.synthetic import (
    DEMO_LEAKAGE_NOTE,
    make_gse278694_like_sheet,
    make_leakage_dataset,
    make_tcr_validation_dataset,
    write_vdj_csvs,
)
from core.tcr_validation import export_audited, parse_vdj_contigs, run_tcr_validation, validate_marker_map
from core.verdict import USE_WARNING


class InputError(ValueError):
    """Errore nei parametri o nei file forniti dall'utente: messaggio senza traceback."""


def _parse_technical(items: list[str]) -> dict[str, str]:
    out = {}
    for item in items:
        if "=" not in item:
            raise InputError(f"--technical vuole ruolo=colonna, ricevuto '{item}'")
        role, col = item.split("=", 1)
        if role not in TECHNICAL_ROLES:
            raise InputError(f"ruolo tecnico '{role}' non riconosciuto; ammessi: {', '.join(TECHNICAL_ROLES)}")
        if role in out:
            raise InputError(f"il ruolo '{role}' è indicato due volte (--technical {role}={out[role]} e "
                             f"--technical {item}): un ruolo, una colonna")
        out[role] = col
    return out


def _parse_comparisons(items: list[str], declared: list[str], patient_col: str) -> list[tuple[str, str, str]]:
    out = []
    for item in items:
        parts = item.split(":")
        if len(parts) != 3 or not all(parts):
            raise InputError(f"--compare vuole colonna:livello_a:livello_b, ricevuto '{item}'")
        if parts[0] == patient_col:
            raise InputError(f"--compare '{item}': il confronto non può essere sulla colonna paziente")
        if parts[0] not in declared:
            raise InputError(f"--compare '{item}': la colonna '{parts[0]}' non è dichiarata. Indicala con "
                             f"--tissue-col, --technical ruolo=colonna oppure --outcome-col. Colonne "
                             f"dichiarate: {', '.join(declared)}")
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
    cols = list(dict.fromkeys([args.patient_col] + ([args.tissue_col] if args.tissue_col else [])
                              + list(technical.values()) + list(args.outcome_col)))
    comparisons = _parse_comparisons(args.compare, cols, args.patient_col)
    notes = []
    if args.meta:
        sheet = read_table(args.meta)
        name = Path(args.meta).stem
    else:
        adata = ad.read_h5ad(args.h5ad, backed="r")
        missing = [c for c in cols if c not in adata.obs.columns]
        if missing:
            raise InputError(f"colonne non trovate in adata.obs: {missing}")
        sheet = sample_sheet_from_obs(adata.obs, cols)
        name = Path(args.h5ad).stem
        if "library" not in technical:
            notes.append(LIBRARY_COLLAPSE_NOTE)
    result = run_design_audit(
        sheet, patient_col=args.patient_col, tissue_col=args.tissue_col,
        technical_cols=technical, outcome_cols=list(args.outcome_col), comparisons=comparisons,
        min_units=args.min_units, extra_notes=notes,
    )
    _print_design(result)
    out = save_report(args.out, design_result=result, dataset_name=name)
    print(f"[ok] report scritto in {out}")


class _ProgressLine:
    """Avanzamento su una sola riga aggiornata (su terminale) o una riga per passo (su file)."""

    def __init__(self) -> None:
        self.tty = sys.stderr.isatty()

    def __call__(self, msg: str) -> None:
        if self.tty and msg.startswith("["):
            sys.stderr.write("\r\033[K" + msg)
        else:
            sys.stderr.write(("\n" if self.tty else "") + msg + "\n")
        sys.stderr.flush()

    def close(self) -> None:
        """Chiude la riga di avanzamento sul terminale."""
        if self.tty:
            sys.stderr.write("\n")


def _cmd_leakage(args: argparse.Namespace) -> None:
    adata = ad.read_h5ad(args.h5ad)
    progress = _ProgressLine()
    result = run_leakage_audit(
        adata, target_col=args.target_col, patient_col=args.patient_col,
        n_folds=args.n_folds, min_patients_for_model_comparison=args.min_patients_for_model_comparison,
        seed=args.seed, benchmark=not args.rapido,
        max_train_cells=args.max_train_cells or None, progress=progress,
    )
    progress.close()
    print(result.narrative)
    print(f"Tempo di esecuzione: {result.elapsed_seconds:.0f} s.")
    out = save_report(args.out, leakage_result=result, dataset_name=Path(args.h5ad).stem)
    print(f"[ok] report scritto in {out}")


def _cmd_serve(args: argparse.Namespace) -> None:
    """Avvia Streamlit legato a localhost, con la telemetria disattivata: nessun dato (né
    statistica d'uso) lascia la macchina."""
    import subprocess
    import threading
    import webbrowser
    app = Path(__file__).resolve().parent / "app.py"
    url = f"http://localhost:{args.port}"
    # headless: niente richiesta interattiva di email al primo avvio; il browser lo apre questo comando
    cmd = [sys.executable, "-m", "streamlit", "run", str(app),
           "--server.address", "localhost", "--server.port", str(args.port),
           "--browser.gatherUsageStats", "false", "--server.headless", "true"]
    print(f"Web app su {url} (Ctrl+C per fermare). {USE_WARNING}", flush=True)
    if not args.no_browser:
        threading.Timer(3.0, lambda: webbrowser.open(url)).start()
    raise SystemExit(subprocess.call(cmd, cwd=str(app.parent)))


def _csv_list(value: str | None) -> list[str] | None:
    return None if value is None else [x.strip() for x in value.split(",") if x.strip()]


def _read_marker_map(path: Path | None) -> dict[str, list[str]] | None:
    if path is None:
        return None
    try:
        mm = json.loads(Path(path).read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as e:
        raise InputError(f"mappa dei marcatori '{path}' illeggibile: {e}") from e
    try:
        return validate_marker_map(mm)
    except ValueError as e:
        raise InputError(str(e)) from e


def _cmd_tcr(args: argparse.Namespace) -> None:
    if args.sperimentale_cd8 and not (args.marker_map and args.reference_compartment):
        raise InputError("--sperimentale-cd8 richiede --marker-map e --reference-compartment")
    if args.export_flags and not (args.marker_map and args.reference_compartment):
        raise InputError("--export-flags richiede --marker-map e --reference-compartment")
    adata = ad.read_h5ad(args.h5ad)
    contigs = parse_vdj_contigs(read_vdj_manifest(args.vdj_manifest))
    marker_map = _read_marker_map(args.marker_map)

    result = run_tcr_validation(
        adata, contigs, patient_col=args.patient_col, compartment_col=args.compartment_col,
        celltype_col=args.celltype_col, barcode_col=args.barcode_col, marker_map=marker_map,
        reference_compartment=args.reference_compartment, gene_col=args.gene_col, convention=args.convention,
        clone_error_labels=_csv_list(args.clone_error_labels),
        clone_marker_priority=_csv_list(args.clone_marker_priority),
        library_col=args.library_col, compartment_order=_csv_list(args.compartment_order),
        tie_rule=args.tie_rule,
    )
    print(result.barcode_match.sentence)
    print(result.narrative)
    if result.flag_coverage is not None:
        print(f"Flag per cellula valutabili (audit_label_vs_reference non NA): "
              f"{result.flag_coverage:.1%} delle cellule.")
    if args.export_flags:
        prefix = Path(args.out).parent / Path(args.h5ad).stem
        h5ad_path, csv_path = export_audited(adata, result, prefix)
        print(f"[ok] copia con i flag scritta in {h5ad_path} (il file originale non è toccato)")
        print(f"[ok] flag per cellula scritti in {csv_path}")
    if args.sperimentale_cd8:
        # flag assegnati per posizione di riga: corretto anche con obs_names ripetuti
        obs = adata.obs.copy()
        for c in result.cell_flags.columns:
            obs[c] = result.cell_flags[c].values
        cd8 = cd8_fraction_intervals(
            obs, patient_col=args.patient_col, compartment_col=args.compartment_col,
            celltype_col=args.celltype_col, reference_col="audit_reference_label",
            target_compartment=args.sperimentale_cd8, cd4_label=args.cd4_label,
            cd8_label=args.cd8_label, n_boot=args.n_boot, seed=args.seed)
        print("\n" + format_cd8_text(cd8))
        print("[nota] la frazione di CD8 è sperimentale: è stampata qui e non entra nel report né "
              "nella sintesi dei controlli.")
    out = save_report(args.out, tcr_result=result, dataset_name=Path(args.h5ad).stem)
    print(f"[ok] report scritto in {out}")


def _cmd_demo(args: argparse.Namespace) -> None:
    out_dir = Path(args.data_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    print(USE_WARNING + "\n")

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
    print(DEMO_LEAKAGE_NOTE)
    print(leakage_result.narrative)

    tcr_adata, contigs = make_tcr_validation_dataset(seed=args.seed)
    write_vdj_csvs(contigs, out_dir / "vdj")
    marker_map = {"CD4T": ["CD4"], "CD8T": ["CD8A", "CD8B"]}
    tcr_result = run_tcr_validation(
        tcr_adata, contigs, patient_col="patient_id", compartment_col="tissue",
        celltype_col="celltype", barcode_col="barcode",
        marker_map=marker_map, reference_compartment="PBMC",
    )
    print("\n=== Modulo B (dataset sintetico) ===")
    print(tcr_result.narrative)

    out = save_report(args.out, leakage_result=leakage_result, tcr_result=tcr_result,
                      dataset_name="demo sintetico", design_result=design_result)
    print(f"\n[ok] report scritto in {out}")
    print(f"[ok] file VDJ di esempio scritti in {out_dir / 'vdj'}")


def main(argv: list[str] | None = None) -> None:
    """Punto d'ingresso del comando ``audit-sc``."""
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="command", required=True)

    p_demo = sub.add_parser("demo", help="genera dati sintetici ed esegue disegno, Modulo A e Modulo B")
    p_demo.add_argument("--out", default="results/demo_report.html", type=Path)
    p_demo.add_argument("--data-dir", default="data/synthetic", type=Path)
    p_demo.add_argument("--seed", type=int, default=0)
    p_demo.set_defaults(func=_cmd_demo)

    p_des = sub.add_parser("design", help="audit del disegno e del confondimento (solo metadati)")
    src = p_des.add_mutually_exclusive_group(required=True)
    src.add_argument("--meta", type=Path, help="CSV/TSV dei metadati, una riga per campione/libreria")
    src.add_argument("--h5ad", type=Path, help="AnnData: si usa solo adata.obs, ridotto ai campioni")
    p_des.add_argument("--patient-col", required=True)
    p_des.add_argument("--tissue-col", default=None)
    p_des.add_argument("--technical", action="append", default=[],
                       help=f"ruolo=colonna, ruoli: {', '.join(TECHNICAL_ROLES)} (ripetibile, un ruolo una volta)")
    p_des.add_argument("--outcome-col", action="append", default=[], help="colonna di esito (ripetibile)")
    p_des.add_argument("--compare", action="append", default=[],
                       help="colonna:livello_a:livello_b (ripetibile; la colonna va dichiarata)")
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
    p_leak.add_argument("--rapido", "--no-benchmark", dest="rapido", action="store_true",
                        help="salta il confronto fra modelli: solo il classificatore di riferimento")
    p_leak.add_argument("--max-train-cells", type=int, default=20_000,
                        help="limite di cellule nel training di ogni fold, con sottocampionamento "
                             "stratificato per paziente (0 = nessun limite)")
    p_leak.add_argument("--out", default="results/leakage_report.html", type=Path)
    p_leak.set_defaults(func=_cmd_leakage)

    p_serve = sub.add_parser("serve", help="avvia la web app locale (solo su questa macchina)")
    p_serve.add_argument("--port", type=int, default=8501)
    p_serve.add_argument("--no-browser", action="store_true", help="non aprire il browser")
    p_serve.set_defaults(func=_cmd_serve)

    p_tcr = sub.add_parser("tcr", help="Modulo B: verifica dell'annotazione via TCR")
    p_tcr.add_argument("--h5ad", required=True, type=Path)
    p_tcr.add_argument("--vdj-manifest", required=True, type=Path,
                       help="CSV con colonne path,patient,compartment[,library]; percorsi relativi al manifest")
    p_tcr.add_argument("--patient-col", required=True)
    p_tcr.add_argument("--compartment-col", required=True)
    p_tcr.add_argument("--celltype-col", required=True)
    p_tcr.add_argument("--barcode-col", required=True)
    p_tcr.add_argument("--library-col", default=None,
                       help="colonna libreria di adata.obs, obbligatoria se il manifest ha la colonna library")
    p_tcr.add_argument("--marker-map", type=Path, default=None,
                       help="JSON {\"etichetta\": [\"gene1\", \"gene2\"]}, opzionale")
    p_tcr.add_argument("--gene-col", default=None,
                       help="colonna di adata.var con i nomi dei geni della mappa dei marcatori "
                            "(default: adata.var_names)")
    p_tcr.add_argument("--reference-compartment", default=None,
                       help="richiesto insieme a --marker-map per il tasso d'errore per compartimento")
    p_tcr.add_argument("--convention", choices=["both", "cell", "clone"], default="both",
                       help="definizione del tasso d'errore: 'cell' (per cellula, etichette fuori "
                            "mappa escluse), 'clone' (per clone-compartimento, etichette fuori mappa "
                            "contate come errore) o 'both' (default: entrambe affiancate)")
    p_tcr.add_argument("--clone-error-labels", default=None,
                       help="convenzione 'clone': etichette fuori mappa ammesse e contate come errore, "
                            "separate da virgola (default: tutte, anche mieloidi e B; la tesi usa NK)")
    p_tcr.add_argument("--clone-marker-priority", default=None,
                       help="convenzione 'clone': ordine di positività esclusiva, es. CD8T,CD4T "
                            "(default: positività non esclusiva)")
    p_tcr.add_argument("--tie-rule", choices=["exclude", "error", "correct"], default="exclude",
                       help="convenzione 'clone': unità con l'etichetta di maggioranza in pareggio "
                            "escluse (default, come la tesi) oppure, come sensibilità, contate come "
                            "errore o come corrette")
    p_tcr.add_argument("--compartment-order", default=None,
                       help="ordine dei compartimenti nelle differenze appaiate A - B, es. "
                            "Tumor,Adjacent_normal (default: alfabetico)")
    p_tcr.add_argument("--export-flags", action="store_true",
                       help="scrive <nome>_audited.h5ad (copia con i flag) e <nome>_audit_flags.csv "
                            "nella cartella di --out; richiede --marker-map")
    p_tcr.add_argument("--sperimentale-cd8", default=None, metavar="COMPARTIMENTO",
                       help="SPERIMENTALE: intervalli della frazione di CD8 per paziente nel compartimento "
                            "indicato, solo stampati qui (non nel report); richiede --marker-map")
    p_tcr.add_argument("--cd4-label", default="CD4T", help="etichetta CD4 per --sperimentale-cd8")
    p_tcr.add_argument("--cd8-label", default="CD8T", help="etichetta CD8 per --sperimentale-cd8")
    p_tcr.add_argument("--n-boot", type=int, default=2000, help="repliche bootstrap di --sperimentale-cd8")
    p_tcr.add_argument("--seed", type=int, default=0, help="seme di --sperimentale-cd8")
    p_tcr.add_argument("--out", default="results/tcr_report.html", type=Path)
    p_tcr.set_defaults(func=_cmd_tcr)

    args = ap.parse_args(argv)
    try:
        args.func(args)
    except (InputError, ValueError, FileNotFoundError, KeyError) as e:
        msg = e.args[0] if isinstance(e, KeyError) and e.args else e
        sys.exit(f"errore: {msg}")


if __name__ == "__main__":
    main()
