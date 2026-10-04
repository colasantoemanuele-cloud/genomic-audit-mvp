"""Test di integrazione sui generatori sintetici: verificano che gli effetti iniettati
deliberatamente (leakage per Modulo A, discordanza cross-compartimento per Modulo B)
vengano effettivamente rilevati dal motore analitico. Non sono test di calibrazione
(quelli sono in test_stats_calibration.py, obbligatori prima dell'interfaccia): qui si
verifica che i moduli si comportino correttamente end-to-end su dati con struttura nota.
"""

from __future__ import annotations

from core.leakage_audit import run_leakage_audit
from core.synthetic import make_leakage_dataset, make_tcr_validation_dataset, write_vdj_csvs
from core.tcr_validation import (
    DiscordanceResult,
    _discordance_narrative,
    parse_vdj_contigs,
    run_tcr_validation,
)


# --------------------------------------------------------------------------- #
# Modulo A
# --------------------------------------------------------------------------- #
def test_leakage_dataset_shows_gap():
    adata = make_leakage_dataset(n_patients=16, cells_per_patient=30, seed=1)
    result = run_leakage_audit(adata, target_col="label", patient_col="patient_id",
                                min_patients_for_model_comparison=8, seed=0)

    assert result.n_patients == 16
    assert result.n_cells == 16 * 30
    # lo split casuale (leakage) deve sovrastimare parecchio l'accuratezza reale
    assert result.gap > 0.15, f"gap troppo piccolo per essere rilevabile: {result.gap:.3f}"
    assert result.random.mean > 0.85, f"split casuale non sfrutta il fingerprint: {result.random.mean:.3f}"
    # lo split per paziente non deve raggiungere le stesse prestazioni (nessuna
    # generalizzazione al fingerprint di un paziente mai visto)
    assert result.grouped.mean < result.random.mean - 0.1
    assert "onesta" in result.narrative and "leakage" not in result.narrative.lower()  # no jargon


def test_leakage_dataset_model_comparison_runs_when_enough_patients():
    adata = make_leakage_dataset(n_patients=10, cells_per_patient=20, seed=2)
    result = run_leakage_audit(adata, target_col="label", patient_col="patient_id",
                                min_patients_for_model_comparison=8, seed=0)
    assert result.model_comparison is not None
    assert set(result.model_comparison.scores) == {"logreg", "random_forest", "hist_gb"}
    assert result.model_comparison.best_model in result.model_comparison.scores
    assert len(result.model_comparison.comparisons) == 2  # gli altri due modelli vs il migliore


def test_leakage_dataset_skips_model_comparison_below_threshold():
    adata = make_leakage_dataset(n_patients=6, cells_per_patient=20, seed=3)
    result = run_leakage_audit(adata, target_col="label", patient_col="patient_id",
                                min_patients_for_model_comparison=8, seed=0)
    assert result.model_comparison is None


# --------------------------------------------------------------------------- #
# Modulo B
# --------------------------------------------------------------------------- #
def test_tcr_validation_detects_injected_discordance():
    adata, contigs = make_tcr_validation_dataset(
        n_patients=12, n_clones_per_patient=15, injected_excess=0.35, seed=0)
    result = run_tcr_validation(
        adata, contigs, patient_col="patient_id", compartment_col="tissue",
        celltype_col="celltype", barcode_col="barcode", n_boot=1000, seed=0,
    )
    assert result.n_cells_with_tcr > 0
    assert result.discordance.sufficient
    assert result.discordance.ci_low > 0, (
        f"l'eccesso iniettato non è rilevato: IC95%=[{result.discordance.ci_low:.3f}, "
        f"{result.discordance.ci_high:.3f}]"
    )
    assert "eccesso di discordanza" in result.narrative


def test_tcr_validation_no_injected_effect_is_not_significant():
    """Controllo negativo: senza eccesso iniettato (stesso mislabel rate ovunque),
    l'eccesso misurato deve restare piccolo in valore assoluto -- di un ordine di
    grandezza sotto l'effetto iniettato nel test precedente (~0.3). Non si richiede che
    l'IC includa esattamente lo zero su un singolo seme: cluster_bootstrap è già
    validato come calibrato (~2-8% di falsi positivi sotto H0, vedi
    test_stats_calibration.py), quindi un'occasionale esclusione dello zero per puro
    rumore campionario è attesa, non un bug."""
    adata, contigs = make_tcr_validation_dataset(
        n_patients=12, n_clones_per_patient=15, injected_excess=0.0, seed=0)
    result = run_tcr_validation(
        adata, contigs, patient_col="patient_id", compartment_col="tissue",
        celltype_col="celltype", barcode_col="barcode", n_boot=1000, seed=0,
    )
    assert result.discordance.sufficient
    assert abs(result.discordance.mean_excess) < 0.05, (
        f"eccesso misurato troppo grande senza alcun effetto iniettato: {result.discordance.mean_excess:+.3f}"
    )


def test_tcr_validation_marker_error_rate():
    adata, contigs = make_tcr_validation_dataset(
        n_patients=12, n_clones_per_patient=15, injected_excess=0.35,
        compartments=("PBMC", "Tumor"), noisy_compartment="Tumor", seed=0,
    )
    marker_map = {"CD4T": ["CD4"], "CD8T": ["CD8A", "CD8B"]}
    result = run_tcr_validation(
        adata, contigs, patient_col="patient_id", compartment_col="tissue",
        celltype_col="celltype", barcode_col="barcode", n_boot=1000, seed=0,
        marker_map=marker_map, reference_compartment="PBMC",
    )
    assert result.marker_error is not None
    assert "Tumor" in result.marker_error.by_compartment
    tumor_res = result.marker_error.by_compartment["Tumor"]
    assert tumor_res.sufficient
    assert tumor_res.mean > 0.15, f"tasso d'errore nel compartimento rumoroso troppo basso: {tumor_res.mean:.3f}"


def test_tcr_validation_insufficient_patients_flagged():
    adata, contigs = make_tcr_validation_dataset(
        n_patients=3, n_clones_per_patient=10, injected_excess=0.35, seed=0)
    result = run_tcr_validation(
        adata, contigs, patient_col="patient_id", compartment_col="tissue",
        celltype_col="celltype", barcode_col="barcode", n_boot=500, seed=0,
        min_patients_for_ci=5,
    )
    assert not result.discordance.sufficient
    assert "insufficiente" in result.narrative


def test_tcr_narrative_handles_significant_negative_excess():
    """Regressione: un IC interamente negativo (es. [-0.010, -0.001]) esclude lo zero
    esattamente quanto uno interamente positivo -- non è "non distinguibile dal rumore".
    Osservato realmente su un dataset con injected_excess=0.0 e pochi pazienti (n=6):
    la vecchia narrativa (if ci_low > 0 ... else "l'intervallo include lo zero") lo
    descriveva in modo scorretto perché controllava solo il lato positivo."""
    d = DiscordanceResult(
        n_pairs=90, n_patients=6, mean_excess=-0.005,
        ci_low=-0.010, ci_high=-0.001, sufficient=True,
        by_compartment_pair=None,
    )
    narrative = _discordance_narrative(d, min_patients_for_ci=5)
    assert "include lo zero" not in narrative
    assert "NEGATIVA" in narrative


def test_vdj_csv_roundtrip(tmp_path):
    """Il percorso di ingestione reale (scrittura CSV Cell Ranger -> parse_vdj_contigs)
    deve produrre lo stesso risultato usato direttamente in memoria dal generatore."""
    _, contigs = make_tcr_validation_dataset(n_patients=3, n_clones_per_patient=5, seed=0)
    files = write_vdj_csvs(contigs, tmp_path)
    assert len(files) == 3 * 2  # pazienti x compartimenti

    parsed = parse_vdj_contigs(files)
    assert set(parsed.patient.unique()) == set(contigs.patient.unique())
    assert set(parsed.chain.unique()) <= {"TRA", "TRB"}
    assert len(parsed) == len(contigs)  # nessun duplicato (barcode,chain) nel sintetico
