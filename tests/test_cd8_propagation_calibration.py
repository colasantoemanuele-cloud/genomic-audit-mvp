"""Calibrazione della propagazione dell'errore sulla frazione di CD8 (Intervento 3).

OBBLIGATORIO prima di esporre core/cd8_propagation.py in CLI, app o report.

Bande dichiarate PRIMA di eseguire i test (non vanno allargate per farli passare):
  - Errore simmetrico iniettato (CD4->CD8 = CD8->CD4 = 0.15, piu' 0.05 verso "altro"):
    copertura dell'IC 95% (scenario 1x) della frazione vera di CD8 di ciascun paziente
    in [90%, 99%], su 200 repliche x 10 pazienti.
  - Errore ASIMMETRICO (CD4->CD8 = 0.24, tre volte CD8->CD4 = 0.08, piu' 0.05 verso
    "altro"): stessa banda [90%, 99%].
  - Controllo di sensibilita' (non una banda): con errore iniettato, l'intervallo della
    frazione RIPORTATA senza correzione (solo Beta sui conteggi) copre il vero meno
    spesso della versione corretta -- dimostra che la correzione fa qualcosa.
  - Rifiuti: con meno di 5 pazienti con cellule di riferimento nel compartimento, o con
    matrice di confusione mal condizionata (J vicino a zero), nessun intervallo.
  - Equivalenza: il ricampionamento dei pazienti usato qui riproduce ESATTAMENTE (1e-12)
    l'IC di core.stats.cluster_bootstrap sullo stesso vettore e con lo stesso seme.

I semi sono fissi: i test sono deterministici.
"""

from __future__ import annotations

import functools

import numpy as np
import pandas as pd

from core.cd8_propagation import cd8_fraction_intervals, patient_resample_draws
from core.stats import cluster_bootstrap
from core.synthetic import make_cd8_fraction_dataset

R = 200
COVERAGE_BAND = (0.90, 0.99)
SYMMETRIC = dict(p_cd4_to_cd8=0.15, p_cd8_to_cd4=0.15, p_to_other=0.05)
ASYMMETRIC = dict(p_cd4_to_cd8=0.24, p_cd8_to_cd4=0.08, p_to_other=0.05)


def _run(errors: dict, seed: int):
    obs, truth = make_cd8_fraction_dataset(n_patients=10, seed=seed, **errors)
    res = cd8_fraction_intervals(obs, patient_col="patient", compartment_col="compartment",
                                 celltype_col="celltype", reference_col="audit_reference_label",
                                 target_compartment="Tumor", n_boot=400, seed=seed)
    return res, truth


@functools.lru_cache(maxsize=4)
def _coverage(kind: str) -> tuple[int, int, int, int]:
    errors = SYMMETRIC if kind == "sym" else ASYMMETRIC
    cov, cov_naive, n, n_refused = 0, 0, 0, 0
    for rep in range(R):
        res, truth = _run(errors, seed=70_000 + rep + (0 if kind == "sym" else 5_000))
        assert res.refused_reason is None, res.refused_reason
        for row in res.patients:
            iv = row.scenarios[1.0]
            if iv.low is None:
                n_refused += 1
                continue
            n += 1
            t = truth[row.patient]
            cov += iv.low <= t <= iv.high
            cov_naive += row.naive_low <= t <= row.naive_high
    return cov, cov_naive, n, n_refused


def _check(kind: str) -> None:
    cov, cov_naive, n, n_refused = _coverage(kind)
    frac = cov / n
    print(f"\n[{kind}] copertura IC95% (scenario 1x): {cov}/{n} = {frac:.3f} (banda {COVERAGE_BAND}); "
          f"frazione riportata senza correzione: {cov_naive}/{n} = {cov_naive / n:.3f}; "
          f"intervalli rifiutati: {n_refused}")
    assert n_refused == 0
    assert COVERAGE_BAND[0] <= frac <= COVERAGE_BAND[1]
    assert cov_naive < cov


def test_coverage_symmetric_error():
    _check("sym")


def test_coverage_asymmetric_error():
    _check("asym")


def test_scenarios_are_ordered_and_all_reported():
    res, _ = _run(ASYMMETRIC, seed=1)
    for row in res.patients:
        assert set(row.scenarios) == {0.5, 1.0, 2.0}
        assert "0.5x" in row.sentence and "2x" in row.sentence


def test_refuses_with_fewer_than_five_reference_patients():
    obs, _ = make_cd8_fraction_dataset(n_patients=10, seed=3, **SYMMETRIC)
    keep_ref = {"P00", "P01", "P02", "P03"}
    obs.loc[~obs.patient.isin(keep_ref), "audit_reference_label"] = pd.NA
    res = cd8_fraction_intervals(obs, "patient", "compartment", "celltype",
                                 "audit_reference_label", "Tumor", n_boot=100, seed=0)
    assert res.refused_reason is not None and "5" in res.refused_reason
    assert all(r.scenarios[1.0].low is None for r in res.patients)


def test_refuses_ill_conditioned_matrix():
    # errore cosi' alto che le due probabilita' di chiamata CD8 quasi coincidono
    obs, _ = make_cd8_fraction_dataset(n_patients=10, seed=4, p_cd4_to_cd8=0.45,
                                       p_cd8_to_cd4=0.45, p_to_other=0.05)
    res = cd8_fraction_intervals(obs, "patient", "compartment", "celltype",
                                 "audit_reference_label", "Tumor", n_boot=200, seed=0)
    assert all(r.scenarios[1.0].low is None for r in res.patients)
    assert "mal condizionata" in res.patients[0].sentence


def test_resampling_matches_cluster_bootstrap_exactly():
    rng = np.random.default_rng(0)
    groups = np.repeat([f"P{i}" for i in range(9)], rng.integers(3, 12, 9))
    values = rng.normal(0.3, 1, len(groups))
    ref = cluster_bootstrap(values, groups, n_boot=300, seed=7)
    means = np.array([values[idx].mean() for idx in patient_resample_draws(groups, n_boot=300, seed=7)])
    lo, hi = np.percentile(means, (2.5, 97.5))
    assert abs(lo - ref.ci_low) <= 1e-12 and abs(hi - ref.ci_high) <= 1e-12


def test_end_to_end_from_module_b_flags():
    """Non e' una calibrazione: la funzione accetta direttamente adata.obs + i flag del
    Modulo B e produce un risultato per ogni paziente del compartimento bersaglio."""
    from core.synthetic import make_tcr_validation_dataset
    from core.tcr_validation import run_tcr_validation

    adata, contigs = make_tcr_validation_dataset(n_patients=12, n_clones_per_patient=15,
                                                 injected_excess=0.2, seed=0)
    res_b = run_tcr_validation(adata, contigs, patient_col="patient_id", compartment_col="tissue",
                               celltype_col="celltype", barcode_col="barcode", n_boot=200, seed=0,
                               marker_map={"CD4T": ["CD4"], "CD8T": ["CD8A", "CD8B"]},
                               reference_compartment="PBMC")
    obs = adata.obs.join(res_b.cell_flags)
    res = cd8_fraction_intervals(obs, "patient_id", "tissue", "celltype", "audit_reference_label",
                                 "Tumor", n_boot=200, seed=0)
    assert res.refused_reason is None
    assert len(res.patients) == 12
    assert res.matrix is not None and res.matrix.shape == (2, 3)
    assert "Assunzioni" in res.narrative
