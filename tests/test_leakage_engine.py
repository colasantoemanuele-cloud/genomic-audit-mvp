"""Test del motore del Modulo A per dati reali: metrica, assenza di leakage nel
preprocessing, sparsita', sottocampionamento, robustezza di formato."""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest
import scipy.sparse as sp

from core.leakage_audit import (
    HVGSelector,
    _pipeline,
    _reference_model,
    cap_training,
    fold_macro_f1,
    run_leakage_audit,
)
from core.synthetic import make_leakage_dataset


def test_macro_f1_ignores_classes_absent_from_test_fold():
    """Regressione: prima un fold perfetto con una classe assente valeva 0.667."""
    y = np.array(["a", "a", "b", "b"])
    score, absent = fold_macro_f1(y, y.copy(), np.array(["a", "b", "c"]))
    assert score == 1.0 and absent == ["c"]
    # una predizione della classe assente resta un errore (abbassa il recall di "a")
    score2, _ = fold_macro_f1(y, np.array(["a", "c", "b", "b"]), np.array(["a", "b", "c"]))
    assert score2 < 1.0


def test_hvg_selection_uses_only_training_cells():
    rng = np.random.default_rng(0)
    X = sp.csr_matrix(rng.poisson(1.0, (100, 20)).astype(float))
    X = X.tolil()
    X[:80, 19] = 0.0     # costante nel training...
    X[80:, 19] = 500.0   # ...variabilissimo SOLO nelle cellule di test
    X = X.tocsr()
    sel = HVGSelector(n_top=5).fit(X[:80])
    assert 19 not in sel.genes_
    assert 19 in HVGSelector(n_top=5).fit(X).genes_


def test_pipeline_never_densifies_the_gene_matrix():
    adata = make_leakage_dataset(n_patients=6, cells_per_patient=20, n_genes=300, seed=0)
    X, y = sp.csr_matrix(adata.X), adata.obs["label"].astype(str).values
    for reduce in (False, True):
        pipe = _pipeline(_reference_model(0), reduce, n_hvg=100, n_svd=10, seed=0)
        Z = X
        for _, step in pipe.steps[:-1]:
            Z = step.fit(Z, y).transform(Z)
            if not sp.issparse(Z):
                assert Z.shape[1] <= 10, "solo l'uscita della SVD puo' essere densa"


def test_cap_training_is_stratified_by_patient_and_deterministic():
    groups = np.repeat(np.array([f"P{i}" for i in range(5)]), [100, 200, 300, 400, 1000])
    tr = np.arange(len(groups))
    a = cap_training(tr, groups, 200, seed=3)
    b = cap_training(tr, groups, 200, seed=3)
    assert np.array_equal(a, b) and len(a) <= 200
    counts = pd.Series(groups[a]).value_counts()
    assert set(counts.index) == {f"P{i}" for i in range(5)}
    assert counts["P4"] > counts["P0"]
    assert np.array_equal(cap_training(tr, groups, 10_000, seed=3), tr)


def _small():
    # seme 6: 4 pazienti per classe (con il seme 4 una classe ha un solo paziente)
    return make_leakage_dataset(n_patients=8, cells_per_patient=15, n_genes=150, seed=6)


def test_single_patient_class_fails_fast_with_clear_message():
    a = make_leakage_dataset(n_patients=8, cells_per_patient=15, n_genes=150, seed=4)  # 7 contro 1
    with pytest.raises(ValueError, match="servono almeno 2 pazienti per classe"):
        run_leakage_audit(a, "label", "patient_id", benchmark=True)


def test_rapido_skips_model_comparison_even_with_enough_patients():
    r = run_leakage_audit(_small(), "label", "patient_id", benchmark=False)
    assert r.model_comparison is None and r.settings["benchmark"] is False


def test_rejects_non_integer_matrix_with_clear_message():
    a = _small()
    a.X = sp.csr_matrix(np.log1p(a.X.toarray()) + 0.1)
    with pytest.raises(ValueError, match="conteggi"):
        run_leakage_audit(a, "label", "patient_id", benchmark=False)


def test_dense_and_categorical_inputs_give_identical_results():
    a = _small()
    ref = run_leakage_audit(a, "label", "patient_id", benchmark=False)
    d = a.copy()
    d.X = np.asarray(d.X.toarray())
    c = a.copy()
    c.obs["label"] = pd.Categorical(c.obs["label"])
    c.obs["patient_id"] = pd.Categorical(c.obs["patient_id"])
    for other in (run_leakage_audit(d, "label", "patient_id", benchmark=False),
                  run_leakage_audit(c, "label", "patient_id", benchmark=False)):
        assert other.grouped.fold_scores == ref.grouped.fold_scores
        assert other.random.fold_scores == ref.random.fold_scores


def test_progress_reports_planned_fits_and_xai_is_descriptive():
    msgs = []
    r = run_leakage_audit(_small(), "label", "patient_id", benchmark=True, progress=msgs.append)
    assert "Addestramenti previsti" in msgs[0]
    total = int(msgs[0].split("Addestramenti previsti: ")[1].split()[0].rstrip("."))
    assert sum(m.startswith("[") for m in msgs) == total
    assert r.model_comparison is not None
    assert r.xai.grouped_matrix.shape == (r.n_folds_used, r.n_folds_used)
    assert 0.0 <= r.xai.grouped_mean <= 1.0


def test_parallel_folds_give_identical_results():
    a = _small()
    seq = run_leakage_audit(a, "label", "patient_id", benchmark=True, n_jobs=1)
    par = run_leakage_audit(a, "label", "patient_id", benchmark=True, n_jobs=3)
    assert par.grouped.fold_scores == seq.grouped.fold_scores
    assert par.random.fold_scores == seq.random.fold_scores
    assert par.grouped.top_genes == seq.grouped.top_genes
    for m in seq.model_comparison.scores:
        assert par.model_comparison.scores[m].fold_scores == seq.model_comparison.scores[m].fold_scores
