"""Nessun errore silenzioso (criterio 3 della revisione, A8-A10) e casi del Modulo B (A4-A5).

Chiavi duplicate, input malformati e scelte d'ordine implicite devono dare un errore chiaro o un
risultato che non dipende dall'ordine delle righe, mai un risultato diverso in silenzio.
"""

from __future__ import annotations

import json
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import pytest
import scipy.sparse as sp

from cli import main as cli_main
from core.design_audit import (
    LIBRARY_COLLAPSE_NOTE,
    parse_comparison_lines,
    run_design_audit,
    sample_sheet_from_obs,
)
from core.io import load_matrix_market, read_table, read_vdj_manifest
from core.synthetic import make_tcr_validation_dataset, write_vdj_csvs
from core.tcr_validation import (
    _majority_error,
    compute_marker_positivity,
    export_audited,
    match_barcodes,
    parse_vdj_contigs,
    run_tcr_validation,
)
from core.verdict import VERDE, tcr_summary

ROOT = Path(__file__).resolve().parents[1]
MARKERS = {"CD4T": ["CD4"], "CD8T": ["CD8A", "CD8B"]}


def _run(adata, contigs, **kw):
    return run_tcr_validation(adata, contigs, patient_col="patient_id", compartment_col="tissue",
                              celltype_col="celltype", barcode_col="barcode", **kw)


# --------------------------------------------------------------------------- #
# A4: filtro del tasso per cellula sulle etichette della marker_map
# --------------------------------------------------------------------------- #
def test_label_never_resolved_as_identity_still_counts_as_error():
    """Il caso della revisione: tutti i cloni sono CD8 e gli errori sono etichette CD4T. Nessun
    clone ha identità CD4T: con il vecchio filtro (etichette fra le identità risolte) le cellule
    CD4T, cioè proprio gli errori, uscivano dal calcolo e il tasso per cellula era 0."""
    adata, contigs = make_tcr_validation_dataset(n_patients=8, n_clones_per_patient=8, identities=("CD8T",),
                                                 distractor_label="CD4T", injected_excess=0.35, seed=4)
    res = _run(adata, contigs, marker_map=MARKERS, reference_compartment="PBMC")
    cell = res.conventions["cell"].by_compartment["Tumor"]
    clone = res.conventions["clone"].by_compartment["Tumor"]
    print(f"\n[A4] per cellula {cell.mean:.3f} su {cell.n_obs} cellule, per clone {clone.mean:.3f}")
    assert cell.mean > 0.25 and clone.mean > 0.25
    tumor_cd4 = (adata.obs["tissue"].values == "Tumor") & (adata.obs["celltype"].values == "CD4T")
    flags = res.cell_flags["audit_label_vs_reference"][tumor_cd4]
    assert flags.notna().all() and flags.astype(bool).all()


# --------------------------------------------------------------------------- #
# A5: compartimenti sbilanciati e sensibilità a pesi uguali
# --------------------------------------------------------------------------- #
def test_unbalanced_compartments_weighted_excess_negative_equal_weights_positive():
    """3 cellule per clone nel sangue, 60 nel tumore con il 40% di errore: l'eccesso pesato ha
    valore atteso negativo anche se le etichette differiscono davvero; quello a pesi uguali no."""
    adata, contigs = make_tcr_validation_dataset(
        n_patients=8, n_clones_per_patient=6, baseline_mislabel_rate=0.0, injected_excess=0.40,
        cells_per_clone_by_compartment={"PBMC": (3, 3), "Tumor": (60, 60)}, seed=2)
    res = _run(adata, contigs)
    d = res.discordance
    print(f"\n[A5] pesato {d.mean_excess:+.3f} [{d.ci_low:+.3f}, {d.ci_high:+.3f}], "
          f"pesi uguali {d.equal_weights.mean:+.3f} [{d.equal_weights.ci_low:+.3f}, {d.equal_weights.ci_high:+.3f}]")
    assert d.mean_excess < 0 < d.equal_weights.mean
    assert d.equal_weights.ci_low > 0
    check = tcr_summary(res).checks[0]
    assert check.state != VERDE
    assert "NEGATIVA" in res.narrative and "≥ 0" not in res.narrative.split("Sensibilità")[0]


def test_non_shared_clones_do_not_create_pairs_and_dropout_ambient_run():
    full, c_full = make_tcr_validation_dataset(n_patients=6, n_clones_per_patient=10, seed=5)
    half, c_half = make_tcr_validation_dataset(n_patients=6, n_clones_per_patient=10, seed=5, shared_fraction=0.5,
                                               marker_dropout=0.3, ambient_rate=0.2)
    r_full, r_half = _run(full, c_full), _run(half, c_half, marker_map=MARKERS, reference_compartment="PBMC")
    assert r_half.discordance.n_pairs < r_full.discordance.n_pairs
    assert r_half.conventions is not None


# --------------------------------------------------------------------------- #
# A8: chiavi duplicate, librerie, scelte d'ordine
# --------------------------------------------------------------------------- #
def test_matrix_market_refuses_duplicate_barcodes_in_metadata(tmp_path):
    import scipy.io
    X = sp.random(5, 4, density=0.5, format="coo", random_state=0)
    scipy.io.mmwrite(tmp_path / "matrix.mtx", X)
    (tmp_path / "barcodes.tsv").write_text("A\nB\nC\nD\n")
    (tmp_path / "features.tsv").write_text("".join(f"G{i}\tG{i}\n" for i in range(5)))
    meta = pd.DataFrame({"bc": ["A", "B", "C", "D", "B"], "patient": list("12345")})
    with pytest.raises(ValueError, match="barcode ripetuto"):
        load_matrix_market(tmp_path, meta, "bc")


def test_read_table_refuses_duplicate_column_names(tmp_path):
    (tmp_path / "t.csv").write_text("patient,patient,tissue\n1,2,T\n")
    with pytest.raises(ValueError, match="stesso nome"):
        read_table(tmp_path / "t.csv")
    with pytest.raises(FileNotFoundError):
        read_table(tmp_path / "manca.csv")


def test_duplicate_patient_compartment_barcode_in_metadata_is_an_error():
    adata, contigs = make_tcr_validation_dataset(n_patients=5, n_clones_per_patient=5, seed=1)
    obs = adata.obs[["patient_id", "tissue", "celltype", "barcode"]].astype(str).copy()
    obs.columns = ["patient", "compartment", "celltype", "barcode"]
    obs.iloc[1, obs.columns.get_loc("barcode")] = obs.iloc[0]["barcode"]
    obs.iloc[1, obs.columns.get_loc("compartment")] = obs.iloc[0]["compartment"]
    obs.iloc[1, obs.columns.get_loc("patient")] = obs.iloc[0]["patient"]
    from core.tcr_validation import build_clonotypes
    with pytest.raises(ValueError, match="--library-col"):
        match_barcodes(obs, build_clonotypes(contigs))


def _two_libraries(tmp_path):
    """Due librerie per ogni coppia paziente-compartimento con gli STESSI barcode e cloni diversi:
    senza la libreria il TCR di una cellula non è attribuibile."""
    a1, c1 = make_tcr_validation_dataset(n_patients=6, n_clones_per_patient=6, seed=7)
    a2 = a1.copy()
    c2 = c1.assign(cdr3_nt="L2_" + c1["cdr3_nt"])
    for a, tag in ((a1, "L1"), (a2, "L2")):
        a.obs["library"] = a.obs["patient_id"].astype(str) + "_" + a.obs["tissue"].astype(str) + "_" + tag
    both = ad.concat([a1, a2], index_unique="-")
    files = []
    for tag, c in (("L1", c1), ("L2", c2)):
        for path, pat, comp in write_vdj_csvs(c, tmp_path / tag):
            files.append((path, pat, comp, f"{pat}_{comp}_{tag}"))
    return both, files


def test_two_vdj_files_on_same_pair_need_a_library(tmp_path):
    both, files = _two_libraries(tmp_path)
    with pytest.raises(ValueError, match="library"):
        parse_vdj_contigs([f[:3] for f in files])
    contigs = parse_vdj_contigs(files)
    with pytest.raises(ValueError, match="--library-col"):
        _run(both, contigs)
    res = _run(both, contigs, library_col="library", marker_map=MARKERS, reference_compartment="PBMC")
    assert res.barcode_match.fraction == 1.0
    assert res.n_cells_with_tcr == both.n_obs


def test_contig_choice_does_not_depend_on_row_order(tmp_path):
    rows = [{"barcode": "AAA-1", "chain": "TRB", "cdr3_nt": "TGTGCC", "umis": 5, "reads": 50},
            {"barcode": "AAA-1", "chain": "TRB", "cdr3_nt": "TGTAAA", "umis": 5, "reads": 50},
            {"barcode": "AAA-1", "chain": "TRB", "cdr3_nt": "TGTCCC", "umis": 5, "reads": 10}]
    picks = set()
    for order in ([0, 1, 2], [2, 1, 0], [1, 2, 0]):
        f = tmp_path / f"v{order[0]}.csv"
        pd.DataFrame([rows[i] for i in order]).to_csv(f, index=False)
        picks.add(parse_vdj_contigs([(f, "P1", "Tumor")])["cdr3_nt"].iloc[0])
    assert picks == {"TGTAAA"}  # stesse UMI e reads: vince la CDR3 minore in ordine alfabetico


def test_marker_positivity_requires_raw_integer_counts():
    adata, _ = make_tcr_validation_dataset(n_patients=3, n_clones_per_patient=3, seed=0)
    adata.X = sp.csr_matrix(np.log1p(adata.X.toarray()))
    with pytest.raises(ValueError, match="conteggi grezzi"):
        compute_marker_positivity(adata, MARKERS)


def test_clone_majority_ties_do_not_depend_on_row_order():
    # regola di pdac-ml: pareggio escluso (NaN); sensibilità: contato come errore o come corretto
    e, tie = _majority_error(pd.Series(["CD8T", "CD4T"]), "CD4T")
    assert tie and np.isnan(e)
    assert _majority_error(pd.Series(["CD4T", "CD8T", "NK"]), "CD4T", "error") == (1.0, True)
    assert _majority_error(pd.Series(["NK", "CD8T"]), "CD4T", "correct") == (0.0, True)
    assert _majority_error(pd.Series(["CD8T", "CD8T", "CD4T"]), "CD4T", "correct") == (1.0, False)
    assert _majority_error(pd.Series(["CD4T", "CD4T", "NK"]), "CD4T", "error") == (0.0, False)
    with pytest.raises(ValueError, match="tie_rule"):
        _majority_error(pd.Series(["CD4T"]), "CD4T", "random")
    adata, contigs = make_tcr_validation_dataset(n_patients=6, n_clones_per_patient=8, injected_excess=0.45, seed=9)
    res = _run(adata, contigs, marker_map=MARKERS, reference_compartment="PBMC")
    perm = np.random.default_rng(0).permutation(adata.n_obs)
    res2 = _run(adata[perm].copy(), contigs.sample(frac=1, random_state=1), marker_map=MARKERS,
                reference_compartment="PBMC")
    a, b = res.conventions["clone"].by_compartment["Tumor"], res2.conventions["clone"].by_compartment["Tumor"]
    assert a.mean == pytest.approx(b.mean, abs=1e-12) and a.n_obs == b.n_obs
    assert res.conventions["clone"].n_ties == res2.conventions["clone"].n_ties


def test_clone_ties_excluded_by_default_and_sensitivities_bracket_it():
    """Le unità in pareggio non entrano nel tasso con la regola di default; le due sensibilità
    usano le stesse unità più i pareggi, contati come errore o come corretti."""
    adata, contigs = make_tcr_validation_dataset(n_patients=6, n_clones_per_patient=8, injected_excess=0.45, seed=9)
    out = {r: _run(adata, contigs, marker_map=MARKERS, reference_compartment="PBMC", tie_rule=r).conventions["clone"]
           for r in ("exclude", "error", "correct")}
    ex, er, co = (out[r].by_compartment["Tumor"] for r in ("exclude", "error", "correct"))
    k = out["exclude"].n_ties["Tumor"]
    assert k > 0 and out["error"].n_ties == out["exclude"].n_ties
    assert er.n_obs == co.n_obs == ex.n_obs + k
    assert er.mean * er.n_obs == pytest.approx(ex.mean * ex.n_obs + k)
    assert co.mean * co.n_obs == pytest.approx(ex.mean * ex.n_obs)
    assert "esclusi" in out["exclude"].definition and "come errore" in out["error"].definition


def test_export_refuses_non_unique_obs_names(tmp_path):
    adata, contigs = make_tcr_validation_dataset(n_patients=5, n_clones_per_patient=5, seed=3)
    adata.obs_names = ["c"] * adata.n_obs
    res = _run(adata, contigs, marker_map=MARKERS, reference_compartment="PBMC")
    with pytest.raises(ValueError, match="unici"):
        export_audited(adata, res, tmp_path / "x")


def _tcr_files(tmp_path, adata, contigs):
    h5 = tmp_path / "dati.h5ad"
    adata.write_h5ad(h5)
    files = write_vdj_csvs(contigs, tmp_path / "vdj")
    man = tmp_path / "manifest.csv"
    pd.DataFrame([{"path": f"vdj/{p.name}", "patient": pat, "compartment": comp} for p, pat, comp in files]).to_csv(
        man, index=False)
    (tmp_path / "markers.json").write_text(json.dumps(MARKERS))
    return h5, man


def test_cli_sperimentale_cd8_assigns_flags_by_position_with_duplicate_obs_names(tmp_path, capsys):
    adata, contigs = make_tcr_validation_dataset(n_patients=6, n_clones_per_patient=10, seed=2)
    adata.obs_names = [f"c{i % 7}" for i in range(adata.n_obs)]  # nomi ripetuti
    h5, man = _tcr_files(tmp_path, adata, contigs)
    cli_main(["tcr", "--h5ad", str(h5), "--vdj-manifest", str(man), "--patient-col", "patient_id",
              "--compartment-col", "tissue", "--celltype-col", "celltype", "--barcode-col", "barcode",
              "--marker-map", str(tmp_path / "markers.json"), "--reference-compartment", "PBMC",
              "--sperimentale-cd8", "Tumor", "--n-boot", "100", "--out", str(tmp_path / "r.html")])
    out = capsys.readouterr().out
    assert "FUNZIONE SPERIMENTALE" in out and "non entra nel report" in out
    assert "Frazione di CD8" not in (tmp_path / "r.html").read_text()


def test_vdj_manifest_paths_relative_to_manifest_and_missing_file(tmp_path):
    adata, contigs = make_tcr_validation_dataset(n_patients=3, n_clones_per_patient=3, seed=2)
    _, man = _tcr_files(tmp_path, adata, contigs)
    assert all(Path(f[0]).is_absolute() and Path(f[0]).exists() for f in read_vdj_manifest(man))
    pd.DataFrame([{"path": "manca.csv", "patient": "P", "compartment": "T"}]).to_csv(tmp_path / "m2.csv", index=False)
    with pytest.raises(ValueError, match="non trovato"):
        read_vdj_manifest(tmp_path / "m2.csv")


@pytest.mark.parametrize("argv,msg", [
    (["--technical", "batch=a", "--technical", "batch=b"], "indicato due volte"),
    (["--compare", "batch:x:y"], "non è dichiarata"),
    (["--compare", "tissue:Tumor"], "colonna:livello_a:livello_b"),
])
def test_cli_design_input_errors_are_clear(tmp_path, argv, msg):
    (tmp_path / "s.csv").write_text("patient,tissue,a,b\nP1,Tumor,x,y\nP1,Normal,x,y\n")
    with pytest.raises(SystemExit) as e:
        cli_main(["design", "--meta", str(tmp_path / "s.csv"), "--patient-col", "patient",
                  "--tissue-col", "tissue", *argv, "--out", str(tmp_path / "r.html")])
    assert msg in str(e.value)


def test_cli_design_h5ad_with_repeated_columns_and_library_note(tmp_path, capsys):
    obs = pd.DataFrame({"patient": ["P1"] * 4 + ["P2"] * 4, "tissue": ["T", "N"] * 4})
    a = ad.AnnData(X=sp.csr_matrix(np.ones((8, 2), dtype=np.float32)), obs=obs)
    a.write_h5ad(tmp_path / "x.h5ad")
    cli_main(["design", "--h5ad", str(tmp_path / "x.h5ad"), "--patient-col", "patient", "--tissue-col", "tissue",
              "--outcome-col", "tissue", "--out", str(tmp_path / "r.html")])
    assert LIBRARY_COLLAPSE_NOTE in capsys.readouterr().out


# --------------------------------------------------------------------------- #
# A9: tabelle per cellula date come foglio campioni
# --------------------------------------------------------------------------- #
def test_cell_level_table_is_refused_as_sample_sheet():
    obs = pd.DataFrame({"patient": np.repeat([f"P{i}" for i in range(6)], 200),
                        "tissue": np.tile(["T", "N"], 600)})
    with pytest.raises(ValueError, match="foglio campioni"):
        run_design_audit(obs, "patient", "tissue")
    sheet = sample_sheet_from_obs(obs, ["patient", "tissue", "tissue"])
    assert len(sheet) == 12 and list(sheet.columns) == ["patient", "tissue"]


def test_few_duplicate_rows_are_declared():
    rows = [{"patient": f"P{p}", "tissue": t} for p in range(6) for t in ("T", "N")] + [{"patient": "P0", "tissue": "T"}]
    res = run_design_audit(pd.DataFrame(rows), "patient", "tissue")
    assert any("ripetono la stessa combinazione" in n for n in res.notes)


def test_low_power_text_does_not_cite_wilcoxon_floor_when_below_alpha():
    rows = [{"patient": f"P{p}", "tissue": t} for p in range(6) for t in ("T", "N")]
    c = run_design_audit(pd.DataFrame(rows), "patient", "tissue", comparisons=[("tissue", "T", "N")],
                         min_units=8).comparisons[0]
    assert c.classification == "stimabile con bassa potenza"
    assert "meno di 8 pazienti appaiati" in c.sentence and "non può scendere" not in c.sentence


def test_app_comparison_lines_report_malformed_rows():
    comps, bad = parse_comparison_lines("tissue:Tumor:Normal\n\nprotocol:sc\nx:y:z:w\n")
    assert comps == [("tissue", "Tumor", "Normal")]
    assert bad == ["protocol:sc", "x:y:z:w"]


# --------------------------------------------------------------------------- #
# A10: l'app non scrive copie dei dati su disco
# --------------------------------------------------------------------------- #
def test_app_reads_uploads_in_memory_without_temporary_files():
    src = (ROOT / "app.py").read_text(encoding="utf-8")
    assert "tempfile" not in src and "NamedTemporaryFile" not in src
    assert "io.BytesIO(up.getvalue())" in src
    a = ad.read_h5ad(__import__("io").BytesIO((ROOT / "data" / "demo" / "GSE125449_demo.h5ad").read_bytes()))
    assert a.n_obs == 1861


# --------------------------------------------------------------------------- #
# Revisione, giro 1: errori d'input del Modulo B che facevano sparire il tasso d'errore
# --------------------------------------------------------------------------- #
def test_reference_compartment_not_in_data_is_an_error():
    adata, contigs = make_tcr_validation_dataset(n_patients=5, n_clones_per_patient=5, seed=3)
    with pytest.raises(ValueError, match="'Blood' non compare.*PBMC"):
        _run(adata, contigs, marker_map=MARKERS, reference_compartment="Blood")


@pytest.mark.parametrize("bad, match", [
    ({"CD4T": ["ENSG00000010610"], "CD8T": ["CD8A"]}, "nessun gene marcatore di 'CD4T'"),
    ({"CD4T": "CD4", "CD8T": ["CD8A"]}, "almeno un gene"),
    (["CD4", "CD8A"], "almeno un gene"),
    ({"CD4T": [], "CD8T": ["CD8A"]}, "almeno un gene"),
])
def test_malformed_marker_map_is_an_error_not_a_vanished_rate(bad, match):
    adata, contigs = make_tcr_validation_dataset(n_patients=5, n_clones_per_patient=5, seed=3)
    with pytest.raises(ValueError, match=match):
        _run(adata, contigs, marker_map=bad, reference_compartment="PBMC")
