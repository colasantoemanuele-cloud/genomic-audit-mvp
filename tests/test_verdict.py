"""Regole della «Sintesi dei controlli» (core/verdict.py), con oggetti risultato veri.

La sintesi descrive i singoli controlli, non lo studio. Ogni regola ha un test: il verde non
compare mai per un controllo con un fatto strutturale, non eseguito, non valutabile, rifiutato o
fuori dal regime calibrato; non esiste un colore complessivo; i testi non contengono giudizi.
I risultati sono quelli prodotti dai moduli (run_design_audit) oppure le stesse dataclass dei
moduli costruite con i valori del caso da provare: nessun oggetto surrogato.
"""

from __future__ import annotations

import re

import numpy as np
import pandas as pd
import pytest

from core.design_audit import MISSING_LEVEL_CLASS, run_design_audit
from core.leakage_audit import (
    GAP_ALERT,
    LeakageAuditResult,
    ModelComparisonResult,
    PairwiseComparison,
    SchemeSummary,
)
from core.report import render_markdown_report, render_report
from core.stats import ClusterInterval
from core.synthetic import make_gse278694_like_sheet
from core.tcr_validation import ConventionResult, DiscordanceResult, TcrValidationResult
from core.verdict import (
    GIALLO,
    GRIGIO,
    ROSSO,
    RULES,
    STATE_LABEL,
    VERDE,
    SectionSummary,
    design_summary,
    leakage_summary,
    standing_limits,
    tcr_summary,
)


def _states(summary: SectionSummary) -> list[str]:
    return [c.state for c in summary.checks]


def _one(summary: SectionSummary, name_part: str):
    return next(c for c in summary.checks if name_part in c.name)


# --------------------------------------------------------------------------- #
# Disegno: risultati veri di run_design_audit
# --------------------------------------------------------------------------- #
def _paired_sheet(n_patients: int, batch_levels: int = 3, seed: int = 0) -> pd.DataFrame:
    rng = np.random.default_rng(seed)
    rows = []
    for p in range(n_patients):
        for t in ("Tumor", "Normal"):
            rows.append({"patient": f"P{p:02d}", "tissue": t, "library": f"L{p}{t}",
                         "batch": f"b{rng.integers(0, batch_levels)}"})
    return pd.DataFrame(rows)


def test_design_estimable_low_power_not_estimable_and_missing_level():
    sheet = make_gse278694_like_sheet()
    d = run_design_audit(sheet, "patient", "tissue", {"protocol": "protocol"},
                         comparisons=[("tissue", "Tumor", "PBMC"), ("tissue", "Tumor", "Adjacent_normal"),
                                      ("protocol", "scRNA", "snRNA"), ("tissue", "Tumor", "Adjacent")])
    s = design_summary(d)
    assert _one(s, "Tumor vs PBMC").state == VERDE  # 10 pazienti appaiati
    assert _one(s, "Tumor vs Adjacent_normal").state == GIALLO  # 5 pazienti: bassa potenza
    assert _one(s, "scRNA vs snRNA").state == ROSSO  # nessun paziente con entrambi
    typo = next(c for c in s.checks if c.name == "Confronto tissue: Tumor vs Adjacent")
    assert d.comparisons[3].classification == MISSING_LEVEL_CLASS
    assert typo.state == GRIGIO and "errore d'input" in typo.text  # refuso: grigio, non rosso


def test_design_without_requested_comparisons_is_grey_not_green():
    s = design_summary(run_design_audit(_paired_sheet(8), "patient", "tissue", {"batch": "batch"}))
    assert _one(s, "Confronti").state == GRIGIO and "non eseguito" in _one(s, "Confronti").text


def test_nested_pair_with_strong_v_is_not_green():
    """Il caso della revisione: batch annidato nel tessuto, 20 righe, V = 0.94. Prima finiva fra le
    coppie «valutabili sotto soglia», cioè in verde."""
    rows = [{"patient": f"P{i}", "tissue": "T" if i < 10 else "N",
             "batch": ("B1", "B2")[i % 2] if i < 10 else ("B3", "B4")[i % 2]} for i in range(20)]
    d = run_design_audit(pd.DataFrame(rows), "patient", "tissue", {"batch": "batch"})
    pair = next(p for p in d.pairs if {p.factor_a, p.factor_b} == {"tissue", "batch"})
    assert pair.cramer_v is not None and pair.cramer_v >= 0.5 and not pair.v_alarm
    s = design_summary(d)
    assert not any(c.state == VERDE and "Associazion" in c.name for c in s.checks)
    nest = _one(s, "batch / tissue")  # batch annidato nel tessuto: esito determinato dal batch
    assert nest.state == GIALLO


@pytest.mark.parametrize("kind,state", [("uno-a-uno", ROSSO), ("annidamento", GIALLO),
                                        ("esito-determinato", GIALLO), ("unità-tecnica", GIALLO)])
def test_every_structural_fact_enters_the_summary(kind, state):
    sheets = {
        "uno-a-uno": pd.DataFrame([{"patient": f"P{p}", "tissue": t, "batch": "r1" if t == "Tumor" else "r2"}
                                   for p in range(6) for t in ("Tumor", "Normal")]),
        "annidamento": pd.DataFrame([{"patient": f"P{p}", "tissue": t, "batch": f"B{p // 2}"}
                                     for p in range(8) for t in ("Tumor", "Normal")]),
        "esito-determinato": pd.DataFrame([{"patient": f"P{p}", "tissue": "Tumor" if p < 4 else "Normal",
                                            "batch": f"B{p % 3}", "rep": r}
                                           for p in range(8) for r in range(2)]),
        "unità-tecnica": make_gse278694_like_sheet(),
    }
    tech = {"unità-tecnica": {"library": "library", "protocol": "protocol"}}.get(kind, {"batch": "batch"})
    if kind == "esito-determinato":
        tech = {"batch": "batch", "date": "rep"}
    d = run_design_audit(sheets[kind], "patient", "tissue", tech)
    found = [f for f in d.findings if f.kind == kind]
    assert found, [f.kind for f in d.findings]
    checks = [c for c in design_summary(d).checks if c.text in {f.sentence for f in found}]
    assert len(checks) == len(found) and all(c.state == state for c in checks)


def test_strong_association_without_structure_is_yellow_and_weak_is_green():
    rng = np.random.default_rng(3)
    rows = []
    for i in range(60):
        a = int(rng.integers(0, 3))
        b = a if rng.random() < 0.9 else int(rng.integers(0, 3))
        rows.append({"patient": f"P{i}", "batch": f"b{a}", "chemistry": f"c{b}", "date": f"d{rng.integers(0, 3)}"})
    d = run_design_audit(pd.DataFrame(rows), "patient",
                         technical_cols={"batch": "batch", "chemistry": "chemistry", "date": "date"})
    s = design_summary(d)
    assert _one(s, "Associazione batch × chemistry").state == GIALLO
    assert _one(s, "Associazioni fra fattori").state == VERDE


def test_not_evaluable_cramer_v_is_grey():
    d = run_design_audit(_paired_sheet(4), "patient", "tissue", {"batch": "batch"},
                         comparisons=[("tissue", "Tumor", "Normal")])
    s = design_summary(d)
    grey = [c for c in s.checks if c.name == "Associazioni fra fattori"]
    assert grey and all(c.state == GRIGIO for c in grey)


# --------------------------------------------------------------------------- #
# Modulo A: dataclass vere di leakage_audit
# --------------------------------------------------------------------------- #
def _scheme(mean=0.8, absent=0, capped=False, nc=0, k=5) -> SchemeSummary:
    return SchemeSummary(scheme="x", fold_scores=[mean + 0.01 * (i - 2) for i in range(k)],
                         n_train=[100] * k, n_test=[25] * k,
                         absent_classes=[["A"] if i < absent else [] for i in range(k)],
                         n_train_before_cap=[200 if capped else 100] * k, top_genes=[],
                         not_converged=[i < nc for i in range(k)])


def _leakage(gap=0.01, model_comparison=None, n_patients=10, **grouped) -> LeakageAuditResult:
    g = _scheme(**grouped)
    r = _scheme(mean=g.mean + gap)
    return LeakageAuditResult(task="t", n_cells=500, n_patients=n_patients, n_classes=2, n_folds_used=5,
                              grouped=g, random=r, gap=r.mean - g.mean, std_ratio=1.0, narrative="",
                              model_comparison=model_comparison)


def _mc() -> ModelComparisonResult:
    ref = _scheme(k=10)
    other = _scheme(mean=0.82, k=10)
    comp = PairwiseComparison("random_forest", other.mean, other.mean - ref.mean, 10, 0, 0,
                              0.01, 0.02, 0.01, 0.02, "esatto", 2 / 2 ** 10, False, "")
    return ModelComparisonResult("fold_lopo", "logreg", "random_forest",
                                 {"logreg": ref, "random_forest": other}, [comp])


def test_leakage_model_comparison_not_run_is_grey_and_run_is_yellow_descriptive():
    c = _one(leakage_summary(_leakage(model_comparison=None)), "Confronto fra modelli")
    assert c.state == GRIGIO and "non eseguito" in c.text
    run = _one(leakage_summary(_leakage(model_comparison=_mc())), "Confronto fra modelli")
    assert run.state == GIALLO and "descrittivo" in run.text  # nessun test calibrato: mai verde


@pytest.mark.parametrize("kw", [{"absent": 2}, {"capped": True}, {"nc": 1}])
def test_leakage_per_patient_estimate_with_limits_is_yellow(kw):
    assert _one(leakage_summary(_leakage(**kw)), "split per paziente").state == GIALLO
    assert _one(leakage_summary(_leakage()), "split per paziente").state == VERDE


def test_leakage_gap_threshold_drives_random_split_state():
    assert _one(leakage_summary(_leakage(gap=GAP_ALERT + 0.001)), "split casuale").state == ROSSO
    assert _one(leakage_summary(_leakage(gap=GAP_ALERT - 0.001)), "split casuale").state == VERDE


def test_leakage_negative_gap_beyond_threshold_is_yellow_not_green():
    """Revisione, giro 1: un divario di -0,108 usciva verde «entro la soglia»."""
    c = _one(leakage_summary(_leakage(gap=-0.108)), "split casuale")
    assert c.state == GIALLO and "più bassa" in c.text and "entro la soglia" not in c.text
    assert _one(leakage_summary(_leakage(gap=-GAP_ALERT + 0.001)), "split casuale").state == VERDE
    html_txt = render_report(_leakage(gap=-0.108), None, "t")
    assert "divario oltre la soglia" in html_txt and "divario entro la soglia" not in html_txt


def test_leakage_with_fewer_than_five_patients_is_never_green():
    s = leakage_summary(_leakage(gap=0.0, n_patients=2))
    assert _one(s, "split per paziente").state == GIALLO and "2 pazienti" in _one(s, "split per paziente").text
    assert _one(s, "split casuale").state == GIALLO


def _partial_sheet(n_pat):
    """Batch associato al tessuto senza separarlo: b1 nel tumore e b3 nel normale per quasi tutti
    i pazienti, b2 in entrambi per gli ultimi 2 (scenario del revisore, giro 1)."""
    rows = []
    for i in range(n_pat):
        late = i >= n_pat - 2
        rows += [(f"P{i}", "Tumor", "b2" if late else "b1"), (f"P{i}", "Normal", "b2" if late else "b3")]
    return pd.DataFrame(rows, columns=["patient", "tissue", "batch"])


@pytest.mark.parametrize("n_pat, motivo", [(12, "fortemente associato"), (7, "non valutabile")])
def test_design_comparison_with_partial_confounder_is_yellow_not_green(n_pat, motivo):
    """Revisione, giro 1: con un batch associato al tessuto (V = 0,78) il confronto restava verde."""
    r = run_design_audit(_partial_sheet(n_pat), "patient", "tissue", technical_cols={"batch": "batch"},
                         comparisons=[("tissue", "Tumor", "Normal")])
    assert r.comparisons[0].classification == "stimabile"
    c = _one(design_summary(r), "Confronto tissue")
    assert c.state == GIALLO and motivo in c.text and "'batch'" in c.text


# --------------------------------------------------------------------------- #
# Modulo B: dataclass vere di tcr_validation e stats
# --------------------------------------------------------------------------- #
def _ci(mean, lo, hi, n=10, ok=True, reason=""):
    return ClusterInterval(mean, lo if ok else float("nan"), hi if ok else float("nan"), n, 100, ok,
                           "t cluster-robusto CR2, gdl Bell-McCaffrey", 5.0, reason)


def _tcr(ci=(-0.01, 0.03), sufficient=True, n_patients=10, equal=(-0.01, 0.02), conventions=None):
    eq = _ci(0.005, *equal) if equal is not None else None
    d = DiscordanceResult(n_pairs=200, n_patients=n_patients, mean_excess=0.01,
                          ci_low=ci[0] if sufficient else float("nan"),
                          ci_high=ci[1] if sufficient else float("nan"), sufficient=sufficient,
                          by_compartment_pair=pd.DataFrame(), equal_weights=eq)
    return TcrValidationResult(n_cells_with_tcr=1000, n_clones_total=100, discordance=d, marker_error=None,
                               narrative="", conventions=conventions)


def test_tcr_insufficient_patients_is_red():
    c = _one(tcr_summary(_tcr(sufficient=False, n_patients=3)), "discordanza")
    assert c.state == ROSSO and "3 pazienti" in c.text


def test_tcr_interval_above_or_below_zero_is_yellow_never_green():
    assert _one(tcr_summary(_tcr(ci=(0.02, 0.09))), "discordanza").state == GIALLO
    below = _one(tcr_summary(_tcr(ci=(-0.09, -0.02))), "discordanza")
    assert below.state == GIALLO and "sotto lo zero" in below.text


def test_tcr_interval_including_zero_says_no_excess_detected():
    c = _one(tcr_summary(_tcr(ci=(-0.01, 0.03))), "discordanza")
    assert c.state == VERDE and "nessun eccesso rilevato" in c.text and "affidabile" not in c.text


def test_tcr_equal_weights_above_zero_makes_it_yellow():
    c = _one(tcr_summary(_tcr(ci=(-0.01, 0.03), equal=(0.01, 0.05))), "discordanza")
    assert c.state == GIALLO and "pesi uguali" in c.text


def test_tcr_error_rate_grey_red_green():
    assert _one(tcr_summary(_tcr(conventions=None)), "Tasso d'errore").state == GRIGIO
    conv = {"cell": ConventionResult("cell", "per cellula", {
        "Tumor": _ci(0.1, 0.05, 0.15), "Adjacent": _ci(0.04, 0, 0, n=4, ok=False, reason="4 pazienti, ne servono almeno 5")},
        {})}
    s = tcr_summary(_tcr(conventions=conv))
    assert _one(s, "in Tumor").state == VERDE
    adj = _one(s, "in Adjacent")
    assert adj.state == ROSSO and "4 pazienti" in adj.text


def test_tcr_rate_refused_for_zero_variance_is_yellow_and_empty_convention_grey():
    """Con 8 pazienti e tasso 0 la stima esiste: giallo, non «stima non possibile». Una
    convenzione senza unità valutabili compare in grigio invece di sparire."""
    conv = {"clone": ConventionResult("clone", "per clone", {
        "Tumor": _ci(0.0, 0, 0, n=8, ok=False, reason="tasso pari a 0 o a 1: intervallo su scala logit non definito")},
        {}), "cell": ConventionResult("cell", "per cellula", {}, {})}
    s = tcr_summary(_tcr(conventions=conv))
    t = _one(s, "[clone] in Tumor")
    assert t.state == GIALLO and "stima 0.000" in t.text
    assert _one(s, "[cell]").state == GRIGIO


def test_tcr_discordance_refused_with_enough_patients_is_yellow():
    d = DiscordanceResult(n_pairs=50, n_patients=8, mean_excess=0.0, ci_low=float("nan"), ci_high=float("nan"),
                          sufficient=False, by_compartment_pair=pd.DataFrame(), equal_weights=None,
                          reason="variabilità nulla fra pazienti")
    r = TcrValidationResult(n_cells_with_tcr=100, n_clones_total=10, discordance=d, marker_error=None, narrative="")
    c = _one(tcr_summary(r), "discordanza")
    assert c.state == GIALLO and "variabilità nulla" in c.text


# --------------------------------------------------------------------------- #
# Nessun colore complessivo, nessun giudizio
# --------------------------------------------------------------------------- #
def test_no_overall_color_only_counts_per_state():
    d = run_design_audit(make_gse278694_like_sheet(), "patient", "tissue", {"protocol": "protocol"},
                         comparisons=[("protocol", "scRNA", "snRNA"), ("tissue", "Tumor", "PBMC")])
    s = design_summary(d)
    assert not hasattr(s, "color")
    assert sum(s.counts.values()) == len(s.checks)
    assert s.counts_text.count("×") == len(s.counts)


JUDGEMENTS = re.compile(r"errat|sbagliat|scars|inaffidabil|pessim|cattiv|mediocr|scorrett|verdett|"
                        r"errore sistematico|bocciat|promoss|affidabile|validat|dimostra|effetto reale",
                        re.IGNORECASE)


def _all_texts() -> list[str]:
    d = run_design_audit(make_gse278694_like_sheet(), "patient", "tissue", {"protocol": "protocol"},
                         comparisons=[("protocol", "scRNA", "snRNA"), ("tissue", "Tumor", "Adjacent_normal"),
                                      ("tissue", "Tumor", "X")])
    sections = [design_summary(d), leakage_summary(_leakage(gap=0.3, absent=1, capped=True, nc=2)),
                leakage_summary(_leakage(model_comparison=_mc())), tcr_summary(_tcr(sufficient=False)),
                tcr_summary(_tcr(ci=(0.02, 0.09))), tcr_summary(_tcr(ci=(-0.09, -0.02)))]
    texts = [t for s in sections for c in s.checks for t in (c.name, c.text, c.label)]
    return texts + list(RULES.values()) + list(STATE_LABEL.values()) + standing_limits()


def test_texts_describe_what_data_allow_without_judgements():
    for t in _all_texts():
        if "Dati sintetici" in t:
            continue
        found = JUDGEMENTS.search(t)
        # "non è un dispositivo medico" e "verificata/verifica" sono ammessi; "validat" no
        assert not found, t


def test_state_labels_describe_the_check():
    assert STATE_LABEL == {
        VERDE: "nessun limite rilevato dalle regole", GIALLO: "stima con limiti",
        ROSSO: "stima non possibile con questi dati", GRIGIO: "controllo non eseguito o non valutabile"}


def test_reports_use_neutral_title_one_row_per_check_and_print_limits():
    d = run_design_audit(make_gse278694_like_sheet(), "patient", "tissue", {"protocol": "protocol"},
                         comparisons=[("protocol", "scRNA", "snRNA"), ("tissue", "Tumor", "Adjacent_normal")])
    md = render_markdown_report(design_result=d, dataset_name="prova")
    page = render_report(design_result=d, dataset_name="prova")
    for text in (md, page):
        assert "Sintesi dei controlli" in text and "erdetto" not in text
        assert "stima non possibile con questi dati" in text and "stima con limiti" in text
        assert "Non è un dispositivo medico" in text and "Limiti dichiarati" in text
    assert "| Disegno | Confronto protocol: scRNA vs snRNA | rosso: stima non possibile con questi dati |" in md
    assert 'class="tag tag-rosso"' in page and 'class="tag tag-giallo"' in page
    for limit in standing_limits():
        assert limit in md
