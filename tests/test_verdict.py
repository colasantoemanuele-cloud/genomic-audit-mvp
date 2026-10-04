"""Regole della «Sintesi dei controlli» (core/verdict.py).

La sintesi descrive i singoli controlli, non lo studio. Questi test fissano le regole di
visualizzazione: il verde non compare mai per un controllo non valutabile, non stimabile,
rifiutato o non eseguito; non esiste un colore complessivo; i testi non contengono giudizi.
I risultati dei moduli sono sostituiti da oggetti minimi con i soli campi letti dalle regole:
nessun calcolo statistico è coinvolto.
"""

from __future__ import annotations

import re
from types import SimpleNamespace as NS

import pytest

from core.report import render_markdown_report, render_report
from core.verdict import (
    GIALLO,
    GRIGIO,
    ROSSO,
    RULES,
    STATE_LABEL,
    VERDE,
    SectionSummary,
    cd8_summary,
    design_summary,
    leakage_summary,
    standing_limits,
    tcr_summary,
)


# --------------------------------------------------------------------------- #
# Oggetti minimi
# --------------------------------------------------------------------------- #
def _comparison(cls: str, n: int = 6):
    return NS(factor="tissue", level_a="Tumor", level_b="Normal", classification=cls, n_units=n)


def _pair(v, alarm=False):
    return NS(factor_a="batch", factor_b="tissue", cramer_v=v, v_alarm=alarm)


def _design(comparisons=(), pairs=()):
    return NS(comparisons=list(comparisons), pairs=list(pairs))


def _scheme(mean=0.8, absent=0, capped=False, nc=0):
    return NS(mean=mean, std=0.05, fold_scores=[mean] * 5, n_folds_with_absent_classes=absent,
              capped=capped, n_not_converged=nc)


def _leakage(gap=0.01, model_comparison=None, **grouped):
    return NS(gap=gap, grouped=_scheme(**grouped), random=_scheme(), model_comparison=model_comparison)


def _mc():
    return NS(best_model="logreg", scores={"logreg": NS(mean=0.8)})


def _boot(sufficient=True, n=10):
    return NS(mean=0.1, ci_low=0.05, ci_high=0.15, n_groups=n, sufficient=sufficient)


def _tcr(ci=(-0.01, 0.03), sufficient=True, conventions=None, n_patients=10):
    return NS(discordance=NS(mean_excess=0.01, ci_low=ci[0], ci_high=ci[1], sufficient=sufficient,
                             n_patients=n_patients), conventions=conventions)


def _conv(**by_comp):
    return {"cell": NS(name="cell", by_compartment=by_comp)}


def _cd8(refused=None, lows=(0.2, None)):
    patients = [NS(scenarios={1.0: NS(low=x, refused_reason=None if x is not None else "J = 0.10")})
                for x in lows]
    return NS(target_compartment="Tumor", refused_reason=refused, patients=patients)


def _states(summary: SectionSummary) -> list[str]:
    return [c.state for c in summary.checks]


def _one(summary: SectionSummary, name_part: str):
    return next(c for c in summary.checks if name_part in c.name)


# --------------------------------------------------------------------------- #
# Il verde non compare mai per un controllo non valutabile / non stimabile / rifiutato / non eseguito
# --------------------------------------------------------------------------- #
def test_design_not_estimable_comparison_is_red_with_its_own_text():
    c = _one(design_summary(_design([_comparison("non stimabile", 0)])), "Confronto")
    assert c.state == ROSSO and c.label == "stima non possibile con questi dati"
    assert "non stimabile" in c.text


def test_design_low_power_is_yellow_and_estimable_is_green():
    s = design_summary(_design([_comparison("stimabile con bassa potenza", 4), _comparison("stimabile", 12)]))
    assert _states(s) == [GIALLO, VERDE]


def test_design_without_requested_comparisons_is_grey_not_green():
    s = design_summary(_design())
    assert _states(s) == [GRIGIO] and "non eseguito" in s.checks[0].text


def test_design_not_evaluable_cramer_v_is_grey_not_green():
    s = design_summary(_design([_comparison("stimabile")], [_pair(None), _pair(None)]))
    assoc = _one(s, "Associazioni")
    assert assoc.state == GRIGIO and "non valutabile per 2 coppie" in assoc.text
    assert sum(c.state == VERDE for c in s.checks) == 1  # solo il confronto stimabile


def test_design_strong_association_is_yellow_and_evaluable_pairs_green():
    s = design_summary(_design([_comparison("stimabile")], [_pair(0.7, alarm=True), _pair(0.2), _pair(None)]))
    assert _one(s, "batch × tissue").state == GIALLO
    assert sorted(c.state for c in s.checks if c.name == "Associazioni fra fattori") == [GRIGIO, VERDE]


def test_leakage_model_comparison_not_run_is_grey_not_green():
    c = _one(leakage_summary(_leakage(model_comparison=None)), "Confronto fra modelli")
    assert c.state == GRIGIO and "non eseguito" in c.text
    assert _one(leakage_summary(_leakage(model_comparison=_mc())), "Confronto fra modelli").state == VERDE


@pytest.mark.parametrize("kw", [{"absent": 2}, {"capped": True}, {"nc": 1}])
def test_leakage_per_patient_estimate_with_limits_is_yellow(kw):
    assert _one(leakage_summary(_leakage(model_comparison=_mc(), **kw)), "split per paziente").state == GIALLO


def test_leakage_gap_threshold_drives_random_split_state():
    assert _one(leakage_summary(_leakage(gap=0.051)), "split casuale").state == ROSSO
    assert _one(leakage_summary(_leakage(gap=0.05)), "split casuale").state == VERDE  # soglia: > 0.05
    assert _states(leakage_summary(_leakage(gap=0.0, model_comparison=_mc()))) == [VERDE, VERDE, VERDE]


def test_tcr_insufficient_patients_is_red_never_green():
    c = _one(tcr_summary(_tcr(sufficient=False, n_patients=3)), "discordanza")
    assert c.state == ROSSO and "3 pazienti" in c.text


def test_tcr_interval_above_zero_is_yellow_and_including_zero_is_green():
    assert _one(tcr_summary(_tcr(ci=(0.02, 0.09))), "discordanza").state == GIALLO
    assert _one(tcr_summary(_tcr(ci=(-0.01, 0.03))), "discordanza").state == VERDE


def test_tcr_error_rate_not_computed_is_grey_and_without_interval_is_red():
    assert _one(tcr_summary(_tcr(conventions=None)), "Tasso d'errore").state == GRIGIO
    s = tcr_summary(_tcr(conventions=_conv(Tumor=_boot(True), Adjacent=_boot(False, n=4))))
    assert _one(s, "in Tumor").state == VERDE
    adj = _one(s, "in Adjacent")
    assert adj.state == ROSSO and "4 pazienti" in adj.text


def test_cd8_is_never_green_refused_is_grey():
    refused = cd8_summary(_cd8(refused="2 pazienti con riferimento, ne servono almeno 5"))
    assert _states(refused) == [GRIGIO] and "rifiutato" in refused.checks[0].text
    produced = cd8_summary(_cd8())
    assert _states(produced) == [GIALLO] and "per 1 su 2 pazienti" in produced.checks[0].text
    none_produced = cd8_summary(_cd8(lows=(None, None)))
    assert _states(none_produced) == [GRIGIO] and "J = 0.10" in none_produced.checks[0].text


# --------------------------------------------------------------------------- #
# Nessun colore complessivo, nessun giudizio
# --------------------------------------------------------------------------- #
def test_no_overall_color_only_counts_per_state():
    s = design_summary(_design([_comparison("non stimabile", 0), _comparison("stimabile")], [_pair(None)]))
    assert not hasattr(s, "color")
    assert s.counts == {VERDE: 1, ROSSO: 1, GRIGIO: 1}
    assert s.counts_text == ("1 × stima affidabile; 1 × stima non possibile con questi dati; "
                             "1 × controllo non eseguito o non valutabile")


JUDGEMENTS = re.compile(r"errat|sbagliat|scars|inaffidabil|pessim|cattiv|mediocr|scorrett|verdett|"
                        r"errore sistematico|bocciat|promoss", re.IGNORECASE)


def _all_texts() -> list[str]:
    sections = [
        design_summary(_design([_comparison("non stimabile", 0), _comparison("stimabile con bassa potenza", 3),
                                _comparison("stimabile")], [_pair(0.8, True), _pair(0.1), _pair(None)])),
        design_summary(_design()),
        leakage_summary(_leakage(gap=0.3, absent=1, capped=True, nc=2)),
        leakage_summary(_leakage(model_comparison=_mc())),
        tcr_summary(_tcr(sufficient=False)), tcr_summary(_tcr(ci=(0.02, 0.09))),
        tcr_summary(_tcr(conventions=_conv(Tumor=_boot(True), Adjacent=_boot(False)))),
        cd8_summary(_cd8()), cd8_summary(_cd8(refused="matrice di confusione mal condizionata")),
    ]
    texts = [t for s in sections for c in s.checks for t in (c.name, c.text, c.label)]
    return texts + list(RULES.values()) + list(STATE_LABEL.values()) + standing_limits(_cd8())


def test_texts_describe_what_data_allow_without_judgements():
    for t in _all_texts():
        assert not JUDGEMENTS.search(t), t


def test_state_labels_describe_the_check():
    assert STATE_LABEL == {
        VERDE: "stima affidabile", GIALLO: "stima con limiti",
        ROSSO: "stima non possibile con questi dati", GRIGIO: "controllo non eseguito o non valutabile"}


def test_reports_use_neutral_title_and_one_row_per_check():
    from core.design_audit import run_design_audit
    from core.synthetic import make_gse278694_like_sheet
    d = run_design_audit(make_gse278694_like_sheet(), "patient", "tissue", {"protocol": "protocol"},
                         comparisons=[("protocol", "scRNA", "snRNA"), ("tissue", "Tumor", "Adjacent_normal")])
    md = render_markdown_report(design_result=d, dataset_name="prova")
    page = render_report(design_result=d, dataset_name="prova")
    for text in (md, page):
        assert "Sintesi dei controlli" in text and "Verdetto" not in text and "verdetto" not in text
        assert "stima non possibile con questi dati" in text and "stima con limiti" in text
    assert "| Disegno | Confronto protocol: scRNA vs snRNA | rosso: stima non possibile con questi dati |" in md
    assert 'class="tag tag-rosso"' in page and 'class="tag tag-giallo"' in page
