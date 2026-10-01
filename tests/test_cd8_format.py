"""B1 -- formato tabellare dell'output della frazione di CD8.

- Regressione: i numeri (matrice, J, frazioni, intervalli, motivi di rifiuto) sono identici
  (tolleranza 1e-12) a quelli prodotti PRIMA della modifica di formato
  (tests/fixtures/cd8_regression_baseline.json, generato sul commit 2f03461).
- Formato: in CLI e nel report HTML avvertenza e assunzioni compaiono una sola volta, c'e'
  una riga per paziente e i rifiuti rimandano a note numerate.
"""

from __future__ import annotations

import json

import pytest

from core.cd8_propagation import ASSUMPTIONS_TEXT, EXPERIMENTAL_NOTE, SCENARIO_WARNING, format_cd8_text
from core.report import render_report
from tests.fixtures.cd8_regression import BASELINE_PATH, CONFIGS, run_config, summarize
from tests.test_tcr_flags import _assert_close


@pytest.mark.parametrize("name", list(CONFIGS))
def test_numbers_identical_to_pre_b1_baseline(name):
    baseline = json.loads(BASELINE_PATH.read_text())[name]
    current = json.loads(json.dumps(summarize(run_config(name)), ensure_ascii=False))
    _assert_close(baseline, current, name)


def test_cli_text_is_a_table_with_warnings_once():
    res = run_config("errore_alto_2x_rifiutato")
    txt = format_cd8_text(res)
    print("\n" + txt)
    assert txt.count(SCENARIO_WARNING) == 1
    assert txt.count(ASSUMPTIONS_TEXT) == 1
    assert txt.count("Nessuno scenario") == 1
    for p in res.patients:
        assert sum(line.startswith(p.patient + " ") for line in txt.splitlines()) == 1
    assert "non prodotto [1]" in txt
    assert "[1] matrice di confusione mal condizionata" in txt


def test_html_section_states_warning_and_reasons_once():
    res = run_config("errore_alto_2x_rifiutato")
    page = render_report(cd8_result=res, dataset_name="test")
    assert page.count("Nessuno scenario") == 1
    assert page.count("matrice di confusione mal condizionata") == 1
    assert page.count("<tr><td>P0") == len(res.patients)


def test_experimental_note_at_top_of_cli_and_html():
    res = run_config("simmetrico")
    txt = format_cd8_text(res)
    assert txt.splitlines()[1] == EXPERIMENTAL_NOTE and txt.count(EXPERIMENTAL_NOTE) == 1
    page = render_report(cd8_result=res, dataset_name="test")
    section = page[page.index("Frazione di CD8 nel compartimento"):]
    assert section.index("FUNZIONE SPERIMENTALE") < section.index("<table>")
