"""Margine dell'identità di riferimento del clone (Modulo B, METODI §4).

Regola documentata: identità = etichetta con la frazione di cellule positive più alta se il margine
sulla seconda è ALMENO 0,20. Prima della correzione il confronto era in virgola mobile
(``best - second < 0.20``): con margine esattamente 0,20 la decisione dipendeva dall'arrotondamento
(3/5 contro 2/5 e 7/10 contro 5/10 rifiutati, 1/5 contro 0/5 e 4/10 contro 2/10 accettati).
Difetto trovato dall'implementazione indipendente della validazione su GSE200996
(validation_esterna_B/GSE200996/REPORT.md).
"""

from __future__ import annotations

import pandas as pd
import pytest

from core.tcr_validation import assign_reference_identity


def _clone(n_pos_cd4: int, n_pos_cd8: int, n: int) -> pd.DataFrame:
    """Un clone nel sangue con n cellule: le prime n_pos_cd4 positive a CD4, le ultime n_pos_cd8 a CD8."""
    return pd.DataFrame({"patient": "P1", "clone_id": "c1", "compartment": "PBMC",
                         "pos_CD4T": [i < n_pos_cd4 for i in range(n)],
                         "pos_CD8T": [i >= n - n_pos_cd8 for i in range(n)]})


@pytest.mark.parametrize("pos4,pos8,n", [(3, 2, 5), (1, 0, 5), (7, 5, 10), (3, 1, 10), (6, 4, 10),
                                         (8, 6, 10), (4, 2, 10), (9, 6, 15)])
def test_margin_exactly_020_is_accepted(pos4, pos8, n):
    ref = assign_reference_identity(_clone(pos4, pos8, n), ["CD4T", "CD8T"])
    assert len(ref) == 1 and ref["reference_identity"].iloc[0] == "CD4T"


@pytest.mark.parametrize("pos4,pos8,n", [(5, 4, 10), (2, 1, 6), (11, 8, 16)])
def test_margin_below_020_is_rejected(pos4, pos8, n):
    assert assign_reference_identity(_clone(pos4, pos8, n), ["CD4T", "CD8T"]).empty
