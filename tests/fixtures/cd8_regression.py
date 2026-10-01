"""Configurazioni e serializzazione per il test di regressione di cd8_fraction_intervals.

``cd8_regression_baseline.json`` e' stato generato con ``python -m
tests.fixtures.cd8_regression`` sul codice del commit 2f03461, PRIMA della modifica di
formato B1 (tabella al posto del testo ripetuto). Non va rigenerato dopo: serve a
verificare che i numeri non cambino.
"""

from __future__ import annotations

import json
import math
from pathlib import Path

from core.cd8_propagation import cd8_fraction_intervals
from core.synthetic import make_cd8_fraction_dataset

BASELINE_PATH = Path(__file__).with_name("cd8_regression_baseline.json")

CONFIGS = {
    "simmetrico": dict(p_cd4_to_cd8=0.15, p_cd8_to_cd4=0.15, p_to_other=0.05, seed=11),
    "asimmetrico": dict(p_cd4_to_cd8=0.24, p_cd8_to_cd4=0.08, p_to_other=0.05, seed=12),
    "errore_alto_2x_rifiutato": dict(p_cd4_to_cd8=0.30, p_cd8_to_cd4=0.25, p_to_other=0.05, seed=13),
}


def _num(x):
    return None if x is None or (isinstance(x, float) and math.isnan(x)) else float(x)


def run_config(name: str):
    obs, _ = make_cd8_fraction_dataset(n_patients=8, **CONFIGS[name])
    return cd8_fraction_intervals(obs, "patient", "compartment", "celltype",
                                  "audit_reference_label", "Tumor", n_boot=300, seed=5)


def summarize(res) -> dict:
    return {
        "J": _num(res.youden_j), "n_ref": res.n_reference_cells, "refused": res.refused_reason,
        "matrix": None if res.matrix is None else [[float(v) for v in r] for r in res.matrix.values],
        "patients": [
            {"patient": p.patient, "n4": p.n_cd4_called, "n8": p.n_cd8_called,
             "reported": _num(p.reported), "naive": [_num(p.naive_low), _num(p.naive_high)],
             "scenarios": {f"{k:g}": [_num(iv.low), _num(iv.high), iv.refused_reason]
                           for k, iv in sorted(p.scenarios.items())}}
            for p in res.patients],
    }


if __name__ == "__main__":
    BASELINE_PATH.write_text(json.dumps({n: summarize(run_config(n)) for n in CONFIGS},
                                        indent=1, ensure_ascii=False))
    print(f"scritto {BASELINE_PATH}")
