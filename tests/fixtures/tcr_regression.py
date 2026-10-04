"""Configurazioni e serializzazione per il test di regressione di run_tcr_validation.

``tcr_regression_baseline.json`` è stato generato eseguendo ``python -m
tests.fixtures.tcr_regression`` sul codice del commit di baseline (cf61e76), PRIMA di
qualunque modifica a core/tcr_validation.py. Non va rigenerato dopo le modifiche: il suo
scopo è proprio verificare che l'output non cambi.
"""

from __future__ import annotations

import json
import math
from pathlib import Path

from core.synthetic import make_tcr_validation_dataset
from core.tcr_validation import run_tcr_validation

BASELINE_PATH = Path(__file__).with_name("tcr_regression_baseline.json")
MARKERS = {"CD4T": ["CD4"], "CD8T": ["CD8A", "CD8B"]}

CONFIGS = {
    "eccesso_035_marcatori": dict(data=dict(n_patients=12, n_clones_per_patient=15, injected_excess=0.35, seed=0),
                                  run=dict(n_boot=500, seed=0, marker_map=MARKERS, reference_compartment="PBMC")),
    "eccesso_0_marcatori": dict(data=dict(n_patients=12, n_clones_per_patient=15, injected_excess=0.0, seed=1),
                                run=dict(n_boot=500, seed=3, marker_map=MARKERS, reference_compartment="PBMC")),
    "tre_compartimenti": dict(data=dict(n_patients=8, n_clones_per_patient=10, injected_excess=0.25,
                                        compartments=("PBMC", "Adjacent", "Tumor"), seed=2),
                              run=dict(n_boot=400, seed=1, marker_map=MARKERS, reference_compartment="PBMC")),
    "senza_marcatori_pochi_pazienti": dict(data=dict(n_patients=4, n_clones_per_patient=10, injected_excess=0.35, seed=4),
                                           run=dict(n_boot=300, seed=0)),
}


def _num(x):
    return None if (x is None or (isinstance(x, float) and math.isnan(x))) else float(x)


def summarize(result) -> dict:
    d = result.discordance
    out = {
        "n_cells_with_tcr": result.n_cells_with_tcr,
        "n_clones_total": result.n_clones_total,
        "narrative": result.narrative,
        "discordance": {"n_pairs": d.n_pairs, "n_patients": d.n_patients, "mean_excess": _num(d.mean_excess),
                        "ci_low": _num(d.ci_low), "ci_high": _num(d.ci_high), "sufficient": bool(d.sufficient)},
        "by_compartment_pair": [
            {k: (_num(v) if isinstance(v, float) else (bool(v) if k == "sufficient" else v))
             for k, v in row.items()}
            for row in d.by_compartment_pair.to_dict(orient="records")
        ],
        "marker_error": None,
    }
    if result.marker_error is not None:
        m = result.marker_error
        out["marker_error"] = {
            "n_resolved_clones": m.n_resolved_clones,
            "by_compartment": {c: {"mean": _num(r.mean), "ci_low": _num(r.ci_low), "ci_high": _num(r.ci_high),
                                   "n_groups": r.n_groups, "n_obs": r.n_obs, "sufficient": bool(r.sufficient)}
                               for c, r in m.by_compartment.items()},
        }
    return out


def run_config(name: str):
    cfg = CONFIGS[name]
    adata, contigs = make_tcr_validation_dataset(**cfg["data"])
    return adata, contigs, run_tcr_validation(
        adata, contigs, patient_col="patient_id", compartment_col="tissue",
        celltype_col="celltype", barcode_col="barcode", **cfg["run"])


if __name__ == "__main__":
    baseline = {name: summarize(run_config(name)[2]) for name in CONFIGS}
    BASELINE_PATH.write_text(json.dumps(baseline, indent=1, ensure_ascii=False, default=int))
    print(f"scritto {BASELINE_PATH}")
