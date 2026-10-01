"""B3 -- robustezza della frazione di CD8 a un errore NON rappresentativo.

Solo una misura, nessun cambio di metodo e nessuna banda nuova. Nella simulazione l'errore di
annotazione delle cellule CON identita' di riferimento (i cloni condivisi con il sangue) e'
``f`` volte quello delle cellule senza riferimento (f = 1 controllo, 2, 4). La matrice
stimata sulle cellule di riferimento sovrastima quindi l'errore delle altre. Due errori di
base (cellule senza riferimento): "basso" CD4->CD8 = CD8->CD4 = 0.05, altro 0.02; "alto"
asimmetrico come nel tumore reale, CD4->CD8 = 0.10, CD8->CD4 = 0.03, altro 0.04. 30% delle
cellule con riferimento; 200 repliche x 10 pazienti = 2000 intervalli per scenario.

Errore medio effettivo sulle cellule di un paziente: 0.3*f*e + 0.7*e. Il fattore di scala che
riporterebbe la matrice al valore giusto e' quindi (0.3 f + 0.7) / f: 1 per f=1, 0.65 per
f=2, 0.475 per f=4. Nessuno dei tre scenari (0.5x, 1x, 2x) coincide esattamente per f=2.
"""

from __future__ import annotations

import functools

import pytest

from core.cd8_propagation import cd8_fraction_intervals
from core.synthetic import make_cd8_fraction_dataset

R = 200
BASES = {"basso": dict(p_cd4_to_cd8=0.05, p_cd8_to_cd4=0.05, p_to_other=0.02),
         "alto": dict(p_cd4_to_cd8=0.10, p_cd8_to_cd4=0.03, p_to_other=0.04)}
NOMINAL_BAND = (0.90, 0.99)  # la banda dichiarata per la calibrazione, solo per confronto


@functools.lru_cache(maxsize=8)
def _coverage(base: str, factor: float) -> dict:
    cov = {k: 0 for k in (0.5, 1.0, 2.0)}
    n = {k: 0 for k in (0.5, 1.0, 2.0)}
    refused = {k: 0 for k in (0.5, 1.0, 2.0)}
    for rep in range(R):
        obs, truth = make_cd8_fraction_dataset(n_patients=10, seed=90_000 + rep,
                                               reference_error_factor=factor, **BASES[base])
        res = cd8_fraction_intervals(obs, "patient", "compartment", "celltype",
                                     "audit_reference_label", "Tumor", n_boot=400, seed=rep)
        for p in res.patients:
            for k, iv in p.scenarios.items():
                if iv.low is None:
                    refused[k] += 1
                    continue
                n[k] += 1
                cov[k] += iv.low <= truth[p.patient] <= iv.high
    return {"cov": cov, "n": n, "refused": refused}


@pytest.mark.parametrize("base", list(BASES))
@pytest.mark.parametrize("factor", [1.0, 2.0, 4.0])
def test_coverage_with_non_representative_error(base, factor):
    r = _coverage(base, factor)
    parts = []
    for k in (0.5, 1.0, 2.0):
        frac = r["cov"][k] / r["n"][k] if r["n"][k] else float("nan")
        parts.append(f"{k:g}x: {r['cov'][k]}/{r['n'][k]} = {frac:.3f} (rifiutati {r['refused'][k]})")
    print(f"\n[B3 {base}] errore cloni con riferimento = {factor:g} x altri -> copertura IC95% | "
          + " | ".join(parts) + f" | banda nominale {NOMINAL_BAND}")
    # Misura, non calibrazione: si verifica solo che ogni scenario sia stato valutato.
    assert all(r["n"][k] + r["refused"][k] == R * 10 for k in (0.5, 1.0, 2.0))
