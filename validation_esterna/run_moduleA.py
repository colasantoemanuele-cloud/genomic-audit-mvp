"""Modulo A sui dati esterni: punteggi fold per fold dello strumento (A2, A4c), controllo
negativo con 20 permutazioni (A3) e varianti di formato (A4a, A4b).

Usa core.leakage_audit.run_leakage_audit, la stessa funzione che chiama `cli.py leakage`,
con gli stessi argomenti (A1 e' verificato a parte con la CLI vera).

Uso: python -m validation_esterna.run_moduleA <passo> [h5ad] [suffisso]
     passo in {base, perm, formati}; default h5ad = data/GSE125449.h5ad
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse as sp

from core.leakage_audit import run_leakage_audit

H5AD = Path(sys.argv[2]) if len(sys.argv) > 2 else Path(__file__).parent / "data" / "GSE125449.h5ad"
SUFFIX = sys.argv[3] if len(sys.argv) > 3 else ""
RES = Path(__file__).parent / "results"
TARGET, PATIENT = "Type", "patient"
NO_LOGO = 10 ** 6  # min_patients_for_model_comparison: salta il confronto LOGO


def _scores(r) -> dict:
    out = {"grouped": r.grouped.fold_scores, "random": r.random.fold_scores,
           "grouped_mean": r.grouped.mean, "random_mean": r.random.mean, "gap": r.gap,
           "tag_leakage_rilevabile": bool(r.gap > 0.05), "narrative": r.narrative}
    if r.model_comparison is not None:
        out["logo"] = {m: s.fold_scores for m, s in r.model_comparison.scores.items()}
        out["logo_mean"] = {m: s.mean for m, s in r.model_comparison.scores.items()}
        out["best_model"] = r.model_comparison.best_model
        out["comparisons"] = [vars(c) for c in r.model_comparison.comparisons]
    return out


def base() -> None:
    adata = ad.read_h5ad(H5AD)
    r = run_leakage_audit(adata, target_col=TARGET, patient_col=PATIENT, seed=0)
    (RES / f"moduleA_base{SUFFIX}.json").write_text(json.dumps(_scores(r), indent=1, ensure_ascii=False, default=float))
    print(r.narrative)


def _perm(seed: int) -> dict:
    adata = ad.read_h5ad(H5AD)
    y = adata.obs[TARGET].astype(str).values
    adata.obs["Type_perm"] = np.random.default_rng(seed).permutation(y)
    r = run_leakage_audit(adata, target_col="Type_perm", patient_col=PATIENT,
                          min_patients_for_model_comparison=NO_LOGO, seed=0)
    s = _scores(r)
    s.pop("narrative")
    return {"seed": seed, **s}


def perm() -> None:
    # sequenziale: ogni esecuzione usa gia' tutti i core per i fold (nessuna sovrascrittura di thread)
    rows = [_perm(seed) for seed in range(20)]
    y = ad.read_h5ad(H5AD, backed="r").obs[TARGET].astype(str)
    k = y.nunique()
    gm = np.array([r["grouped_mean"] for r in rows])
    n_tag = sum(r["tag_leakage_rilevabile"] for r in rows)
    summary = {"K": k, "atteso_1_su_K": 1 / k, "media_raggruppata": float(gm.mean()),
               "min_max_raggruppata": [float(gm.min()), float(gm.max())],
               "media_divario": float(np.mean([r["gap"] for r in rows])),
               "permutazioni_con_tag_leakage": int(n_tag),
               "A3_i": bool(abs(gm.mean() - 1 / k) <= 0.05), "A3_ii": bool(n_tag <= 1),
               "permutazioni": rows}
    (RES / f"moduleA_perm{SUFFIX}.json").write_text(json.dumps(summary, indent=1, default=float))
    print(json.dumps({k: v for k, v in summary.items() if k != "permutazioni"}, indent=1))


def formati() -> None:
    adata = ad.read_h5ad(H5AD)
    run = lambda a: run_leakage_audit(a, target_col=TARGET, patient_col=PATIENT,  # noqa: E731
                                      min_patients_for_model_comparison=NO_LOGO, seed=0)
    ref = run(adata)
    out = {"obs_dtype_originale": {c: str(adata.obs[c].dtype) for c in (TARGET, PATIENT)},
           "X_originale_sparsa": bool(sp.issparse(adata.X))}
    dense = adata.copy()
    dense.X = np.asarray(dense.X.toarray())
    rd = run(dense)
    out["A4a_densa_identica"] = (rd.grouped.fold_scores == ref.grouped.fold_scores
                                 and rd.random.fold_scores == ref.random.fold_scores)
    strs = adata.copy()
    strs.obs[TARGET] = strs.obs[TARGET].astype(str).astype(object)
    strs.obs[PATIENT] = strs.obs[PATIENT].astype(str).astype(object)
    cats = adata.copy()
    cats.obs[TARGET] = pd.Categorical(cats.obs[TARGET].astype(str))
    cats.obs[PATIENT] = pd.Categorical(cats.obs[PATIENT].astype(str))
    rs, rc = run(strs), run(cats)
    out["A4b_stringhe_identica"] = rs.grouped.fold_scores == ref.grouped.fold_scores and rs.random.fold_scores == ref.random.fold_scores
    out["A4b_categoriche_identica"] = rc.grouped.fold_scores == ref.grouped.fold_scores and rc.random.fold_scores == ref.random.fold_scores
    out["fold_raggruppati"] = {"riferimento": ref.grouped.fold_scores, "densa": rd.grouped.fold_scores,
                               "stringhe": rs.grouped.fold_scores, "categoriche": rc.grouped.fold_scores}
    (RES / f"moduleA_formati{SUFFIX}.json").write_text(json.dumps(out, indent=1, default=float))
    print(json.dumps({k: v for k, v in out.items() if k != "fold_raggruppati"}, indent=1))


if __name__ == "__main__":
    {"base": base, "perm": perm, "formati": formati}[sys.argv[1]]()
