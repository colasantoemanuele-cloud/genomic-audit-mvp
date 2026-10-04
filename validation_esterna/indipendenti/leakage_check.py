"""Verifica indipendente del Modulo A (criterio A2) -- NESSUN import da core/.

Riesegue con scikit-learn la valutazione raggruppata per paziente secondo le regole dichiarate
(CRITERI.md, sezione 6): StratifiedGroupKFold a min(5, n pazienti) fold con shuffle e seme 0;
CPM a 1e4 + log1p; 2000 geni a varianza piu' alta stimati sul solo training (parita' ->
indice piu' basso); StandardScaler senza centratura; LogisticRegression(C=1,
class_weight="balanced", max_iter=1000, seme 0). Riesegue anche lo split casuale
(StratifiedKFold, controllo negativo) e, con almeno 8 pazienti, LeaveOneGroupOut per la
regressione logistica.

Per ogni fold registra tre macro-F1:
- "classi_del_test": labels = classi presenti in y_true del fold (definizione dichiarata dal
  nuovo motore, sezione 6);
- "standard": f1_score(average="macro") con le etichette presenti nel fold (vere o predette),
  il default di scikit-learn;
- "tutte_le_classi": labels = tutte le classi del dataset, zero_division=0. Una classe assente
  dal fold e mai predetta conta come F1 = 0.
Il confronto fra le due mostra se una classe assente viene contata come zero.

Con `--pipeline originale` lo script applica invece la pipeline dichiarata nei criteri
ORIGINALI (CRITERI.md, sezione 3, commit 0b5a847): tutti i geni, senza selezione dei 2000 a
varianza piu' alta. Serve a valutare lo strumento attuale anche con i criteri scritti prima di
vedere i dati (REPORT.md, sezione "Criteri originali e modifiche").
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import anndata as ad
import numpy as np
import scipy.sparse as sp
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import f1_score
from sklearn.model_selection import LeaveOneGroupOut, StratifiedGroupKFold, StratifiedKFold
from sklearn.preprocessing import StandardScaler


def log_cpm(X):
    X = sp.csr_matrix(X, dtype=np.float64)
    tot = np.asarray(X.sum(axis=1)).ravel()
    tot[tot == 0] = 1.0
    X = sp.csr_matrix(sp.diags(1e4 / tot) @ X)
    X.data = np.log1p(X.data)
    return X


def hvg_columns(Xtr, n_top=2000):
    """Indici delle colonne a varianza piu' alta, calcolati sul solo training."""
    m = np.asarray(Xtr.mean(axis=0)).ravel()
    v = np.asarray(Xtr.multiply(Xtr).mean(axis=0)).ravel() - m ** 2
    return np.sort(np.argsort(-v, kind="stable")[: min(n_top, Xtr.shape[1])])


N_TOP = 2000  # None = tutti i geni (pipeline dei criteri originali)


def fit_predict(X, y, tr, te, seed=0):
    if N_TOP is not None:
        X = X[:, hvg_columns(X[tr], N_TOP)]
    sc = StandardScaler(with_mean=False).fit(X[tr])
    clf = LogisticRegression(C=1.0, class_weight="balanced", max_iter=1000, random_state=seed)
    clf.fit(sc.transform(X[tr]), y[tr])
    return clf.predict(sc.transform(X[te]))


def fold_scores(X, y, splits, classes):
    test_cls, std, allc, absent = [], [], [], []
    for tr, te in splits:
        pred = fit_predict(X, y, tr, te)
        test_cls.append(float(f1_score(y[te], pred, average="macro", labels=np.unique(y[te]), zero_division=0)))
        std.append(float(f1_score(y[te], pred, average="macro", zero_division=0)))
        allc.append(float(f1_score(y[te], pred, average="macro", labels=classes, zero_division=0)))
        absent.append(sorted(set(classes) - set(y[te]) - set(pred)))
    return {"classi_del_test": test_cls, "standard": std, "tutte_le_classi": allc,
            "classi_assenti_e_non_predette": absent}


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--h5ad", type=Path, required=True)
    ap.add_argument("--target-col", required=True)
    ap.add_argument("--patient-col", required=True)
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--no-logo", action="store_true")
    ap.add_argument("--pipeline", choices=["sezione6", "originale"], default="sezione6")
    a = ap.parse_args()
    global N_TOP
    N_TOP = 2000 if a.pipeline == "sezione6" else None
    adata = ad.read_h5ad(a.h5ad)
    obs = adata.obs
    keep = obs[a.target_col].notna().values
    y = obs[a.target_col].astype(str).values[keep]
    g = obs[a.patient_col].astype(str).values[keep]
    X = log_cpm(adata.X[keep])
    classes = np.array(sorted(np.unique(y)))
    n_pat = len(np.unique(g))
    res = {"pipeline": a.pipeline, "n_cellule": int(len(y)), "n_pazienti": n_pat, "classi": classes.tolist()}

    sgkf = StratifiedGroupKFold(n_splits=min(5, n_pat), shuffle=True, random_state=0)
    res["raggruppato"] = fold_scores(X, y, list(sgkf.split(X, y, groups=g)), classes)
    n_rand = min(5, int(min(np.unique(y, return_counts=True)[1])))
    skf = StratifiedKFold(n_splits=n_rand, shuffle=True, random_state=0)
    res["casuale"] = fold_scores(X, y, list(skf.split(X, y)), classes)
    if n_pat >= 8 and not a.no_logo:
        logo = list(LeaveOneGroupOut().split(X, y, groups=g))
        res["logo_logreg"] = fold_scores(X, y, logo, classes)
        res["logo_pazienti"] = [str(g[te][0]) for _, te in logo]
    for k in ("raggruppato", "casuale", "logo_logreg"):
        if k in res:
            for m in ("classi_del_test", "standard", "tutte_le_classi"):
                res[k][f"media_{m}"] = float(np.mean(res[k][m]))
    a.out.write_text(json.dumps(res, indent=1, ensure_ascii=False))
    print(json.dumps({k: {m: v.get(f"media_{m}") for m in ("classi_del_test", "standard", "tutte_le_classi")}
                      for k, v in res.items() if isinstance(v, dict)}, indent=1))


if __name__ == "__main__":
    main()
