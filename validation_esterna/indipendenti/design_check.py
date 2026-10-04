"""Verifica indipendente dell'audit del disegno (criteri D2, D3, D4).

NON importa nulla da core/: reimplementa, a partire dalle regole dichiarate in CRITERI.md, i
fatti strutturali, le classi dei confronti e la valutabilità del Cramér V. Lettura dei
metadati identica a quella documentata della CLI: CSV letto come stringhe e valori mancanti
trattati come il livello "nan".
"""

from __future__ import annotations

from math import comb

import numpy as np
import pandas as pd

TECHNICAL = {"library", "batch", "chemistry", "protocol", "date"}


def read_sheet(path, cols: list[str]) -> pd.DataFrame:
    # valori mancanti -> livello esplicito "nan" (con pandas 3, astype(str) li lascia mancanti)
    return pd.read_csv(path, dtype=str)[cols].fillna("nan").astype(str).reset_index(drop=True)


def nested(df: pd.DataFrame, a: str, b: str) -> bool:
    return all(len(set(g)) == 1 for _, g in df.groupby(a)[b])


def structural_facts(df: pd.DataFrame, roles: dict[str, str], outcome_like: list[str],
                     patient: str, tissue: str | None) -> set[tuple]:
    n = len(df)
    constant = {c for c in roles if df[c].nunique() < 2}
    ident = {c for c in roles if c not in constant and df[c].nunique() == n}
    active = [c for c in roles if c not in constant and c not in ident]
    facts: set[tuple] = set()
    for i in range(len(active)):
        for j in range(i + 1, len(active)):
            a, b = active[i], active[j]
            ab, ba = nested(df, a, b), nested(df, b, a)
            if ab and ba:
                facts.add(("uno-a-uno", a, b))
                continue
            for inner, outer, ok in ((a, b, ab), (b, a, ba)):
                if not ok:
                    continue
                if outer in outcome_like and roles[inner] != "outcome":
                    facts.add(("esito-determinato", inner, outer))
                else:
                    facts.add(("annidamento", inner, outer))
    if tissue is not None:
        combo = df[patient] + "|" + df[tissue]
        for c, r in roles.items():
            if r in TECHNICAL and df[c].nunique() >= 2:
                t = pd.DataFrame({"t": df[c], "u": combo})
                if nested(t, "t", "u") and nested(t, "u", "t"):
                    facts.add(("unità-tecnica", c))
    return facts


def classify(df: pd.DataFrame, roles: dict[str, str], patient: str, factor: str, a: str, b: str,
             min_units: int = 5) -> tuple[str, int]:
    rows = df[df[factor].isin([a, b])]
    pa = set(rows.loc[rows[factor] == a, patient])
    pb = set(rows.loc[rows[factor] == b, patient])
    if not pa or not pb:
        return "non stimabile", 0
    both = pa & pb
    role = roles.get(factor, "outcome")
    if both:
        units = len(both)
        sub = rows[rows[patient].isin(both)]
        low = units < min_units or 2 / 2 ** units > 0.05
    elif role == "outcome":
        units = min(len(pa), len(pb))
        sub = rows
        low = units < min_units or min(1.0, 2 / comb(len(pa) + len(pb), len(pa))) > 0.05
    else:
        return "non stimabile", 0
    for c, r in roles.items():
        if r in TECHNICAL and c != factor:
            k = sub[c].nunique()
            if 2 <= k < len(sub) and nested(sub, c, factor):
                return "non stimabile", units
    return ("stimabile con bassa potenza" if low else "stimabile"), units


def cramer_v_expected(df: pd.DataFrame, roles: dict[str, str]) -> dict[frozenset, float | None]:
    """Per ogni coppia di fattori attivi: None se 'non valutabile' secondo le regole
    (n < 10 o n/(r*c) < 2), altrimenti il V di Bergsma calcolato qui da zero."""
    n = len(df)
    constant = {c for c in roles if df[c].nunique() < 2}
    ident = {c for c in roles if c not in constant and df[c].nunique() == n}
    active = [c for c in roles if c not in constant and c not in ident]
    out = {}
    for i in range(len(active)):
        for j in range(i + 1, len(active)):
            a, b = active[i], active[j]
            t = pd.crosstab(df[a], df[b]).to_numpy(dtype=float)
            r, k = t.shape
            if n < 10 or n / (r * k) < 2:
                out[frozenset((a, b))] = None
                continue
            exp = t.sum(1, keepdims=True) * t.sum(0, keepdims=True) / n
            chi2 = float(((t - exp) ** 2 / exp).sum())
            phi2c = max(0.0, chi2 / n - (k - 1) * (r - 1) / (n - 1))
            rc, kc = r - (r - 1) ** 2 / (n - 1), k - (k - 1) ** 2 / (n - 1)
            den = min(kc - 1, rc - 1)
            out[frozenset((a, b))] = float(np.sqrt(phi2c / den)) if den > 0 else None
    return out
