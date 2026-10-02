"""Tabelle dei campioni dai GEO series matrix (criterio E3: parsing versionato, nessuna
interpretazione oltre a quella dichiarata qui).

Per ogni campione: titolo, source, piattaforma e ogni campo "chiave: valore" di
!Sample_characteristics_ch1, con la chiave come nome di colonna e il valore senza modifiche.
Unica derivazione: per GSE125449 il paziente non e' fra le caratteristiche ma solo nel
titolo (formato "S<campione>_P<n>_LCP<id>"); la colonna `patient_from_title` ne estrae l'ultimo
token "LCP<id>" con una regex esplicita, e lo script verifica che il formato valga per tutti i
titoli.

Uso: python -m validation_esterna.parse_geo
"""

from __future__ import annotations

import gzip
import re
from pathlib import Path

import pandas as pd

DATA = Path(__file__).parent / "data"
OUT = Path(__file__).parent / "results"


def parse_series_matrix(path: Path) -> pd.DataFrame:
    fields: dict[str, list[list[str]]] = {}
    with gzip.open(path, "rt") as fh:
        for line in fh:
            if not line.startswith("!Sample_"):
                continue
            key, *vals = line.rstrip("\n").split("\t")
            fields.setdefault(key, []).append([v.strip('"') for v in vals])
    n = len(fields["!Sample_geo_accession"][0])
    df = pd.DataFrame({"gsm": fields["!Sample_geo_accession"][0],
                       "title": fields["!Sample_title"][0],
                       "source": fields["!Sample_source_name_ch1"][0],
                       "platform": fields["!Sample_platform_id"][0]})
    # Le righe di caratteristiche possono mescolare chiavi diverse (un campione senza un campo
    # fa "scorrere" i successivi): le coppie chiave:valore si raccolgono campione per campione.
    per_sample: list[dict[str, str]] = [{} for _ in range(n)]
    for row in fields.get("!Sample_characteristics_ch1", []):
        assert len(row) == n
        for i, v in enumerate(row):
            if ":" in v:
                k, val = v.split(":", 1)
                assert k.strip() not in per_sample[i], f"chiave ripetuta {k} nel campione {i}"
                per_sample[i][k.strip()] = val.strip()
    chars = pd.DataFrame(per_sample)
    return pd.concat([df, chars], axis=1)


def main() -> None:
    OUT.mkdir(exist_ok=True)
    for gse in ("GSE132465", "GSE131907"):
        df = parse_series_matrix(DATA / f"{gse}_series_matrix.txt.gz")
        df.to_csv(OUT / f"{gse}_samples.csv", index=False)
        print(f"{gse}: {len(df)} campioni, colonne {list(df.columns)}")
    parts = [parse_series_matrix(p) for p in sorted(DATA.glob("GSE125449-GPL*_series_matrix.txt.gz"))]
    df = pd.concat(parts, ignore_index=True)
    pat = df["title"].str.extract(r"^S\d+_P\d+_(LCP\d+)$")[0]
    assert pat.notna().all(), f"titoli fuori formato: {df.title[pat.isna()].tolist()}"
    df["patient_from_title"] = pat
    df.to_csv(OUT / "GSE125449_samples.csv", index=False)
    print(f"GSE125449: {len(df)} campioni, colonne {list(df.columns)}")


if __name__ == "__main__":
    main()
