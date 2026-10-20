"""Verifica della condizione di eleggibilità E1 su GSE139555 (validation_esterna_B/CRITERI.md).

Legge solo tre file piccoli scaricati da GEO, senza l'archivio GSE139555_RAW.tar:
- ``filelist.txt``: elenco dei file contenuti nell'archivio (matrici e file VDJ per campione);
- ``GSE139555_series_matrix.txt.gz``: paziente e regione di ogni campione GSM;
- ``GSE139555_tcell_metadata.txt.gz``: metadati delle cellule T degli autori (controllo incrociato).

E1: per almeno 5 pazienti esistono, sia per il sangue sia per il tumore, la matrice dei conteggi e il
file VDJ filtrato. Scrive ``results/eleggibilita.json``.

Uso: python -I validation_esterna_B/eleggibilita.py ../dati_esterni/GSE139555/meta
"""

from __future__ import annotations

import gzip
import json
import re
import sys
from collections import Counter, defaultdict
from pathlib import Path

HERE = Path(__file__).resolve().parent
MIN_PATIENTS = 5
# nome dei file nell'archivio: GSM<id>_SAM<id>-<codice campione>.<tipo>
FILE_RE = re.compile(r"^(GSM\d+)_SAM\d+-([a-z]+\d+)\.(barcodes\.tsv|genes\.tsv|matrix\.mtx|"
                     r"filtered_contig_annotations\.csv)\.gz$")
MATRIX_FILES = {"barcodes.tsv", "genes.tsv", "matrix.mtx"}
VDJ_FILE = "filtered_contig_annotations.csv"


def series_matrix(path: Path) -> dict[str, dict[str, str]]:
    """GSM -> {title, patient, region} dal series matrix."""
    rows: dict[str, list[list[str]]] = defaultdict(list)
    with gzip.open(path, "rt", encoding="utf-8") as f:
        for line in f:
            if line.startswith(("!Sample_geo_accession", "!Sample_title", "!Sample_characteristics_ch1")):
                key, *vals = line.rstrip("\n").split("\t")
                rows[key].append([v.strip('"') for v in vals])
    gsm, title = rows["!Sample_geo_accession"][0], rows["!Sample_title"][0]
    out = {g: {"title": t} for g, t in zip(gsm, title)}
    for vals in rows["!Sample_characteristics_ch1"]:
        for g, v in zip(gsm, vals):
            k, _, val = v.partition(": ")
            if k in ("patient", "region"):
                out[g][k] = val
    return out


def file_list(path: Path) -> dict[str, set[str]]:
    """GSM -> tipi di file presenti nell'archivio."""
    out: dict[str, set[str]] = defaultdict(set)
    for line in path.read_text(encoding="utf-8").splitlines()[1:]:
        parts = line.split("\t")
        m = FILE_RE.match(parts[1]) if len(parts) > 1 and parts[0] == "File" else None
        if m:
            out[m.group(1)].add(m.group(3))
    return out


def tcell_metadata(path: Path) -> dict[str, dict[str, int]]:
    """Paziente -> {sorgente: cellule T} dai metadati degli autori."""
    counts: Counter = Counter()
    with gzip.open(path, "rt", encoding="utf-8") as f:
        header = f.readline().rstrip("\n").split("\t")
        ip, isrc = header.index("patient"), header.index("source")  # l'intestazione inizia con il campo vuoto dell'indice
        for line in f:
            p = line.rstrip("\n").split("\t")
            counts[(p[ip], p[isrc])] += 1
    out: dict[str, dict[str, int]] = defaultdict(dict)
    for (pat, src), n in sorted(counts.items()):
        out[pat][src] = n
    return dict(out)


def _complete(regions: dict, region: str) -> bool:
    """Matrice e VDJ filtrato presenti per quella regione del paziente."""
    return region in regions and regions[region]["matrice"] and regions[region]["vdj_filtrato"]


def main(meta_dir: Path) -> None:
    sm = series_matrix(meta_dir / "GSE139555_series_matrix.txt.gz")
    fl = file_list(meta_dir / "filelist.txt")
    tc = tcell_metadata(meta_dir / "GSE139555_tcell_metadata.txt.gz")
    samples = []
    for gsm, info in sorted(sm.items()):
        kinds = fl.get(gsm, set())
        samples.append({"gsm": gsm, **info, "matrice": MATRIX_FILES <= kinds, "vdj_filtrato": VDJ_FILE in kinds})
    patients: dict[str, dict] = {}
    for s in samples:
        p = patients.setdefault(s["patient"], {"regioni": {}})
        p["regioni"][s["region"]] = {"campione": s["title"], "matrice": s["matrice"], "vdj_filtrato": s["vdj_filtrato"]}
    for pat, p in patients.items():
        p["sangue_completo"], p["tumore_completo"] = _complete(p["regioni"], "Blood"), _complete(p["regioni"], "Tumor")
        p["e1"] = p["sangue_completo"] and p["tumore_completo"]
        p["cellule_T_metadati"] = tc.get(pat, {})
    eligible = sorted(p for p, v in patients.items() if v["e1"])
    out = {
        "accessione": "GSE139555",
        "file_letti": ["filelist.txt", "GSE139555_series_matrix.txt.gz", "GSE139555_tcell_metadata.txt.gz"],
        "n_campioni": len(samples),
        "n_pazienti": len(patients),
        "campioni_per_regione": dict(Counter(s["region"] for s in samples)),
        "campioni": samples,
        "pazienti": patients,
        "E1": {"richiesti": MIN_PATIENTS, "pazienti_eleggibili": eligible, "n": len(eligible),
               "soddisfatta": len(eligible) >= MIN_PATIENTS},
        "pazienti_con_sangue_nei_metadati": sorted(p for p, v in tc.items() if v.get("Blood", 0) > 0),
    }
    (HERE / "results").mkdir(exist_ok=True)
    (HERE / "results" / "eleggibilita.json").write_text(json.dumps(out, indent=2, ensure_ascii=False) + "\n")
    print(f"campioni {len(samples)}, pazienti {len(patients)}, per regione {out['campioni_per_regione']}")
    for pat in sorted(patients):
        v = patients[pat]
        print(f"  {pat:7s} sangue {'sì' if v['sangue_completo'] else 'no':2s}  tumore "
              f"{'sì' if v['tumore_completo'] else 'no':2s}  cellule T {v['cellule_T_metadati']}")
    print(f"E1: {len(eligible)} pazienti con sangue e tumore completi ({', '.join(eligible)}); "
          f"richiesti {MIN_PATIENTS}: {'SODDISFATTA' if out['E1']['soddisfatta'] else 'NON SODDISFATTA'}")


if __name__ == "__main__":
    main(Path(sys.argv[1]))
