"""Assembla validation_esterna_B/REPORT.md dai file di risultato (nessun numero trascritto a mano).

Uso: python -I validation_esterna_B/make_report.py "<riga della suite di test>"
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
R = HERE / "results"
REGIONI = {"Blood": "sangue", "Tumor": "tumore", "Normal adjacent tissue": "adiacente"}


def _cell(reg: dict | None) -> str:
    """Contenuto di una cella della tabella dei pazienti: campione e file presenti."""
    if reg is None:
        return "—"
    files = [n for n, ok in (("matrice", reg["matrice"]), ("VDJ", reg["vdj_filtrato"])) if ok]
    return f"{reg['campione']} ({', '.join(files) or 'nessun file'})"


def main() -> None:
    """Compila REPORT.md da results/eleggibilita.json."""
    suite = sys.argv[1] if len(sys.argv) > 1 else "(non disponibile)"
    el = json.loads((R / "eleggibilita.json").read_text(encoding="utf-8"))
    e1 = el["E1"]
    per_reg = ", ".join(f"{n} {REGIONI[r]}" for r, n in el["campioni_per_regione"].items())
    rows = []
    for pat in sorted(el["pazienti"]):
        v = el["pazienti"][pat]
        r = v["regioni"]
        ct = v["cellule_T_metadati"]
        n = " / ".join(f"{ct.get(k, 0):,}".replace(",", ".") for k in ("Blood", "Tumor", "NAT"))
        rows.append(f"| {pat} | {_cell(r.get('Blood'))} | {_cell(r.get('Tumor'))} | "
                    f"{_cell(r.get('Normal adjacent tissue'))} | {n} | {'sì' if v['e1'] else 'no'} |")
    eleg = ", ".join(e1["pazienti_eleggibili"])
    meta_blood = ", ".join(el["pazienti_con_sangue_nei_metadati"])
    stop = "ci si ferma alla prima condizione non soddisfatta"
    be5 = "non valutabile su questo dataset: GSE139555 non ha proteine di superficie; rimandato alla voce 2 della roadmap"
    nv = "non valutabile: dataset non eleggibile (E1), criterio non eseguito"
    text = f"""# Validazione tecnica esterna del Modulo B su GSE139555

Criteri: `validation_esterna_B/CRITERI.md`, scritti e committati prima di qualsiasi download (commit
`117b9a1`) e non modificati dopo. Numeri e tabelle di questo documento sono assemblati da
`make_report.py` a partire da `results/eleggibilita.json`, prodotto da `eleggibilita.py`.

**Esito: GSE139555 non è eleggibile. La condizione E1 non è soddisfatta: solo {e1['n']} pazienti
({eleg}) hanno sia il sangue sia il tumore, ne servono almeno {e1['richiesti']}. Come previsto da
`CRITERI.md`, mi sono fermato senza eseguire i criteri: il Modulo B non è stato eseguito su
GSE139555, e la validazione tecnica esterna del Modulo B resta da fare.**

Versione dello strumento: 0.9.0 (master = `cb9c1ec`), branch `validazione-b`. Nessun file di
`core/`, `cli.py` o `app.py` è stato modificato.

## 1. Dataset ed eleggibilità

[GSE139555](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE139555) (Wu et al., *Nature*
2020): scRNA-seq + scTCR-seq di tumore, tessuto adiacente normale e sangue periferico, prima del
trattamento, in {el['n_pazienti']} pazienti con tumori di polmone, endometrio, colon-retto e rene.

**File verificati prima del download.** L'elenco dei file supplementari della serie, letto sul
server FTP di GEO, contiene l'archivio `GSE139555_RAW.tar` (609 MB), gli oggetti Seurat integrati
(`.rds.gz`, da 2,8 a 3,8 GB) e tre tabelle di metadati per cellula, fra cui
`GSE139555_tcell_metadata.txt.gz` (3,8 MB; colonne `ident`, cioè il cluster, `patient`, `sample`,
`source`, `clonotype`). `filelist.txt` elenca il contenuto dell'archivio: per ogni campione
`barcodes.tsv.gz`, `genes.tsv.gz`, `matrix.mtx.gz` e `filtered_contig_annotations.csv.gz`.

**Scaricato** in `../dati_esterni/GSE139555/meta/`, fuori dal repository: {', '.join(f'`{f}`' for f in el['file_letti'])}
(circa 4 MB in tutto, meno di 2 s; spazio libero sul disco prima del download: 582 GB). L'archivio
`GSE139555_RAW.tar` non è stato scaricato: la verifica di E1 non ne ha bisogno e, con E1 non
soddisfatta, le regole vietano di proseguire.

**Campioni.** Il series matrix elenca {el['n_campioni']} campioni GSM: {per_reg}. Paziente e regione
sono nelle caratteristiche del campione (`patient: ...`, `region: ...`).

| Paziente | Sangue | Tumore | Adiacente | Cellule T nei metadati (sangue / tumore / adiacente) | E1 |
|---|---|---|---|---|---|
{chr(10).join(rows)}

Le tre fonti concordano: il series matrix, l'elenco dei file dell'archivio e i metadati delle
cellule T degli autori hanno il sangue solo per {meta_blood}. Per questi 4 pazienti matrice e file
VDJ filtrato esistono sia per il sangue sia per il tumore; che le matrici contengano conteggi interi
non è stato verificato, perché con 4 pazienti E1 fallisce comunque.

| Condizione | Esito | Numeri |
|---|---|---|
| E1 (≥ {e1['richiesti']} pazienti con matrice e VDJ filtrato per sangue e tumore) | **non soddisfatta** | {e1['n']} pazienti ({eleg}) |
| E2 (etichetta di cluster da cui derivare la linea CD4/CD8) | non valutata | {stop}; `MAPPA_ETICHETTE.md` non è stata scritta |
| E3 (≥ 50% delle cellule T etichettate con TCR produttivo abbinato) | non valutata | come sopra |

## 2. Criteri × esito

| Criterio | Esito | Numeri |
|---|---|---|
| B-E1 esecuzione | {nv} | — |
| B-E2 implementazione indipendente | {nv}; `indipendente.py` non è stato scritto | — |
| B-E3 controllo negativo | {nv} | — |
| B-E4 controllo positivo | {nv} | — |
| B-E5 riferimento proteico | {be5} | — |
| B-E6 sensibilità alle scelte | {nv} | — |

## 3. Risultati descrittivi

Nessuno: il modulo non è stato eseguito. Una nota per chi riprenderà il lavoro su un altro dataset
o su questo: `GSE139555_tcell_metadata.txt.gz` contiene solo cellule T, quindi gli errori verso NK
osservati nella tesi non potrebbero comparire, e il confronto con la tesi riguarderebbe solo
CD4↔CD8.

## 4. Difetti trovati e correzioni

Nessun difetto dello strumento: lo strumento non è stato eseguito. Nello script di verifica
`eleggibilita.py` la prima esecuzione contava zero cellule T per paziente (indice della colonna
sbagliato di uno, perché l'intestazione dei metadati inizia con il campo vuoto della colonna
indice); corretto prima di produrre `results/eleggibilita.json`. Non cambiava l'esito di E1, che
si basa sul series matrix e sull'elenco dei file.

## 5. Verdetto con i criteri originali

**Non eseguibile: GSE139555 non soddisfa E1 (4 pazienti con sangue e tumore su 5 richiesti).**
Nessun criterio B-E è stato valutato; B-E5 sarebbe stato comunque non valutabile. Nessun criterio è
stato modificato dopo aver visto i dati, quindi non c'è una sezione «Criteri modificati dopo aver
visto i dati».

## 6. Che cosa questa validazione dimostra e che cosa no

Non dimostra nulla sul Modulo B: non si sa ancora se si esegue su dati di un altro laboratorio, se
i suoi calcoli coincidono con un'implementazione indipendente, se non segnala errori dove non ce ne
sono, se rileva errori iniettati e se è stabile alle scelte dichiarate. Resta vero quanto scritto
nei limiti della 0.9.0: il Modulo B è stato eseguito su dati reali solo su GSE278694, il dataset su
cui è nato (verifica di fedeltà, non validazione).

Dimostra soltanto che GSE139555, indicato nella roadmap come dataset per la validazione esterna,
non ha abbastanza pazienti con il sangue: la roadmap lo aveva scelto senza averlo verificato. Anche
quando si troverà un dataset eleggibile, l'esito potrà essere al massimo una validazione tecnica
esterna: per sapere se le etichette segnalate come errate lo sono davvero serve un riferimento
proteico (B-E5, voce 2 della roadmap).

## 7. Strade possibili (da decidere, non eseguite)

1. **Un altro dataset** con scRNA-seq + scTCR-seq 10x, sangue e tumore per almeno 5 pazienti,
   conteggi grezzi e file VDJ filtrati per campione, ed etichette di cluster degli autori. I criteri
   B-E1…B-E6 si possono riusare, ma vanno riscritti con il nome del nuovo dataset e committati prima
   del suo download.
2. **Criteri modificati su GSE139555 con 4 pazienti**, da dichiarare come tali in una sezione
   «Criteri modificati dopo aver visto i dati»; il verdetto principale resterebbe «non eseguibile».
   Senza abbassare il minimo di 5 pazienti (vietato dalle regole) nessun intervallo verrebbe
   prodotto: resterebbero valutabili B-E1, i conteggi e le stime puntuali di B-E2, la media di B-E3
   (non la copertura degli intervalli), B-E4 e B-E6 senza intervalli.
3. **Usare il tessuto adiacente come riferimento** non è previsto dai criteri (riferimento = sangue)
   e cambierebbe il significato del tasso d'errore: sconsigliato.

## 8. Test

Nessun file di `core/`, `cli.py` o `app.py` è cambiato: eseguito solo `pytest -m "not slow"` nella
cartella di lavoro. Esito: {suite}.
"""
    (HERE / "REPORT.md").write_text(text, encoding="utf-8")
    print(f"[ok] scritto {HERE / 'REPORT.md'}")


if __name__ == "__main__":
    main()
