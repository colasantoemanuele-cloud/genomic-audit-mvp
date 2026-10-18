# Registro della revisione

Memoria del lavoro sul ramo `revisione`. Aggiornato alla fine di ogni fase e di ogni giro.

## Fase 0 — Criteri e velocità dei test (2026-10-18)

- Verifica iniziale: `git remote -v` → `github.com/colasantoemanuele-cloud/genomic-audit-mvp` (corretto).
  Dati di GSE278694 presenti in `../pdac-ml` (`data/interim/sc_raw.h5ad`, `data/raw/tcr/`);
  `GSE278694_DIR` non impostata. Dati della validazione esterna presenti in `validation_esterna/data/`
  (GSE125449 completo e ridotto, annotazioni e series matrix di GSE131907 e GSE132465).
- Commit `1b35b99`: `revisione/CRITERI.md` copiato dal prompt prima di ogni altra modifica.
- `pytest-xdist` e `ruff` aggiunti alle dipendenze `dev`.
- Suite di partenza (codice di `master`, copia separata, con i dati reali):
  `-n auto`: **115 passed, 2 xfailed in 267 s** (`revisione/log/baseline_xdist.txt`).

## Fase 1 — Correzioni

### A1, scelta del metodo (esplorazione, prima della calibrazione)

Simulazioni esplorative in scratchpad (non versionate), solo per scegliere il metodo; le bande della
calibrazione vera sono dichiarate dopo, nel docstring del test, e non dipendono da questi numeri.
Falsi positivi sotto H0 a soglia 0,05 (400 repliche, dimensioni dei gruppi 3-14, valori normali):

| gruppi | percentile (attuale) | BCa | t CR1, G−1 gdl | t CR2, gdl Bell-McCaffrey |
|---|---|---|---|---|
| 5 | 0,180 | 0,170 | 0,058 | 0,045-0,050 |
| 10 | 0,087 | 0,087 | 0,045 | 0,040-0,043 |

Con tassi binari molto eterogenei fra pazienti (come i tassi d'errore per paziente di GSE278694:
media 0,09, d.s. 0,06 per cellula; media 0,18, d.s. 0,125 per clone) la versione su scala naturale
arriva all'8-10%; su scala logit resta fra 3,6% e 6,9% per 4-15 gruppi (1000 repliche).
Scelta: intervallo t cluster-robusto con varianza CR2 e gradi di libertà di Bell-McCaffrey
(Bell e McCaffrey 2002; Imbens e Kolesár 2016), sulla media pesata per cellule (stessa stima
puntuale di prima); per i tassi (valori in [0, 1]) calcolato su scala logit e riportato indietro.
