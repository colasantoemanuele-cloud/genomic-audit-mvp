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

### Fase 1 — correzioni applicate (blocchi A, B, C)

Commit: `10fd706` (calibrazione A1 dichiarata prima dell'esecuzione), `58f0733` (calibrazione A2
dichiarata prima), `3bd512b` (blocco A parte 1), `d072810` (blocco A parte 2 e B), `b4581ba`
(calibrazione del Cramér V su fogli piccoli, dichiarata e committata prima dell'esecuzione).

- **A1.** `core.stats.cluster_t_interval`: t cluster-robusto CR2, gdl di Bell-McCaffrey, logit per i
  tassi. Calibrazione (1000 repliche, 5/8/10/15 gruppi, 5 scenari): 3,6%-6,3% (banda 2-8%);
  percentile a 5 gruppi 14,5% (`revisione/log/a1_calibrazione.txt`). Calibrato a 5 gruppi:
  `MIN_PATIENTS_FOR_CI_DEFAULT` resta 5.
- **A2.** Calibrazione in LOGO (8/10/15/20 gruppi, 400 repliche): Nadeau-Bengio 0,018/0,007/0,007/
  0,005 (fuori banda, conservativo); Wilcoxon 0,060/0,083/0,072/0,045 (fuori banda a 10); t non
  corretto 0,065/0,077/0,080/0,043 (in banda, non usato). NB e Wilcoxon `xfail(strict=True)` e non
  esposti (`revisione/log/a2_calibrazione.txt`). Test tautologico rimosso; spiegazione dello xfail a
  k=5 corretta; README corretto.
- **A3.** Riferimento fissato (regressione logistica), Holm sui confronti (calcolato ma non esposto,
  per A2); confronto esposto descrittivo (Δ e pazienti meglio/peggio/uguale). Metodo di Wilcoxon
  scelto esplicitamente (esatto / permutazione / normale) e dichiarato.
- **A4.** Filtro sulle chiavi della marker_map in `_judged_cells` (tasso per cellula, unità per
  cellula, flag). Caso sintetico: master 0,000 su 267 cellule, ora 0,386 su 435
  (`tests/test_input_errors.py`). GSE278694: per cellula invariato, 0,0930.
- **A5.** Testo e docstring corretti; sensibilità a pesi uguali (`excess_equal`, valore atteso
  ½‖pA − pB‖²); generatore con compartimenti sbilanciati, cloni non condivisi, dropout, RNA
  ambientale. Caso 3 contro 60 cellule con 40% di errore: pesato −0,131, pesi uguali +0,125.
- **A6.** Fatti strutturali nella sintesi (uno-a-uno rosso, gli altri gialli), coppie strutturali
  escluse dalle valutabili, livello assente in grigio, verde = «nessun limite rilevato dalle
  regole»; `tests/test_verdict.py` riscritto con oggetti veri.
- **A7.** Frazione di CD8 fuori da app, report e sintesi; solo `--sperimentale-cd8`; rifiuto fuori
  dall'intervallo invertibile (prima: troncamento, con intervalli «0,00-0,00»).
- **A8.** Barcode duplicati, chiavi (paziente, compartimento, libreria, barcode) duplicate, più file
  VDJ sulla stessa coppia senza libreria: errore; colonna libreria supportata. Contig a pari UMI:
  più reads, poi CDR3 alfabetica (346 cellule su GSE278694). Pareggi della maggioranza per clone:
  probabilità di errore con scelta a caso (11 unità su 288 nel tumore). Flag per posizione in CLI;
  export rifiuta obs_names non unici. Conteggi interi richiesti per i marcatori. App: righe di
  confronto malformate segnalate, file illeggibili con messaggio chiaro. CLI: `--compare` su
  colonna non dichiarata, ruolo ripetuto, colonne ripetute in `--h5ad`.
- **A9.** Tabelle con più della metà di righe ripetute rifiutate; nota sulle librerie che
  collassano con `--h5ad` senza colonna libreria; testo del pavimento di Wilcoxon corretto.
- **A10.** Upload letto in memoria (nessun file temporaneo); README: workstation monoutente;
  Streamlit minimo 1.51 (1.45-1.50 falliscono, `revisione/log/streamlit_versione_minima.txt`).
- **A11.** Limiti dichiarati nel report HTML, `GAP_ALERT`, soglia stampata, piè di pagina, colori
  della tabella dei confronti, «fuori dalla marker_map».
- **B1-B8.** Ipotesi biologiche in README, report e docstring; B2 misurato: le NK spiegano il 71%
  della differenza fra convenzioni, non il 43% indicato nel prompt (`validation/results/
  b2_quota_nk.json`); la cifra «72-79% di doppietti» è superata nella tesi stessa (README di
  pdac-ml: maggioranza su tutti i tessuti, non nel tumore) ed è riportata così, senza numeri.
  «Validazione» → «verifica di fedeltà del porting»; note sugli intervalli con 4 pazienti nei
  risultati di `validation/`; Modulo A riformulato; citazioni complete; docstring aggiunte.
- **C1.** `docs/GUIDA_UTENTE.md`. **C2.** Prova su GSE144735: `revisione/C2_prova_guida.md`.

Numeri reali cambiati (GSE278694, codice nuovo; vecchio → nuovo):
- eccesso di discordanza +0,059 [+0,034; +0,098] → +0,059 [+0,020; +0,098]; nuovo: pesi uguali
  +0,043 [+0,012; +0,074];
- per cellula, Tumor 0,093 [0,050; 0,144] → 0,093 [0,043; 0,188];
- per clone (opzioni della tesi), Tumor 0,195 [0,123; 0,274] → 0,186 [0,113; 0,290]; Adjacent
  0,059 → 0,064; differenza appaiata +0,080 → +0,064 (88 cloni, IC non prodotto);
- per clone (default), Tumor 0,197 [0,124; 0,276] → 0,187 [0,114; 0,293]; Adjacent 0,051 → 0,055;
  differenza appaiata +0,079 → +0,064 (89 cloni, IC non prodotto);
- fedeltà del porting (S5): 0,1951 [0,1230; 0,2737] → 0,1944 [0,1225; 0,2731] (contig a pari UMI);
- cellule di riferimento della frazione di CD8: 4851 → 4853.

Incidente: `revisione/esegui_validazioni.sh` modificato mentre bash lo eseguiva; fermato lo script
(il processo Python in corso è proseguito) e lanciati i passi restanti da una copia fissa. D'ora in
poi lo script si lancia da una copia.
