# Validazione esterna di genomic-audit-mvp

Criteri: `validation_esterna/CRITERI.md`.
- Sezioni 1-4: scritte e committate prima di toccare i dati esterni (commit `0b5a847`).
- Sezione 5: coorte ridotta, decisa dall'utente.
- Sezione 6: nuovo motore del Modulo A, deciso dall'utente.

Le sezioni 5 e 6 sono state aggiunte DOPO aver visto dei risultati, e sono dichiarate come tali.

Ambiente: `.venv` del progetto gestito con uv (anndata 0.13.2, scikit-learn 1.9.0, pandas 3.0.5).
Workstation locale con 12 core e 31 GB di RAM. Numeri e output di questo documento sono
assemblati dai file in `results/` da `make_report.py`.

## 1. Dataset valutati

| Accessione | Tumore | Pazienti | Campioni / cellule | Uso | Eleggibilita' |
|---|---|---|---|---|---|
| [GSE132465](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE132465) | carcinoma colorettale | 23 | 33 campioni (10 pazienti con tumore + mucosa normale) | audit del disegno | eleggibile: E1-E3 (`patient_id`, `tissue type` nel series matrix); nomi di colonna con spazi; stadio e regione mancanti per i normali |
| [GSE131907](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE131907) | adenocarcinoma polmonare | 44 | 58 campioni, 7 tessuti | audit del disegno | eleggibile: E1-E3 (`patient id`, `tissue origin abbrevation`) |
| [GSE125449](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE125449) | tumori primitivi del fegato (HCC, iCCA) | 19 | 19 campioni (uno per paziente) su 2 piattaforme; 9.946 cellule, 18.372 geni comuni | audit del disegno e Modulo A | eleggibile: E1-E3 (paziente ricavato dal titolo `LCP<id>`, regex dichiarata in `parse_geo.py`) e M1-M4 (conteggi interi in mtx, 19 pazienti, etichetta `Type` degli autori, meno di 50.000 cellule) |

Tutte le accessioni sono state verificate leggendo le pagine GEO e scaricando i file (series
matrix e file supplementari). Il Modulo B non e' testato qui: resta validato su un solo
dataset (GSE278694).

## 2. Criteri × dataset (esecuzione finale)

### Audit del disegno

| Dataset | D1 (CLI) | D2 (fatti strutturali) | D3 (classi dei confronti) | D4 (Cramér V) |
|---|---|---|---|---|
{{DROW_GSE132465}}
{{DROW_GSE131907}}
{{DROW_GSE125449}}

**Prima esecuzione**, prima della correzione dei valori mancanti (commit `7bb30ca`):
- GSE132465 **falliva D2 e D4**. L'audit riportava due fatti falsi, {{SOLO_AUDIT_PRIMA}}, e non
  produceva il V per due tabelle valutabili (tissue type × tumor stage, tissue type × region).
- GSE131907 e GSE125449 passavano.

Dopo la correzione ho rieseguito tutto da zero, e poi una seconda volta dopo la correzione del
testo del Cramér V (sezione 4).

### Modulo A (GSE125449 completo: 9.946 cellule, 19 pazienti, 8 classi)

| Criterio | Esito | Numeri misurati |
|---|---|---|
| A1 CLI | OK | `audit-sc leakage --h5ad GSE125449.h5ad --target-col Type --patient-col patient`: exit 0, 324 s reali, 1,2 GB, confronto fra modelli incluso |
| A2 verifica indipendente | OK | {{A2_TXT}} |
| A3 controllo negativo | OK | {{A3_TXT}} |
| A4a matrice sparsa / densa | OK | punteggi dei fold identici (uguaglianza esatta) |
| A4b colonne categoriche / stringhe | OK | punteggi dei fold identici (uguaglianza esatta) |
| A4c classe assente da un fold | OK | {{A4C_TXT}} |

**Prima esecuzione del Modulo A** (motore precedente, criteri delle sezioni 1-5):
- **A1:** NON terminato dopo 5 h 26 min reali e 18 h 52 min di CPU. Interrotto per decisione
  dell'utente; era rallentato anche dalle mie esecuzioni in parallelo.
- **A4c FALLITO, per un difetto di calcolo:** un fold perfetto con una classe assente valeva
  0,667 invece di 1,0. Sui dati LeaveOneGroupOut la media era 0,705, contro 0,759 con la
  definizione standard.

Come previsto dai criteri, mi ero fermato senza correggere. La correzione e' stata richiesta
in seguito dall'utente.

## 3. Output reali

### Audit del disegno — GSE132465
```
{{CLI_GSE132465}}
```

### Audit del disegno — GSE131907
```
{{CLI_GSE131907}}
```

### Audit del disegno — GSE125449
```
{{CLI_GSE125449}}
```

### Modulo A — GSE125449 (stdout di `audit-sc leakage`)
```
{{A1_STDOUT}}
```
Confronto fra modelli (LeaveOneGroupOut): {{LOGO_TXT}}.

## 4. Difetti trovati e corretti

| Difetto | Tipo | Correzione | Commit |
|---|---|---|---|
| Valori mancanti nei metadati scartati in silenzio da crosstab/groupby ma contati in n (pandas >= 3): fatti strutturali falsi e V mancanti | lettura dati | livello esplicito "NA", dichiarato nelle note; test di regressione, che fallisce sul codice precedente | `7bb30ca` |
| La macro-F1 contava come zero le classi assenti dal fold di test | calcolo: riportato senza correggerlo, poi corretto su richiesta dell'utente | macro-F1 sulle classi presenti, classi assenti dichiarate | `bbe7bd0` |
| Modulo A non praticabile su dati reali (oltre 5,5 h senza terminare) | architettura | matrice sparsa, HVG e SVD stimati nel fold, training limitato, fold in parallelo: 403 s sul dataset completo | `bbe7bd0` |
| Classe presente in un solo paziente: traceback di scikit-learn | robustezza | errore esplicito con il motivo | commit finale |
| Mancata convergenza del classificatore non dichiarata (solo un avviso su stderr) | onesta' dell'output | conteggio dei fold, frase nella narrativa, verdetto giallo | commit finale |
| Frase sbagliata: "V = 0,95 [...] sotto la soglia di 0,5" quando la coppia ha un fatto strutturale | testo | frase specifica per le coppie strutturali; test di regressione | commit finale |
| Classi assenti non dichiarate per i fold del confronto fra modelli | onesta' dell'output | elenco nell'app e nel report Markdown | commit finale |
| Telemetria di Streamlit attiva per default; richiesta interattiva di email al primo avvio | privacy e usabilita' | `gatherUsageStats=false`, solo localhost, avvio headless; verificato che il server ascolta solo su 127.0.0.1 e rifiuta le connessioni dall'IP di rete | `6938623` e commit finale |

Suite completa: 74 passed, 2 xfailed prima di questa fase; esecuzione finale: {{SUITE_FINALE}}.

## 5. Risultati descrittivi (non criteri)

Su GSE125449 lo split casuale sovrastima la macro-F1 dello split per paziente di
**{{GAP_TXT}}. La deviazione standard fra fold e' 15,4 volte piu' ampia nello split onesto.
Stabilita' delle spiegazioni (Jaccard fra i 50 geni principali): {{XAI_TXT}}.

Tempi:
- GSE125449 completo: 403 s nel benchmark e 324 s in A1;
- dati sintetici con 25.000 cellule, 20.000 geni e 20 pazienti: {{SINT_S}} s, misurati prima
  della parallelizzazione dei fold. Su questi dati il segnale e' troppo facile (macro-F1 = 1,0):
  la misura vale per i tempi, non per l'accuratezza.

## 6. Verdetto

**PRONTO**, secondo i criteri: dopo le correzioni dichiarate, tutti i criteri D1-D4 sono
soddisfatti sui 3 dataset del disegno, e tutti i criteri A1-A4 sul dataset del Modulo A.

Il verdetto vale con queste precisazioni, che ne fanno parte:
- La prima esecuzione era **NON PRONTO** (GSE132465: D2 e D4; Modulo A: A1 non terminato, A4c
  fallito).
- I criteri sono stati modificati due volte DOPO aver visto dei risultati, per decisione
  dell'utente (sezioni 5 e 6 di CRITERI.md).
- A3 e A4 sono stati eseguiti prima delle ultime modifiche, che aggiungono soltanto messaggi
  e controlli: dichiarazione della non convergenza, errore esplicito per una classe in un solo
  paziente, classi assenti nel confronto fra modelli. I calcoli non sono cambiati, ma A3 e A4
  non sono stati rieseguiti dopo queste modifiche.

## Cosa questa validazione non dimostra

- Il **Modulo B** resta validato su un solo dataset reale (GSE278694, PDAC) e non e' testato qui.
- La generalita' del codice non implica la sua **utilita'** per un utente clinico: lo strumento
  non e' stato provato da nessun utente esterno.
- I dataset sono **pochi**: 3 per il disegno e uno solo per il Modulo A, con una sola etichetta
  (tipo cellulare). Il Modulo A non e' stato provato su un'etichetta di condizione o di esito,
  ne' su dataset reali con piu' di 25.000 cellule.
- Il limite di 20.000 cellule per fold non e' mai entrato in funzione sui dati reali testati.
- Sul dataset reale l'obiettivo dei 10 minuti e' stato misurato su 9.946 cellule; il caso da
  25.000 cellule e' stato misurato solo su dati sintetici.
