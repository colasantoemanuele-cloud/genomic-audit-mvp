# Validazione esterna di genomic-audit-mvp

Criteri: `validation_esterna/CRITERI.md`. (Il file cita `REPORT_INTERVENTI.md` e i `logs/`: oggi sono
in `docs/sviluppo/`. `CRITERI.md` non è stato aggiornato di proposito, per lasciarlo com'era.)
- Sezioni 1-4: scritte e committate prima di toccare i dati esterni (commit `0b5a847`).
- Sezione 5: coorte ridotta, decisa dall'utente (commit `b08b2d2`).
- Sezione 6: nuovo motore del Modulo A, deciso dall'utente (commit `0662a40`).

Le sezioni 5 e 6 sono state aggiunte DOPO aver visto dei risultati, e sono dichiarate come tali.
**Verdetto con i criteri originali: NON PRONTO. Verdetto con i criteri modificati: PRONTO.**
Testo dei criteri prima e dopo, motivi ed esiti con entrambe le versioni: sezioni 6 e 7.

Ambiente: `.venv` del progetto gestito con uv (anndata 0.13.2, scikit-learn 1.9.0, pandas 3.0.5).
Workstation locale con 12 core e 31 GB di RAM. Numeri e output di questo documento sono
assemblati dai file in `results/` da `make_report.py`.

## 1. Dataset valutati

| Accessione | Tumore | Pazienti | Campioni / cellule | Uso | Eleggibilità |
|---|---|---|---|---|---|
| [GSE132465](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE132465) | carcinoma colorettale | 23 | 33 campioni (10 pazienti con tumore + mucosa normale) | audit del disegno | eleggibile: E1-E3 (`patient_id`, `tissue type` nel series matrix); nomi di colonna con spazi; stadio e regione mancanti per i normali |
| [GSE131907](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE131907) | adenocarcinoma polmonare | 44 | 58 campioni, 7 tessuti | audit del disegno | eleggibile: E1-E3 (`patient id`, `tissue origin abbrevation`) |
| [GSE125449](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE125449) | tumori primitivi del fegato (HCC, iCCA) | 19 | 19 campioni (uno per paziente) su 2 piattaforme; 9.946 cellule, 18.372 geni comuni | audit del disegno e Modulo A | eleggibile: E1-E3 (paziente ricavato dal titolo `LCP<id>`, regex dichiarata in `parse_geo.py`) e M1-M4 (conteggi interi in mtx, 19 pazienti, etichetta `Type` degli autori, meno di 50.000 cellule) |

Tutte le accessioni sono state verificate leggendo le pagine GEO e scaricando i file (series
matrix e file supplementari). Il Modulo B non è testato qui: resta validato su un solo
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
| A1 CLI | OK | `audit-sc leakage --h5ad GSE125449.h5ad --target-col Type --patient-col patient`: {{A1_TXT}}, confronto fra modelli incluso |
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

Come previsto dai criteri, mi ero fermato senza correggere. La correzione è stata richiesta
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
| Mancata convergenza del classificatore non dichiarata (solo un avviso su stderr) | onestà dell'output | conteggio dei fold, frase nella narrativa, verdetto giallo | commit finale |
| Frase sbagliata: "V = 0,95 [...] sotto la soglia di 0,5" quando la coppia ha un fatto strutturale | testo | frase specifica per le coppie strutturali; test di regressione | commit finale |
| Classi assenti non dichiarate per i fold del confronto fra modelli | onestà dell'output | elenco nell'app e nel report Markdown | commit finale |
| Telemetria di Streamlit attiva per default; richiesta interattiva di email al primo avvio | privacy e usabilità | `gatherUsageStats=false`, solo localhost, avvio headless; verificato che il server ascolta solo su 127.0.0.1 e rifiuta le connessioni dall'IP di rete | `6938623` e commit finale |

Suite completa: 74 passed, 2 xfailed prima di questa fase; esecuzione finale: {{SUITE_FINALE}}.

## 5. Risultati descrittivi (non criteri)

Su GSE125449 lo split casuale sovrastima la macro-F1 dello split per paziente di
**{{GAP_TXT}}. La deviazione standard fra fold è 15,4 volte più ampia nello split onesto.
Stabilità delle spiegazioni (Jaccard fra i 50 geni principali): {{XAI_TXT}}.

Tempi:
- GSE125449 completo: 403 s nel benchmark; in A1, 324 s alla prima misura e {{A1_S}} s nella
  riesecuzione con il codice finale (2026-10-04, macchina condivisa con altri processi);
- dati sintetici con 25.000 cellule, 20.000 geni e 20 pazienti: {{SINT_S}} s, misurati prima
  della parallelizzazione dei fold. Su questi dati il segnale è troppo facile (macro-F1 = 1,0):
  la misura vale per i tempi, non per l'accuratezza.

## 6. Criteri originali e modifiche

I criteri sono stati scritti prima di toccare i dati esterni (sezioni 1-4 di `CRITERI.md`, commit
`0b5a847`) e poi modificati **due volte dopo aver visto dei risultati**. Criteri cambiati dopo
aver visto i risultati tolgono valore a una validazione: per questo qui sotto ci sono, per
ciascuna modifica, il testo prima e dopo, il motivo, e l'esito di ogni dataset con entrambe le
versioni. `CRITERI.md` non è stato ritoccato in questa revisione, neppure nella forma.

Tre versioni dei criteri:
- **O** = sezioni 1-4, commit `0b5a847` (2026-10-01): gli originali;
- **M1** = O + sezione 5, commit `b08b2d2` (2026-10-02): prima modifica;
- **M2** = O + sezione 6, commit `0662a40` (2026-10-02): seconda modifica, quella usata per il
  verdetto finale. La sezione 6 annulla la 5 (si torna al dataset completo).

I criteri D1-D4 dell'audit del disegno **non sono mai stati modificati**: le due modifiche
riguardano solo il Modulo A.

### Modifica 1 — coorte ridotta per il Modulo A (sezione 5, commit `b08b2d2`)

- **Prima** (sezione 1, commit `0b5a847`): il Modulo A si valuta su un dataset con «M3. Almeno 8
  pazienti» e «M5. Se ci sono più di 50.000 cellule: sottocampionamento stratificato per
  paziente». GSE125449 ha 9.946 cellule: andava quindi valutato **intero** (19 pazienti), con
  A1-A4.
- **Dopo** (sezione 5): «10 pazienti estratti a caso fra i 19 di GSE125449 (numpy
  default_rng(0), scelta senza reinserimento); per ciascun paziente al massimo 200 cellule
  estratte a caso [...]. Su questa coorte si valutano A1-A4 con gli stessi criteri della
  sezione 3.» Risultato: 10 pazienti, 1.861 cellule.
- **Motivo:** sul dataset completo `cli.py leakage` non era terminato dopo 5 h 26 min reali
  (18 h 52 min di CPU).
- **Che cosa era già noto quando è stata scritta:** che A1 non terminava sul dataset completo.
  La modifica rende il criterio più facile da soddisfare (meno dati) e va letta così.

### Modifica 2 — pipeline dichiarata del nuovo motore (sezione 6, commit `0662a40`)

- **Prima** (sezione 3, A2, commit `0b5a847`): lo script indipendente riesegue «CPM a 1e4 +
  log1p, StandardScaler senza centratura; LogisticRegression(C=1, class_weight="balanced",
  max_iter=1000, seme 0); macro-F1», cioè **tutti i geni** e la macro-F1 standard di
  scikit-learn, e confronta le medie (split per paziente e LeaveOneGroupOut) con tolleranza
  0,02. A4(c): se rivela un difetto nel calcolo della metrica, «mi fermo e lo riporto senza
  correggerlo».
- **Dopo** (sezione 6): la pipeline che lo script indipendente deve riprodurre diventa quella
  del nuovo motore: «i 2000 geni con varianza più alta dei valori log-CPM, stimata sul SOLO
  training del fold» e «macro-F1 calcolata sulle classi presenti nel fold di test». A2
  confronta i punteggi per fold (anche dello split casuale), stessa tolleranza 0,02. A4(c) è
  «ora atteso come gestito». Il dataset torna a essere GSE125449 completo.
- **Motivo:** l'utente ha chiesto di correggere l'architettura del Modulo A, non praticabile su
  dati reali, e di correggere la macro-F1 che contava come zero le classi assenti.
- **Che cosa era già noto quando è stata scritta:** tutti gli esiti della prima esecuzione,
  compreso il fallimento di A4(c). Inoltre la sezione 6 è stata committata **nello stesso
  commit** dei risultati finali (`0662a40`): la cronologia di git non permette di dimostrare
  che sia stata scritta prima di ottenerli. Con questa modifica lo strumento e lo script che lo
  verifica sono stati cambiati insieme: A2 secondo M2 verifica che due implementazioni della
  stessa pipeline coincidano, non che la pipeline sia quella dichiarata all'inizio.

### Esiti con ciascuna versione dei criteri

**Audit del disegno** (criteri D1-D4 identici in O, M1 e M2):

| Dataset | Prima esecuzione | Dopo la correzione dei valori mancanti (`7bb30ca`) |
|---|---|---|
| GSE132465 | **NON soddisfatti: D2 e D4 falliti** (due fatti strutturali falsi, due V mancanti) | D1-D4 soddisfatti |
| GSE131907 | D1-D4 soddisfatti | D1-D4 soddisfatti |
| GSE125449 | D1-D4 soddisfatti | D1-D4 soddisfatti |

Qui non è cambiato il criterio: è stato corretto un difetto dello strumento, trovato proprio dai
criteri originali.

**Modulo A** (GSE125449). «Prima esecuzione» = motore di allora; «codice finale» = motore
attuale, rivalutato in questa revisione senza modificarlo: è cambiata solo la regola di
valutazione nello script di verifica (`criteri_originali.py`, `leakage_check.py --pipeline
originale`; numeri in `results/criteri_originali.json`).

| Criteri | Dataset | Codice | A1 | A2 (tolleranza 0,02) | A3 | A4a, A4b | A4c | Tutti soddisfatti |
|---|---|---|---|---|---|---|---|---|
| O | completo, 19 pazienti | prima esecuzione | **NO**: non terminato dopo 5 h 26 min | non valutabile (nessun output) | non valutabile | non valutabili | **NO**: fold perfetto con una classe assente = 0,667 invece di 1,0 | **NO** |
| M1 | ridotto, 10 pazienti | motore precedente | {{M1_PRIMA}} |
{{CO_ROW_O}}
{{CO_ROW_M1}}
{{CO_ROW_M2}}

Due precisazioni sulla riga «O, codice finale», perché non venga letta come più favorevole di
quanto è:
- A2 è soddisfatto alla lettera, con poco margine. Il criterio originale confronta la media
  dello strumento con la macro-F1 **standard** dello script indipendente su tutti i geni:
  LeaveOneGroupOut {{O_LOGO_TXT}}. Ma lo strumento finale calcola la macro-F1 sulle sole classi
  presenti nel fold: con la stessa definizione dai due lati la differenza è {{O_LOGO_LIKE}},
  **oltre** la tolleranza. Le due differenze (selezione dei geni e definizione della metrica)
  hanno segno opposto e in parte si compensano.
- È una rivalutazione fatta **dopo** aver visto i risultati e dopo aver riscritto lo strumento:
  non ha il valore di un test scritto in anticipo.

Sulla riga «M1, codice finale»: A2 **fallisce** sul confronto LeaveOneGroupOut
({{M1_LOGO_TXT}}). In 7 fold su 10 il paziente lasciato fuori non ha tutte le classi, e lì la
macro-F1 standard dello script e quella sulle classi presenti dello strumento divergono; con la
stessa definizione dai due lati la differenza è {{M1_LOGO_LIKE}}. Il criterio, così come era
scritto, non è soddisfatto.

## 7. Verdetto, con i criteri originali e con quelli modificati

**Con i criteri originali (O, commit `0b5a847`): NON PRONTO.** È l'esito della prima esecuzione,
l'unica fatta senza conoscere i risultati: GSE132465 falliva D2 e D4; il Modulo A non terminava
(A1) e contava come zero le classi assenti (A4c). Secondo i criteri stessi («Non esistono
verdetti parziali favorevoli») questo è il verdetto della validazione pre-dichiarata.

**Con i criteri modificati (M2, sezione 6, commit `0662a40`): PRONTO.** Dopo le correzioni
dichiarate nella sezione 4, tutti i criteri D1-D4 sono soddisfatti sui 3 dataset del disegno e
tutti i criteri A1-A4 sul dataset completo del Modulo A.

Fra i due:
- il codice finale, rivalutato con i criteri originali O, li soddisfa tutti (A2 con le riserve
  scritte sopra);
- il codice finale, rivalutato con la prima modifica M1 (coorte ridotta), **non** soddisfa A2.

Che cosa si può concludere: i difetti trovati dai criteri originali sono stati corretti e il
codice corretto supera gli stessi controlli. Che cosa **non** si può concludere: che lo strumento
abbia superato una validazione indipendente al primo tentativo. Non l'ha superata.

Altre precisazioni che fanno parte del verdetto:
- A3 e A4 sul dataset completo sono stati eseguiti prima delle ultime modifiche, che aggiungono
  soltanto messaggi e controlli: dichiarazione della non convergenza, errore esplicito per una
  classe in un solo paziente, classi assenti nel confronto fra modelli. I calcoli non sono
  cambiati, ma A3 e A4 sul dataset completo non sono stati rieseguiti dopo queste modifiche.
  Sulla coorte ridotta A1, A3 e A4 sono stati eseguiti con il codice attuale (riga M1).
- Il verdetto riguarda la correttezza e la generalità del codice sull'audit del disegno e sul
  Modulo A. Non riguarda il Modulo B.

## Cosa questa validazione non dimostra

- Il **Modulo B** resta validato su un solo dataset reale (GSE278694, PDAC) e non è testato qui.
- La generalità del codice non implica la sua **utilità** per un utente clinico: lo strumento
  non è stato provato da nessun utente esterno.
- I dataset sono **pochi**: 3 per il disegno e uno solo per il Modulo A, con una sola etichetta
  (tipo cellulare). Il Modulo A non è stato provato su un'etichetta di condizione o di esito,
  né su dataset reali con più di 25.000 cellule.
- Il limite di 20.000 cellule per fold non è mai entrato in funzione sui dati reali testati.
- Sul dataset reale l'obiettivo dei 10 minuti è stato misurato su 9.946 cellule; il caso da
  25.000 cellule è stato misurato solo su dati sintetici.
