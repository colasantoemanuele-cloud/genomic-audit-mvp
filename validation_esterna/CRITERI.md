# Criteri della validazione esterna (scritti PRIMA di toccare i dati esterni)

Scopo: verificare che **l'audit del disegno** e il **Modulo A** funzionino correttamente su
dataset diversi da GSE278694, con layout di metadati e formati diversi. E' una validazione
di correttezza e generalita' del codice, non di utilita' ne' di risultati scientifici. Il
Modulo B non e' testato qui e resta validato su un solo dataset (GSE278694).

Ambiente: `.venv` del progetto gestito con uv (`uv run ...`). E' l'ambiente indicato dal
README e dal `uv.lock`; gli import reali sono verificati (anndata 0.13.2, scikit-learn 1.9.0,
pandas 3.0.5). L'ambiente conda `base` non contiene anndata.

Questi criteri non vengono modificati dopo aver visto i risultati. Un dataset viene scartato
solo per i motivi di eleggibilita' elencati qui sotto, e il motivo viene documentato.

## 1. Eleggibilita' dei dataset

Per ogni candidato valutato si documenta: accessione (verificata leggendo davvero la pagina
GEO o scaricando i file), tumore, eleggibile si'/no e perche'.

**Audit del disegno**: almeno 2 dataset, ciascuno con tutti questi requisiti:
- E1. Pubblico su GEO, accessione verificata.
- E2. Tumore diverso dal PDAC, e i due dataset con tumori diversi fra loro.
- E3. Metadati per campione disponibili (GEO series matrix, SOFT/MINiML o tabelle
  supplementari), con almeno una colonna paziente/soggetto e una tessuto/condizione
  ricavabili senza interpretazioni arbitrarie. Se servono passaggi di parsing (es. "patient:
  P01" nelle caratteristiche GEO), sono scritti in uno script versionato.
- E4. Preferibilmente un layout di colonne diverso da GSE278694. Non e' un requisito
  escludente.

**Modulo A**: almeno 1 dataset con tutti questi requisiti:
- M1. scRNA-seq pubblico su GEO, accessione verificata, diverso da GSE278694.
- M2. Conteggi grezzi (interi non normalizzati) a livello di cellula.
- M3. Almeno 8 pazienti.
- M4. Un'etichetta di classe per cellula (tipo cellulare o condizione) fornita dagli autori,
  non ricalcolata da me.
- M5. Se ci sono piu' di 50.000 cellule: sottocampionamento stratificato per paziente
  (stessa frazione per paziente, seme fisso 0), dichiarato con i numeri prima e dopo.

Se dopo una ricerca ragionevole nessun dataset soddisfa M1-M4, l'esito e' "nessun dataset
eleggibile per il Modulo A": i requisiti non vengono abbassati.

## 2. Criteri di successo — Audit del disegno (tutti obbligatori per ciascun dataset)

- **D1.** `cli.py design` termina senza errori (exit code 0, nessun traceback) usando solo i
  flag documentati nel README (`--meta`/`--h5ad`, `--patient-col`, `--tissue-col`,
  `--technical ruolo=colonna`, `--outcome-col`, `--compare colonna:a:b`), con i nomi di
  colonna originali del dataset.
- **D2.** Uno script indipendente (`validation_esterna/indipendenti/design_check.py`, solo
  pandas, NESSUN import da `core/`) ricalcola, per le stesse colonne e la stessa tabella dei
  campioni:
  - (a) annidamenti (A annidato in B se ogni livello di A compare in un solo livello di B,
    con fattori non costanti e non identificativi di riga);
  - (b) coincidenze uno-a-uno;
  - (c) esiti determinati da un fattore;
  - (d) fattori tecnici in corrispondenza 1:1 con la coppia paziente-tessuto;
  - (e) per ogni confronto, il numero di pazienti con entrambi i livelli.

  Successo solo se l'insieme dei fatti strutturali dell'audit e quello dello script
  coincidono esattamente: nessun fatto in piu' ne' in meno da nessuna delle due parti, e
  stesso numero di pazienti per ogni confronto.
- **D3.** Le classi dei confronti (stimabile / stimabile con bassa potenza / non stimabile)
  coincidono con quelle ottenute applicando a mano, nello script indipendente, le regole
  dichiarate nel README e in `REPORT_INTERVENTI.md`:
  - almeno un paziente con entrambi i livelli: disegno appaiato, unita' = quei pazienti;
    "bassa potenza" se le unita' sono meno di 5 OPPURE se 2/2^unita' > 0.05;
  - nessun paziente con entrambi i livelli e fattore = tessuto o tecnico: "non stimabile";
  - fattore = esito, fra pazienti: unita' = pazienti del gruppo piu' piccolo; "bassa
    potenza" se sono meno di 5 OPPURE se 2/C(n_a+n_b, n_a) > 0.05;
  - in ogni caso, un fattore tecnico non identificativo di riga, con almeno 2 livelli, che
    separa perfettamente i due livelli: "non stimabile".
- **D4.** Nessun valore di Cramér V e' prodotto per tabelle che le regole dichiarano "non
  valutabili", cioe' con meno di 10 righe o meno di 2 righe attese per cella (n/(r·c) < 2).
  Per tutte le altre coppie viene prodotto un numero.

## 3. Criteri di successo — Modulo A

- **A1.** `cli.py leakage` termina senza errori sul dataset esterno, con i nomi di colonna di
  quel dataset.
- **A2.** Uno script indipendente (`validation_esterna/indipendenti/leakage_check.py`,
  scikit-learn, NESSUN import da `core/`) riesegue la valutazione raggruppata per paziente
  con le stesse regole dichiarate dal modulo:
  - StratifiedGroupKFold a min(5, n pazienti) fold, shuffle, seme 0;
  - CPM a 1e4 + log1p, StandardScaler senza centratura;
  - LogisticRegression(C=1, class_weight="balanced", max_iter=1000, seme 0);
  - macro-F1.

  Confronta la macro-F1 media sui fold raggruppati per paziente e, se il modulo esegue il
  confronto LeaveOneGroupOut (≥ 8 pazienti), la macro-F1 media per paziente della
  regressione logistica. Successo se ciascuna differenza assoluta e' ≤ 0.02. Se una
  differenza supera 0.02, la causa viene individuata prima di concludere.
- **A3.** Controllo negativo su dati reali: 20 permutazioni globali delle etichette (semi
  0-19). Valore atteso per etichette casuali: 1/K, cioe' la macro-F1 attesa di un predittore
  indipendente dalla verita' che predice con le proporzioni di classe, dove F1 della classe
  k = p_k.

  Successo se entrambe le condizioni valgono:
  - (i) la media sulle permutazioni della macro-F1 raggruppata per paziente sta entro ±0.05
    da 1/K;
  - (ii) il divario casuale-vs-raggruppato non risulta "leakage rilevabile" secondo il
    criterio del modulo (soglia del report: divario > 0.05) in piu' di 1 permutazione su 20.
- **A4.** Robustezza di formato. Il modulo deve gestire correttamente, oppure rifiutare con
  un messaggio chiaro (ValueError con testo esplicito, mai un traceback non gestito ne'
  numeri silenziosamente sbagliati), queste varianti dello stesso dataset:
  - (a) matrice sparsa contro densa: risultati identici;
  - (b) colonne `obs` categoriche contro stringhe: risultati identici;
  - (c) una classe assente da almeno un fold di test: la macro-F1 del fold non deve contare
    come F1 = 0 una classe che nel fold non e' presente e non e' predetta, senza
    dichiararlo. In alternativa il modulo deve rifiutare o avvisare esplicitamente.

  Se (c) rivela un difetto nel CALCOLO della metrica (e non di lettura o formato), mi fermo
  e lo riporto senza correggerlo, come richiesto.

## 4. Verdetto

"PRONTO" solo se tutti i criteri (D1-D4 per ogni dataset dell'audit, A1-A4 per ogni dataset
del Modulo A) sono soddisfatti su tutti i dataset eleggibili testati. Altrimenti "NON
PRONTO", con l'elenco dei criteri falliti. Non esistono verdetti parziali favorevoli.

## 5. Modifica del piano (2026-10-02), decisa dall'utente DOPO l'esecuzione sul dataset completo

Questa sezione e' stata aggiunta dopo aver visto dei risultati, ed e' quindi dichiarata come
tale. I criteri delle sezioni 1-4 restano invariati.

**Motivo.** `cli.py leakage` su GSE125449 completo (9.946 cellule, 18.372 geni, 19
pazienti) non era terminato dopo oltre 5 ore e 18 ore di CPU. Era rallentato anche dalle mie
esecuzioni in parallelo. Il confronto LeaveOneGroupOut addestra 3 modelli (tra cui random
forest e gradient boosting su matrice densa) su tutti i geni, per ciascuno dei 19 pazienti.
L'utente ha chiesto di testare il Modulo A su una piccola coorte.

**Coorte ridotta del Modulo A** (`validation_esterna/build_small.py`), regole fissate prima
di costruirla:
- 10 pazienti estratti a caso fra i 19 di GSE125449 (numpy default_rng(0),
  scelta senza reinserimento);
- per ciascun paziente al massimo 200 cellule estratte a caso (default_rng(0)), tutte se
  sono meno di 200;
- nessun filtro sui geni ne' sulle classi, etichetta `Type` degli autori invariata.

Su questa coorte si valutano A1-A4 con gli stessi criteri della sezione 3. Il tempo
dell'esecuzione sul dataset completo e' riportato come risultato. L'esecuzione completa resta
attiva a priorita' minima: se termina, i suoi numeri saranno riportati. L'audit del disegno
non cambia: e' gia' stato valutato sui dataset completi.
