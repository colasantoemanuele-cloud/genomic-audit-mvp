# genomic-audit

Prototipo di strumento di audit per dati genomici a coorte piccola (10-100 pazienti,
trascrittomica a singola cellula). Pensato per un pilota con un singolo laboratorio
partner: applicazione locale, non multi-tenant, senza autenticazione ne' pagamenti.

## Cosa fa

**Audit del disegno e del confondimento (solo metadati).** Legge una tabella con una riga
per campione o libreria (oppure `adata.obs`, ridotto alle combinazioni distinte dei fattori)
e l'indicazione del ruolo delle colonne: paziente, tessuto, fattori tecnici (libreria,
batch, chimica, protocollo, data), esiti. Riporta:
(1) i fatti strutturali del disegno, rilevati in modo deterministico: annidamenti, fattori
coincidenti, esiti determinati da un fattore, fattori tecnici in corrispondenza 1:1 con la
coppia paziente-tessuto; per ogni coppia di fattori riporta anche il Cramér V con la
correzione di Bergsma, "non valutabile" su tabelle troppo piccole, con allarme a V >= 0.5;
(2) per ogni confronto richiesto, il numero di unita' indipendenti (pazienti) e la classe
`stimabile` / `stimabile con bassa potenza` / `non stimabile`.
Funziona anche prima di sequenziare, sul foglio di disegno dei campioni.

**Modulo A -- Audit del leakage per paziente.** Per un task di classificazione a scelta
dell'utente, confronta una valutazione onesta (split a 5 fold raggruppato per paziente,
`StratifiedGroupKFold`) con un controllo negativo (split casuale sulle singole cellule,
che ignora il paziente e quindi lascia "gemelle" della stessa cellula sia in train sia in
test). Il divario fra le due misura quanto una pipeline che non raggruppa per paziente
sovrastimerebbe l'accuratezza e sottostimerebbe l'incertezza. Con almeno 8 pazienti,
ripete il confronto con `LeaveOneGroupOut` su tre modelli (regressione logistica, random
forest, gradient boosting), con la correzione di Nadeau-Bengio (anticonservativa la
versione non corretta, perche' i training set dei fold si sovrappongono) e il test di
Wilcoxon appaiato, segnalando quando le due conclusioni divergono.

Motore pensato per dati reali: la matrice resta sparsa; in ogni fold, e stimati sul solo
training, si calcolano log-CPM, i 2.000 geni piu' variabili e (per random forest e gradient
boosting) 50 componenti SVD. Il training e' limitato a 20.000 cellule per fold con
sottocampionamento stratificato per paziente, dichiarato. La macro-F1 e' calcolata sulle sole
classi presenti nel fold di test, e le classi assenti sono dichiarate. I fold della
regressione logistica girano in parallelo, con risultati identici a quelli sequenziali.
Tempi misurati su questa workstation (12 core): GSE125449 reale (9.946 cellule, 18.372 geni,
19 pazienti) 6,7 minuti; dati sintetici con 25.000 cellule, 20.000 geni e 20 pazienti 6,2
minuti, misurati prima della parallelizzazione dei fold. In entrambi i casi il confronto fra
modelli e' incluso. `--rapido` lo esclude. Mostra anche la stabilita' delle spiegazioni
(Jaccard fra i 50 geni principali dei fold), una misura descrittiva.

**Modulo B -- Validazione dell'annotazione via TCR.** Usa il repertorio T-cell receptor
come ancora indipendente dal clustering: se le cellule dello stesso clone T (stessa
sequenza CDR3 della catena TRB) ricevono etichette di tipo cellulare diverse a seconda
del compartimento tissutale in cui si trovano, quell'eccesso di discordanza (oltre il
rumore di base entro compartimento) e' un segnale di errore sistematico
nell'annotazione -- tipicamente piu' frequente nel tessuto tumorale che nel sangue.
L'intervallo di confidenza usa il cluster bootstrap sui PAZIENTI (non sui cloni, che sono
annidati nei pazienti e non sono osservazioni indipendenti). Se fornita una mappa
etichetta -> geni marcatori canonici, calcola anche il tasso d'errore per compartimento
usando come riferimento un solo compartimento (tipicamente il sangue), per evitare la
circolarita' di stimare il riferimento sugli stessi dati che poi si giudicano.

**Modulo B -- Convenzioni del tasso d'errore.** Il tasso d'errore rispetto al riferimento
dal sangue e' riportato con due convenzioni affiancate, ciascuna con nome e definizione
(`--convention both|cell|clone`, default `both`):
- *per cellula* (etichette fuori mappa escluse);
- *per clone* (un'unita' per coppia clone-compartimento con almeno 3 cellule, etichetta di
  maggioranza, etichette fuori mappa contate come errore).

Con `--clone-error-labels NK --clone-marker-priority CD8T,CD4T` la convenzione per clone
riproduce esattamente i numeri della tesi su GSE278694 (tumore 0,1951 contro 0,093 per
cellula). Prima di ogni calcolo il Modulo B riporta la frazione di barcode VDJ ritrovati nei
metadati e normalizza, dichiarandolo, il suffisso "-N" quando serve. Sotto il 50% di match
si ferma con un errore.

**Modulo B -- Flag per cellula.** Con la mappa dei marcatori e il compartimento di
riferimento, il Modulo B scrive in una COPIA dell'AnnData (`<nome>_audited.h5ad`) e in un
CSV due colonne per cellula: `audit_reference_label`, l'identita' del clone stimata nel
sangue, e `audit_label_vs_reference`, True/False se l'etichetta assegnata discorda o
concorda. Il valore e' `NA` dove la cellula non era verificabile. La media dei flag
valutabili di un compartimento coincide esattamente con il tasso d'errore del Modulo B.

**Modulo B -- Frazione di CD8 con l'errore propagato (SPERIMENTALE).** Per ogni paziente, nel
compartimento scelto: frazione di cellule etichettate CD8 sul totale delle CD4+CD8, e
intervallo plausibile al 95% dopo l'inversione di una matrice di confusione con direzione
(vera CD4/CD8 -> chiamata CD4/CD8/altro), stimata sulle cellule con identita' di
riferimento e pooled fra pazienti. Tre scenari di errore sono sempre riportati (0.5x, 1x,
2x), senza indicarne uno come "il risultato". Nessun intervallo con meno di 5 pazienti con
riferimento o con una matrice mal condizionata (J < 0.2, oppure un intervallo di J che
include lo zero).

## Cosa NON fa

- Nessuna autenticazione, nessun utente multiplo, nessuna fatturazione, nessun deployment
  cloud, nessun database persistente. Gira in locale, i dati non lasciano la macchina.
- Non e' una validazione biologica indipendente: il Modulo B usa il TCR come ancora, ma
  resta un segnale statistico, non una conferma sperimentale (es. citometria, sorting).
- Con pochi fold (tipico di coorti piccole) la potenza statistica e' bassa: "non
  significativo" non vuol dire "equivalente". Lo strumento lo segnala esplicitamente
  invece di nasconderlo.
- Il Modulo B richiede che il dataset abbia una colonna di compartimento tissutale e file
  VDJ Cell Ranger; senza quelli, restano disponibili solo l'audit del disegno e il Modulo A.
- **Audit del disegno:** non guarda l'espressione genica e non stima effetti biologici.
  La scomposizione della varianza su pseudobulk (paziente/tessuto/batch) esiste in
  `core/design_audit.py` ma NON e' esposta: la calibrazione del suo intervallo e' fuori
  banda per la quota del tessuto (vedi `NOTE.md`). Il Cramér V a soglia 0.5 ha potenza
  0.46 a V = 0.5: un'associazione moderata spesso non viene segnalata. Non verifica che una
  colonna dichiarata come esito sia davvero a livello di paziente.
- **Flag per cellula:** non correggono nulla, e le etichette originali non vengono mai
  modificate. Non esiste un flag dei doppietti: lo fanno gia' scDblFinder e Scrublet. Le
  cellule senza TCR, i cloni senza cellule sufficienti nel sangue, le etichette fuori dalla
  mappa dei marcatori e le cellule del compartimento di riferimento restano `NA`, e `NA`
  non significa "corretta".
- **Frazione di CD8:** assume che l'identita' dal sangue sia corretta e che i cloni
  condivisi con il sangue (quelli espansi) siano rappresentativi dei non condivisi.
  Quest'ultima assunzione non e' verificabile: nei dati PDAC originali l'errore cresce con
  la dimensione del clone. La matrice e' unica per tutti i pazienti. Doppietti e cellule
  non-T etichettate CD4/CD8 non sono modellati.
  **Nessuno dei tre scenari garantisce la copertura nominale in tutti i casi** (misura B3,
  `tests/test_cd8_robustness.py`, 2000 intervalli per caso, errore asimmetrico CD4->CD8 0.10,
  CD8->CD4 0.03, altro 0.04). Copertura dell'IC 95%:
  - errore dei cloni condivisi uguale a quello degli altri: 0.5x 0.878; 1x 0.944; 2x 0.755;
  - errore dei cloni condivisi doppio: 0.5x 0.923; 1x 0.839; 2x 0.270;
  - errore dei cloni condivisi quadruplo: 0.5x 0.945; 1x 0.285; 2x mai calcolabile.

  Lo scenario che copre il valore vero dipende quindi da quanto i cloni condivisi sono
  rappresentativi, cosa che sui dati non si puo' verificare. Sui dati reali di GSE278694
  la matrice stimata nel tumore ha P(chiamata CD8 | vera CD4) = 0.45, un errore molto
  piu' alto di quelli simulati: li' lo scenario 2x non e' calcolabile e l'intervallo 1x
  arriva a 0 in 8 pazienti su 13.

## Installazione

Richiede [uv](https://docs.astral.sh/uv/).

```bash
uv sync --extra dev
```

## Web app locale (demo con un clic)

```bash
uv run audit-sc serve          # apre http://localhost:8501
```

L'app ascolta solo su `localhost`, con la telemetria di Streamlit disattivata: nessun dato e
nessuna statistica d'uso lasciano la macchina. Nella barra laterale ci sono due modalita':
- **Demo immediata:** un clic carica metadati reali GEO per il disegno e un sottoinsieme
  reale di GSE125449 per il Modulo A (`data/demo/`); il Modulo B usa dati sintetici
  dichiarati;
- **Carica studio:** percorsi locali o upload di CSV/TSV, `.h5ad`, cartelle 10x Matrix
  Market e manifest VDJ.

Le sezioni sono Disegno, Modulo A, Modulo B e Verdetto. Il verdetto a semaforo usa regole
dichiarate in `core/verdict.py` ed e' esportabile in Markdown o HTML.

## Demo da riga di comando

```bash
uv run python cli.py demo --out results/demo_report.html
```

La demo include l'audit del disegno su metadati sintetici con la struttura di GSE278694.

La demo genera due dataset sintetici (uno per modulo, per mostrare chiaramente l'effetto
che ciascuno misura): per il Modulo A, un "fingerprint" genico casuale specifico per
paziente che uno split casuale sulle cellule puo' sfruttare ma uno split per paziente no;
per il Modulo B, cloni T con discordanza cross-compartimento iniettata deliberatamente in
un compartimento "rumoroso", con l'espressione dei geni marcatori legata all'identita'
VERA del clone (non all'etichetta assegnata, che puo' essere sbagliata).

## Uso con dati propri

```bash
# Audit del disegno (solo metadati: CSV con una riga per campione, oppure --h5ad)
uv run python cli.py design --meta campioni.csv --patient-col patient --tissue-col tissue \
    --technical batch=run --technical protocol=protocol \
    --compare tissue:Tumor:Adjacent_normal --out results/design_report.html

# Modulo A (--rapido: senza confronto fra modelli; --max-train-cells: limite per fold)
uv run python cli.py leakage --h5ad dati.h5ad \
    --target-col tissue --patient-col patient_id --out results/leakage_report.html

# Modulo B: serve un manifest CSV (path,patient,compartment) per i file VDJ Cell Ranger,
# perche' quei CSV non contengono queste informazioni e non c'e' una convenzione di nome
# file universale per dedurle.
uv run python cli.py tcr --h5ad dati.h5ad --vdj-manifest vdj_manifest.csv \
    --patient-col patient_id --compartment-col tissue --celltype-col celltype \
    --barcode-col barcode --marker-map markers.json --reference-compartment PBMC \
    --export-flags --cd8-compartment Tumor \
    --out results/tcr_report.html
# --export-flags  -> results/dati_audited.h5ad (copia) e results/dati_audit_flags.csv
# --cd8-compartment -> intervalli sulla frazione di CD8 per paziente nel report
```

Oppure `uv run audit-sc serve`, modalita' "Carica studio", per fare lo stesso
dall'interfaccia web. Tutti i comandi sono disponibili anche come `uv run audit-sc <comando>`.

## Test

```bash
uv run pytest
```

`tests/test_stats_calibration.py` e' **obbligatoria e va eseguita per prima**: verifica
su dati simulati (>=200 repliche, semi fissi -- deterministica, non flaky) che gli
stimatori statistici (correzione di Nadeau-Bengio, cluster bootstrap) siano calibrati
sotto l'ipotesi nulla (2-8% di falsi positivi a soglia 0.05) e abbiano potenza sotto un
effetto reale iniettato. Se un test di calibrazione fallisce, lo stimatore va corretto
prima di usare l'interfaccia -- non e' stato aggirato per far passare la build in nessun
punto di questo repository. Un test (`test_nadeau_bengio_calibration_null_k5_default`) e'
marcato `xfail(strict=True)`: a k=5 (il default reale usato da `leakage_audit.py`) la
correzione e' conservativa (falsi positivi ~1.8-1.9%, misurato, non dedotto), sotto la
banda teorica [2%, 8%] ma sul lato sicuro. La banda resta quella teorica, non e' stata
allargata per far passare il numero osservato: se in futuro questo test tornasse a
passare inaspettatamente, la suite fallirebbe (segnale da investigare).

Altre calibrazioni obbligatorie, con le bande dichiarate nel docstring di ciascun file:
`tests/test_design_audit_calibration.py` (falsi allarmi del Cramér V, annidamento e
confondimento rilevati al 100%, caso GSE278694; la copertura dell'intervallo della quota
del tessuto e' `xfail(strict=True)`) e `tests/test_cd8_propagation_calibration.py`
(copertura dell'IC 95% della frazione di CD8 con errore simmetrico e asimmetrico).
`tests/test_tcr_flags.py` verifica che i flag coincidano esattamente con il tasso d'errore
e che l'output del Modulo B sia identico a quello precedente alla modifica
(`tests/fixtures/tcr_regression_baseline.json`). Il resoconto degli interventi, con tutti i
numeri misurati, e' in `REPORT_INTERVENTI.md`.

`tests/test_synthetic_data.py` verifica invece, end-to-end, che gli effetti iniettati nei
generatori sintetici (`core/synthetic.py`) vengano effettivamente rilevati dai due moduli.

## Struttura

```
core/               motore analitico, installabile e testabile senza Streamlit
  design_audit.py   audit del disegno e del confondimento
  leakage_audit.py  Modulo A
  tcr_validation.py Modulo B (inclusi i flag per cellula)
  cd8_propagation.py frazione di CD8 con propagazione dell'errore
  stats.py          Nadeau-Bengio, cluster bootstrap (validati con calibrazione)
  synthetic.py       generatori di dati sintetici (demo + test)
  report.py         report HTML autocontenuto
tests/
  test_stats_calibration.py   OBBLIGATORIA, eseguita per prima
  test_synthetic_data.py
app.py              interfaccia Streamlit (file sottile)
cli.py              esecuzione da riga di comando
data/synthetic/     output della demo (file VDJ di esempio)
```

## Provenienza della logica statistica

La correzione di Nadeau-Bengio e il cluster bootstrap sui pazienti, cosi' come lo
stimatore di discordanza a coppie del Modulo B, sono portati (generalizzati: nomi di
colonna configurabili, non piu' legati a un dataset specifico) da un progetto precedente
di ricerca su un dataset scRNA-seq + TCR-seq di adenocarcinoma duttale pancreatico, dove
erano gia' stati validati. Non sono stati reimplementati da zero.
