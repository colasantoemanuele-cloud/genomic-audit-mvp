# genomic-audit

Prototipo di strumento di audit per dati genomici a coorte piccola (10-100 pazienti,
trascrittomica a singola cellula). Pensato per un pilota con un singolo laboratorio
partner: applicazione locale, non multi-tenant, senza autenticazione ne' pagamenti.

## Cosa fa

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

## Cosa NON fa

- Nessuna autenticazione, nessun utente multiplo, nessuna fatturazione, nessun deployment
  cloud, nessun database persistente. Gira in locale, i dati non lasciano la macchina.
- Non e' una validazione biologica indipendente: il Modulo B usa il TCR come ancora, ma
  resta un segnale statistico, non una conferma sperimentale (es. citometria, sorting).
- Con pochi fold (tipico di coorti piccole) la potenza statistica e' bassa: "non
  significativo" non vuol dire "equivalente". Lo strumento lo segnala esplicitamente
  invece di nasconderlo.
- Il Modulo B richiede che il dataset abbia una colonna di compartimento tissutale e file
  VDJ Cell Ranger; senza quelli, resta disponibile solo il Modulo A.

## Installazione

Richiede [uv](https://docs.astral.sh/uv/).

```bash
uv sync --extra dev
```

## Demo (senza dati propri)

```bash
# interfaccia web
uv run streamlit run app.py
# spunta "Usa dati sintetici di esempio" per una demo senza upload

# oppure da riga di comando
uv run python cli.py demo --out results/demo_report.html
```

La demo genera due dataset sintetici (uno per modulo, per mostrare chiaramente l'effetto
che ciascuno misura): per il Modulo A, un "fingerprint" genico casuale specifico per
paziente che uno split casuale sulle cellule puo' sfruttare ma uno split per paziente no;
per il Modulo B, cloni T con discordanza cross-compartimento iniettata deliberatamente in
un compartimento "rumoroso", con l'espressione dei geni marcatori legata all'identita'
VERA del clone (non all'etichetta assegnata, che puo' essere sbagliata).

## Uso con dati propri

```bash
# Modulo A
uv run python cli.py leakage --h5ad dati.h5ad \
    --target-col tissue --patient-col patient_id --out results/leakage_report.html

# Modulo B: serve un manifest CSV (path,patient,compartment) per i file VDJ Cell Ranger,
# perche' quei CSV non contengono queste informazioni e non c'e' una convenzione di nome
# file universale per dedurle.
uv run python cli.py tcr --h5ad dati.h5ad --vdj-manifest vdj_manifest.csv \
    --patient-col patient_id --compartment-col tissue --celltype-col celltype \
    --barcode-col barcode --marker-map markers.json --reference-compartment PBMC \
    --out results/tcr_report.html
```

Oppure `uv run streamlit run app.py` senza la spunta "dati sintetici", per caricare i
file da interfaccia web (upload guidato, selezione colonne, editor per l'associazione
file VDJ -> paziente/compartimento).

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

`tests/test_synthetic_data.py` verifica invece, end-to-end, che gli effetti iniettati nei
generatori sintetici (`core/synthetic.py`) vengano effettivamente rilevati dai due moduli.

## Struttura

```
core/               motore analitico, installabile e testabile senza Streamlit
  leakage_audit.py  Modulo A
  tcr_validation.py Modulo B
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
