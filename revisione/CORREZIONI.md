# Correzioni richieste dalla revisione esterna

Copia dell'elenco «Correzioni» del prompt di revisione, passato al revisore indipendente a ogni giro.

### Blocco A — Errori e calibrazioni

**A1. Cluster bootstrap con pochi pazienti** (`core/stats.py:102-142`). Il test di calibrazione usa 25 gruppi, ma lo strumento lo usa con 5-11 pazienti. Una simulazione sotto H0 (400 repliche, dimensioni dei gruppi casuali fra 3 e 14) dà il 16% di falsi positivi con 5 gruppi e l'11% con 10, fuori dalla banda 2-8%. Il bootstrap percentile con pochi cluster è noto per essere anticonservativo.
- Aggiungi calibrazioni a 5, 8, 10 e 15 gruppi, con dimensioni diseguali.
- Sostituisci il metodo con uno che entra nella banda in questo regime: intervallo t sulle medie per paziente con G − 1 gradi di libertà (attento alla pesatura per cellule), bootstrap-t, BCa o wild cluster bootstrap. Sceglilo in base alla calibrazione.
- Se nessuno è calibrato a 5 gruppi, alza `MIN_PATIENTS_FOR_CI_DEFAULT` alla soglia più bassa calibrata e scrivilo.
- Riesegui i risultati reali e riporta i nuovi intervalli (per esempio 0,195 [0,123; 0,274]).

**A2. Nadeau-Bengio testato fuori dal suo regime d'uso** (`tests/test_stats_calibration.py`, `core/leakage_audit.py:288-305`, `core/stats.py:42-44`).
- Lo strumento usa NB solo su LeaveOneGroupOut, con k ≥ 8 fold diseguali e `mean(n_test)/mean(n_train)`; la calibrazione invece è a k = 3 e k = 5.
- Lo scenario «nullo» usa due `DummyClassifier`, le cui differenze fra fold non sono correlate: il fenomeno che NB corregge non c'è.
- Scrivi una calibrazione con due classificatori veri e ugualmente buoni, addestrati sugli stessi training sovrapposti, in LOGO con 8-20 gruppi diseguali.
- Correggi o rimuovi il test tautologico `test_nadeau_bengio_less_anticonservative_than_uncorrected_ttest`.
- Correggi la spiegazione dello xfail a k = 5: la causa non sono i gradi di libertà.
- Correggi il README, che dice che k = 5 è «il default reale usato da leakage_audit.py»: è falso.
- Documenta in `stats.py` che il rapporto è fra medie.

**A3. Confronto fra modelli** (`leakage_audit.py:288`). Il migliore è scelto per media e poi testato contro gli altri senza correzione per confronti multipli. Usa un modello di riferimento fissato prima (la regressione logistica) contro gli altri, oppure applica Holm a tutti i confronti, e dichiaralo nel testo. Il fallback all'approssimazione normale di scipy con ties o zeri va dichiarato nella nota sul p minimo di Wilcoxon.

**A4. Modulo B, filtro del tasso per cellula** (`core/tcr_validation.py:336`, `:394`, `:511`). Oggi tiene solo le etichette fra le identità di riferimento *risolte*, non quelle della `marker_map`. Se nessun clone si risolve in una delle identità, tutte le cellule con quell'etichetta, cioè proprio gli errori, escono dal calcolo. Esempio sintetico: per clone 0,40, per cellula 0,000 su 277 cellule.
- Filtra su `list(marker_map)`.
- Aggiungi un test con questo caso.
- Aggiorna `tests/fixtures/tcr_regression_baseline.json`, documentando nel commit perché cambia.
- Riesegui la validazione su GSE278694 e riporta il nuovo valore per cellula (oggi 0,093).
- Allinea `FLAG_DEFINITIONS` al codice.

**A5. Eccesso di discordanza** (`tcr_validation.py:574-576`). Il testo dice che un errore di annotazione «produrrebbe eccesso ≥ 0»: è falso quando n_A ≠ n_B, perché lo stimatore pesato ha media zero sotto H0 ma non è una divergenza non negativa (esempio: 3 cellule CD8 nel sangue contro 60 nel tumore con il 40% di errore danno −0,065).
- Correggi il testo e la docstring del modulo.
- Aggiungi al generatore sintetico casi sbilanciati e cloni non condivisi.
- In `verdict.py`, un IC interamente sotto lo zero non va in verde, e un IC che include lo zero mostra «nessun eccesso rilevato», non «stima affidabile».
- Valuta, come sensibilità, la versione a pesi uguali ½‖pA − pB‖², che è non negativa e non distorta.

**A6. Sintesi dei controlli** (`core/verdict.py:117-127`). Una coppia annidata o 1:1 con V ≥ 0,5 non fa scattare l'allarme e finisce fra le «valutabili sotto soglia», cioè in verde (riprodotto: batch annidato nel tessuto, n = 20, V = 0,94 → verde).
- Escludi dalle valutabili le coppie con fatti strutturali.
- Fai entrare ogni fatto strutturale nella sintesi come controllo giallo o rosso.
- Un livello mancante (refuso) va in grigio, come errore d'input, non in rosso.
- Riscrivi `tests/test_verdict.py` con oggetti `DesignAuditResult` veri.
- Cambia l'etichetta verde da «stima affidabile» a «nessun limite rilevato dalle regole», in codice, report, app e README.

**A7. Frazione di CD8.** Sui dati reali è inutilizzabile: P(chiamata CD8 | vera CD4) = 0,45, 8 pazienti su 13 vengono troncati a 0, gli scenari 0,5× e 1× danno intervalli disgiunti, e la protezione su J (0,48) non intercetta il problema.
- Toglila dall'app e dal report HTML/Markdown.
- Toglila dalla Sintesi dei controlli.
- Lasciala solo in CLI dietro un flag esplicito `--sperimentale-cd8`.
- Quando la frazione osservata cade fuori dall'intervallo invertibile, rifiuta invece di troncare.
- Togli dall'app il titolo «Frazione di CD8 validata dal repertorio».
- Aggiorna README e test.

**A8. Errori silenziosi d'input.**
- `core/io.py:57`: `drop_duplicates(barcode_col)` tiene la prima riga per barcode duplicato. Deve dare errore.
- `match_barcodes` e `parse_vdj_contigs:69`: chiavi (paziente, compartimento, barcode) duplicate, oppure più file VDJ sulla stessa coppia paziente-compartimento senza dimensione di libreria. Devono dare errore, oppure gestire una colonna libreria.
- `cli.py:186` e `app.py:466`: `obs.join(cell_flags)` con obs_names non unici. Deve dare errore, oppure assegnare per posizione.
- `app.py:468-469`: le etichette CD4/CD8 sono prese dall'ordine delle chiavi della mappa; servono selettori espliciti, se la funzione resta nell'app.
- `compute_marker_positivity`: verifica che X contenga conteggi interi non negativi.
- `clone_level_units:440`: pareggi della maggioranza decisi dall'ordine delle righe. Rendi la regola deterministica e documentala.
- App:
  - `app.py:206`: righe di confronto malformate scartate in silenzio; vanno segnalate.
  - `app.py:401-403`: JSON e tabelle illeggibili danno un traceback grezzo; serve un messaggio chiaro.
- CLI:
  - un fattore in `--compare` non dichiarato deve dare un messaggio chiaro;
  - un ruolo ripetuto (due `batch=`) non deve sovrascrivere in silenzio;
  - colonne ripetute in `--h5ad` non devono andare in crash.

**A9. Audit del disegno su input per cellula.** Con `--meta` o nell'app un CSV con una riga per cellula viene analizzato come se ogni cellula fosse un campione. Avvisa, o rifiuta, quando le righe sono molte più dei pazienti. Con `--h5ad` senza colonna libreria, le librerie replicate collassano: avvisa. In `design_audit.py:400-402` la spiegazione cita il pavimento di Wilcoxon anche quando è sotto 0,05; correggila.

**A10. App e privacy.**
- `app.py:135-137`: l'`.h5ad` caricato viene scritto in un file temporaneo mai cancellato, e riscritto a ogni clic. Va scritto una volta per caricamento e cancellato dopo la lettura.
- Scrivi nel README che l'app è per una workstation con un solo utente: «localhost» non protegge da altri utenti della stessa macchina o di un server condiviso.
- Verifica la versione minima di Streamlit che supporta `width="stretch"` e aggiorna `pyproject.toml`.

**A11. Report.**
- Il report HTML deve stampare i limiti dichiarati (`standing_limits`).
- `report.py:191`: usa `GAP_ALERT`, non 0,05 fisso.
- `report.py:218`: stampa la soglia, non il numero di pazienti.
- Il piè di pagina parla di «due moduli»: le sezioni sono quattro.
- Allinea i colori della tabella dei confronti a quelli della sintesi.
- `report.py:402`: «fuori dalla marker_map», non «non CD4/CD8».

### Blocco B — Affermazioni, biologia, README

**B1. Ipotesi biologiche del Modulo B.** Dichiarale nel README, nel report e nella docstring di `tcr_validation.py`. Il metodo assume che lo stesso clone abbia la stessa linea CD4/CD8; eccezioni note:
- cellule CD4 con CD8αα;
- CD4 citotossiche, probabile causa reale del 45%;
- MAIT, NKT e cellule γδ, non catturate da un clonotipo basato solo su TRB;
- doppietti T-NK;
- RNA ambientale specifico del compartimento.

L'eccesso di discordanza è interpretabile solo con etichette di linea (CD4/CD8), non di stato (per esempio effettrici contro esauste). Togli «un effetto reale» e «segnale di errore sistematico» dove non è dimostrato. «Tipicamente più frequente nel tumore» diventa «nel dataset GSE278694».

**B2. Etichette NK e convenzione per clone.** Scrivi che il default per clone conta come errore *ogni* etichetta fuori mappa, anche mieloide e B. Su GSE278694 le etichette NK sono circa il 43% della differenza fra la convenzione per cellula e quella per clone, e la tesi attribuisce il 72-79% delle NK con TCR a doppietti.

**B3. «Validazione».** La riproduzione dei numeri della tesi sullo stesso dataset è una verifica di fedeltà del porting, non una validazione: rinominala così nel README. Correggi «validato su un solo dataset» e «già validati» (README:107, :115, :409-412, `verdict.standing_limits`). Gli script «indipendenti» verificano che il codice segua le sue regole, non che le regole siano giuste: scrivilo.

**B4. Intervalli con 4 pazienti.** Negli output di validazione e in `docs/sviluppo/REPORT_INTERVENTI.md`, gli intervalli con 4 pazienti (per esempio [+0,021; +0,170]) vanno annotati: lo strumento li rifiuta, e con 4 pazienti esistono solo 35 ricampionamenti distinti.

**B5. Modulo A, formulazioni.**
- «Lo split casuale sottostima l'incertezza reale di quel fattore»: la deviazione standard fra 5 fold non è una stima valida dell'incertezza. Riformula.
- Su GSE125449 la piattaforma coincide con il paziente: il divario misura anche il passaggio fra piattaforme. Scrivilo.
- Non indicare una direzione («sovrastima») per divari sotto una soglia di rumore dichiarata.
- La demo da CLI genera un'etichetta casuale a livello di paziente con un'impronta per paziente, quindi il divario è garantito per costruzione: dillo nell'output della demo e nel README.
- «Jaccard descrittivo»: togli la regola di lettura inferenziale dalla docstring di `XaiStability`.

**B6. Audit del disegno, formulazioni.**
- «Potenza 0,46 a V = 0,5»: indica le condizioni (tabella 3×3, n = 60, una sola configurazione simulata).
- La regola di valutabilità è sulla media delle righe attese per cella, non sul minimo: correggi il testo.
- «Stimabile con bassa potenza» con 5 coppie significa potenza nulla al 5%: dillo.
- Scrivi che i fatti strutturali sono controlli deterministici, utili come lista di verifica; l'unico elemento inferenziale è il Cramér V.

**B7. README.**
- Citazioni complete: GSE278694 = Chen M-M et al., *Cancer Cell* 2025;43(9):1656-1676.e10, doi:10.1016/j.ccell.2025.06.020, sia nell'estratto sia nell'elenco; GSE125449 = Ma L. et al., *Cancer Cell* 2019;36(4):418-430.e6.
- Esiti dei test misurati (Fase 3).
- Sezione «Ipotesi e limiti», che raccolga B1-B6.
- Tabella della Sintesi dei controlli allineata ad A5, A6 e A7.

**B8. Commenti.**
- Docstring per `_compare`, `_pair_sentence`, `_partial_confounders`, le funzioni `*_summary` e `cap_training` (con l'arrotondamento che può dare qualche cellula in più o in meno).
- Togli il parametro `name` inutilizzato di `read_table`.
- `convention_definition` e `_paired_differences` non devono avere «3 cellule» e «Tumor» fissi.
- `leakage_check.py:18`: le metriche elencate sono tre, non due.
- `make_tcr_validation_dataset`: documenta che tutti i cloni sono condivisi e che i marcatori separano perfettamente, e aggiungi una variante con dropout e RNA ambientale.
- `verdict.py`: aggiorna la docstring dopo A6.
- `sample_sheet_from_obs`: aggiungi l'avviso sulle librerie che collassano.

### Blocco C — Materiale per chi lo userà

**C1. Guida per l'utente** (`docs/GUIDA_UTENTE.md`, in italiano, per un biologo o un clinico, non per uno sviluppatore):
- cosa serve in ingresso per ciascuna parte: per l'audit del disegno un foglio campioni; per il Modulo A un `.h5ad` con conteggi grezzi, paziente ed etichetta; per il Modulo B anche i file VDJ di Cell Ranger, un manifest e un compartimento di riferimento;
- come convertire un oggetto Seurat (`.rds`) in `.h5ad`, con i comandi R di zellkonverter o SeuratDisk, e quali colonne controllare dopo la conversione. Molti gruppi clinici lavorano in R;
- come leggere ogni riga della «Sintesi dei controlli», con un esempio per stato, e cosa **non** significa il verde;
- i limiti di B1-B6 in forma breve.

**C2. Prova su dati non visti.** Esegui la CLI dall'inizio alla fine, seguendo solo la guida, su un dataset che lo sviluppo non ha mai usato, per esempio un sottoinsieme di GSE131907 o di GSE132465 per audit del disegno e Modulo A. Riporta nel resoconto ogni punto in cui la guida non bastava, e correggilo.
