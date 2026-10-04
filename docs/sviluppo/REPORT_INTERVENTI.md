# Resoconto degli interventi su genomic-audit-mvp

Punto di partenza (commit `cf61e76`, baseline): suite **18 passed, 1 xfailed** in 309 s (`logs/pytest_baseline.txt`).
Tutti i numeri qui sotto sono misurati; i log completi di pytest sono in `logs/`.

## Intervento 1 — Audit del disegno e del confondimento

### 1. File modificati e creati

- `core/design_audit.py` (NUOVO): modulo dell'audit del disegno. Contiene il Cramér V con la correzione di Bergsma, i fatti strutturali (annidamento, uno-a-uno, esito determinato, unità tecnica, costante, identificatore), la classificazione dei confronti e la scomposizione della varianza su pseudobulk (implementata e testata, NON esposta).
- `core/synthetic.py` (ESISTENTE, solo aggiunte in coda): generatori di metadati sintetici. Ci sono gli scenari indipendente, annidato, confuso con l'esito e con associazione a V noto, la struttura di GSE278694, le unità con quote di varianza note e un AnnData per la pseudobulk. Le funzioni esistenti non sono state toccate.
- `core/report.py` (ESISTENTE): nuova sezione `_design_section`; `render_report` e `save_report` accettano `design_result` come parametro opzionale in coda, quindi le chiamate esistenti restano compatibili. Una frase del footer è stata estesa.
- `cli.py` (ESISTENTE): nuovo comando `design` (`--meta` CSV oppure `--h5ad`, da cui si usa solo `obs`; `--technical ruolo=colonna`, `--outcome-col`, `--compare colonna:a:b`). La demo esegue anche l'audit del disegno sulla struttura di GSE278694.
- `app.py` (ESISTENTE): nuova scheda "Audit del disegno" e report scaricabile che la include.
- `tests/test_design_audit_calibration.py` (NUOVO): calibrazioni, con le bande dichiarate nel docstring prima dell'esecuzione.
- `tests/test_design_audit.py` (NUOVO): test funzionali sui casi limite.
- `logs/pytest_baseline.txt`, `logs/pytest_intervento1.txt` (NUOVI): output reale della suite.
- `core/stats.py`: NON toccato. `core/tcr_validation.py`: NON toccato in questo intervento.

### 2. Test

Suite completa dopo l'Intervento 1: **36 passed, 2 xfailed** in 476 s (`logs/pytest_intervento1.txt`). Prima: 18 passed, 1 xfailed. I nuovi xfail sono 1: `test_variance_share_bootstrap_coverage_tissue`, strict.

Output reale di `pytest -v -s` per i test nuovi (le righe tra parentesi quadre sono stampate dai test):

```
tests/test_design_audit.py::test_cramers_v_perfect_and_null PASSED
tests/test_design_audit.py::test_small_table_is_not_evaluable PASSED
tests/test_design_audit.py::test_batch_perfectly_separating_tissue_makes_comparison_not_estimable PASSED
tests/test_design_audit.py::test_paired_comparison_with_enough_patients_is_estimable PASSED
tests/test_design_audit.py::test_patient_level_outcome_is_between_patient_comparison PASSED
tests/test_design_audit.py::test_missing_level_is_not_estimable PASSED
tests/test_design_audit.py::test_sample_sheet_from_cell_level_obs PASSED
tests/test_design_audit.py::test_protocol_note_present_only_when_protocol_declared PASSED
[a] falsi allarmi Cramér V: 1/200 = 0.005 (banda (0.0, 0.05)); repliche con fatti strutturali: 0/200
PASSED
[b] annidamento paziente-in-batch rilevato: 200/200
PASSED
[c] esito determinato dal paziente rilevato: 200/200; confronto tessuto classificato 'non stimabile': 200/200
PASSED
[potenza] detection rate (soglia V >= 0.5, n=60, 3x3): V=0.3: 0.020, V=0.5: 0.460, V=0.7: 0.980, V=0.9: 1.000
PASSED
[quote] patient: vero 0.5, media stime 0.500 (bias -0.000), entro +-0.1: 1.000
[quote] tissue: vero 0.3, media stime 0.300 (bias +0.000), entro +-0.1: 1.000
PASSED
[copertura IC95%] patient: 193/200 = 0.965 (banda (0.9, 0.99))
PASSED
[copertura IC95%] tissue: 200/200 = 1.000 (banda (0.9, 0.99))
XFAIL
[pseudobulk] quote stimate: {'patient_id': 0.49802726409817816, 'tissue': 0.317798575636533}
PASSED
tests/test_design_audit_calibration.py::test_variance_decomposition_refuses_confounded_factors PASSED
tests/test_design_audit_calibration.py::test_variance_decomposition_refuses_too_few_patients PASSED
[GSE278694] Confronto 'tissue': Tumor vs Adjacent_normal: stimabile con bassa potenza, disegno appaiato (entro paziente), 5 pazienti con entrambi i livelli (soglia dichiarata: 5 pazienti, la stessa del Modulo B; inoltre con 5 pazienti appaiati anche un test di Wilcoxon esatto non può scendere sotto p = 0.062). Un risultato non significativo qui non indica assenza di effetto. Attenzione: 16 pazienti con un solo livello esclusi dal confronto appaiato.
[GSE278694] Confronto 'protocol': scRNA vs snRNA: non stimabile -- nessun paziente ha entrambi i livelli: la differenza fra 'scRNA' e 'snRNA' coincide con la differenza fra due gruppi di pazienti diversi (fattore confuso con il paziente).
PASSED
```

Numeri delle calibrazioni, con la banda dichiarata nel docstring del file di test prima dell'esecuzione:

| Test | Repliche | Misurato | Banda / requisito | Esito |
|---|---|---|---|---|
| (a) fattori indipendenti, falsi allarmi sul Cramér V (V >= 0.5) | 200 | 1/200 = 0.005 | [0, 0.05] | passa |
| (a) fatti strutturali sotto indipendenza | 200 | 0/200 | 0 | passa |
| (b) annidamento paziente-in-batch iniettato | 200 | 200/200 | 100% | passa |
| (c) esito determinato dal paziente | 200 | 200/200; confronto "non stimabile" 200/200 | 100% | passa |
| Potenza del V (3x3, n=60) | 200 per punto | V=0.3: 0.020; V=0.5: 0.460; V=0.7: 0.980; V=0.9: 1.000 | curva, nessuna banda (solo monotonia) | passa |
| Quote di varianza, paziente (vero 0.5) | 200 | media 0.500 (distorsione -0.000); entro ±0.10: 1.000 | ≥ 0.95 entro ±0.10; distorsione ≤ 0.03 | passa (dopo correzione 1) |
| Quote di varianza, tessuto (vero 0.3) | 200 | media 0.300 (distorsione +0.000); entro ±0.10: 1.000 | idem | passa (dopo correzione 1) |
| Copertura IC95% quota paziente | 200 | 193/200 = 0.965 | [0.90, 0.99] | passa (dopo correzione 2) |
| Copertura IC95% quota tessuto | 200 | 200/200 = 1.000 | [0.90, 0.99] | **fuori banda → xfail(strict=True), stima NON esposta** |
| Rifiuto con fattori confusi / < 5 pazienti | — | quote non prodotte | rifiuto | passa |
| Caso GSE278694 | — | tumore/adiacente: 5 unità, "stimabile con bassa potenza"; scRNA/snRNA: 0 unità, "non stimabile" | come da prompt | passa |

Storia della calibrazione della scomposizione della varianza (due correzioni del metodo, banda mai modificata):
1. Prima esecuzione: distorsione della quota paziente -0.044 (fuori dalla tolleranza ±0.03); copertura 0/200 per entrambe le quote. Le cause: ω² con divisore N anche per il fattore casuale, e un bootstrap percentile spostato verso il basso dai pazienti duplicati.
2. Correzione 1: componenti con il metodo dei momenti (divisore N·(L-1)/L per i fattori casuali, N per quelli fissi), quota come rapporto fra le medie sui geni, intervallo "basic". Risultato: distorsione 0.000, ma il bootstrap restava spostato di più della sua ampiezza.
3. Correzione 2: intervallo normale, stima ± 1.96 · SE bootstrap. Risultato: paziente 0.965 (in banda), tessuto 1.000 (sopra la banda: intervallo conservativo). Il test è marcato xfail(strict=True) e la scomposizione non è esposta in CLI, app o report.

### 3. Output reale

`python cli.py design --meta gse_like.csv --patient-col patient --tissue-col tissue --technical protocol=protocol --technical library=library --compare tissue:Tumor:Adjacent_normal --compare protocol:scRNA:snRNA`, dove `gse_like.csv` è generato da `make_gse278694_like_sheet()`:

```
Audit del disegno su 37 unità (righe dei metadati) e 4 fattori. 2 fatti strutturali rilevati (annidamenti, coincidenze, esiti determinati da un fattore): Nel disegno attuale il fattore 'patient' è annidato in 'protocol': ogni livello di 'patient' compare in un solo livello di 'protocol'. Un confronto fra livelli di 'protocol' non distingue l'effetto di 'protocol' da quello dei livelli di 'patient' che contiene. Ogni livello di 'library' corrisponde a una sola coppia paziente-tessuto e viceversa: l'effetto tecnico di 'library' non è separabile da quello della coppia paziente-tessuto. È la norma senza multiplexing o librerie replicate, ma significa che ogni differenza fra tessuti dello stesso paziente include anche la differenza fra due librerie. 0 coppie di fattori con associazione forte (Cramér V >= 0.5). Confronto 'tissue': Tumor vs Adjacent_normal: stimabile con bassa potenza, disegno appaiato (entro paziente), 5 pazienti con entrambi i livelli (soglia dichiarata: 5 pazienti, la stessa del Modulo B; inoltre con 5 pazienti appaiati anche un test di Wilcoxon esatto non può scendere sotto p = 0.062). Un risultato non significativo qui non indica assenza di effetto. Attenzione: 16 pazienti con un solo livello esclusi dal confronto appaiato. Confronto 'protocol': scRNA vs snRNA: non stimabile -- nessun paziente ha entrambi i livelli: la differenza fra 'scRNA' e 'snRNA' coincide con la differenza fra due gruppi di pazienti diversi (fattore confuso con il paziente).
 - 'patient' e 'tissue': Cramér V non valutabile (37 unità per 22x3 celle; servono almeno 10 unità e 2 unità attese per cella). Un numero su una tabella così piccola sarebbe instabile.
 - 'patient' e 'protocol': Cramér V non valutabile (37 unità per 22x2 celle; servono almeno 10 unità e 2 unità attese per cella). Un numero su una tabella così piccola sarebbe instabile.
 - 'tissue' e 'protocol': Cramér V corretto = 0.40 su 37 unità, sotto la soglia di 0.5: nessuna segnalazione (non significa che il disegno sia perfettamente bilanciato).
 * Nota fissa sui protocolli: lo snRNA-seq (nuclei) sottorappresenta le cellule immunitarie e cambia la composizione misurata rispetto allo scRNA-seq; i protocolli per tessuto fissato (FFPE, 10x Flex) misurano l'RNA tramite sonde e danno profili non direttamente sovrapponibili a quelli da tessuto fresco. Un confronto fra protocolli è interpretabile solo se gli stessi pazienti sono misurati con entrambi.
 * Le unità di questo audit sono le righe della tabella dei metadati (campioni o librerie), non le cellule: le cellule dello stesso campione non sono osservazioni indipendenti del disegno.
[ok] report scritto in results/design_report_gse_like.html
```

### 4. Decisioni e assunzioni

- **Unità dell'audit = righe dei metadati.** Un CSV è preso così com'è (una riga per campione o libreria). Un `adata.obs` a livello di cellula viene ridotto alle combinazioni distinte dei fattori scelti, perché le cellule non sono unità indipendenti del disegno.
- **Quando il V è valutabile:** almeno 10 righe e almeno 2 righe attese per cella in media (n/(r·c) ≥ 2). Sotto questa soglia l'output è "non valutabile", mai un numero. La soglia è stata fissata prima del test.
- **Soglia di allarme sul V = 0.5**, come da prompt. L'allarme non scatta se la coppia ha già un fatto strutturale, che viene riportato al suo posto.
- **Fattori costanti o identificativi di riga** (un livello per riga, es. `library` in un foglio per libreria) vengono esclusi dalle coppie. Il controllo "unità tecnica" verifica se un fattore tecnico corrisponde 1:1 alla coppia paziente-tessuto.
- **"Esito determinato"** si applica quando il fattore esterno è il tessuto o una colonna di esito. Per gli altri fattori il fatto si chiama "annidamento".
- **Classificazione dei confronti.** Con almeno un paziente che ha entrambi i livelli, il disegno è appaiato e le unità sono quei pazienti. Senza pazienti con entrambi i livelli: un confronto su tessuto o su un fattore tecnico è "non stimabile" (confuso con il paziente); una colonna dichiarata come esito diventa un confronto fra pazienti, con unità pari ai pazienti del gruppo più piccolo. Qualunque fattore tecnico che separi perfettamente i due livelli rende il confronto "non stimabile".
- **DEVIAZIONE sulla soglia di bassa potenza.** Il prompt definisce "bassa potenza" come meno di 5 unità, ma si aspetta che GSE278694 (5 pazienti) risulti a bassa potenza, e le due richieste si contraddicono. Ho applicato la lettura più prudente: bassa potenza se le unità sono meno di 5 (la soglia del Modulo B, dichiarata) OPPURE se il p-value minimo raggiungibile da un test esatto supera 0.05. Nel disegno appaiato quel minimo è quello di Wilcoxon, `core.stats.wilcoxon_min_pvalue`, pari a 0.0625 con 5 coppie; fra pazienti è quello di Mann-Whitney, 2/C(n_a+n_b, n_a). Con 5 coppie risulta quindi "bassa potenza", con 6 "stimabile".
- **Scomposizione della varianza:** il paziente (e il batch) è trattato come fattore casuale, il tessuto come fisso. Per l'intervallo servono almeno 5 pazienti. La scomposizione viene rifiutata se più del 10% delle repliche bootstrap non è identificabile. Le unità con meno di 20 cellule vengono escluse; si usano 1000 geni ad alta varianza, con CPM+log1p tramite `LogCPM(target_sum=1e6)` del Modulo A.
- **Scomposizione della varianza NON esposta** (regola 2 di "Quando fermarti"): è presente in `core/`, tolta dall'orchestratore `run_design_audit` e annotata in `NOTE.md`.
- La nota fissa sui protocolli compare solo se viene dichiarata una colonna con ruolo `protocol`.
- **Baseline committata:** il repository non aveva commit. Ho creato il commit `cf61e76` con lo stato iniziale, così che ogni intervento sia un diff revisionabile.

### 5. Bug o anomalie

- Nel primo test funzionale cercavo la parola "assente" nella frase. Il codice usa "non compare": ho corretto il test, non il codice.
- Nell'app, la colonna "Cramér V" mescolava numeri e testo ("non valutabile") e faceva fallire la conversione pyarrow (trovato con `streamlit.testing.AppTest`). Ora è formattata come stringa.
- Un errore di sintassi (apice sbagliato) nel footer di `core/report.py`, introdotto da me, è stato corretto prima di qualunque esecuzione dei test.
- Nessuna modifica al codice esistente di calcolo: le calibrazioni già approvate (`test_stats_calibration.py`) sono state comunque rieseguite nella suite completa, senza variazioni (stesse righe PASSED/XFAIL della baseline).

### 6. Limiti noti

- La scomposizione della varianza non è disponibile per l'utente (vedi sopra). La sua calibrazione copre un solo disegno: bilanciato, 12 pazienti x 2 tessuti, 300 geni.
- Il V di Cramér è calibrato solo su uno scenario indipendente (10 pazienti x 8 campioni, fattori con 2-4 livelli) e su tabelle 3x3 per la potenza; la potenza a V = 0.5 è 0.46, quindi un'associazione moderata viene spesso non segnalata.
- I disegni sbilanciati o parzialmente annidati sono gestiti in modo approssimato (le componenti di varianza non sono più esattamente non distorte).
- Un confronto fra pazienti su una colonna di esito presume che l'esito sia davvero a livello di paziente; lo strumento non lo verifica.
- Nel foglio di GSE278694 la struttura della coorte snRNA (8 pazienti, un campione tumorale ciascuno) è un'approssimazione: i dettagli reali non sono stati verificati.
- L'app è stata verificata solo in modalità demo con `AppTest` (nessun upload reale).

---

## Intervento 2 — Flag per cellula dal Modulo B

### 1. File modificati e creati

- `core/tcr_validation.py` (ESISTENTE, **toccato** come consentito dal prompt; solo aggiunte, +91 righe e nessuna riga rimossa):
  - `TcrValidationResult` ha due nuovi campi opzionali in coda, con default: `cell_flags` e `flag_coverage`. I costruttori esistenti restano compatibili.
  - `run_tcr_validation` aggiunge la colonna interna `_row` (posizione di riga originale) prima della merge e calcola i flag quando sono dati `marker_map` e `reference_compartment`. `_row` non entra in nessun calcolo.
  - Nuove funzioni `cell_flags` (stessi filtri di `marker_error_rate`) ed `export_audited` (copia `*_audited.h5ad` + CSV), più `FLAG_COLUMNS` e `FLAG_DEFINITIONS`.
  - `pairwise_excess_discordance`, `assign_reference_identity`, `marker_error_rate` e i narrativi NON sono stati modificati.
- `cli.py` (ESISTENTE): opzione `tcr --export-flags` e stampa della copertura dei flag.
- `core/report.py` (ESISTENTE): sottosezione "Flag per cellula" nel Modulo B, con la copertura.
- `app.py` (ESISTENTE): copertura, download del CSV dei flag e della copia `.h5ad`.
- `tests/test_tcr_flags.py` (NUOVO), `tests/fixtures/tcr_regression.py` e `tests/fixtures/__init__.py` (NUOVI).
- `tests/fixtures/tcr_regression_baseline.json` (NUOVO): output di `run_tcr_validation` su 4 configurazioni sintetiche, generato dal codice del commit di baseline PRIMA di modificare `tcr_validation.py`. Prima di generarlo ho verificato con `git diff cf61e76 -- core/tcr_validation.py core/stats.py` che non ci fossero differenze.
- `core/stats.py`: NON toccato.

### 2. Test

Suite completa dopo l'Intervento 2: **46 passed, 2 xfailed** in 475 s (`logs/pytest_intervento2.txt`), contro 36 passed e 2 xfailed dopo l'Intervento 1. Nessun test preesistente è cambiato di esito. Le calibrazioni già approvate (`test_stats_calibration.py`) e i test del Modulo B (`test_synthetic_data.py`) sono stati rieseguiti dopo la modifica a `tcr_validation.py` e passano tutti.

Output reale di `pytest -v`, test nuovi:

```
tests/test_tcr_flags.py::test_regression_identical_to_baseline[eccesso_035_marcatori] PASSED
tests/test_tcr_flags.py::test_regression_identical_to_baseline[eccesso_0_marcatori] PASSED
tests/test_tcr_flags.py::test_regression_identical_to_baseline[tre_compartimenti] PASSED
tests/test_tcr_flags.py::test_regression_identical_to_baseline[senza_marcatori_pochi_pazienti] PASSED
tests/test_tcr_flags.py::test_flags_mean_equals_marker_error_rate_exactly PASSED
tests/test_tcr_flags.py::test_reference_compartment_and_unverifiable_cells_are_na PASSED
tests/test_tcr_flags.py::test_flags_follow_rows_when_order_is_shuffled PASSED
tests/test_tcr_flags.py::test_export_writes_copy_and_leaves_original_untouched PASSED
tests/test_tcr_flags.py::test_export_refuses_existing_flag_columns PASSED
tests/test_tcr_flags.py::test_no_flags_without_marker_map PASSED
```

Questo intervento non introduce nessuna nuova stima statistica: i flag riportano, cellula per cellula, la stessa quantità del tasso d'errore esistente. Il requisito di calibrazione è quindi sostituito dai due test richiesti.
- **Coerenza:** su 5 dataset sintetici con 3 compartimenti, per ogni compartimento non di riferimento la media di `audit_label_vs_reference` sulle cellule valutabili è uguale con `==` (uguaglianza esatta in virgola mobile) alla stima puntuale di `marker_error_rate`. Anche il numero di cellule valutabili coincide con `n_obs` del bootstrap.
- **Regressione:** 4 configurazioni (eccesso 0.35 con marcatori; eccesso 0 con marcatori; 3 compartimenti; 4 pazienti senza marcatori). Tutti i numeri, le tabelle per coppia di compartimenti e le frasi narrative sono identici al riferimento generato prima della modifica, con tolleranza 1e-12 (di fatto identici). Come prova di sensibilità, cambiando il solo seme del bootstrap il test fallisce, quindi il confronto non è banale.

### 3. Output reale

Dataset sintetico (`make_tcr_validation_dataset(n_patients=12, n_clones_per_patient=15, injected_excess=0.35, seed=0)`, salvato come `sintetico.h5ad` con i CSV VDJ e il manifest), comando `python cli.py tcr ... --marker-map markers.json --reference-compartment PBMC --export-flags --n-boot 1000`:

```
Le cellule dello stesso clone T (stessa sequenza CDR3 della catena TRB) cambiano etichetta di tipo cellulare fra compartimenti tissutali con un eccesso di discordanza di +0.107 (IC95% [+0.088, +0.127], su 180 coppie clone-compartimenti da 12 pazienti) rispetto al rumore di base entro lo stesso compartimento: un effetto reale (l'intervallo esclude lo zero). Tasso d'errore dell'annotazione rispetto all'identità clonale dai marcatori (riferimento: compartimento 'PBMC', 180 cloni risolti): 'Tumor': 0.283 (IC95% [0.250, 0.312], 12 pazienti).
Flag per cellula valutabili (audit_label_vs_reference non NA): 39.6% delle cellule.
[ok] copia con i flag scritta in out/sintetico_audited.h5ad (il file originale non è toccato)
[ok] flag per cellula scritti in out/sintetico_audit_flags.csv
[ok] report scritto in out/tcr_report.html
```

Prime righe di `sintetico_audit_flags.csv` (le prime cellule sono nel sangue, quindi `audit_label_vs_reference` è vuoto, cioè NA):

```
obs_name,audit_reference_label,audit_label_vs_reference
cell0,CD8T,
cell1,CD8T,
cell2,CD8T,
cell3,CD8T,
```

Nella copia `sintetico_audited.h5ad`, `audit_label_vs_reference` è di tipo `BooleanDtype`, con NA 1381, False 648 e True 256 su 2285 cellule. Le colonne originali hanno gli stessi valori. Frase generata nel report HTML (sezione "Flag per cellula"): "Copertura: 2,285 cellule su 2,285 (100.0%) hanno un'identità di riferimento; 904 (39.6%) sono valutabili, e di queste 256 sono discordanti. Tutte le altre sono NA: lo strumento non poteva verificarle (…). NA non significa "corretta"."

### 4. Decisioni e assunzioni

- **Cellule del compartimento di riferimento:** `audit_label_vs_reference` = NA, anche quando il clone ha un riferimento. `marker_error_rate` le esclude perché il riferimento è stimato proprio da loro, e il prompt chiede la stessa logica. `audit_reference_label` invece è valorizzato anche per loro.
- `audit_reference_label` è valorizzato per ogni cellula il cui clone ha un riferimento, compresa un'etichetta come "NK": si tratta dell'identità del clone, non di un giudizio sulla cellula.
- **Copertura riportata** = frazione di TUTTE le cellule dell'AnnData con `audit_label_vs_reference` non NA. Sui dati reali, dove il sangue e le cellule senza TCR sono molte, sarà bassa: è voluto.
- **Export:** si rifiuta se `adata.obs` ha già colonne con quei nomi (nessuna sovrascrittura) o se l'indice dei flag non coincide con `obs_names`. Le definizioni vanno in `uns['genomic_audit_flags']`; il CSV contiene solo `obs_name` e i due flag, senza righe di commento, così resta leggibile da qualunque parser. Da CLI i file vanno nella cartella di `--out`, con il nome del file `.h5ad` di input.
- Nessun flag dei doppietti, come da prompt.

### 5. Bug o anomalie

- Nessun bug trovato nel codice esistente. Un'anomalia di formato: scrivendo la copia, anndata converte le colonne di testo in categoriali; i valori non cambiano. Il test verifica l'uguaglianza dei valori come stringhe, non del tipo.
- I 10 test nuovi sono passati al primo tentativo. Per non fidarmi di un test di regressione potenzialmente banale, ho fatto la verifica di sensibilità descritta sopra.

### 6. Limiti noti

- I flag hanno la stessa definizione di identità di riferimento del Modulo B (≥ 3 cellule nel sangue, margine ≥ 0.20, marcatori "conta > 0"). Le cellule di cloni piccoli restano NA.
- Non ho testato AnnData con `obs_names` duplicati, né file `.h5ad` aperti in modalità `backed`.
- Nell'app, la copia `.h5ad` viene costruita in memoria: su dataset grandi potrebbe essere lenta o pesante. Non l'ho misurato.

---

## Intervento 3 — Propagazione dell'errore sulla frazione di CD8

### 1. File modificati e creati

- `core/cd8_propagation.py` (NUOVO): matrice di confusione con direzione, inversione, ricampionamento dei pazienti, Monte Carlo sui conteggi, tre scenari, rifiuti espliciti, testo delle assunzioni.
- `core/synthetic.py` (ESISTENTE, solo aggiunta in coda): `make_cd8_fraction_dataset`, con frazione vera per paziente nota ed errore con direzione iniettato.
- `core/report.py` (ESISTENTE): sezione `_cd8_section`; `render_report` e `save_report` hanno un nuovo parametro opzionale in coda, `cd8_result`.
- `cli.py` (ESISTENTE): opzioni `tcr --cd8-compartment`, `--cd4-label` e `--cd8-label`.
- `app.py` (ESISTENTE): sottosezione "Frazione di CD8" nella scheda del Modulo B.
- `tests/test_cd8_propagation_calibration.py` (NUOVO).
- `README.md` e `NOTE.md`: aggiornati alla fine del lavoro, con cosa fa e cosa non fa ogni modulo e con le idee scartate o rimandate.
- `core/stats.py` e `core/tcr_validation.py`: NON toccati in questo intervento. `marker_error_rate` resta invariata; la categoria "altro" esiste solo in `cd8_propagation.py`.

### Formula di inversione (come richiesto)

Notazione: per una cellula vera CD4, a4 = P(chiamata CD4), b4 = P(chiamata CD8), o4 = P(chiamata altro); per una vera CD8, a8, b8, o8. Le probabilità sono stimate sulle cellule del compartimento con identità di riferimento, pooled fra pazienti. Per il paziente, r è la frazione CD8 osservata, cioè n_CD8 / (n_CD4 + n_CD8). Con x e y le cellule vere CD4 e CD8:

    n_CD4 = a4·x + a8·y        n_CD8 = b4·x + b8·y

Da cui, con (n_CD4, n_CD8) ∝ (1 − r, r):

    x ∝ b8·(1 − r) − a8·r      y ∝ a4·r − b4·(1 − r)      frazione vera = y / (x + y), troncata a [0, 1]

Trattamento di "altro": le cellule chiamate "altro" non entrano né nel numeratore né nel denominatore della frazione riportata. La loro perdita, diversa fra CD4 e CD8, è già contenuta nelle a e b, che per riga sommano a 1 − o. Non serve quindi una categoria "altro" nell'inversione, e le cellule non-T chiamate CD4/CD8 non sono modellate.

Condizionamento: J = b8/(a8+b8) − b4/(a4+b4). Rifiuto se J < 0.2, oppure se il 2.5° percentile di J nel bootstrap è ≤ 0.

Scenari: l'errore di riga, 1 − chiamata corretta, è moltiplicato per k ∈ {0.5, 1, 2}, mantenendo la proporzione fra "chiamata sbagliata" e "altro".

Incertezza: le repliche della matrice ricampionano i pazienti con lo schema di `core.stats.cluster_bootstrap`; per ciascuna si estrae r* da Beta(n_CD8 + ½, n_CD4 + ½) e si inverte. L'intervallo è dato dai percentili 2.5–97.5.

### 2. Test

Suite completa finale: **53 passed, 2 xfailed** in 524 s (`logs/pytest_intervento3.txt`), contro 46 passed e 2 xfailed dopo l'Intervento 2 e 18 passed e 1 xfailed alla baseline. Nessun FAILED.

Output reale di `pytest -v -s`, test nuovi:

```
[sym] copertura IC95% (scenario 1x): 1910/2000 = 0.955 (banda (0.9, 0.99)); frazione riportata senza correzione: 1368/2000 = 0.684; intervalli rifiutati: 0
[asym] copertura IC95% (scenario 1x): 1916/2000 = 0.958 (banda (0.9, 0.99)); frazione riportata senza correzione: 597/2000 = 0.298; intervalli rifiutati: 0
tests/test_cd8_propagation_calibration.py::test_scenarios_are_ordered_and_all_reported PASSED
tests/test_cd8_propagation_calibration.py::test_refuses_with_fewer_than_five_reference_patients PASSED
tests/test_cd8_propagation_calibration.py::test_refuses_ill_conditioned_matrix PASSED
tests/test_cd8_propagation_calibration.py::test_resampling_matches_cluster_bootstrap_exactly PASSED
tests/test_cd8_propagation_calibration.py::test_end_to_end_from_module_b_flags PASSED
```

(Le due righe `[sym]` e `[asym]` sono stampate da `test_coverage_symmetric_error` e `test_coverage_asymmetric_error`, entrambe PASSED; pytest con `-s` stampa "PASSED" su una riga a parte, qui filtrata.)

| Calibrazione | Repliche | Copertura IC95% (scenario 1x) | Banda dichiarata | Senza correzione (solo conteggi) | Rifiuti |
|---|---|---|---|---|---|
| Errore simmetrico (CD4→CD8 = CD8→CD4 = 0.15, altro 0.05) | 200 × 10 pazienti = 2000 intervalli | 1910/2000 = **0.955** | [0.90, 0.99] | 1368/2000 = 0.684 | 0 |
| Errore asimmetrico (CD4→CD8 = 0.24 = 3 × CD8→CD4 = 0.08, altro 0.05) | 2000 | 1916/2000 = **0.958** | [0.90, 0.99] | 597/2000 = 0.298 | 0 |

Entrambe sono in banda al primo tentativo, senza correzioni del metodo. La colonna "senza correzione" mostra che l'inversione conta: l'intervallo dei soli conteggi manca il valore vero nel 32% dei casi con errore simmetrico e nel 70% con errore asimmetrico. Altri test: equivalenza esatta (1e-12) con `cluster_bootstrap`; rifiuto con 4 pazienti di riferimento; rifiuto con J ≈ 0.05 (errore 0.45/0.45); tre scenari sempre presenti; test end-to-end dai flag del Modulo B.

### 3. Output reale

Stesso dataset sintetico dell'Intervento 2 (errore nel tumore circa 42%), comando `python cli.py tcr ... --marker-map markers.json --reference-compartment PBMC --cd8-compartment Tumor --n-boot 1000`, output completo:

```
Le cellule dello stesso clone T (stessa sequenza CDR3 della catena TRB) cambiano etichetta di tipo cellulare fra compartimenti tissutali con un eccesso di discordanza di +0.107 (IC95% [+0.088, +0.127], su 180 coppie clone-compartimenti da 12 pazienti) rispetto al rumore di base entro lo stesso compartimento: un effetto reale (l'intervallo esclude lo zero). Tasso d'errore dell'annotazione rispetto all'identità clonale dai marcatori (riferimento: compartimento 'PBMC', 180 cloni risolti): 'Tumor': 0.283 (IC95% [0.250, 0.312], 12 pazienti).
Flag per cellula valutabili (audit_label_vs_reference non NA): 39.6% delle cellule.
Matrice di confusione stimata su 1131 cellule con identità di riferimento da 12 pazienti (pooled fra pazienti), J = 0.43. Frazione CD8 nel compartimento 'Tumor' del paziente PT000 (CD8T sul totale delle cellule etichettate CD4T o CD8T, n = 71): riportata 0.56 (IC95% dei soli conteggi 0.44-0.67, senza correzione). Intervallo plausibile al 95% tenendo conto dell'errore di annotazione misurato, per scenario di errore sui cloni non condivisi con il sangue -- 0.5x: fra 0.43 e 0.74; 1x: fra 0.37 e 0.92; 2x: non prodotto (matrice di confusione mal condizionata (J = -0.51; servono J >= 0.2 e un intervallo di J che escluda lo zero)). Nessuno scenario è "il risultato": la loro distanza mostra quanto la conclusione dipende dall'assunzione. Frazione CD8 nel compartimento 'Tumor' del paziente PT001 (CD8T sul totale delle cellule etichettate CD4T o CD8T, n = 74): riportata 0.43 (IC95% dei soli conteggi 0.31-0.55, senza correzione). Intervallo plausibile al 95% tenendo conto dell'errore di annotazione misurato, per scenario di errore sui cloni non condivisi con il sangue -- 0.5x: fra 0.25 e 0.57; 1x: fra 0.05 e 0.61; 2x: non prodotto (matrice di confusione mal condizionata (J = -0.51; servono J >= 0.2 e un intervallo di J che escluda lo zero)). Nessuno scenario è "il risultato": la loro distanza mostra quanto la conclusione dipende dall'assunzione. Frazione CD8 nel compartimento 'Tumor' del paziente PT002 (CD8T sul totale delle cellule etichettate CD4T o CD8T, n = 70): riportata 0.47 (IC95% dei soli conteggi 0.36-0.59, senza correzione). Intervallo plausibile al 95% tenendo conto dell'errore di annotazione misurato, per scenario di errore sui cloni non condivisi con il sangue -- 0.5x: fra 0.31 e 0.62; 1x: fra 0.18 e 0.72; 2x: non prodotto (matrice di confusione mal condizionata (J = -0.51; servono J >= 0.2 e un intervallo di J che escluda lo zero)). Nessuno scenario è "il risultato": la loro distanza mostra quanto la conclusione dipende dall'assunzione. Frazione CD8 nel compartimento 'Tumor' del paziente PT003 (CD8T sul totale delle cellule etichettate CD4T o CD8T, n = 74): riportata 0.65 (IC95% dei soli conteggi 0.54-0.75, senza correzione). Intervallo plausibile al 95% tenendo conto dell'errore di annotazione misurato, per scenario di errore sui cloni non condivisi con il sangue -- 0.5x: fra 0.56 e 0.83; 1x: fra 0.60 e 1.00; 2x: non prodotto (matrice di confusione mal condizionata (J = -0.51; servono J >= 0.2 e un intervallo di J che escluda lo zero)). Nessuno scenario è "il risultato": la loro distanza mostra quanto la conclusione dipende dall'assunzione. Frazione CD8 nel compartimento 'Tumor' del paziente PT004 (CD8T sul totale delle cellule etichettate CD4T o CD8T, n = 59): riportata 0.49 (IC95% dei soli conteggi 0.37-0.62, senza correzione). Intervallo plausibile al 95% tenendo conto dell'errore di annotazione misurato, per scenario di errore sui cloni non condivisi con il sangue -- 0.5x: fra 0.33 e 0.66; 1x: fra 0.20 e 0.79; 2x: non prodotto (matrice di confusione mal condizionata (J = -0.51; servono J >= 0.2 e un intervallo di J che escluda lo zero)). Nessuno scenario è "il risultato": la loro distanza mostra quanto la conclusione dipende dall'assunzione. Frazione CD8 nel compartimento 'Tumor' del paziente PT005 (CD8T sul totale delle cellule etichettate CD4T o CD8T, n = 90): riportata 0.51 (IC95% dei soli conteggi 0.41-0.60, senza correzione). Intervallo plausibile al 95% tenendo conto dell'errore di annotazione misurato, per scenario di errore sui cloni non condivisi con il sangue -- 0.5x: fra 0.38 e 0.65; 1x: fra 0.29 e 0.76; 2x: non prodotto (matrice di confusione mal condizionata (J = -0.51; servono J >= 0.2 e un intervallo di J che escluda lo zero)). Nessuno scenario è "il risultato": la loro distanza mostra quanto la conclusione dipende dall'assunzione. Frazione CD8 nel compartimento 'Tumor' del paziente PT006 (CD8T sul totale delle cellule etichettate CD4T o CD8T, n = 80): riportata 0.56 (IC95% dei soli conteggi 0.46-0.66, senza correzione). Intervallo plausibile al 95% tenendo conto dell'errore di annotazione misurato, per scenario di errore sui cloni non condivisi con il sangue -- 0.5x: fra 0.45 e 0.72; 1x: fra 0.41 e 0.90; 2x: non prodotto (matrice di confusione mal condizionata (J = -0.51; servono J >= 0.2 e un intervallo di J che escluda lo zero)). Nessuno scenario è "il risultato": la loro distanza mostra quanto la conclusione dipende dall'assunzione. Frazione CD8 nel compartimento 'Tumor' del paziente PT007 (CD8T sul totale delle cellule etichettate CD4T o CD8T, n = 74): riportata 0.46 (IC95% dei soli conteggi 0.35-0.57, senza correzione). Intervallo plausibile al 95% tenendo conto dell'errore di annotazione misurato, per scenario di errore sui cloni non condivisi con il sangue -- 0.5x: fra 0.30 e 0.60; 1x: fra 0.14 e 0.68; 2x: non prodotto (matrice di confusione mal condizionata (J = -0.51; servono J >= 0.2 e un intervallo di J che escluda lo zero)). Nessuno scenario è "il risultato": la loro distanza mostra quanto la conclusione dipende dall'assunzione. Frazione CD8 nel compartimento 'Tumor' del paziente PT008 (CD8T sul totale delle cellule etichettate CD4T o CD8T, n = 89): riportata 0.45 (IC95% dei soli conteggi 0.35-0.56, senza correzione). Intervallo plausibile al 95% tenendo conto dell'errore di annotazione misurato, per scenario di errore sui cloni non condivisi con il sangue -- 0.5x: fra 0.30 e 0.58; 1x: fra 0.15 e 0.64; 2x: non prodotto (matrice di confusione mal condizionata (J = -0.51; servono J >= 0.2 e un intervallo di J che escluda lo zero)). Nessuno scenario è "il risultato": la loro distanza mostra quanto la conclusione dipende dall'assunzione. Frazione CD8 nel compartimento 'Tumor' del paziente PT009 (CD8T sul totale delle cellule etichettate CD4T o CD8T, n = 77): riportata 0.48 (IC95% dei soli conteggi 0.37-0.60, senza correzione). Intervallo plausibile al 95% tenendo conto dell'errore di annotazione misurato, per scenario di errore sui cloni non condivisi con il sangue -- 0.5x: fra 0.33 e 0.63; 1x: fra 0.21 e 0.75; 2x: non prodotto (matrice di confusione mal condizionata (J = -0.51; servono J >= 0.2 e un intervallo di J che escluda lo zero)). Nessuno scenario è "il risultato": la loro distanza mostra quanto la conclusione dipende dall'assunzione. Frazione CD8 nel compartimento 'Tumor' del paziente PT010 (CD8T sul totale delle cellule etichettate CD4T o CD8T, n = 70): riportata 0.40 (IC95% dei soli conteggi 0.29-0.51, senza correzione). Intervallo plausibile al 95% tenendo conto dell'errore di annotazione misurato, per scenario di errore sui cloni non condivisi con il sangue -- 0.5x: fra 0.22 e 0.52; 1x: fra 0.01 e 0.53; 2x: non prodotto (matrice di confusione mal condizionata (J = -0.51; servono J >= 0.2 e un intervallo di J che escluda lo zero)). Nessuno scenario è "il risultato": la loro distanza mostra quanto la conclusione dipende dall'assunzione. Frazione CD8 nel compartimento 'Tumor' del paziente PT011 (CD8T sul totale delle cellule etichettate CD4T o CD8T, n = 76): riportata 0.36 (IC95% dei soli conteggi 0.25-0.46, senza correzione). Intervallo plausibile al 95% tenendo conto dell'errore di annotazione misurato, per scenario di errore sui cloni non condivisi con il sangue -- 0.5x: fra 0.17 e 0.45; 1x: fra 0.00 e 0.43; 2x: non prodotto (matrice di confusione mal condizionata (J = -0.51; servono J >= 0.2 e un intervallo di J che escluda lo zero)). Nessuno scenario è "il risultato": la loro distanza mostra quanto la conclusione dipende dall'assunzione. Assunzioni. (1) L'identità di riferimento stimata dal sangue è considerata corretta. (2) La matrice di confusione è stimata solo sui cloni condivisi con il sangue, che per costruzione sono quelli espansi (almeno 3 cellule nel sangue e identità netta), ed è unica per tutti i pazienti (pooled), non per paziente. Nei dati originali (PDAC, GSE278694) l'errore di annotazione cresce con la dimensione del clone (rho = +0.177): la direzione plausibile del bias è quindi una SOVRASTIMA dell'errore sui cloni non condivisi, ma non è stata verificata, perché quelle cellule non hanno riferimento. Per questo sono riportati sempre tre scenari (errore 0.5x, 1x, 2x), nessuno dei quali è "il risultato". (3) Le cellule etichettate CD4 o CD8 sono vere cellule T; doppietti e altre cellule non-T etichettate CD4/CD8 non sono modellati.
[ok] report scritto in out3/tcr_cd8_report.html
```


### 4. Decisioni e assunzioni

- **La matrice è stimata nel compartimento bersaglio**, non in tutti i compartimenti: l'errore dipende dal compartimento (0.195 nel tumore contro 0.059 nel tessuto adiacente in `pdac-ml`).
- **Denominatore:** cellule etichettate CD4T o CD8T nel compartimento bersaglio di quel paziente. È dichiarato nella frase e nel report; le etichette sono configurabili (`--cd4-label`, `--cd8-label`).
- **Ricampionamento dei pazienti.** `core.stats.cluster_bootstrap` restituisce solo media e IC, non le repliche, e non va modificata. Ho quindi scritto `patient_resample_draws`, che ne riproduce lo schema chiamata per chiamata. Un test verifica che la media sulle stesse repliche dia esattamente (1e-12) lo stesso IC di `cluster_bootstrap`. È una **deviazione dalla lettera** ("importali"), ma non dalla sostanza: lo stimatore validato non è stato modificato né reimplementato in modo diverso.
- **Monte Carlo:** una frazione estratta da Beta(n_CD8 + ½, n_CD4 + ½), cioè con prior di Jeffreys, per ogni replica della matrice; 1000 repliche di default (`--n-boot`).
- **Soglie dichiarate:** J minimo 0.2; almeno 20 cellule CD4+CD8 per paziente; rifiuto se oltre il 10% delle repliche è non valido (matrice senza una delle due identità, errore scalato ≥ 1, inversione non definita); almeno 5 pazienti con cellule di riferimento nel compartimento (la stessa soglia del Modulo B).
- **Troncamento a [0, 1]** delle frazioni invertite. Vicino ai bordi l'intervallo può toccare 0 o 1.
- **Il rifiuto per condizionamento vale per singolo scenario:** con errore 2x lo scenario può essere rifiutato mentre 0.5x e 1x vengono prodotti. Succede nel dataset di esempio, dove l'errore nel tumore è circa il 42% e raddoppiato rende J negativo.
- **Nessuno scenario è indicato come "il risultato"** e il testo delle assunzioni compare sempre accanto agli intervalli, come chiesto.

### 5. Bug o anomalie

- **Motivo di rifiuto sbagliato (corretto).** Nel primo output reale lo scenario 2x veniva rifiutato con "le probabilità superano 1", mentre la causa vera era J = −0.51 (inversione mal condizionata): bastava una sola replica bootstrap con errore raddoppiato ≥ 1 per far scattare il rifiuto meno informativo. Ho riordinato i controlli: stima puntuale, poi J, poi repliche non valide (rifiuto solo oltre il 10%). Ho anche allineato il numero di estrazioni Beta alle repliche rimaste. Le calibrazioni sono state rieseguite dopo la correzione e danno numeri identici (1x non è toccato).
- Nessuna correzione a codice esistente in questo intervento.

### 6. Limiti noti

- La calibrazione copre un solo disegno: 10 pazienti, 150–300 cellule T ciascuno, 30% con riferimento scelto a caso, frazione vera fra 0.2 e 0.7. In simulazione l'assunzione chiave (cloni con riferimento rappresentativi) è vera per costruzione: la copertura con cloni NON rappresentativi (errore che dipende dalla dimensione del clone) non è stata misurata. È proprio il caso che gli scenari 0.5x/2x dovrebbero coprire, ma senza garanzia di copertura.
- Matrice pooled: se l'errore varia molto fra pazienti, l'intervallo del singolo paziente può essere troppo stretto. Non misurato.
- Le cellule non-T etichettate CD4/CD8 (es. doppietti) non sono modellate.
- Con errori grandi (come nel dataset di esempio, circa 42%) gli intervalli 1x sono molto larghi (es. 0.37–0.92): è l'informazione corretta, ma va detto a chi legge.
- L'app ricalcola gli intervalli a ogni clic, senza cache; non misurato su dataset reali grandi.

---

## Riepilogo finale

- **Commit:** `cf61e76` baseline (stato iniziale, prima di ogni modifica); `4079895` Intervento 1; `e56d1b9` Intervento 2; il commit "Intervento 3" contiene questo resoconto (vedi `git log`).
- **Suite:** 18 passed / 1 xfailed alla baseline → 36 / 2 → 46 / 2 → **53 passed / 2 xfailed**. Gli xfail sono quello preesistente (Nadeau-Bengio a k=5, conservativo) e quello nuovo (`test_variance_share_bootstrap_coverage_tissue`, conservativo: copertura 1.000 > 0.99).
- **Stime NON esposte per calibrazione fuori banda:** la scomposizione della varianza dell'audit del disegno.
- **`core/stats.py` non è mai stato modificato** (`git diff cf61e76 -- core/stats.py` è vuoto). `core/tcr_validation.py` ha solo aggiunte (Intervento 2); l'output del Modulo B è identico a prima della modifica (test di regressione).
- `README.md` aggiornato con cosa fa e cosa non fa ogni modulo; `NOTE.md` con le idee scartate o rimandate.

---

## Prompt finale — Parte A (PARZIALE: interrotta per limite di utilizzo)

Fatti: A1 (tabella + scomposizione), A2. Da fare: A3 in forma tabellare, A4 (copertura per compartimento), Parte B (B1, B2, B3), suite prima/dopo.

### A1 — Modulo B su GSE278694 reale (`validation/gse278694_moduleB.py`, output in `validation/results/`)

Esecuzione `cli.py tcr` sui dati completi (199.184 cellule, nessun sottocampionamento): 69 s, picco di memoria 6,1 GB. Input preparati in scratch, senza modificare `pdac-ml`: copie dei CSV VDJ con il suffisso "-1" rimosso dal barcode (295.451 righe), perché in `sc_raw.h5ad` il barcode è di 16 nt; manifest PA_xx→PAxx, Normal→Adjacent_normal. Marcatori: CD4T=[CD4], CD8T=[CD8A, CD8B], come in `10_loco.py`.

| Quantità | MVP | Tesi | Differenza | Causa |
|---|---|---|---|---|
| Discordanza entro compartimento (baseline) | 0,124 | 0,107 | +0,017 (+16%) | L'MVP include tutte le cellule con TCR (105.809), la tesi solo quelle etichettate CD4T/CD8T/NK (104.964). Con lo stesso filtro l'MVP dà 0,1070 |
| Eccesso di discordanza | +0,059 [+0,034; +0,098], 11 pz, 1.124 coppie | +0,059 [+0,034; +0,100], 11 pz, 1.114 coppie | < 1% | Coincide; con il filtro della tesi: +0,0590 [+0,0342; +0,1004], 1.114 coppie |
| Errore solo sangue, Tumor | 0,093 [0,050; 0,144], 10 pz | 0,195 [0,123; 0,274], 10 pz | −52% | Definizioni diverse, vedi la scomposizione |
| Errore solo sangue, Adjacent_normal | 0,041, IC NON prodotto (4 pz < 5) | 0,059 [0,000; 0,084], 4 pz | −30% | Stesse definizioni del Tumor; in più la soglia di 5 pazienti dell'MVP rifiuta l'IC |
| Tumor − Adjacent (stessi cloni) | non calcolato dall'MVP; con le definizioni MVP: +0,021 [−0,035; +0,075], 132 cloni, 4 pz | +0,080 [+0,021; +0,170], 88 cloni | — | Stesse cause |

Scomposizione del tasso d'errore Tumor, passando dalle definizioni dell'MVP a quelle della tesi un cambiamento alla volta (stessi dati, stesso bootstrap):

| Passo | Tumor | Adjacent | Tumor − Adj |
|---|---|---|---|
| S0 MVP (per cellula, tutte le cellule con TCR, CD4>0 non esclusivo, solo etichette CD4T/CD8T) | 0,093 | 0,041 | +0,021 |
| S1 + solo cellule CD4T/CD8T/NK | 0,093 | 0,042 | +0,021 |
| S2 + CD4 positivo solo se CD8A e CD8B negativi | 0,092 | 0,043 | +0,021 |
| S3 + etichetta **NK contata come errore** | **0,163** | 0,053 | +0,081 |
| S4 + ≥ 3 cellule del clone anche nel compartimento giudicato | 0,156 | 0,052 | +0,062 |
| S5 + unità = coppia clone-compartimento con etichetta maggioritaria (= tesi) | **0,1951 [0,1230; 0,2737]** | **0,0593 [0,000; 0,0839]** | **+0,0795 [+0,021; +0,170], 88 cloni** |

S5 riproduce esattamente i numeri della tesi, quindi la scomposizione è completa. Le due cause principali sono:
1. `marker_error_rate` dell'MVP esclude le cellule etichettate NK, mentre la tesi le conta come errore: +0,070.
2. La tesi conta un errore per coppia clone-compartimento, pesando i cloni e non le cellule, e richiede ≥ 3 cellule anche nel compartimento giudicato: +0,033 in totale.

La regola di positività di CD4 e il filtro sulle cellule linfoidi hanno effetto trascurabile (≤ 0,001). Nessun parametro dell'MVP è stato modificato.

### A2 — Audit del disegno sui metadati reali

Output completo in `validation/results/cli_design_real.txt` (sc + sn, 37 unità) e in `cli_design_real_sc.txt` (solo sc, via `--h5ad`). L'output coincide con quello del generatore "GSE278694-like" riga per riga, a parte il nome della colonna (`patients`). Verifica strutturale: stessi tessuti per ciascuno dei 14 pazienti sc, 37 righe, 8 pazienti sn disgiunti, tutti tumore. Il generatore NON va corretto.

Una segnalazione: le librerie della coorte che GEO chiama "snRNA" si chiamano `*_FFPE`, e 4 su 8 hanno barcode di sonda (probabile protocollo su FFPE con sonde). Il protocollo va quindi descritto come "snRNA da FFPE", non come snRNA da tessuto congelato.

### A3/A4 — risultati grezzi già ottenuti (da rendere in tabella)

- J pooled = 0,48, su 4.851 cellule di riferimento da 10 pazienti.
- Scenario 2x: non prodotto per i pazienti visti ("probabilità di errore stimata supera 1").
- Scenario 1x: sposta fortemente verso il basso la frazione di CD8 (es. PA01: riportata 0,61, 1x 0,00–0,35).
- Flag valutabili: 4,2% di tutte le cellule; manca il dettaglio per compartimento.

### A3 — Frazione di CD8 sui dati reali (completata dopo l'interruzione)

`validation/gse278694_cd8_flags.py` ricalcola sui flag esportati dalla CLI con gli stessi parametri (2000 repliche, seme 0). I numeri coincidono con l'output CLI (`validation/results/cli_tcr_real.txt`).

Matrice di confusione pooled nel tumore: 4.851 cellule con identità di riferimento, da 10 pazienti.

| | chiamata CD4T | chiamata CD8T | chiamata altro |
|---|---|---|---|
| vera CD4T | 0,425 | 0,447 | 0,127 |
| vera CD8T | 0,006 | 0,913 | 0,081 |

J = 0,48, sopra la soglia di 0,2.

Intervalli prodotti per scenario, su 13 pazienti con tumore (PA15 ha solo PBMC):

| Scenario | Intervalli prodotti | Motivo del rifiuto |
|---|---|---|
| 0.5x | 13/13 | — |
| 1x | 13/13 | — |
| 2x | 0/13 | "con errore 2x la probabilità di errore stimata supera 1": l'errore stimato sulle vere CD4 è 1 − 0,425 = 0,575, e raddoppiato supera 1 |

| Paziente | n CD4T+CD8T | Riportata | IC95% soli conteggi | 0.5x | 1x | 2x |
|---|---|---|---|---|---|---|
| PA01 | 2802 | 0,61 | 0,59–0,63 | 0,44–0,53 | 0,00–0,35 | non prodotto |
| PA02 | 3311 | 0,65 | 0,64–0,67 | 0,50–0,58 | 0,06–0,41 | non prodotto |
| PA04 | 3485 | 0,48 | 0,46–0,50 | 0,25–0,37 | 0,00–0,14 | non prodotto |
| PA05 | 1575 | 0,72 | 0,70–0,74 | 0,59–0,67 | 0,23–0,53 | non prodotto |
| PA06 | 1827 | 0,51 | 0,49–0,54 | 0,29–0,41 | 0,00–0,19 | non prodotto |
| PA07 | 3164 | 0,67 | 0,65–0,69 | 0,52–0,60 | 0,10–0,44 | non prodotto |
| PA08 | 4255 | 0,69 | 0,68–0,70 | 0,55–0,62 | 0,15–0,47 | non prodotto |
| PA09 | 2990 | 0,46 | 0,45–0,48 | 0,22–0,35 | 0,00–0,11 | non prodotto |
| PA10 | 2114 | 0,48 | 0,46–0,50 | 0,25–0,37 | 0,00–0,15 | non prodotto |
| PA11 | 5745 | 0,57 | 0,55–0,58 | 0,37–0,47 | 0,00–0,27 | non prodotto |
| PA12 | 3537 | 0,67 | 0,65–0,68 | 0,52–0,59 | 0,10–0,43 | non prodotto |
| PA13 | 3501 | 0,46 | 0,44–0,47 | 0,21–0,34 | 0,00–0,11 | non prodotto |
| PA14 | 571 | 0,54 | 0,50–0,58 | 0,31–0,46 | 0,00–0,24 | non prodotto |

Cosa produce il metodo, senza interpretazione:
- In tutti i 13 pazienti gli intervalli 0.5x e 1x stanno sotto la frazione riportata.
- In 8 pazienti su 13 l'intervallo 1x arriva a 0,00, cioè viene troncato al bordo.
- Il motivo è che la matrice attribuisce il 45% delle vere CD4 (per marcatori nel sangue) all'etichetta CD8T nel tumore.
- Se questa matrice, stimata sui cloni condivisi con il sangue, valga anche per gli altri cloni non è verificabile su questi dati: è il tema di B3.

### A4 — Copertura dei flag (`validation/results/a4_flag_coverage.csv`)

| Compartimento | Cellule | Con riferimento | Valutabili (non NA) | Discordanti | % valutabili | % valutabili fra CD4T+CD8T |
|---|---|---|---|---|---|---|
| Adjacent_normal | 38.264 | 3.921 | 3.858 | 160 | 10,08 | 15,61 |
| PBMC | 89.506 | 16.460 | 0 | 0 | 0,00 | 0,00 |
| Tumor | 71.414 | 4.851 | 4.417 | 409 | 6,19 | 11,36 |
| **Totale** | 199.184 | 25.232 | 8.275 | 569 | **4,15** | 7,15 |

Il PBMC è sempre NA per costruzione, perché è il compartimento di riferimento. Nel tumore il 93,8% delle cellule non è verificabile con questo schema.

---

## Prompt finale — Parte B: correzioni minori

Suite completa **prima** della Parte B (dopo la Parte A): **53 passed, 2 xfailed** (`logs/pytest_prima_parteB.txt`). **Dopo**: **66 passed, 2 xfailed** in 555 s (`logs/pytest_dopo_parteB.txt`). I 13 test nuovi sono 5 di B1, 2 di B2 e 6 di B3; nessun FAILED.

Dopo quell'esecuzione ho cambiato solo il testo delle assunzioni (vedi B3). Ho rieseguito i test che leggono quel testo (`test_cd8_format.py` e la parte non di copertura di `test_cd8_propagation_calibration.py`): 10 passed.

### File toccati
- `core/cd8_propagation.py`:
  - B1: nuovi `SCENARIO_WARNING`, `refusal_notes` e `format_cd8_text`; `narrative` non ripete più le frasi per paziente; la frase per paziente non ripete l'avvertenza.
  - B3: attenuata la frase sulla direzione del bias in `ASSUMPTIONS_TEXT`.
  - Nessun calcolo modificato.
- `core/report.py` (B1): l'avvertenza compare una volta, i rifiuti sono note numerate sotto la tabella.
- `cli.py` (B1): stampa `format_cd8_text`.
- `core/design_audit.py` (B2): `ComparisonAssessment` ha due campi opzionali in coda, `min_pvalue` e `min_pvalue_test`; nuova costante `MIN_PVALUE_NOTE`; le frasi dei confronti sono estese. Soglie e classificazione sono invariate.
- `core/synthetic.py` (B3): nuovo parametro `reference_error_factor` (default 1.0) in `make_cd8_fraction_dataset`. Le estrazioni casuali sono le stesse, quindi con fattore 1 i dati sono identici; lo verifica il fixture di regressione di B1, che passa.
- Test nuovi: `tests/test_cd8_format.py`, `tests/test_cd8_robustness.py`, `tests/fixtures/cd8_regression.py` e `cd8_regression_baseline.json` (generato prima della modifica B1). In `tests/test_design_audit.py` due test nuovi.
- `README.md` (limite di B3), `NOTE.md` (note dalla validazione reale).
- **`core/stats.py` e `core/tcr_validation.py`: NON toccati** nella Parte A né nella Parte B.

### B1 — output tabellare della frazione di CD8
I numeri sono identici al fixture pre-modifica (3 configurazioni, tolleranza 1e-12), compresa una configurazione con scenario 2x rifiutato. Anche i numeri sui dati reali coincidono con quelli dell'esecuzione precedente. Output reale CLI su GSE278694 (`validation/results/cli_tcr_real_B1.txt`):

```
Frazione di CD8 nel compartimento 'Tumor'
Matrice di confusione stimata su 4851 cellule con identità di riferimento da 10 pazienti (pooled fra pazienti), J = 0.48.
Per ogni paziente: frazione di cellule etichettate CD8 sul totale delle cellule etichettate CD4 o CD8 nel compartimento (denominatore dichiarato), intervallo dei soli conteggi (senza correzione) e intervalli plausibili al 95% dopo la correzione per l'errore di annotazione, in tre scenari di errore sui cloni non condivisi con il sangue (0.5x, 1x, 2x). Nessuno scenario è "il risultato": la loro distanza mostra quanto la conclusione dipende dall'assunzione.

paziente  n CD4+CD8  riportata  IC95% conteggi  0.5x       1x         2x              
--------  ---------  ---------  --------------  ---------  ---------  ----------------
PA01      2802       0.61       0.59-0.63       0.44-0.53  0.00-0.35  non prodotto [1]
PA02      3311       0.65       0.64-0.67       0.50-0.58  0.06-0.41  non prodotto [1]
PA04      3485       0.48       0.46-0.50       0.25-0.37  0.00-0.14  non prodotto [1]
PA05      1575       0.72       0.70-0.74       0.59-0.67  0.23-0.53  non prodotto [1]
PA06      1827       0.51       0.49-0.54       0.29-0.41  0.00-0.19  non prodotto [1]
PA07      3164       0.67       0.65-0.69       0.52-0.60  0.10-0.44  non prodotto [1]
PA08      4255       0.69       0.68-0.70       0.55-0.62  0.15-0.47  non prodotto [1]
PA09      2990       0.46       0.45-0.48       0.22-0.35  0.00-0.11  non prodotto [1]
PA10      2114       0.48       0.46-0.50       0.25-0.37  0.00-0.15  non prodotto [1]
PA11      5745       0.57       0.55-0.58       0.37-0.47  0.00-0.27  non prodotto [1]
PA12      3537       0.67       0.65-0.68       0.52-0.59  0.10-0.43  non prodotto [1]
PA13      3501       0.46       0.44-0.47       0.21-0.34  0.00-0.11  non prodotto [1]
PA14      571        0.54       0.50-0.58       0.31-0.46  0.00-0.24  non prodotto [1]

[1] con errore 2x la probabilità di errore stimata supera 1
```
(Segue, una sola volta, il paragrafo "Assunzioni". Il testo stampato in questa esecuzione è quello precedente all'attenuazione di B3.)

### B2 — dicitura dei confronti
La classificazione è invariata; ogni confronto riporta ora unità, p-value minimo e la frase fissa. Output reale sul foglio GSE278694:

> Confronto 'tissue': Tumor vs Adjacent_normal: stimabile con bassa potenza, […] Unità indipendenti: 5; p-value minimo raggiungibile con un test esatto di Wilcoxon appaiato: 0.062. Il p-value minimo non misura la potenza: con questa numerosità solo effetti grandi sono rilevabili.
> Confronto 'protocol': scRNA vs snRNA: non stimabile -- […] Unità indipendenti: 0; p-value minimo non definito (nessun test possibile).

Decisione: con 0 unità il p-value minimo non esiste, e il testo lo dice al posto della frase fissa.

### B3 — robustezza a un errore non rappresentativo (`tests/test_cd8_robustness.py`)
Configurazione: 200 repliche × 10 pazienti = 2000 intervalli per caso. L'errore delle cellule con riferimento è f volte quello delle altre. Copertura dell'IC 95%:

| Errore di base | f | 0.5x | 1x | 2x |
|---|---|---|---|---|
| basso (0.05/0.05, altro 0.02) | 1 | 0.938 | 0.944 | 0.929 |
| basso | 2 | 0.939 | 0.931 | 0.799 |
| basso | 4 | 0.949 | 0.835 | 0.500 (solo 20 prodotti, 1980 rifiutati) |
| alto (CD4→CD8 0.10, CD8→CD4 0.03, altro 0.04) | 1 | 0.878 | 0.944 | 0.755 |
| alto | 2 | 0.923 | 0.839 | 0.270 |
| alto | 4 | 0.945 | 0.285 | non calcolabile (2000 rifiutati) |

Quale scenario copre il valore vero:
- con f = 1 (errore rappresentativo), lo scenario 1x;
- con f = 2 e f = 4, lo scenario 0.5x, che è il più vicino al fattore corretto (0.3f+0.7)/f = 0.65 e 0.475;
- **nessuno scenario garantisce la copertura nominale in tutti i casi.** Lo scenario 1x scende fino a 0.285; lo 0.5x scende a 0.878 quando l'errore è rappresentativo.

L'ho scritto nel README come limite, senza toccare la banda. Prima di aggiungere la seconda base d'errore, la sola base "bassa" dava una misura troppo favorevole: con errori piccoli la distorsione è piccola rispetto all'ampiezza degli intervalli.

Frase sulla "direzione plausibile del bias": riformulata in modo meno sicuro ("è quindi possibile che la matrice sovrastimi…, non verificabile"), con i numeri della simulazione. Lo 0.5x non è molto sotto la copertura dichiarata quando il bias va davvero in quella direzione, ma lo è quando non ci va: ho applicato la lettura prudente.

### Limiti della Parte B
- B3 misura solo errori fino a CD4→CD8 = 0.40 (base "alta" × 4). Nel tumore reale la matrice stimata ha 0.45 già a 1x; i casi intermedi (es. f = 3) non sono misurati.
- Nel report HTML la colonna del p-value minimo non è separata: sta nella frase del confronto.
- Nell'app la tabella CD8 mostra "non prodotto" senza la nota con il motivo (B1 chiedeva CLI e report).

---

## Riepilogo del prompt finale

- **Commit:**
  - `8bb3d1a` Parte A (parziale), prima dell'interruzione per limite di utilizzo;
  - `2f03461` Parte A: validazione su dati reali;
  - il commit "Parte B: correzioni" contiene questa sezione.
- **Discrepanze con la tesi:**
  - la discordanza entro compartimento (0.124 contro 0.107) dipende dal filtro sulle sole cellule linfoidi;
  - il tasso d'errore nel tumore (0.093 contro 0.195) dipende da come sono trattate le cellule NK (+0.070) e dall'unità di conteggio, clone-compartimento invece di cellula (+0.033);
  - l'Adjacent_normal non ha IC perché 4 pazienti < 5.

  Tutte le differenze sono spiegate da definizioni o soglie, riprodotte esattamente e senza tarare parametri.
- **Dati reali:** nessun sottocampionamento (picco 6,1 GB). Serve un passo manuale di normalizzazione dei barcode (annotato in NOTE.md).

---

## Pulizia finale (4 correzioni)

Suite **prima**: 66 passed, 2 xfailed (`logs/pytest_prima_pulizia.txt`). **Dopo**: **74 passed, 2 xfailed**, 0 FAILED, 0 SKIPPED, 580 s (`logs/pytest_dopo_pulizia.txt`). I 3 test sui dati reali sono stati eseguiti: i dati di `pdac-ml` sono presenti su questa macchina; dove mancano, il test viene saltato e ne dichiara il motivo.

**File toccati:**
- `core/tcr_validation.py` (punti 1, 2, 3: solo aggiunte, più la compressione nell'export);
- `cli.py`, `core/report.py`, `app.py`, `core/cd8_propagation.py` (solo testo: avvertenza sperimentale);
- `.gitignore`, `README.md`, `NOTE.md`;
- test: `tests/test_tcr_flags.py` (regressione: i numeri restano a 1e-12, il testo narrativo è cambiato di proposito), `tests/test_real_data_gse278694.py`, `tests/test_barcode_match.py`, `tests/test_cd8_format.py`;
- `validation/gse278694_cd8_diagnostic.py` e `validation/results/*`.

**`core/stats.py` NON toccato** (`git diff cf61e76 -- core/stats.py` è vuoto).

### 1. Convenzioni affiancate
Output reale (`cli.py tcr` su GSE278694 con `--clone-error-labels NK --clone-marker-priority CD8T,CD4T`, file `validation/results/cli_tcr_real_conventions_tesi.txt`):

> [cell] per cellula: ogni cellula giudicata è un'unità; etichette fuori dalla mappa dei marcatori escluse -- 'Adjacent_normal': 0.041, IC non prodotto (4 pazienti, ne servono almeno 5). 'Tumor': 0.093 (IC95% [0.050, 0.144], 10 pazienti). Differenza 'Tumor' - 'Adjacent_normal' sugli stessi cloni: +0.021, IC non prodotto (4 pazienti), 132 cloni. [clone] per clone: un'unità per coppia clone-compartimento con almeno 3 cellule, etichetta di maggioranza; etichette fuori mappa contate come errore (etichette ammesse fuori mappa: NK); positività esclusiva in ordine CD8T > CD4T -- 'Adjacent_normal': 0.059, IC non prodotto (4 pazienti, ne servono almeno 5). 'Tumor': 0.195 (IC95% [0.123, 0.274], 10 pazienti). Differenza 'Tumor' - 'Adjacent_normal' sugli stessi cloni: +0.080, IC non prodotto (4 pazienti), 88 cloni.

Il test sui dati reali verifica la riproduzione esatta (tolleranza 1e-12): Tumor 0.1951 [0.1230; 0.2737], Adjacent 0.0593, Tumor − Adj +0.0795 su 88 cloni. Con i default della convenzione per clone (tutte le etichette fuori mappa contate come errore, positività non esclusiva) il tumore è 0.197 [0.124; 0.276].

Decisioni:
- servono due parametri dichiarati (`--clone-error-labels`, `--clone-marker-priority`), perché senza uno dei due la riproduzione non è esatta (0.1979 senza l'esclusività, 0.1938 senza il filtro sulle etichette);
- default `--convention both`;
- l'IC della differenza non è prodotto perché ci sono 4 pazienti, sotto la soglia di 5 dell'MVP;
- i flag per cellula restano nella convenzione per cellula.

### 2. Match dei barcode
Output reale con i file VDJ originali di `pdac-ml`, senza preparazione manuale:

> Barcode VDJ (cellule con TRB) ritrovati nei metadati per paziente, compartimento e barcode: 105,809/123,280 (85.8%). Normalizzato il suffisso -1: rimosso da 123,280 barcode VDJ e da 0 barcode dei metadati (prima della normalizzazione: 0.0%).

I numeri sono identici a quelli ottenuti prima con i barcode preparati a mano.

Decisioni:
- la chiave di match resta (paziente, compartimento, barcode);
- la normalizzazione si applica solo se aumenta il match e non fonde barcode diversi (es. "-1" e "-2" dello stesso campione);
- sotto il 50% viene sollevato un errore con 3 esempi per lato.

Test sintetici: stesso suffisso sui due lati, nessuna normalizzazione (100%); suffisso solo sul lato VDJ, normalizzato e dichiarato, con numeri e flag identici; collisione, normalizzazione rifiutata ed errore; nessun match, errore "solo il 0.0%…" con esempi.

### 3. Export compresso e motivi di rifiuto nell'app
- `sc_raw_audited.h5ad` su GSE278694: **2,904,147,974 byte senza compressione → 844,778,382 con gzip** (−71%; l'input è 844,714,396 byte). I flag riletti sono identici. Il tempo della CLI con export passa da 69 s a 95 s.
- Nell'app la colonna dello scenario rifiutato mostra ora il motivo, verificato con AppTest sulla demo: "non prodotto: matrice di confusione mal condizionata (J = -0.51; …)".

### 4. Frazione di CD8: "sperimentale" + diagnostica
L'avvertenza fissa `EXPERIMENTAL_NOTE` compare in cima alla sezione in CLI (seconda riga), nel report HTML (prima della tabella, verificato da test) e nell'app (`st.warning`). Nessun calcolo è cambiato.

Diagnostica (`validation/gse278694_cd8_diagnostic.py`, tumore, 10 pazienti):

| Margine | Min cellule nel sangue | Cloni di riferimento (CD4/CD8) | Cloni usati nel tumore | P(chiamata CD8 \| vera CD4) | P(chiamata CD4 \| vera CD8) | J | Cloni CD4 rif. con ≥1 cellula CD8A+ nel sangue |
|---|---|---|---|---|---|---|---|
| 0.20 | 3 | 869 (264/605) | 525 | 0.447 | 0.006 | 0.481 | 22.0% (58/264) |
| 0.20 | 5 | 496 (143/353) | 361 | 0.443 | 0.005 | 0.476 | 35.7% (51/143) |
| 0.40 | 3 | 774 (181/593) | 476 | 0.429 | 0.006 | 0.496 | 20.4% (37/181) |
| 0.40 | 5 | 452 (106/346) | 327 | 0.446 | 0.005 | 0.469 | 31.1% (33/106) |

Interpretazione prudente:
- Il 45% **non scende** con soglie più stringenti (0.429–0.447): non è spiegato da margine e numero minimo di cellule.
- Questo non lo rende un risultato biologico. La diagnostica non tocca altri aspetti del riferimento, cioè la positività definita come "conta > 0" (sensibile al dropout di CD4 e all'RNA ambientale) e il fatto che dal 20 al 36% dei cloni "CD4" di riferimento ha almeno una cellula CD8A+ nel sangue.
- La quota con almeno una cellula CD8A+ cresce con il minimo di cellule, come ci si aspetta da un "almeno una" su più cellule.
- Annotato in NOTE.md come punto aperto.

### Anomalie
- **Errore mio nella Parte A.** `.gitignore` conteneva `results/`, che escludeva anche `validation/results/`. Nei commit della Parte A gli output reali NON erano versionati, anche se il resoconto li indicava così. Corretto in Pulizia 1, con un'eccezione in `.gitignore` e i file aggiunti.
- Nel test della collisione ho usato `DataFrame._append`, che non esiste in questa versione di pandas: corretto in `pd.concat`.

### Limiti
- La convenzione per clone è riprodotta esattamente solo con i parametri della tesi. I suoi default (tutte le etichette fuori mappa, positività non esclusiva) danno numeri vicini ma diversi.
- La soglia del 50% sul match dei barcode è una scelta mia, non calibrata.
- La normalizzazione tratta solo il suffisso "-N", non altri formati (prefissi di campione, ecc.).
- La compressione allunga l'export di circa 26 s su questo dataset.
- La diagnostica del punto 4 non spiega il 45%: esclude solo due cause.

---

## Motore del Modulo A per dati reali, web app e validazione esterna

Il report completo, con criteri × dataset, output reali, difetti e verdetto, è in
`validation_esterna/REPORT.md`. In sintesi:

- **Modulo A riscritto per dati reali:**
  - matrice sparsa; dentro ogni fold, e sul solo training, log-CPM, 2.000 HVG e SVD a 50
    componenti;
  - macro-F1 sulle classi presenti nel fold di test;
  - training limitato a 20.000 cellule per fold e fold in parallelo, con risultati identici;
  - fail-fast su conteggi non interi e su una classe presente in un solo paziente;
  - non convergenza dichiarata.

  Su GSE125449 completo (9.946 cellule, 19 pazienti) impiega 403 s; prima non terminava dopo
  5,5 h.
- **Web app locale** (`uv run audit-sc serve`): demo con un clic, caricamento di studi reali,
  semafori con regole dichiarate, export Markdown e HTML; solo localhost, telemetria spenta.
- **Validazione esterna:** audit del disegno D1-D4 superati su 3 dataset GEO (GSE132465,
  GSE131907, GSE125449), dopo la correzione dei valori mancanti. Modulo A A1-A4 superati su
  GSE125449 completo. La prima esecuzione era NON PRONTO, e i criteri sono stati modificati
  due volte su decisione dell'utente: entrambe le cose sono dichiarate.
- **Suite finale:** 97 passed, 2 xfailed (calibrazioni conservative note), 0 failed
  (`logs/pytest_finale.txt`). `core/stats.py` non è stato modificato.

---

## Rifinitura finale della repo pubblica (2026-10-04)

Trasparenza, documentazione e igiene. Nessuna nuova funzionalità e nessuna modifica ai calcoli:
`core/stats.py` e `validation_esterna/CRITERI.md` sono identici byte per byte a prima; negli
altri moduli di calcolo la differenza è la sola conversione dei testi in lettere accentate
(verificato file per file). Cinque commit, `Rifinitura 0` … `Rifinitura 4`, sul ramo
`rifinitura`.

### 0. Prova di installazione pulita

Clone dell'URL pubblico (commit `0662a40`) in una cartella vuota, senza modificare nulla. Output
in `docs/sviluppo/logs/installazione_pulita.txt`.

| Strada | Installazione | `cli.py demo` | `pytest` |
|---|---|---|---|
| uv, seguendo il README | `uv sync --extra dev`: 0,9 s (cache di uv già presente) | exit 0, 2 min 59 s, report di 126.588 byte | 97 passed, 2 xfailed in 938 s |
| venv + `pip install -e ".[dev]"` | exit 0, 3 min 01 s | exit 0, 3 min 02 s | 97 passed, 2 xfailed in 970 s |

Entrambe funzionano. Mancanze del README, corrette al punto 3: la strada con pip non era
documentata; il requisito Python >= 3.12 non era scritto; non era detto che i 3 test sui dati
reali di GSE278694 girano solo dove quei dati esistono (usavano un percorso assoluto della
macchina di sviluppo). I tempi sono un limite superiore: la macchina eseguiva anche altro.

### 1. Trasparenza della validazione esterna

File: `validation_esterna/REPORT.md` e `REPORT_template.md` (nuove sezioni 6 e 7),
`make_report.py`, `criteri_originali.py` (nuovo), `indipendenti/leakage_check.py` (opzione
`--pipeline originale`), nuovi file in `results/`.

Il codice attuale è stato rivalutato con le tre versioni dei criteri senza modificarlo: cambia
solo la regola nello script di verifica.

| Criteri | Dataset del Modulo A | Codice | A1 | A2 (tolleranza 0,02) | A3 | A4 | Tutti soddisfatti |
|---|---|---|---|---|---|---|---|
| O, originali (`0b5a847`) | completo | prima esecuzione | NO: non terminato in 5 h 26 min | non valutabile | non valutabile | A4c NO (0,667 invece di 1,0) | **NO** |
| M1, sezione 5 (`b08b2d2`) | ridotto | motore precedente, rieseguito oggi | non concluso: interrotto dopo 45 min | — | — | — | non determinato |
| O, originali | completo | finale | OK | OK: 0.7776 contro 0.7880 (0.0104); LeaveOneGroupOut 0.7712 contro 0.7593 (0.0119) | OK | OK | sì |
| M1, sezione 5 | ridotto | finale | OK, 156 s | **NO**: LeaveOneGroupOut 0.7530 contro 0.7127 (0.0403) | OK | OK | **NO** |
| M2, sezione 6 (`0662a40`) | completo | finale | OK | OK: 0.7776 contro 0.7779; 0.7712 contro 0.7695 | OK | OK | sì |

Audit del disegno: i criteri D1-D4 non sono mai cambiati. Prima esecuzione: GSE132465 falliva D2
e D4; dopo la correzione `7bb30ca`: soddisfatti su tutti e tre i dataset (rieseguito oggi con i
testi attuali: fatti strutturali, classi e V identici).

**Verdetto con i criteri originali: NON PRONTO** (prima esecuzione, l'unica fatta senza conoscere
i risultati). **Verdetto con i criteri modificati: PRONTO.** Il codice finale soddisfa anche i
criteri originali, ma solo alla lettera per A2 (a parità di definizione della macro-F1 la
differenza LeaveOneGroupOut è 0.0212, oltre la tolleranza), e questa rivalutazione è successiva
ai risultati. Con la prima modifica (coorte ridotta) il codice finale non soddisfa A2. Inoltre la
sezione 6 dei criteri è stata committata insieme ai risultati finali: git non dimostra che sia
stata scritta prima.

### 2. Sintesi dei controlli (ex «verdetto a semaforo»)

File: `core/verdict.py`, `core/report.py`, `app.py`, `tests/test_verdict.py` (nuovo, 18 test),
`tests/test_app_and_io.py`, `tests/test_leakage_engine.py`, README. Nessun calcolo e nessuna
soglia è cambiata.

Prima: un colore per sezione (verde, giallo, rosso).
- Disegno: rosso se un confronto è non stimabile; giallo con bassa potenza o Cramér V >= 0.5;
  altrimenti verde, **anche senza confronti richiesti e con V non valutabili**.
- Modulo A: rosso se il divario supera 0.05 («qualunque valutazione che non separa i pazienti è
  inaffidabile»); giallo per classi assenti, training sottocampionato, non convergenza o
  confronto fra modelli non eseguito; altrimenti verde.
- Modulo B: rosso se l'intervallo dell'eccesso di discordanza è sopra lo zero («errore
  sistematico di annotazione»); giallo con numerosità insufficiente o tasso senza intervallo;
  verde altrimenti.

Dopo: nessun colore per sezione o per studio; una riga per controllo, con quattro stati: verde
«stima affidabile», giallo «stima con limiti», rosso «stima non possibile con questi dati»,
grigio «controllo non eseguito o non valutabile». Le soglie sono le stesse (0.05 sul divario,
0.5 sul V, 5 pazienti, classi dell'audit del disegno). Che cosa cambia nella visualizzazione:
- confronto fra modelli non eseguito, V non valutabile, nessun confronto richiesto, tasso
  d'errore non calcolato, frazione di CD8 rifiutata: grigio, mai verde;
- numerosità insufficiente nel Modulo B: da giallo a rosso («stima non possibile»);
- intervallo dell'eccesso di discordanza sopra lo zero: da rosso a giallo, con un testo che
  dice che cosa si può stimare invece di «errore sistematico»;
- frazione di CD8: entra nella sintesi, mai verde;
- report HTML: sintesi in apertura; il bordo verde di una scheda compare solo se tutti i
  controlli della sezione sono verdi.

La tabella completa delle regole, con le soglie, è nel README.

### 3. README e lettere accentate

README riscritto nell'ordine richiesto, con tutti i contenuti tecnici di prima. Estratti di
output reali: disegno su GSE278694 (confronto scRNA/snRNA non stimabile), Modulo A su GSE125449
(0.778 contro 0.895, divario +0.117), Modulo B su GSE278694 (0.093 per cellula, 0.195 per
clone). Screenshot in `docs/img/report_demo.png`.

Accenti: `docs/sviluppo/accenti.py` converte un elenco esplicito di parole (non ogni «vocale +
apostrofo», perché l'apostrofo chiude anche le citazioni come 'tissue type'). Convertiti i
file `.md` e `.py`. Non convertiti di proposito: `CRITERI.md` (resta com'era), `core/stats.py`
(13 righe di commenti con l'apostrofo: il file non è stato toccato), gli output grezzi storici
in `results/` e in `logs/`. La baseline `tests/fixtures/tcr_regression_baseline.json` cambia
solo in 4 righe di testo narrativo; i numeri sono identici.

Output rigenerati con i testi attuali: audit del disegno sui tre dataset esterni, A1 sul
dataset completo (stessi numeri: 0.778, 0.895, +0.117; 419 s), disegno su GSE278694.

### 4. Igiene

- **Percorsi assoluti nei file attuali:** erano in 10 file (4 script o test con il percorso dei
  dati di GSE278694, 6 output versionati, di cui uno con 200 righe di avvisi). Rimossi: gli
  script leggono ora `GSE278694_DIR` oppure una cartella `pdac-ml` accanto alla repo; negli
  output i percorsi sono relativi o sostituiti da `<tmp>`. Ricerca finale della cartella home,
  della cartella temporanea e del nome utente nei file tracciati: 0 risultati.
- **Cronologia** (non riscritta): il percorso della home, con il nome utente, compare in 6
  commit (`8bb3d1a`, `2f03461`, `61e9ab9`, `aa3af26`, `bbe7bd0`, `0662a40`); un percorso
  temporaneo in `61e9ab9`; il nome di un istituto in `NOTE.md`, dal commit `2ba6b89` fino a
  questa rifinitura. L'indirizzo email dell'autore non compare in nessun file, ma è nei metadati
  di tutti i commit. Nessun token, chiave o password, né nei file né nella cronologia.
- **Nomi di istituti:** in `NOTE.md` la «valutazione strategica» citava un istituto per nome;
  ora dice «per un laboratorio partner». Nessun altro nome di istituto o di persona nei file.
- **File più pesanti:** `data/demo/GSE125449_demo.h5ad` 8,4 MB; `uv.lock` 378 kB;
  `docs/img/report_demo.png` 231 kB; `validation_esterna/results/A3_perm.txt` 117 kB;
  `docs/sviluppo/REPORT_INTERVENTI.md` 65 kB; `validation_esterna/results/leakage_GSE125449.html`
  64 kB; `core/tcr_validation.py` 39 kB; `core/design_audit.py` 34 kB; `core/report.py` 32 kB;
  `app.py` 31 kB.
- `data/demo/README.md`: accessioni, articoli, riduzione, dicitura sui dati pubblici
  ridistribuiti.
- `REPORT_INTERVENTI.md`, `NOTE.md` e `logs/` spostati in `docs/sviluppo/`.
- `LICENSE` MIT; `pyproject.toml` con autore e licenza.

### Suite

- Prima (clone pulito del commit `0662a40`): 97 passed, 2 xfailed.
- Dopo: 115 passed, 2 xfailed, 0 failed, 0 skipped in 666 s (`docs/sviluppo/logs/pytest_rifinitura.txt`). I test in più sono i 18 di
  `tests/test_verdict.py`. I 2 xfail sono le calibrazioni conservative già note. Su una macchina
  senza i dati di GSE278694 i 3 test di `test_real_data_gse278694.py` vengono saltati.

### Limiti e punti aperti

- Il riferimento completo dell'articolo di GSE278694 (Chen et al., *Cancer Cell* 2025) è
  riportato come nel progetto di tesi: titolo e pagine sono da completare.
- A3 e A4 sul dataset completo non sono stati rieseguiti dopo le ultime modifiche ai messaggi
  (sono stati eseguiti sulla coorte ridotta).
- Il motore precedente sulla coorte ridotta è stato interrotto dopo 45 minuti: l'esito dei
  criteri M1 con il codice di allora resta non determinato.
- La cronologia pubblica contiene ancora i percorsi e il nome dell'istituto indicati sopra.
