# Resoconto degli interventi su genomic-audit-mvp

Punto di partenza (commit `cf61e76`, baseline): suite **18 passed, 1 xfailed** in 309 s (`logs/pytest_baseline.txt`).
Tutti i numeri qui sotto sono misurati; i log completi di pytest sono in `logs/`.

## Intervento 1 — Audit del disegno e del confondimento

### 1. File modificati e creati

- `core/design_audit.py` (NUOVO): modulo dell'audit del disegno. Contiene il Cramér V con la correzione di Bergsma, i fatti strutturali (annidamento, uno-a-uno, esito determinato, unita' tecnica, costante, identificatore), la classificazione dei confronti e la scomposizione della varianza su pseudobulk (implementata e testata, NON esposta).
- `core/synthetic.py` (ESISTENTE, solo aggiunte in coda): generatori di metadati sintetici. Ci sono gli scenari indipendente, annidato, confuso con l'esito e con associazione a V noto, la struttura di GSE278694, le unita' con quote di varianza note e un AnnData per la pseudobulk. Le funzioni esistenti non sono state toccate.
- `core/report.py` (ESISTENTE): nuova sezione `_design_section`; `render_report` e `save_report` accettano `design_result` come parametro opzionale in coda, quindi le chiamate esistenti restano compatibili. Una frase del footer e' stata estesa.
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
[GSE278694] Confronto 'tissue': Tumor vs Adjacent_normal: stimabile con bassa potenza, disegno appaiato (entro paziente), 5 pazienti con entrambi i livelli (soglia dichiarata: 5 pazienti, la stessa del Modulo B; inoltre con 5 pazienti appaiati anche un test di Wilcoxon esatto non puo' scendere sotto p = 0.062). Un risultato non significativo qui non indica assenza di effetto. Attenzione: 16 pazienti con un solo livello esclusi dal confronto appaiato.
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
| Caso GSE278694 | — | tumore/adiacente: 5 unita', "stimabile con bassa potenza"; scRNA/snRNA: 0 unita', "non stimabile" | come da prompt | passa |

Storia della calibrazione della scomposizione della varianza (due correzioni del metodo, banda mai modificata):
1. Prima esecuzione: distorsione della quota paziente -0.044 (fuori dalla tolleranza ±0.03); copertura 0/200 per entrambe le quote. Le cause: ω² con divisore N anche per il fattore casuale, e un bootstrap percentile spostato verso il basso dai pazienti duplicati.
2. Correzione 1: componenti con il metodo dei momenti (divisore N·(L-1)/L per i fattori casuali, N per quelli fissi), quota come rapporto fra le medie sui geni, intervallo "basic". Risultato: distorsione 0.000, ma il bootstrap restava spostato di piu' della sua ampiezza.
3. Correzione 2: intervallo normale, stima ± 1.96 · SE bootstrap. Risultato: paziente 0.965 (in banda), tessuto 1.000 (sopra la banda: intervallo conservativo). Il test e' marcato xfail(strict=True) e la scomposizione non e' esposta in CLI, app o report.

### 3. Output reale

`python cli.py design --meta gse_like.csv --patient-col patient --tissue-col tissue --technical protocol=protocol --technical library=library --compare tissue:Tumor:Adjacent_normal --compare protocol:scRNA:snRNA`, dove `gse_like.csv` e' generato da `make_gse278694_like_sheet()`:

```
Audit del disegno su 37 unita' (righe dei metadati) e 4 fattori. 2 fatti strutturali rilevati (annidamenti, coincidenze, esiti determinati da un fattore): Nel disegno attuale il fattore 'patient' e' annidato in 'protocol': ogni livello di 'patient' compare in un solo livello di 'protocol'. Un confronto fra livelli di 'protocol' non distingue l'effetto di 'protocol' da quello dei livelli di 'patient' che contiene. Ogni livello di 'library' corrisponde a una sola coppia paziente-tessuto e viceversa: l'effetto tecnico di 'library' non e' separabile da quello della coppia paziente-tessuto. E' la norma senza multiplexing o librerie replicate, ma significa che ogni differenza fra tessuti dello stesso paziente include anche la differenza fra due librerie. 0 coppie di fattori con associazione forte (Cramér V >= 0.5). Confronto 'tissue': Tumor vs Adjacent_normal: stimabile con bassa potenza, disegno appaiato (entro paziente), 5 pazienti con entrambi i livelli (soglia dichiarata: 5 pazienti, la stessa del Modulo B; inoltre con 5 pazienti appaiati anche un test di Wilcoxon esatto non puo' scendere sotto p = 0.062). Un risultato non significativo qui non indica assenza di effetto. Attenzione: 16 pazienti con un solo livello esclusi dal confronto appaiato. Confronto 'protocol': scRNA vs snRNA: non stimabile -- nessun paziente ha entrambi i livelli: la differenza fra 'scRNA' e 'snRNA' coincide con la differenza fra due gruppi di pazienti diversi (fattore confuso con il paziente).
 - 'patient' e 'tissue': Cramér V non valutabile (37 unita' per 22x3 celle; servono almeno 10 unita' e 2 unita' attese per cella). Un numero su una tabella cosi' piccola sarebbe instabile.
 - 'patient' e 'protocol': Cramér V non valutabile (37 unita' per 22x2 celle; servono almeno 10 unita' e 2 unita' attese per cella). Un numero su una tabella cosi' piccola sarebbe instabile.
 - 'tissue' e 'protocol': Cramér V corretto = 0.40 su 37 unita', sotto la soglia di 0.5: nessuna segnalazione (non significa che il disegno sia perfettamente bilanciato).
 * Nota fissa sui protocolli: lo snRNA-seq (nuclei) sottorappresenta le cellule immunitarie e cambia la composizione misurata rispetto allo scRNA-seq; i protocolli per tessuto fissato (FFPE, 10x Flex) misurano l'RNA tramite sonde e danno profili non direttamente sovrapponibili a quelli da tessuto fresco. Un confronto fra protocolli e' interpretabile solo se gli stessi pazienti sono misurati con entrambi.
 * Le unita' di questo audit sono le righe della tabella dei metadati (campioni o librerie), non le cellule: le cellule dello stesso campione non sono osservazioni indipendenti del disegno.
[ok] report scritto in results/design_report_gse_like.html
```

### 4. Decisioni e assunzioni

- **Unita' dell'audit = righe dei metadati.** Un CSV e' preso cosi' com'e' (una riga per campione o libreria). Un `adata.obs` a livello di cellula viene ridotto alle combinazioni distinte dei fattori scelti, perche' le cellule non sono unita' indipendenti del disegno.
- **Quando il V e' valutabile:** almeno 10 righe e almeno 2 righe attese per cella in media (n/(r·c) ≥ 2). Sotto questa soglia l'output e' "non valutabile", mai un numero. La soglia e' stata fissata prima del test.
- **Soglia di allarme sul V = 0.5**, come da prompt. L'allarme non scatta se la coppia ha gia' un fatto strutturale, che viene riportato al suo posto.
- **Fattori costanti o identificativi di riga** (un livello per riga, es. `library` in un foglio per libreria) vengono esclusi dalle coppie. Il controllo "unita' tecnica" verifica se un fattore tecnico corrisponde 1:1 alla coppia paziente-tessuto.
- **"Esito determinato"** si applica quando il fattore esterno e' il tessuto o una colonna di esito. Per gli altri fattori il fatto si chiama "annidamento".
- **Classificazione dei confronti.** Con almeno un paziente che ha entrambi i livelli, il disegno e' appaiato e le unita' sono quei pazienti. Senza pazienti con entrambi i livelli: un confronto su tessuto o su un fattore tecnico e' "non stimabile" (confuso con il paziente); una colonna dichiarata come esito diventa un confronto fra pazienti, con unita' pari ai pazienti del gruppo piu' piccolo. Qualunque fattore tecnico che separi perfettamente i due livelli rende il confronto "non stimabile".
- **DEVIAZIONE sulla soglia di bassa potenza.** Il prompt definisce "bassa potenza" come meno di 5 unita', ma si aspetta che GSE278694 (5 pazienti) risulti a bassa potenza, e le due richieste si contraddicono. Ho applicato la lettura piu' prudente: bassa potenza se le unita' sono meno di 5 (la soglia del Modulo B, dichiarata) OPPURE se il p-value minimo raggiungibile da un test esatto supera 0.05. Nel disegno appaiato quel minimo e' quello di Wilcoxon, `core.stats.wilcoxon_min_pvalue`, pari a 0.0625 con 5 coppie; fra pazienti e' quello di Mann-Whitney, 2/C(n_a+n_b, n_a). Con 5 coppie risulta quindi "bassa potenza", con 6 "stimabile".
- **Scomposizione della varianza:** il paziente (e il batch) e' trattato come fattore casuale, il tessuto come fisso. Per l'intervallo servono almeno 5 pazienti. La scomposizione viene rifiutata se piu' del 10% delle repliche bootstrap non e' identificabile. Le unita' con meno di 20 cellule vengono escluse; si usano 1000 geni ad alta varianza, con CPM+log1p tramite `LogCPM(target_sum=1e6)` del Modulo A.
- **Scomposizione della varianza NON esposta** (regola 2 di "Quando fermarti"): e' presente in `core/`, tolta dall'orchestratore `run_design_audit` e annotata in `NOTE.md`.
- La nota fissa sui protocolli compare solo se viene dichiarata una colonna con ruolo `protocol`.
- **Baseline committata:** il repository non aveva commit. Ho creato il commit `cf61e76` con lo stato iniziale, cosi' che ogni intervento sia un diff revisionabile.

### 5. Bug o anomalie

- Nel primo test funzionale cercavo la parola "assente" nella frase. Il codice usa "non compare": ho corretto il test, non il codice.
- Nell'app, la colonna "Cramér V" mescolava numeri e testo ("non valutabile") e faceva fallire la conversione pyarrow (trovato con `streamlit.testing.AppTest`). Ora e' formattata come stringa.
- Un errore di sintassi (apice sbagliato) nel footer di `core/report.py`, introdotto da me, e' stato corretto prima di qualunque esecuzione dei test.
- Nessuna modifica al codice esistente di calcolo: le calibrazioni gia' approvate (`test_stats_calibration.py`) sono state comunque rieseguite nella suite completa, senza variazioni (stesse righe PASSED/XFAIL della baseline).

### 6. Limiti noti

- La scomposizione della varianza non e' disponibile per l'utente (vedi sopra). La sua calibrazione copre un solo disegno: bilanciato, 12 pazienti x 2 tessuti, 300 geni.
- Il V di Cramér e' calibrato solo su uno scenario indipendente (10 pazienti x 8 campioni, fattori con 2-4 livelli) e su tabelle 3x3 per la potenza; la potenza a V = 0.5 e' 0.46, quindi un'associazione moderata viene spesso non segnalata.
- I disegni sbilanciati o parzialmente annidati sono gestiti in modo approssimato (le componenti di varianza non sono piu' esattamente non distorte).
- Un confronto fra pazienti su una colonna di esito presume che l'esito sia davvero a livello di paziente; lo strumento non lo verifica.
- Nel foglio di GSE278694 la struttura della coorte snRNA (8 pazienti, un campione tumorale ciascuno) e' un'approssimazione: i dettagli reali non sono stati verificati.
- L'app e' stata verificata solo in modalita' demo con `AppTest` (nessun upload reale).

---

## Intervento 2 — Flag per cellula dal Modulo B

### 1. File modificati e creati

- `core/tcr_validation.py` (ESISTENTE, **toccato** come consentito dal prompt; solo aggiunte, +91 righe e nessuna riga rimossa):
  - `TcrValidationResult` ha due nuovi campi opzionali in coda, con default: `cell_flags` e `flag_coverage`. I costruttori esistenti restano compatibili.
  - `run_tcr_validation` aggiunge la colonna interna `_row` (posizione di riga originale) prima della merge e calcola i flag quando sono dati `marker_map` e `reference_compartment`. `_row` non entra in nessun calcolo.
  - Nuove funzioni `cell_flags` (stessi filtri di `marker_error_rate`) ed `export_audited` (copia `*_audited.h5ad` + CSV), piu' `FLAG_COLUMNS` e `FLAG_DEFINITIONS`.
  - `pairwise_excess_discordance`, `assign_reference_identity`, `marker_error_rate` e i narrativi NON sono stati modificati.
- `cli.py` (ESISTENTE): opzione `tcr --export-flags` e stampa della copertura dei flag.
- `core/report.py` (ESISTENTE): sottosezione "Flag per cellula" nel Modulo B, con la copertura.
- `app.py` (ESISTENTE): copertura, download del CSV dei flag e della copia `.h5ad`.
- `tests/test_tcr_flags.py` (NUOVO), `tests/fixtures/tcr_regression.py` e `tests/fixtures/__init__.py` (NUOVI).
- `tests/fixtures/tcr_regression_baseline.json` (NUOVO): output di `run_tcr_validation` su 4 configurazioni sintetiche, generato dal codice del commit di baseline PRIMA di modificare `tcr_validation.py`. Prima di generarlo ho verificato con `git diff cf61e76 -- core/tcr_validation.py core/stats.py` che non ci fossero differenze.
- `core/stats.py`: NON toccato.

### 2. Test

Suite completa dopo l'Intervento 2: **46 passed, 2 xfailed** in 475 s (`logs/pytest_intervento2.txt`), contro 36 passed e 2 xfailed dopo l'Intervento 1. Nessun test preesistente e' cambiato di esito. Le calibrazioni gia' approvate (`test_stats_calibration.py`) e i test del Modulo B (`test_synthetic_data.py`) sono stati rieseguiti dopo la modifica a `tcr_validation.py` e passano tutti.

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

Questo intervento non introduce nessuna nuova stima statistica: i flag riportano, cellula per cellula, la stessa quantita' del tasso d'errore esistente. Il requisito di calibrazione e' quindi sostituito dai due test richiesti.
- **Coerenza:** su 5 dataset sintetici con 3 compartimenti, per ogni compartimento non di riferimento la media di `audit_label_vs_reference` sulle cellule valutabili e' uguale con `==` (uguaglianza esatta in virgola mobile) alla stima puntuale di `marker_error_rate`. Anche il numero di cellule valutabili coincide con `n_obs` del bootstrap.
- **Regressione:** 4 configurazioni (eccesso 0.35 con marcatori; eccesso 0 con marcatori; 3 compartimenti; 4 pazienti senza marcatori). Tutti i numeri, le tabelle per coppia di compartimenti e le frasi narrative sono identici al riferimento generato prima della modifica, con tolleranza 1e-12 (di fatto identici). Come prova di sensibilita', cambiando il solo seme del bootstrap il test fallisce, quindi il confronto non e' banale.

### 3. Output reale

Dataset sintetico (`make_tcr_validation_dataset(n_patients=12, n_clones_per_patient=15, injected_excess=0.35, seed=0)`, salvato come `sintetico.h5ad` con i CSV VDJ e il manifest), comando `python cli.py tcr ... --marker-map markers.json --reference-compartment PBMC --export-flags --n-boot 1000`:

```
Le cellule dello stesso clone T (stessa sequenza CDR3 della catena TRB) cambiano etichetta di tipo cellulare fra compartimenti tissutali con un eccesso di discordanza di +0.107 (IC95% [+0.088, +0.127], su 180 coppie clone-compartimenti da 12 pazienti) rispetto al rumore di base entro lo stesso compartimento: un effetto reale (l'intervallo esclude lo zero). Tasso d'errore dell'annotazione rispetto all'identita' clonale dai marcatori (riferimento: compartimento 'PBMC', 180 cloni risolti): 'Tumor': 0.283 (IC95% [0.250, 0.312], 12 pazienti).
Flag per cellula valutabili (audit_label_vs_reference non NA): 39.6% delle cellule.
[ok] copia con i flag scritta in out/sintetico_audited.h5ad (il file originale non e' toccato)
[ok] flag per cellula scritti in out/sintetico_audit_flags.csv
[ok] report scritto in out/tcr_report.html
```

Prime righe di `sintetico_audit_flags.csv` (le prime cellule sono nel sangue, quindi `audit_label_vs_reference` e' vuoto, cioe' NA):

```
obs_name,audit_reference_label,audit_label_vs_reference
cell0,CD8T,
cell1,CD8T,
cell2,CD8T,
cell3,CD8T,
```

Nella copia `sintetico_audited.h5ad`, `audit_label_vs_reference` e' di tipo `BooleanDtype`, con NA 1381, False 648 e True 256 su 2285 cellule. Le colonne originali hanno gli stessi valori. Frase generata nel report HTML (sezione "Flag per cellula"): "Copertura: 2,285 cellule su 2,285 (100.0%) hanno un'identita' di riferimento; 904 (39.6%) sono valutabili, e di queste 256 sono discordanti. Tutte le altre sono NA: lo strumento non poteva verificarle (…). NA non significa "corretta"."

### 4. Decisioni e assunzioni

- **Cellule del compartimento di riferimento:** `audit_label_vs_reference` = NA, anche quando il clone ha un riferimento. `marker_error_rate` le esclude perche' il riferimento e' stimato proprio da loro, e il prompt chiede la stessa logica. `audit_reference_label` invece e' valorizzato anche per loro.
- `audit_reference_label` e' valorizzato per ogni cellula il cui clone ha un riferimento, compresa un'etichetta come "NK": si tratta dell'identita' del clone, non di un giudizio sulla cellula.
- **Copertura riportata** = frazione di TUTTE le cellule dell'AnnData con `audit_label_vs_reference` non NA. Sui dati reali, dove il sangue e le cellule senza TCR sono molte, sara' bassa: e' voluto.
- **Export:** si rifiuta se `adata.obs` ha gia' colonne con quei nomi (nessuna sovrascrittura) o se l'indice dei flag non coincide con `obs_names`. Le definizioni vanno in `uns['genomic_audit_flags']`; il CSV contiene solo `obs_name` e i due flag, senza righe di commento, cosi' resta leggibile da qualunque parser. Da CLI i file vanno nella cartella di `--out`, con il nome del file `.h5ad` di input.
- Nessun flag dei doppietti, come da prompt.

### 5. Bug o anomalie

- Nessun bug trovato nel codice esistente. Un'anomalia di formato: scrivendo la copia, anndata converte le colonne di testo in categoriali; i valori non cambiano. Il test verifica l'uguaglianza dei valori come stringhe, non del tipo.
- I 10 test nuovi sono passati al primo tentativo. Per non fidarmi di un test di regressione potenzialmente banale, ho fatto la verifica di sensibilita' descritta sopra.

### 6. Limiti noti

- I flag hanno la stessa definizione di identita' di riferimento del Modulo B (≥ 3 cellule nel sangue, margine ≥ 0.20, marcatori "conta > 0"). Le cellule di cloni piccoli restano NA.
- Non ho testato AnnData con `obs_names` duplicati, ne' file `.h5ad` aperti in modalita' `backed`.
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

Notazione: per una cellula vera CD4, a4 = P(chiamata CD4), b4 = P(chiamata CD8), o4 = P(chiamata altro); per una vera CD8, a8, b8, o8. Le probabilita' sono stimate sulle cellule del compartimento con identita' di riferimento, pooled fra pazienti. Per il paziente, r e' la frazione CD8 osservata, cioe' n_CD8 / (n_CD4 + n_CD8). Con x e y le cellule vere CD4 e CD8:

    n_CD4 = a4·x + a8·y        n_CD8 = b4·x + b8·y

Da cui, con (n_CD4, n_CD8) ∝ (1 − r, r):

    x ∝ b8·(1 − r) − a8·r      y ∝ a4·r − b4·(1 − r)      frazione vera = y / (x + y), troncata a [0, 1]

Trattamento di "altro": le cellule chiamate "altro" non entrano ne' nel numeratore ne' nel denominatore della frazione riportata. La loro perdita, diversa fra CD4 e CD8, e' gia' contenuta nelle a e b, che per riga sommano a 1 − o. Non serve quindi una categoria "altro" nell'inversione, e le cellule non-T chiamate CD4/CD8 non sono modellate.

Condizionamento: J = b8/(a8+b8) − b4/(a4+b4). Rifiuto se J < 0.2, oppure se il 2.5° percentile di J nel bootstrap e' ≤ 0.

Scenari: l'errore di riga, 1 − chiamata corretta, e' moltiplicato per k ∈ {0.5, 1, 2}, mantenendo la proporzione fra "chiamata sbagliata" e "altro".

Incertezza: le repliche della matrice ricampionano i pazienti con lo schema di `core.stats.cluster_bootstrap`; per ciascuna si estrae r* da Beta(n_CD8 + ½, n_CD4 + ½) e si inverte. L'intervallo e' dato dai percentili 2.5–97.5.

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
Le cellule dello stesso clone T (stessa sequenza CDR3 della catena TRB) cambiano etichetta di tipo cellulare fra compartimenti tissutali con un eccesso di discordanza di +0.107 (IC95% [+0.088, +0.127], su 180 coppie clone-compartimenti da 12 pazienti) rispetto al rumore di base entro lo stesso compartimento: un effetto reale (l'intervallo esclude lo zero). Tasso d'errore dell'annotazione rispetto all'identita' clonale dai marcatori (riferimento: compartimento 'PBMC', 180 cloni risolti): 'Tumor': 0.283 (IC95% [0.250, 0.312], 12 pazienti).
Flag per cellula valutabili (audit_label_vs_reference non NA): 39.6% delle cellule.
Matrice di confusione stimata su 1131 cellule con identita' di riferimento da 12 pazienti (pooled fra pazienti), J = 0.43. Frazione CD8 nel compartimento 'Tumor' del paziente PT000 (CD8T sul totale delle cellule etichettate CD4T o CD8T, n = 71): riportata 0.56 (IC95% dei soli conteggi 0.44-0.67, senza correzione). Intervallo plausibile al 95% tenendo conto dell'errore di annotazione misurato, per scenario di errore sui cloni non condivisi con il sangue -- 0.5x: fra 0.43 e 0.74; 1x: fra 0.37 e 0.92; 2x: non prodotto (matrice di confusione mal condizionata (J = -0.51; servono J >= 0.2 e un intervallo di J che escluda lo zero)). Nessuno scenario e' "il risultato": la loro distanza mostra quanto la conclusione dipende dall'assunzione. Frazione CD8 nel compartimento 'Tumor' del paziente PT001 (CD8T sul totale delle cellule etichettate CD4T o CD8T, n = 74): riportata 0.43 (IC95% dei soli conteggi 0.31-0.55, senza correzione). Intervallo plausibile al 95% tenendo conto dell'errore di annotazione misurato, per scenario di errore sui cloni non condivisi con il sangue -- 0.5x: fra 0.25 e 0.57; 1x: fra 0.05 e 0.61; 2x: non prodotto (matrice di confusione mal condizionata (J = -0.51; servono J >= 0.2 e un intervallo di J che escluda lo zero)). Nessuno scenario e' "il risultato": la loro distanza mostra quanto la conclusione dipende dall'assunzione. Frazione CD8 nel compartimento 'Tumor' del paziente PT002 (CD8T sul totale delle cellule etichettate CD4T o CD8T, n = 70): riportata 0.47 (IC95% dei soli conteggi 0.36-0.59, senza correzione). Intervallo plausibile al 95% tenendo conto dell'errore di annotazione misurato, per scenario di errore sui cloni non condivisi con il sangue -- 0.5x: fra 0.31 e 0.62; 1x: fra 0.18 e 0.72; 2x: non prodotto (matrice di confusione mal condizionata (J = -0.51; servono J >= 0.2 e un intervallo di J che escluda lo zero)). Nessuno scenario e' "il risultato": la loro distanza mostra quanto la conclusione dipende dall'assunzione. Frazione CD8 nel compartimento 'Tumor' del paziente PT003 (CD8T sul totale delle cellule etichettate CD4T o CD8T, n = 74): riportata 0.65 (IC95% dei soli conteggi 0.54-0.75, senza correzione). Intervallo plausibile al 95% tenendo conto dell'errore di annotazione misurato, per scenario di errore sui cloni non condivisi con il sangue -- 0.5x: fra 0.56 e 0.83; 1x: fra 0.60 e 1.00; 2x: non prodotto (matrice di confusione mal condizionata (J = -0.51; servono J >= 0.2 e un intervallo di J che escluda lo zero)). Nessuno scenario e' "il risultato": la loro distanza mostra quanto la conclusione dipende dall'assunzione. Frazione CD8 nel compartimento 'Tumor' del paziente PT004 (CD8T sul totale delle cellule etichettate CD4T o CD8T, n = 59): riportata 0.49 (IC95% dei soli conteggi 0.37-0.62, senza correzione). Intervallo plausibile al 95% tenendo conto dell'errore di annotazione misurato, per scenario di errore sui cloni non condivisi con il sangue -- 0.5x: fra 0.33 e 0.66; 1x: fra 0.20 e 0.79; 2x: non prodotto (matrice di confusione mal condizionata (J = -0.51; servono J >= 0.2 e un intervallo di J che escluda lo zero)). Nessuno scenario e' "il risultato": la loro distanza mostra quanto la conclusione dipende dall'assunzione. Frazione CD8 nel compartimento 'Tumor' del paziente PT005 (CD8T sul totale delle cellule etichettate CD4T o CD8T, n = 90): riportata 0.51 (IC95% dei soli conteggi 0.41-0.60, senza correzione). Intervallo plausibile al 95% tenendo conto dell'errore di annotazione misurato, per scenario di errore sui cloni non condivisi con il sangue -- 0.5x: fra 0.38 e 0.65; 1x: fra 0.29 e 0.76; 2x: non prodotto (matrice di confusione mal condizionata (J = -0.51; servono J >= 0.2 e un intervallo di J che escluda lo zero)). Nessuno scenario e' "il risultato": la loro distanza mostra quanto la conclusione dipende dall'assunzione. Frazione CD8 nel compartimento 'Tumor' del paziente PT006 (CD8T sul totale delle cellule etichettate CD4T o CD8T, n = 80): riportata 0.56 (IC95% dei soli conteggi 0.46-0.66, senza correzione). Intervallo plausibile al 95% tenendo conto dell'errore di annotazione misurato, per scenario di errore sui cloni non condivisi con il sangue -- 0.5x: fra 0.45 e 0.72; 1x: fra 0.41 e 0.90; 2x: non prodotto (matrice di confusione mal condizionata (J = -0.51; servono J >= 0.2 e un intervallo di J che escluda lo zero)). Nessuno scenario e' "il risultato": la loro distanza mostra quanto la conclusione dipende dall'assunzione. Frazione CD8 nel compartimento 'Tumor' del paziente PT007 (CD8T sul totale delle cellule etichettate CD4T o CD8T, n = 74): riportata 0.46 (IC95% dei soli conteggi 0.35-0.57, senza correzione). Intervallo plausibile al 95% tenendo conto dell'errore di annotazione misurato, per scenario di errore sui cloni non condivisi con il sangue -- 0.5x: fra 0.30 e 0.60; 1x: fra 0.14 e 0.68; 2x: non prodotto (matrice di confusione mal condizionata (J = -0.51; servono J >= 0.2 e un intervallo di J che escluda lo zero)). Nessuno scenario e' "il risultato": la loro distanza mostra quanto la conclusione dipende dall'assunzione. Frazione CD8 nel compartimento 'Tumor' del paziente PT008 (CD8T sul totale delle cellule etichettate CD4T o CD8T, n = 89): riportata 0.45 (IC95% dei soli conteggi 0.35-0.56, senza correzione). Intervallo plausibile al 95% tenendo conto dell'errore di annotazione misurato, per scenario di errore sui cloni non condivisi con il sangue -- 0.5x: fra 0.30 e 0.58; 1x: fra 0.15 e 0.64; 2x: non prodotto (matrice di confusione mal condizionata (J = -0.51; servono J >= 0.2 e un intervallo di J che escluda lo zero)). Nessuno scenario e' "il risultato": la loro distanza mostra quanto la conclusione dipende dall'assunzione. Frazione CD8 nel compartimento 'Tumor' del paziente PT009 (CD8T sul totale delle cellule etichettate CD4T o CD8T, n = 77): riportata 0.48 (IC95% dei soli conteggi 0.37-0.60, senza correzione). Intervallo plausibile al 95% tenendo conto dell'errore di annotazione misurato, per scenario di errore sui cloni non condivisi con il sangue -- 0.5x: fra 0.33 e 0.63; 1x: fra 0.21 e 0.75; 2x: non prodotto (matrice di confusione mal condizionata (J = -0.51; servono J >= 0.2 e un intervallo di J che escluda lo zero)). Nessuno scenario e' "il risultato": la loro distanza mostra quanto la conclusione dipende dall'assunzione. Frazione CD8 nel compartimento 'Tumor' del paziente PT010 (CD8T sul totale delle cellule etichettate CD4T o CD8T, n = 70): riportata 0.40 (IC95% dei soli conteggi 0.29-0.51, senza correzione). Intervallo plausibile al 95% tenendo conto dell'errore di annotazione misurato, per scenario di errore sui cloni non condivisi con il sangue -- 0.5x: fra 0.22 e 0.52; 1x: fra 0.01 e 0.53; 2x: non prodotto (matrice di confusione mal condizionata (J = -0.51; servono J >= 0.2 e un intervallo di J che escluda lo zero)). Nessuno scenario e' "il risultato": la loro distanza mostra quanto la conclusione dipende dall'assunzione. Frazione CD8 nel compartimento 'Tumor' del paziente PT011 (CD8T sul totale delle cellule etichettate CD4T o CD8T, n = 76): riportata 0.36 (IC95% dei soli conteggi 0.25-0.46, senza correzione). Intervallo plausibile al 95% tenendo conto dell'errore di annotazione misurato, per scenario di errore sui cloni non condivisi con il sangue -- 0.5x: fra 0.17 e 0.45; 1x: fra 0.00 e 0.43; 2x: non prodotto (matrice di confusione mal condizionata (J = -0.51; servono J >= 0.2 e un intervallo di J che escluda lo zero)). Nessuno scenario e' "il risultato": la loro distanza mostra quanto la conclusione dipende dall'assunzione. Assunzioni. (1) L'identita' di riferimento stimata dal sangue e' considerata corretta. (2) La matrice di confusione e' stimata solo sui cloni condivisi con il sangue, che per costruzione sono quelli espansi (almeno 3 cellule nel sangue e identita' netta), ed e' unica per tutti i pazienti (pooled), non per paziente. Nei dati originali (PDAC, GSE278694) l'errore di annotazione cresce con la dimensione del clone (rho = +0.177): la direzione plausibile del bias e' quindi una SOVRASTIMA dell'errore sui cloni non condivisi, ma non e' stata verificata, perche' quelle cellule non hanno riferimento. Per questo sono riportati sempre tre scenari (errore 0.5x, 1x, 2x), nessuno dei quali e' "il risultato". (3) Le cellule etichettate CD4 o CD8 sono vere cellule T; doppietti e altre cellule non-T etichettate CD4/CD8 non sono modellati.
[ok] report scritto in out3/tcr_cd8_report.html
```


### 4. Decisioni e assunzioni

- **La matrice e' stimata nel compartimento bersaglio**, non in tutti i compartimenti: l'errore dipende dal compartimento (0.195 nel tumore contro 0.059 nel tessuto adiacente in `pdac-ml`).
- **Denominatore:** cellule etichettate CD4T o CD8T nel compartimento bersaglio di quel paziente. E' dichiarato nella frase e nel report; le etichette sono configurabili (`--cd4-label`, `--cd8-label`).
- **Ricampionamento dei pazienti.** `core.stats.cluster_bootstrap` restituisce solo media e IC, non le repliche, e non va modificata. Ho quindi scritto `patient_resample_draws`, che ne riproduce lo schema chiamata per chiamata. Un test verifica che la media sulle stesse repliche dia esattamente (1e-12) lo stesso IC di `cluster_bootstrap`. E' una **deviazione dalla lettera** ("importali"), ma non dalla sostanza: lo stimatore validato non e' stato modificato ne' reimplementato in modo diverso.
- **Monte Carlo:** una frazione estratta da Beta(n_CD8 + ½, n_CD4 + ½), cioe' con prior di Jeffreys, per ogni replica della matrice; 1000 repliche di default (`--n-boot`).
- **Soglie dichiarate:** J minimo 0.2; almeno 20 cellule CD4+CD8 per paziente; rifiuto se oltre il 10% delle repliche e' non valido (matrice senza una delle due identita', errore scalato ≥ 1, inversione non definita); almeno 5 pazienti con cellule di riferimento nel compartimento (la stessa soglia del Modulo B).
- **Troncamento a [0, 1]** delle frazioni invertite. Vicino ai bordi l'intervallo puo' toccare 0 o 1.
- **Il rifiuto per condizionamento vale per singolo scenario:** con errore 2x lo scenario puo' essere rifiutato mentre 0.5x e 1x vengono prodotti. Succede nel dataset di esempio, dove l'errore nel tumore e' circa il 42% e raddoppiato rende J negativo.
- **Nessuno scenario e' indicato come "il risultato"** e il testo delle assunzioni compare sempre accanto agli intervalli, come chiesto.

### 5. Bug o anomalie

- **Motivo di rifiuto sbagliato (corretto).** Nel primo output reale lo scenario 2x veniva rifiutato con "le probabilita' superano 1", mentre la causa vera era J = −0.51 (inversione mal condizionata): bastava una sola replica bootstrap con errore raddoppiato ≥ 1 per far scattare il rifiuto meno informativo. Ho riordinato i controlli: stima puntuale, poi J, poi repliche non valide (rifiuto solo oltre il 10%). Ho anche allineato il numero di estrazioni Beta alle repliche rimaste. Le calibrazioni sono state rieseguite dopo la correzione e danno numeri identici (1x non e' toccato).
- Nessuna correzione a codice esistente in questo intervento.

### 6. Limiti noti

- La calibrazione copre un solo disegno: 10 pazienti, 150–300 cellule T ciascuno, 30% con riferimento scelto a caso, frazione vera fra 0.2 e 0.7. In simulazione l'assunzione chiave (cloni con riferimento rappresentativi) e' vera per costruzione: la copertura con cloni NON rappresentativi (errore che dipende dalla dimensione del clone) non e' stata misurata. E' proprio il caso che gli scenari 0.5x/2x dovrebbero coprire, ma senza garanzia di copertura.
- Matrice pooled: se l'errore varia molto fra pazienti, l'intervallo del singolo paziente puo' essere troppo stretto. Non misurato.
- Le cellule non-T etichettate CD4/CD8 (es. doppietti) non sono modellate.
- Con errori grandi (come nel dataset di esempio, circa 42%) gli intervalli 1x sono molto larghi (es. 0.37–0.92): e' l'informazione corretta, ma va detto a chi legge.
- L'app ricalcola gli intervalli a ogni clic, senza cache; non misurato su dataset reali grandi.

---

## Riepilogo finale

- **Commit:** `cf61e76` baseline (stato iniziale, prima di ogni modifica); `4079895` Intervento 1; `e56d1b9` Intervento 2; il commit "Intervento 3" contiene questo resoconto (vedi `git log`).
- **Suite:** 18 passed / 1 xfailed alla baseline → 36 / 2 → 46 / 2 → **53 passed / 2 xfailed**. Gli xfail sono quello preesistente (Nadeau-Bengio a k=5, conservativo) e quello nuovo (`test_variance_share_bootstrap_coverage_tissue`, conservativo: copertura 1.000 > 0.99).
- **Stime NON esposte per calibrazione fuori banda:** la scomposizione della varianza dell'audit del disegno.
- **`core/stats.py` non e' mai stato modificato** (`git diff cf61e76 -- core/stats.py` e' vuoto). `core/tcr_validation.py` ha solo aggiunte (Intervento 2); l'output del Modulo B e' identico a prima della modifica (test di regressione).
- `README.md` aggiornato con cosa fa e cosa non fa ogni modulo; `NOTE.md` con le idee scartate o rimandate.
