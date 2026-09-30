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
