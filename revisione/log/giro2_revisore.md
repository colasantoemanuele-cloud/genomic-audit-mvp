# Resoconto del revisore indipendente — giro 2 (chiusura)

Repo `genomic-audit-mvp`, ramo `revisione`, commit `df4f539`. L'albero era pulito prima della revisione ed è
pulito dopo (`git status --short` vuoto). Non ho modificato nessun file della repo né di `../pdac-ml`.
Script e output delle prove sono in
`/tmp/claude-1000/-home-nemo-Uni-PROGETTO-ML-genomic-audit-mvp/22405c93-e5d2-4170-9aed-2c2f9b779264/scratchpad/revisore2/`:
`partial_conf_report.py`, `keys_mismatch.py`, `labels_typo.py`, `disc_zero_var.py`, `sim_v_shapes.py`,
`sim_v_unbal.py`, `sim_v_edge.py`, `cd8_few.py` e `gap_null.py`, ciascuno con il suo `.txt`, più
`pytest_fast.txt`. I comandi che usano tutti i core sono stati lanciati uno alla volta.

## Esiti

- `uv run pytest -n auto -m "not slow" -q`: **130 passed, 62 warnings in 126,37 s**, exit 0. Nessun test
  saltato: i dati di GSE278694 sono presenti, quindi anche i 3 test sui dati reali sono girati.
  Raccolta: 199 test, di cui 130 «not slow» e 69 «slow».
- `uv run ruff check .`: **All checks passed!**
- Log della suite completa sul clone pulito
  (`.../scratchpad/clone_finale_log/clone_uv_suite_completa.txt`: clone uv del commit df4f539, con
  `GSE278694_DIR`): demo con exit 0 in 97 s; suite completa **192 passed, 7 xfailed, 64 warnings in
  1520,89 s (25 min 20 s)**, nessun fallimento e nessun test saltato. I 7 xfail sono quelli dichiarati:
  NB a k=5 con i DummyClassifier, NB in LOGO a 8/10/15/20, Wilcoxon in LOGO a 10, copertura della quota
  del tessuto. Una nota sul criterio 6, che non conto come problema (MEDIA 5 non va riproposto): il log
  copre solo la suite **con** i dati reali. Le righe del README «senza i dati reali» (pip, Python 3.12)
  vanno o rimisurate sul commit finale o tolte.

**Conteggi: ALTA 1 · MEDIA 2 · BASSA 8.**

### Verifiche passate

Le riporto per completezza, non sono problemi.

- **Correzione 5 (pareggi).** Il default esclude i pareggi (`core/tcr_validation.py:540-552`, `TIE_RULE_DEFAULT
  = "exclude"`) e le sensibilità sono disponibili con `--tie-rule error|correct` (`cli.py:352`). La tesi in
  `../pdac-ml` (ramo `revisione`, HEAD `05947d0`, verificato con `git log`) dà
  `results/tables/error_rate_blood.csv: Tumor,277,10,10,0.16967509025270758`. Lo strumento dà
  0,16967509025270758 su 277 unità (`validation/results/cli_tcr_real_conventions_tesi.txt`,
  `a1_moduleB.json` S5, `tests/test_real_data_gse278694.py:61`, passato nel mio run). L'intervallo è diverso
  ([0,100; 0,274] contro [0,095; 0,252]) ed è dichiarato. Il vecchio 0,1951 è indicato come superato in
  README:109 e :124, in `validation_esterna/REPORT.md:30`, nel test e nello script di validazione.
  `docs/sviluppo/REPORT_INTERVENTI.md` ha la nota «superato» in testa.
- **Numeri del README** ritrovati nei file di risultati prodotti dal codice attuale:
  - disegno: 37 unità (`cli_design_real.txt`, rigenerato in df4f539);
  - Modulo A: 0.778 ± 0.060, 0.895 ± 0.004, +0.117 (`A1_cli_stdout.txt`); 9.946 cellule, 19 pazienti,
    8 classi (`A1_cli_stderr.txt`);
  - Modulo B: 0.093 [0.043, 0.188], 0.170 [0.100, 0.274], 277 unità, 11 pareggi
    (`cli_tcr_real_conventions_tesi.txt`);
  - quote NK: 72% e 0,115 (`b2_quota_nk.json`); 0,093 → 0,163 al passo S2 → S3 (0,0927 → 0,1632, in
    `a1_moduleB.json`);
  - fedeltà del porting: 0,201/0,163 su 288, 0,199/0,164 su 287, 0,152 su 231 e 0,154 su 227
    (`a1_moduleB.json` `fedelta_S5`); +0,059 [0,000; 0,141] su 85 cloni (S5);
  - frazione di CD8: 0,45, 7 pazienti su 13 rifiutati nello scenario 1x, 0.5x e 1x disgiunti (per esempio
    PA01 0.44-0.53 contro 0.03-0.35) (`a3_summary.json`, `cli_tcr_real.txt`);
  - calibrazioni: 3,6%-6,3% e 14,5% (`a1_calibrazione.txt`); 0,5-1,8% e 8,3% (`a2_calibrazione.txt`); 0,46
    e 0,0-0,7% (`design_calibrazione.txt`);
  - validazione esterna: 0,7530/0,7127 e 7 fold (`A2_*_ridotta.txt`, `criteri_originali.json`);
  - Streamlit 1.51 (`streamlit_versione_minima.txt`).

  Unica eccezione è la formulazione di MEDIA 1 qui sotto: il numero è giusto, la frase no.
- **Correzioni del giro 1 applicate bene:** ALTA 1, 3, 4; MEDIA 1, 2, 3, 4, 7, 8; BASSA 1, 2, 3, 4 (con la
  regola scelta), 6, 7, 9, 10, 11, 12, 13. Fanno eccezione ALTA 2 (applicata solo nella sintesi: ALTA 1 qui
  sotto), MEDIA 6 (applicata in un solo punto: MEDIA 1 qui sotto) e BASSA 14 (applicata in parte: BASSA 3).
- **Cramér V da 20 righe (ALTA 1 del giro 1)**, con simulazioni mie su forme di tabella che la calibrazione
  non copre (fattori indipendenti):
  - tabelle più grandi (`sim_v_shapes.txt`, 400 repliche per caso): 2×5 a 20 righe 0,050; 3×4 0,007;
    2×10 a 40 righe 0,010; 3×6 a 37 righe 0,003;
  - margini sbilanciati (`sim_v_edge.txt`, 4000 repliche): 2×2 a 20 righe con livelli al 10% 0,051
    (IC 0,044-0,059); 2×5 a 20 righe 0,042.

  I valori stanno in banda o al suo bordo. Non lo segnalo come problema.
- **Frazione di CD8 sotto i 10 pazienti** (`cd8_few.txt`, 200 repliche, stesso generatore e stessa regola
  di conteggio del test): copertura 0,960/0,960 con 5 pazienti e 0,958/0,963 con 6, in banda. Il solo
  difetto è formale (BASSA 7).
- **Divario del Modulo A senza leakage** (`gap_null.txt`): con 5, 6 e 10 pazienti nessuna replica oltre
  +0,05; una sola sotto −0,05 su 40 con 10 pazienti, oggi gialla.

---

## ALTA

### ALTA 1 — Confronto «stimabile» con un fattore tecnico fortemente associato, o con l'associazione non valutabile: la sintesi è gialla, ma la tabella dei confronti del report HTML e dell'app lo mostra ancora in verde

- **file:riga**:
  - `core/report.py:38-40`: `COMPARISON_STATE = {"stimabile": "verde", ...}`, con il commento «lo stesso
    della riga corrispondente nella sintesi»;
  - `core/report.py:133-139`: tag della tabella «Confronti richiesti»;
  - `app.py:82-83`: `COLOR_OF = {"stimabile": "verde", ...}`;
  - `app.py:269`: badge della tabella dei confronti.
- **gravità**: ALTA.
- **criterio violato**:
  - 2, perché il verde compare per un controllo in cui una regola segnala un limite;
  - 7, perché il commento a `report.py:38` contraddice il codice.

  È la correzione ALTA 2 del giro 1 applicata solo a `core/verdict.py`. Il resoconto del giro 1 diceva in
  modo esplicito «Lo stesso verde compare nella tabella dei confronti del report HTML e dell'app». La
  correzione A11 chiedeva «Allinea i colori della tabella dei confronti a quelli della sintesi».
- **motivo**: in `design_summary` (`core/verdict.py:141-142`) un «stimabile» con `partial_confounders` o
  `unassessed_technical` diventa giallo. Le due tabelle colorano invece per `c.classification` e restano
  verdi.
  - Prova costruita (`partial_conf_report.py`, 12 pazienti appaiati, batch associato al tessuto, V = 0,65):
    ```
    CLASSE: stimabile | partial: (('batch', 0.6474503217945908),) | unassessed: ()
    SINTESI: giallo | Confronto tissue: Tumor vs Normal | stimabile con 12 unità indipendenti; 'batch' fortemente associato al confronto (Cramér V = 0.65)
    REPORT HTML, tabella dei confronti: [('tag-verde', 'stimabile')]
    ```
  - Stesso caso su dati reali, in un file di risultati versionato (`validation_esterna/results/design_GSE125449.html`,
    rigenerato in df4f539):
    ```
    SINTESI giallo | Confronto cancer type: Hepatocellular carcinoma vs Intrahepatic cholangiocarcinoma | stimabile con 9 unità indipendenti; associazione con 'platform' non valutabile (tabella troppo picco…
    TABELLA ('cancer type: Hepatocellular carcinoma vs Intrahepatic cholangiocarcinoma', 'tag-verde', 'stimabile')
    ```
  - Il nuovo test (`tests/test_verdict.py:215-222`) controlla solo `design_summary`, non le due tabelle,
    quindi il problema passa la suite.
- **correzione proposta**: derivare il colore della tabella dei confronti dal `Check` della sintesi, invece
  che dalla classe:
  - nel report, `design_summary(result)` e il `Check` con nome `Confronto {factor}: {a} vs {b}`;
  - nell'app, lo stesso mapping;
  - in alternativa, una funzione unica `comparison_state(c)` in `core/verdict.py`, usata da sintesi, report
    e app.

  Va poi corretto il commento a `report.py:38`. Al test di `test_verdict.py:216` vanno aggiunte due
  asserzioni: il tag della tabella in `render_report` non è `tag-verde`, e il colore dell'app (stesso helper)
  è giallo. Infine vanno rigenerati `validation_esterna/results/design_GSE125449.html` e il relativo
  REPORT.

---

## MEDIA

### MEDIA 1 — README:94-95, «Il 72% della differenza è dovuto alle etichette NK»: la correzione MEDIA 6 del giro 1 è applicata solo nella sezione dei limiti

- **file:riga**: `README.md:94-95`; `validation/gse278694_b2_nk.py:1-2`.
- **gravità**: MEDIA.
- **criterio violato**: 4. È un'affermazione causale che dipende dall'ordine della scomposizione.
- **motivo**: il REGISTRO (riga «MEDIA 6: quota NK descritta con entrambi gli ordini») e README:239-248
  correggono il testo: «La scomposizione non è additiva e la quota dipende dall'ordine: descrive il peso
  delle etichette NK, non attribuisce una causa». Nell'estratto in apertura del README, che è la parte
  più letta, è rimasta la frase vecchia: `grep -n "72%" README.md` →
  `94:definizione. Il 72% della differenza è dovuto alle etichette NK contate come errore`. La stessa
  formulazione causale è nella docstring di `validation/gse278694_b2_nk.py:1-2` («quanto della
  differenza ... è dovuto alle etichette NK»), che era nella lista file:riga di MEDIA 6 del giro 1.
- **correzione proposta**: in README:94 usare la stessa formula della sezione dei limiti. Per esempio:
  «Togliendo le NK dalla convenzione per clone il tasso scende da 0,170 a 0,115 (72% della differenza);
  la quota dipende dall'ordine della scomposizione (vedi Ipotesi e limiti)». Allineare la docstring dello
  script.

### MEDIA 2 — Etichette indicate dall'utente che non esistono nei dati: il tasso d'errore cambia in silenzio e resta verde

- **file:riga**:
  - `core/tcr_validation.py:445`: filtro di `_judged_cells` sulle chiavi della mappa;
  - `:581`: filtro di `clone_level_units` su `marker_labels + error_labels`;
  - `:822-823` e `:351-361`: `validate_marker_map` controlla solo la forma;
  - `cli.py:207-208`: `--clone-error-labels` e `--clone-marker-priority` passati senza controllo.
- **gravità**: MEDIA.
- **criterio violato**: 3. È un input malformato che produce un risultato diverso invece di un errore
  chiaro.
- **motivo**: dopo MEDIA 1 del giro 1 lo strumento si ferma se un gene marcatore non esiste. Non si ferma se
  un'**etichetta** della mappa o di `--clone-error-labels` non compare nella colonna del tipo cellulare. È
  un errore facile: la guida stessa avverte, a `docs/GUIDA_UTENTE.md:207-208`, «`CD4T` e non `CD4+ T
  cells`».
  - `keys_mismatch.py`: stessi dati sintetici, con `CD8T` rinominata `CD8+ T` in `obs` e mappa invariata.
    ```
    --- etichette coerenti
       verde | Tasso d'errore [cell] in Tumor | 0.283 [0.246, 0.324], 12 pazienti
       verde | Tasso d'errore [clone] in Tumor | 0.181 [0.119, 0.266], 12 pazienti
    --- CD8T rinominata 'CD8+ T' in obs
       verde | Tasso d'errore [cell] in Tumor | 0.258 [0.187, 0.344], 12 pazienti
       verde | Tasso d'errore [clone] in Tumor | 0.577 [0.505, 0.646], 12 pazienti
    ```
    Nessun errore e nessun avviso. Il tasso per clone triplica, perché tutte le cellule CD8 diventano
    «fuori mappa» e quindi errori. Il tasso per cellula giudica le sole CD4. Entrambi restano verdi.
  - `labels_typo.py`, con un refuso in `--clone-error-labels`:
    ```
    {'clone_error_labels': ['NK']} -> 0.181 su 149 unità
    {'clone_error_labels': ['Nk']} -> 0.132 su 144 unità
    {'clone_marker_priority': ['CD8', 'CD4T']} -> KeyError 'pos_CD8'
    ```
    Le cellule NK escono dal calcolo in silenzio. Su GSE278694 è il passaggio da 0,170 a 0,115. Il refuso
    in `--clone-marker-priority` arriva invece all'utente come `errore: pos_CD8` (`cli.py:375-377`), un
    messaggio che non dice qual è il problema.
- **correzione proposta**: in `run_tcr_validation`, dare un errore esplicito con l'elenco delle etichette
  presenti in tre casi:
  - una chiave di `marker_map` non compare in `celltype_col`;
  - un'etichetta di `clone_error_labels` non compare in `celltype_col`;
  - un'etichetta di `clone_marker_priority` non è una chiave della mappa.

  Aggiungere un test per ciascun caso in `tests/test_input_errors.py`.

---

## BASSA

1. **`core/verdict.py:266-267` (con `:11-12` e README:285): un rosso «stima non possibile con questi dati»
   accanto a una stima.** Sotto 5 pazienti la riga del tasso d'errore è rossa con l'etichetta «stima non
   possibile con questi dati», ma il testo stampa la stima, per esempio «stima 0.060, intervallo non
   prodotto» per Adjacent_normal su GSE278694. La regola ora è dichiarata (BASSA 4 del giro 1), ma
   l'etichetta contraddice la riga. Correzione: un testo del rosso che dica «intervallo non possibile con
   questi dati», oppure grigio per i rifiuti, come diceva prima la definizione.
2. **`core/report.py:302-303` e `:359-360`: tag «numerosità insufficiente» anche quando l'intervallo
   dell'eccesso non è prodotto per variabilità nulla con abbastanza pazienti.** `disc_zero_var.py`, 8
   pazienti:
   ```
   SINTESI: ('giallo', 'stima +0.000, intervallo non prodotto: nessuna variabilità fra pazienti: ...')
   TAG REPORT: ['numerosità insufficiente']
   ```
   La sintesi e la narrativa sono state aggiornate nel giro 1, il tag del report e la tabella per coppia
   di compartimenti no. Correzione: usare `d.reason` (e `r.reason`) come in `_discordance_narrative`.
3. **`app.py:446-447`: BASSA 14 del giro 1 applicata in parte.** Il riferimento è predefinito solo se un
   compartimento si chiama come il sangue, ma altrimenti resta `index=0`, cioè il primo in ordine
   alfabetico. Con compartimenti `Normal`, `PB`, `Tumor` il riferimento predefinito è `Normal` (criterio
   3, scelta d'ordine implicita). Correzione: `index=None` con selezione obbligatoria quando nessun nome
   corrisponde.
4. **`app.py:197-198` e `:210-212`: la matrice dei fattori colora in rosso (`#b3261e`, lo stesso del
   badge rosso) ogni coppia strutturale, anche un annidamento.** Nella sintesi lo stesso annidamento è
   giallo (`STRUCTURAL_STATE`): due colori diversi per lo stesso fatto nella stessa pagina. Correzione:
   colorare le coppie strutturali con `STRUCTURAL_STATE` (rosso solo per uno-a-uno).
5. **`docs/METODI.md:122-124`: la sezione sul divario non riporta le regole nuove (ALTA 4 e MEDIA 3 del
   giro 1).** Il testo dice «sopra 0,05 la sintesi lo segnala in rosso» e tace sul giallo per un divario
   sotto −0,05 e per meno di 5 pazienti. README:294 e `RULES` sono invece aggiornati. Correzione: una
   frase con le due regole gialle.
6. **`core/leakage_audit.py:53-54`: commento non aggiornato.** Dice che oltre la soglia «il report segnala
   "leakage rilevabile"», ma nel codice quella stringa non esiste più: il tag è «divario oltre la soglia
   di ±0.05» (`report.py:200`) e la sintesi parla di «non stima la prestazione su pazienti nuovi».
   Correzione: aggiornare il commento.
7. **`tests/test_cd8_propagation_calibration.py:43`: la frazione di CD8 è calibrata solo con 10
   pazienti**, mentre lo strumento la produce da 5 (`core/cd8_propagation.py:223`, `:247`). Il difetto è
   solo formale, rispetto al criterio 1: la mia simulazione (`cd8_few.txt`) la trova in banda a 5 e 6
   pazienti (0,958-0,963). Correzione: parametrizzare il test con 5, 8 e 10 pazienti.
8. **`core/report.py:448-449` e `:505`: «Strumento diagnostico» nell'intestazione di ogni report**, una
   riga sotto «Non è un dispositivo medico». Per un istituto clinico «diagnostico» suggerisce proprio un
   uso clinico, e la frase è visibile anche nell'immagine del README. Correzione: «Strumento di controllo
   per la ricerca: non è una certificazione, è un supporto per chi analizza i dati».

Nessun altro problema ALTA o MEDIA trovato nelle parti rilette: `core/`, `cli.py`, `app.py`, i test
modificati, gli script di validazione modificati, README, GUIDA, METODI, ROADMAP, `pyproject.toml` e
REGISTRO. ROADMAP e README non dichiarano come fatte voci che sono rimandate alla 1.0.
