# Resoconto del revisore indipendente — giro 1

Repo `genomic-audit-mvp`, ramo `revisione`, commit `70279da` (albero pulito prima e dopo la revisione).
Nessun file della repo né di `../pdac-ml` è stato modificato. Script e output delle prove sono in
`/tmp/claude-1000/-home-nemo-Uni-PROGETTO-ML-genomic-audit-mvp/22405c93-e5d2-4170-9aed-2c2f9b779264/scratchpad/revisore1/`
(`sim_v.py`, `partial_conf.py`, `neg_gap.py`, `neg_gap2.py`, `few_pat.py`, `mk.py`, `mk2.py`, `reftypo.py`,
`sim_lowrate.py`, `sim_paired.py`, `gap_null.py`, `demo_stdout.txt`).

## Esiti

- `uv run pytest -n auto -m "not slow" -q`: **119 passed, 57 warnings in 124,48 s** (nessun test saltato:
  i dati di GSE278694 sono presenti accanto alla repo, quindi i 3 test sui dati reali sono girati e passano).
- `uv run ruff check .`: **All checks passed!**
- Raccolta: 188 test in tutto, 119 «not slow» e 69 «slow» (i lenti non sono stati eseguiti, come da mandato).
- `uv run audit-sc demo --out <scratch>/demo_report.html --data-dir <scratch>/synthetic`: exit 0, 1 min 34 s.

**Conteggi: ALTA 4 · MEDIA 8 · BASSA 14.**

Verifiche passate (non sono problemi, le riporto per completezza):
- Correzione 5 (pareggi): il codice esclude i pareggi (`core/tcr_validation.py:515-527`) e
  `tests/test_real_data_gse278694.py:61` verifica, sui dati reali e con il codice attuale, 0,16967509025270758
  su 277 unità e 10 pazienti, uguale a `../pdac-ml/results/tables/error_rate_blood.csv` (Tumor, 277, 10,
  0.16967509025270758). Sensibilità 0,201/0,163 su 288 unità in
  `validation/results/cli_tcr_real_conventions_tesi_pareggi_{error,correct}.txt`. Il vecchio 0,1951 è
  indicato come superato in README, METODI, REPORT esterno, REPORT_INTERVENTI e nello script di
  validazione (eccezione: un file di risultati, BASSA 2).
- Numeri del README ritrovati nei file di risultati: 37 librerie/22 pazienti
  (`validation/results/gse278694_real_sheet.csv`, `cli_design_real.txt`); 0.778±0.060, 0.895±0.004, +0.117,
  9.946 cellule, 19 pazienti, 8 classi (`validation_esterna/results/A1_cli_stdout.txt`, `A1_cli_stderr.txt`);
  0.093 [0.043, 0.188], 0.170 [0.100, 0.274], 277 unità, 11 pareggi (`cli_tcr_real_conventions_tesi.txt`);
  72% (`b2_quota_nk.json`); [0,095; 0,252], 0,199/0,164 su 287, 0,154 su 227, 0,152 su 231 (`a1_moduleB.json`,
  `fedelta_S5`); 0,45 e 7 pazienti su 13 (`a3_summary.json`, `cli_tcr_real.txt`); 3,6%-6,3%, 14,5%
  (`revisione/log/a1_calibrazione.txt`); 0,5-1,8% e 8,3% (`a2_calibrazione.txt`); 0,46
  (`design_calibrazione.txt`); 5 h 26 min (`validation_esterna/results/leakage_GSE125449_completo_interrotto.txt`);
  20 permutazioni (`moduleA_perm.json`); Streamlit 1.51 (`revisione/log/streamlit_versione_minima.txt`).
  Eccezioni: MEDIA 4 e MEDIA 5.
- Calibrazione dell'intervallo t cluster-robusto anche fuori dagli scenari dei test (mie simulazioni,
  1000 repliche, `sim_lowrate.py`): tassi bassi per cellula (0,02-0,04) 4,0%-6,7%, per clone (0,03-0,05)
  2,4%-3,8% di mancati, in banda. Nessun problema lì.

---

## ALTA

### ALTA 1 — Cramér V esposto (allarmi e verde) sotto le 20 righe, dove non è calibrato
- **file:riga**: `core/design_audit.py:46` (`MIN_ROWS_FOR_V = 10`), `:281`; `core/verdict.py:149-152`;
  README:146-148; `docs/METODI.md:160-166`.
- **gravità**: ALTA.
- **criterio violato**: 1 (calibrazione nel regime d'uso); 2 (verde fuori dal regime calibrato).
- **motivo**: il V è «valutabile» da 10 righe (con almeno 2 righe per cella in media), ma la
  calibrazione dei falsi allarmi parte da 20 righe (`tests/test_design_audit_calibration.py:90-109`: 20, 30,
  40 righe, tabelle 2×2 e 3×3). Il regime 10-19 righe è proprio quello delle coorti pilota (5-9 pazienti × 2
  tessuti) e dei sotto-fogli usati da `_partial_confounders`. Simulazione con lo stesso generatore del test
  (fattori indipendenti, 400 repliche, `sim_v.py`), frazione di coppie valutabili con allarme, banda
  dichiarata (0; 0,05):
  ```
  10 righe, 2x2: allarmi 38/396 = 0.096
  12 righe, 2x2: allarmi 23/399 = 0.058
  12 righe, 2x3: allarmi 28/400 = 0.070
  15 righe, 2x2: 0.033   16 righe, 2x3: 0.028   20 righe, 2x5: 0.043
  ```
  Fuori banda a 10-12 righe. Anche il verde «Associazioni fra fattori» (`verdict.py:149-152`) compare in
  quel regime, dove la potenza non è mai stata misurata (l'unica misura è 3×3 con 60 righe).
- **correzione proposta**: portare `MIN_ROWS_FOR_V` a 20 (la soglia più bassa calibrata) oppure calibrare
  10-19 righe e tabelle non quadrate con banda dichiarata prima, rifiutando («non valutabile», grigio) dove
  esce dalla banda; aggiornare README, METODI e il test.

### ALTA 2 — Confronto in verde anche quando la regola del confondimento parziale lo segnala
- **file:riga**: `core/design_audit.py:486-501` (confondente parziale solo aggiunto ai `reasons`, la classe
  resta «stimabile»); `core/verdict.py:126-127`; `core/report.py:39`; `app.py:76`.
- **gravità**: ALTA.
- **criterio violato**: 2 (verde dove una regola segnala un limite); 4 (l'etichetta verde «nessun limite
  rilevato dalle regole» è falsa in quel caso).
- **motivo**: `partial_conf.py`, 7 pazienti appaiati Tumor/Normal, batch fortemente associato al tessuto ma
  senza separarlo del tutto, 3 tessuti in tutto (16 righe):
  ```
  CLASSE: stimabile
  FRASE: ... stimabile, disegno appaiato (entro paziente), 7 pazienti con entrambi i livelli. ...
         Attenzione: 'batch' è fortemente associato al confronto (Cramér V = 0.78). ...
  PAIR tissue batch None False
  SINTESI: verde | Confronto tissue: Tumor vs Normal | stimabile con 7 unità indipendenti
  SINTESI: grigio | Associazioni fra fattori | Cramér V non valutabile per 3 coppie
  ```
  La regola del modulo (V ≥ 0,5 fra fattore tecnico e confronto) trova il limite, ma la sintesi mostra
  verde e il testo della riga non lo riporta; nessun'altra riga è gialla, perché sulla tabella intera il V
  non è valutabile. Lo stesso verde compare nella tabella dei confronti del report HTML e dell'app.
- **correzione proposta**: con un confondente parziale (o, in generale, con `reasons` diversi dal solo
  testo di potenza) il confronto va in giallo, con il motivo nel testo della riga; test con un
  `DesignAuditResult` vero come quello sopra.

### ALTA 3 — Il README mostra come «report della demo» un'immagine del report vecchio, con i verdi e le diciture eliminati
- **file:riga**: README:113-115; `docs/img/report_demo.png` (ultimo commit `818ba98`, 2026-10-04,
  prima della revisione).
- **gravità**: ALTA.
- **criterio violato**: 4 («affidabile» non letteralmente vero); 2 (verdi mostrati dove le regole attuali non
  li danno). Correzioni A6 e A7 non applicate al README.
- **motivo**: l'immagine (letta con il visualizzatore) mostra la «Sintesi dei controlli» con lo stato
  **«stima affidabile»** per «Associazioni fra fattori», «Valutazione con split per paziente», «Confronto fra
  modelli» (in verde, con «macro-F1 più alta: logreg») e per i tassi d'errore, una riga «Frazione di CD8 —
  Intervalli della frazione di CD8 in Tumor — stima con limiti», nessuna riga per i fatti strutturali, e
  l'intestazione «è un supporto alla decisione per chi analizza i dati». Il codice attuale produce invece «nessun limite rilevato dalle
  regole», confronto fra modelli giallo, nessuna riga CD8 (verificato con la demo rieseguita). Il README la
  presenta come «Il report HTML della demo».
- **correzione proposta**: rigenerare l'immagine dal report della demo prodotto dal codice attuale (o
  toglierla) e aggiungere un controllo (anche manuale, nel registro) che la rigeneri a ogni cambio di regole.

### ALTA 4 — Divario negativo oltre la soglia: verde con il testo falso «entro la soglia»
- **file:riga**: `core/verdict.py:174-181`; `core/report.py:200-201`; README:284; `docs/METODI.md:121-122`.
- **gravità**: ALTA.
- **criterio violato**: 4 (affermazione falsa nella sintesi e nel report); 2.
- **motivo**: la regola usa `r.gap > GAP_ALERT`; un divario di −0,108 finisce nel ramo verde con il testo
  «entro la soglia di 0.05», mentre la narrativa dello stesso risultato dice «oltre la soglia di rumore
  dichiarata». Caso naturale, non costruito (`neg_gap2.py`, dati sintetici senza segnale, 6 pazienti, seme 0):
  ```
  seed 0 gap -0.108 grouped 0.552 random 0.444 | verde divario rispetto allo split per paziente -0.108, entro la soglia di 0.05
  ```
  e nel report HTML il tag «divario entro la soglia di 0.05» (`neg_gap.py`). Un divario negativo oltre la
  soglia è un segnale anomalo (per esempio classi assenti dai fold per paziente), non «nessun limite».
  `tests/test_verdict.py:185-187` prova solo divari positivi.
- **correzione proposta**: usare `abs(gap) > GAP_ALERT` per uscire dal verde; per il divario negativo uno
  stato giallo con un testo che dica la direzione e che va letto con le classi assenti; test con divario
  −0,1; allineare README:284 e METODI.

---

## MEDIA

### MEDIA 1 — Il controllo del tasso d'errore sparisce dalla sintesi con un compartimento di riferimento sbagliato o con marcatori assenti
- **file:riga**: `core/tcr_validation.py:361-365` (geni mancanti: solo `warnings.warn`), `:387` (riferimento
  non presente: nessun controllo), `core/verdict.py:221-232` (nessuna riga se `by_compartment` è vuoto);
  `cli.py:203-211`.
- **gravità**: MEDIA.
- **criterio violato**: 3 (input malformato che produce in silenzio un risultato diverso); 2 (un controllo
  non eseguito deve comparire in grigio).
- **motivo**: `reftypo.py` con `reference_compartment="Blood"` (i dati hanno `PBMC`) e `mk2.py` con
  marcatori indicati per ID Ensembl su `var_names` simbolici: nessun errore, narrativa «nessuna unità
  valutabile», e la sintesi contiene solo
  ```
  SINTESI: giallo | Eccesso di discordanza fra compartimenti
  ```
  cioè la riga del tasso d'errore non c'è, né grigia né rossa; copertura dei flag 0,0. Nell'audit del
  disegno un refuso equivalente (livello assente) è invece un errore d'input in grigio.
- **correzione proposta**: errore esplicito se `reference_compartment` non è fra i compartimenti presenti e
  se per un'etichetta della mappa non si trova nessun gene; in `tcr_summary`, una riga grigia «nessuna unità
  valutabile» per ogni convenzione senza compartimenti.

### MEDIA 2 — App: la mappa dei marcatori non è validata (traceback grezzo o risultato vuoto in silenzio)
- **file:riga**: `app.py:444` (`json.loads(markers)` senza controllo della forma), `app.py:463`
  (`except ValueError` soltanto).
- **gravità**: MEDIA.
- **criterio violato**: 3. Correzione A8 («JSON e tabelle illeggibili danno un traceback grezzo») applicata
  solo in parte: la CLI valida (`cli.py:181-191`), l'app no.
- **motivo**: `mk.py`: con il JSON `["CD4", "CD8A"]` → `AttributeError 'list' object has no attribute
  'items'`, non intercettato dall'app (traceback di Streamlit); con `{"CD4T": "CD4", "CD8T": "CD8A"}`
  (stringa invece di lista, refuso facile nella casella di testo) nessun errore: i geni cercati diventano
  `['C', 'D', '4']`, solo avvisi in console, e il tasso d'errore sparisce dalla sintesi (MEDIA 1).
- **correzione proposta**: riusare nell'app la validazione di `_read_marker_map` (spostandola in `core/`),
  con messaggio `st.error`.

### MEDIA 3 — Modulo A in verde con 2 pazienti
- **file:riga**: `core/verdict.py:171-181`; `core/leakage_audit.py:517`.
- **gravità**: MEDIA.
- **criterio violato**: 2 (verde fuori dal regime in cui la regola è stata provata).
- **motivo**: `few_pat.py`, 2 pazienti, 2 fold:
  ```
  verde | Valutazione con split per paziente | macro-F1 = 0.534 ± 0.026 su 2 fold
  verde | Valutazione con split casuale delle cellule | divario ... -0.004, entro la soglia di 0.05
  ```
  Lo strumento rifiuta ogni intervallo sotto 5 pazienti e l'audit del disegno chiama «bassa potenza» meno di
  5 unità; la soglia ±0,05 è stata provata solo su GSE125449 con 19 pazienti (A3, 20 permutazioni). Con 2-4
  pazienti la stima per paziente poggia su 2-4 fold, ma nessuna regola lo segnala.
- **correzione proposta**: giallo (o grigio) per entrambe le righe del Modulo A sotto 5 pazienti, con il
  motivo nel testo; dichiararlo nella tabella del README.

### MEDIA 4 — README:243 cita come intervallo della tesi un valore superato e non tracciabile
- **file:riga**: README:242-245.
- **gravità**: MEDIA.
- **criterio violato**: 5.
- **motivo**: «Gli intervalli della tesi con 4 pazienti (per esempio la differenza Tumor − Adjacent_normal,
  [+0,021; +0,170]) sono riportati negli output di verifica». Il valore viene dalla vecchia regola dei
  pareggi (+0,0795 su 88 cloni, `docs/sviluppo/REPORT_INTERVENTI.md:348`, documento «superato»). Con la
  regola attuale la tesi dà +0,059 [0,000; 0,160] su 85 cloni (`../pdac-ml/results/tables/error_rate_tie_sensitivity.csv`,
  riga `blood,ties_excluded,diff_Tumor_minus_Adjacent_normal`) e gli output di verifica riportano
  [0,0; 0,141] (`validation/results/a1_moduleB.json`, S5, `Tumor-Adjacent`). `grep` di `0,170]`/`+0,021;`
  nei file di risultati: nessuna occorrenza.
- **correzione proposta**: sostituire con l'esempio attuale (+0,059, 85 cloni, 4 pazienti, con la fonte)
  o indicare il vecchio valore come superato.

### MEDIA 5 — Esiti dei test del README non misurati sul commit consegnato
- **file:riga**: README:356-367.
- **gravità**: MEDIA.
- **criterio violato**: 6. Correzione B7 («Esiti dei test misurati») incompleta.
- **motivo**: il README stesso dice «commit 69cf239, prima dell'allineamento dei pareggi; da ripetere sul
  commit finale». I conteggi sono già superati: oggi la raccolta dà 188 test (README: 177+3+7 = 187) e 119
  «not slow» (README: 115+3 = 118); il commit `70279da` ha cambiato `core/tcr_validation.py` e
  `tests/test_input_errors.py`. Anche `validation_esterna/REPORT.md:124` porta la stessa nota.
- **correzione proposta**: rieseguire `revisione/prova_clone_pulito.sh` e la suite completa sul commit
  finale, salvare i log in `docs/sviluppo/logs/` e aggiornare la tabella.

### MEDIA 6 — «Il 72% della differenza è dovuto alle etichette NK»: la quota dipende dall'ordine della scomposizione
- **file:riga**: README:94-95 e :234-236; `validation/gse278694_b2_nk.py:4-5`.
- **gravità**: MEDIA.
- **criterio violato**: 4 (affermazione causale più forte di quanto i dati sostengano).
- **motivo**: 72% = (per clone con NK errore − per clone senza fuori mappa)/(per clone con NK − per
  cellula), cioè l'effetto NK misurato **dopo** il passaggio all'unità clone. Partendo dalla convenzione per
  cellula, la scomposizione dello stesso progetto (`validation/results/a1_moduleB.json`, S2 → S3: 0,0927 →
  0,1632, cioè +0,0705 su una differenza totale 0,1697 − 0,0930 = 0,0767) attribuisce alle NK il **92%**. La
  scomposizione non è additiva: «dovuto al 72%» non è una proprietà dei dati.
- **correzione proposta**: riportare entrambe le quote (72% e 92%, secondo l'ordine) o descrivere il
  calcolo («passando, nella convenzione per clone, dalle NK escluse alle NK contate come errore il tasso
  sale da 0,115 a 0,170»).

### MEDIA 7 — CD4 citotossiche «causa probabile di gran parte» delle chiamate CD8: più forte della fonte
- **file:riga**: README:222-225.
- **gravità**: MEDIA.
- **criterio violato**: 4.
- **motivo**: la fonte citata (README di `../pdac-ml`, riga 35) scrive che la maggiore espressione di GNLY
  «suggerisce che **una parte** siano CD4 citotossiche [...]: è una **possibile** causa dell'errore, non
  un'alternativa all'errore», e che le cellule discordanti «sono quindi della linea CD4». Il README dello
  strumento la trasforma in «causa probabile di gran parte» e la elenca fra le eccezioni che il modulo «non
  distingue da un errore di annotazione», mentre per la tesi una CD4 citotossica etichettata CD8T è proprio
  un errore di linea (che il modulo rileva correttamente). La docstring (`core/tcr_validation.py:12-13`,
  «probabile causa reale di parte») è più vicina alla fonte di quanto lo sia il README. Nota: la formulazione nasce dalla
  correzione B1 («probabile causa reale del 45%»), che già andava oltre la fonte.
- **correzione proposta**: «possibile causa di una parte delle chiamate CD8 sulle vere CD4»; spiegare che
  una CD4 citotossica etichettata CD8 è un errore di linea, e che l'eccezione all'ipotesi riguarda le
  etichette funzionali («citotossica») usate come se fossero di linea.

### MEDIA 8 — Docstring di `stats.py` con numeri attribuiti al test che il test non produce
- **file:riga**: `core/stats.py:8-10` e `:254-256`.
- **gravità**: MEDIA.
- **criterio violato**: 7.
- **motivo**: «sotto H0, 16-18% di falsi positivi con 5 gruppi e 9-11% con 10 [...]
  (tests/test_stats_calibration.py)». Il test (`test_percentile_bootstrap_is_anticonservative_with_five_groups`)
  misura 14,5% a 5 gruppi (`revisione/log/a1_calibrazione.txt`: `58/400 = 0.145`) e non misura i 10 gruppi;
  «banda rispettata solo da circa 25 gruppi»: c'è un test a 25 gruppi, ma nessuno fra 10 e 25, quindi «solo da»
  non è misurato. METODI:52-56 dà la versione corretta (14,5% dal test; 18% e 9% da
  simulazioni esplorative nel REGISTRO).
- **correzione proposta**: allineare le due docstring a METODI (14,5% misurato dal test; gli altri numeri
  come simulazioni esplorative, con la fonte).

---

## BASSA

1. **README:53 e :65 — citazioni incomplete negli estratti** (correzione B7). «Chen M-M et al., *Cancer
   Cell* 2025» e «Ma L. et al., *Cancer Cell* 2019» senza volume, pagine e doi, mentre B7 chiedeva la
   citazione completa «sia nell'estratto sia nell'elenco». Correzione: riportare la citazione completa
   anche negli estratti (o rimandare esplicitamente all'elenco).
2. **`validation/results/cli_tcr_real_vdj_originali.txt:2` — risultato del codice vecchio non annotato**
   (correzione 5 e B1). File del 2026-10-01 (`07b1311`) con «'Tumor': 0.195 (IC95% [0.123, 0.274]» senza
   «superato» e con «un effetto reale (l'intervallo esclude lo zero)», frase tolta da B1. Correzione:
   rigenerarlo con il codice attuale o spostarlo fra i documenti storici con una nota «superato».
3. **`tests/fixtures/tcr_regression_baseline.json:83` e `tests/test_tcr_flags.py:54-69` — regressione
   cieca sulla convenzione per clone** (criterio 9). La narrativa salvata dice ancora «pareggi: probabilità
   di errore con scelta a caso fra le etichette a pari merito», e il test confronta solo discordanza e tasso
   per cellula: sul caso `eccesso_035_marcatori` il tasso per clone è passato da 0,246 (fixture) a 0,181 con
   31 pareggi esclusi su 180 unità senza che nessun test lo segnali (verificato con
   `tests.fixtures.tcr_regression.run_config`). Correzione: aggiungere alla fixture i numeri della
   convenzione per clone (con `n_ties`), rigenerarla documentando il motivo, aggiornare la docstring di
   `tests/fixtures/tcr_regression.py:3-9`.
4. **`core/verdict.py:201-202` e `:231-232` contro `:13-14` — stati incoerenti con le definizioni.** Le
   definizioni dicono grigio per un controllo «rifiutato dallo strumento», ma l'intervallo rifiutato sotto 5
   pazienti è rosso; un tasso pari a 0 con 8 pazienti (fixture `tre_compartimenti`, Adjacent 0,000) diventa
   rosso «stima non possibile con questi dati», anche se la stima esiste. Correzione: decidere una regola e
   scriverla in modo coerente (per esempio grigio per i rifiuti; per il tasso nullo un testo che dica che la
   stima è 0 e l'intervallo non è definito su scala logit).
5. **`core/tcr_validation.py:480-501` — differenze appaiate senza un test di calibrazione dedicato**
   (criterio 1, solo formale). Gli intervalli delle differenze Tumor − Adjacent (valori in {−1, 0, 1} per
   clone) sono esposti ma nessuno scenario del test li simula. Mia simulazione (`sim_paired.py`, 1000
   repliche, 5-15 pazienti): 2,0%-4,6% di mancati, in banda. Correzione: aggiungere lo scenario al test con
   banda dichiarata.
6. **`docs/METODI.md:74-75` contro `core/tcr_validation.py:30-31`** — lo stesso esempio (3 cellule nel
   sangue, 60 nel tumore, 40% di errore) vale «circa −0,13» in METODI e «circa −0,06» nella docstring; la
   differenza dipende da quante etichette sbagliate ci sono (2 nel test citato, 1 nella docstring), che
   nessuno dei due testi dice. Correzione: dichiarare il numero di etichette in entrambi.
7. **`core/io.py:9-10` contro `:85`** — «nessuna riga viene scartata in silenzio», ma `meta.loc[bc.values]`
   scarta senza avviso le righe dei metadati senza barcode nella matrice. Correzione: dichiarare quante righe
   restano fuori (o correggere la docstring).
8. **Test più deboli del loro nome** (criterio 9). `tests/test_input_errors.py:229`
   (`..._assigns_flags_by_position_with_duplicate_obs_names`) verifica solo che la CLI non vada in errore e
   che il report non contenga la frazione di CD8, non che i flag stiano sulle righe giuste;
   `tests/test_input_errors.py:308` verifica l'assenza di scritture su disco cercando la stringa
   «tempfile» nel sorgente. Correzione: confrontare i flag usati da `cd8_fraction_intervals` con quelli
   per posizione; per l'app, un AppTest che controlli che la cartella temporanea resti vuota.
9. **`tests/test_stats_calibration.py:6-11` e `:17`** — l'intestazione elenca ancora «Nadeau-Bengio,
   cluster bootstrap» come stimatori e dà il regime «4-11 pazienti» (lo strumento rifiuta sotto 5 e la
   calibrazione è a 5-15). Correzione: aggiornare l'intestazione.
10. **README:231-232** — «Nel dataset GSE278694 l'errore è più frequente nel tumore che nel sangue»: con lo
    schema dello strumento il sangue è il riferimento e non viene mai giudicato; nessun file di risultati
    dello strumento contiene un tasso per il sangue (il confronto viene dagli schemi naive/LOCO della tesi).
    Correzione: citare la fonte (tesi, `error_rate_tie_sensitivity.csv`, schemi naive e loco) o dire «nel
    tumore che nel tessuto adiacente (4 pazienti, intervallo non prodotto)».
11. **`docs/GUIDA_UTENTE.md:96` contro `:264`** — lo snippet per accorpare le classi rare fa
    `a.write_h5ad("studio.h5ad")`, cioè sovrascrive il file di partenza, mentre la sezione 8 dice «I dati di
    partenza non vengono mai modificati». Correzione: scrivere `studio_classi.h5ad`.
12. **`core/leakage_audit.py:5-7`** — «Il divario [...] misura quanto della performance [...] sarebbe un
    artefatto»: più forte del README, che giustamente dice che il divario include tutto ciò che distingue i
    pazienti (su GSE125449 anche la piattaforma). Correzione: allineare la docstring.
13. **`app.py:494`** — titolo «Cellule con TCR per etichetta (contaminazioni cross-lineage)», mentre la
    didascalia sotto dice che la tabella non distingue doppietti, RNA ambientale ed errori. Correzione:
    togliere «contaminazioni».
14. **`app.py:440`** — il compartimento di riferimento ha come valore predefinito il primo in ordine
    alfabetico (su GSE278694 `Adjacent_normal`, non `PBMC`): un clic su «Esegui» senza cambiare il selettore
    stima le identità dal tessuto adiacente. Correzione: nessun valore predefinito (selezione obbligatoria)
    o predefinito solo se un compartimento si chiama come il sangue.

Nota non contata: `validation_esterna/REPORT.md:220` riporta che sulla coorte ridotta di GSE125449 (10
pazienti) la verifica indipendente A2 **fallisce** (LeaveOneGroupOut
0,7530 contro 0,7127); il README (riga 123) riporta solo l'esito sul dataset completo. Il REPORT è
trasparente e il README vi rimanda, quindi non lo conto come problema, ma una riga nel README sarebbe più
onesta.
