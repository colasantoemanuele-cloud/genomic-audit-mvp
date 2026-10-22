# Validazione tecnica esterna del Modulo B su GSE200996

Criteri: `validation_esterna_B/GSE200996/CRITERI.md`, scritti e committati prima di qualsiasi download
(commit `8b10bb0`) e non modificati dopo. Mappa delle etichette: `MAPPA_ETICHETTE.md` (commit
`2e67f4d`), prima di eseguire il Modulo B. Secondo tentativo della voce 1 della roadmap, dopo
GSE139555 (non eleggibile: `validation_esterna_B/REPORT.md`). Numeri e tabelle di questo documento
sono assemblati da `make_report.py` a partire dai file in `results/`.

Versione dello strumento: 0.9.0, branch `validazione-b-gse200996`, creato da `validazione-b` perché
master non ne contiene ancora il merge. Ambiente: `.venv` del progetto (uv); workstation con 12 core
e 31 GB di RAM.

## 1. Dataset ed eleggibilità

[GSE200996](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE200996) (Luoma et al., *Cell* 2022,
PMC9508682): carcinoma squamoso del cavo orale, terapia neoadiuvante anti-PD-1 ± anti-CTLA-4
(NCT02919683); scRNA-seq 10x 5′ + scTCR-seq di tumore (prima e dopo la terapia) e di sangue
(prelievi B1-B3, corsie con più pazienti marcati con hashtag). Cell Ranger 4.0.0, GRCh38-2020-A.

**File e download.** Elenco dei file supplementari letto sul server FTP di GEO. Scaricati in
`../dati_esterni/GSE200996/`, fuori dal repository (spazio libero prima del download: 582 GB):
- prima, solo metadati per cellula (CD4, CD8, CD45, Myeloid × PBMC/tumor), le due chiavi degli
  hashtag, `filelist.txt` e le 5 series matrix (9,3 MB, 43 s con pause di 1 s fra i file; il nome
  `GSE200996_series_matrix.txt.gz` dà 404, la serie ha una matrix per piattaforma);
- poi, dopo il commit di `MAPPA_ETICHETTE.md`, `GSE200996_RAW.tar` (1.591.296.000 byte, uguale a
  `filelist.txt`; 53 s, al primo tentativo), da cui sono estratti solo i 61 file single-cell di tumore
  e sangue (corsie non selezionate) dei pazienti eleggibili (814 MB): niente TCR bulk, niente corsie
  selezionate.
Nessuna limitazione delle richieste (reCAPTCHA, HTTP 429). Testo dell'articolo da NCBI efetch (PMC
XML; Europe PMC rispondeva 500).

**Campioni GEO single-cell.** 25 matrici di tumore (6 prima e 19 dopo la terapia) e 24 file TCR
single-cell di tumore (5 prima, 19 dopo); 12 corsie di sangue «total PBMC» (matrice + TCR); 4
corsie di sangue selezionate CD38+ HLA-DR+ CD4 o CD8 (P13, P19, P26, P32; prelievi B1 e B2).

### E5 — origine del sangue

Metodi, «Preparation of total PBMC»: le cellule di ogni prelievo, marcate con un hashtag TotalSeq-C,
sono state selezionate al FACS come «Live CD45+ cells»: tutti i leucociti vivi, non una linea. Le
corsie «sorted» (GSM6048338-345) sono un esperimento separato su cellule T attivate CD38+ HLA-DR+
(Metodi, «Sorting of activated blood T cells»), e non sono usate. I metadati del sangue contengono
145.584 cellule, lo stesso numero del «blood dataset» dei Metodi; ogni coppia
paziente-prelievo dei metadati compare nella chiave delle corsie non selezionate (fuori chiave:
0). Tutti i
18 pazienti eleggibili hanno nel
sangue sia cellule CD4 sia CD8. **E5 soddisfatta.**

### E2 — origine delle etichette

Etichetta = appartenenza al file CD4 o CD8 degli autori (mappa in `MAPPA_ETICHETTE.md`, commit
`2e67f4d`, prima di eseguire il modulo). Metodi: trasferimento di etichette da un riferimento CITE-seq
di PBMC (Seurat `MapQuery`, Hao et al. 2021) con punteggio di predizione > 0,6, cioè **annotazione
computazionale sull'espressione**, non sorting; inoltre gli autori hanno tolto le cellule CD4 con
CD8A/CD8B rilevati e le CD8 con CD4 rilevato. Sangue: 33.446 CD4 e 14.207 CD8;
tumore: 24.032 CD4 e 18.126 CD8; nessuna cellula in entrambi i file. Cellule T del
file CD45 senza etichetta (escluse): 6.801 nel sangue e
6.309 nel tumore (tutti i pazienti). Cellule etichettate con tipo diverso
da T nel file CD45 (tenute, come deciso nella mappa): sangue {'T cell': 47060, 'Cycling': 458, 'NK cell': 135};
tumore {'T cell': 39527, 'Cycling': 2232, 'NK cell': 294, 'B cell': 52, 'Myeloid': 50, 'Mast': 3}. **E2 soddisfatta.**

### E1, E3, E4 — demultiplazione, abbinamento, formati

Demultiplazione degli autori (Metodi): GMM-demux sulle conte degli hashtag, singoletti con
confidenza ≥ 0,8; in media 83,46% singoletti, 15,53% doppietti o multipletti, 0,28% negativi. I
metadati contengono solo i singoletti assegnati; l'ID di cellula porta paziente e prelievo. Barcode
VDJ di ciascuna corsia non assegnati dagli autori (doppietti, negativi e cellule tolte dal controllo
di qualità, non separabili senza rifare la demultiplazione), esclusi:

| Corsia | Pazienti | Barcode VDJ (is_cell) | Assegnati a un paziente | Non assegnati, esclusi |
|---|---|---|---|---|
| GSM6048239 | P13-P14 | 6.755 | 5.100 | 1.655 |
| GSM6048240 | P15-P17 | 6.547 | 4.949 | 1.598 |
| GSM6048241 | P16-P24-P26 | 5.556 | 3.234 | 2.322 |
| GSM6048242 | P18-P29-P30 | 3.976 | 3.074 | 902 |
| GSM6048243 | P19-P20 | 5.212 | 3.806 | 1.406 |
| GSM6048244 | P21-P25-P28 | 7.709 | 5.693 | 2.016 |
| GSM6048245 | P22-P23 | 5.549 | 4.466 | 1.083 |
| GSM6048246 | P31-P32 | 3.812 | 2.957 | 855 |

Matrici: `raw_feature_bc_matrix.h5` di Cell Ranger, conteggi int32 interi non negativi in tutti i
30 campioni usati: sì. Feature diverse da
Gene Expression: Hashtag1_TotalSeqC, Hashtag2_TotalSeqC, Hashtag3_TotalSeqC, Hashtag4_TotalSeqC, Hashtag5_TotalSeqC, Hashtag6_TotalSeqC, Hashtag7_TotalSeqC (solo hashtag, nessuna proteina di superficie). File VDJ:
`filtered_contig_annotations.csv` di Cell Ranger, con le colonne attese dal Modulo B; i file del
sangue (uno per corsia) sono divisi per paziente copiando le righe senza modificarle
(`build_gse200996.py`). **E4 soddisfatta.**

Per campione (analisi principale e descrittiva): cellule CD45 assegnate dagli autori, cellule T
etichettate CD4/CD8, cellule T non etichettate (escluse dalla mappa), cellule etichettate con almeno
un contig TRB produttivo (is_cell, high_confidence, full_length, productive, CDR3 presente):

| Campione | Compartimento | Cellule CD45 | T etichettate | T escluse | Con TRB produttiva |
|---|---|---|---|---|---|
| P13_Post-Tx | Tumor | 1.871 | 442 | 141 | 316 (71,5%) |
| P14_Post-Tx | Tumor | 1.881 | 1.162 | 115 | 1.074 (92,4%) |
| P15_Post-Tx | Tumor | 3.480 | 2.354 | 549 | 1.994 (84,7%) |
| P16_Post-Tx | Tumor | 3.170 | 2.109 | 353 | 1.932 (91,6%) |
| P18_Post-Tx | Tumor | 3.386 | 2.069 | 290 | 1.738 (84,0%) |
| P19_Post-Tx | Tumor | 2.569 | 1.707 | 254 | 1.578 (92,4%) |
| P20_Post-Tx | Tumor | 4.069 | 2.690 | 371 | 2.397 (89,1%) |
| P21_Post-Tx | Tumor | 3.169 | 2.035 | 179 | 1.914 (94,1%) |
| P22_Post-Tx | Tumor | 3.049 | 2.023 | 243 | 1.902 (94,0%) |
| P23_Post-Tx | Tumor | 3.009 | 1.698 | 194 | 1.571 (92,5%) |
| P24_Post-Tx | Tumor | 1.369 | 476 | 88 | 447 (93,9%) |
| P25_Post-Tx | Tumor | 1.812 | 716 | 203 | 594 (83,0%) |
| P26_Post-Tx | Tumor | 1.888 | 854 | 148 | 756 (88,5%) |
| P28_Post-Tx | Tumor | 3.190 | 1.953 | 336 | 1.790 (91,7%) |
| P29_Post-Tx | Tumor | 4.108 | 2.896 | 201 | 2.731 (94,3%) |
| P30_Post-Tx | Tumor | 2.712 | 693 | 80 | 586 (84,6%) |
| P31_Post-Tx | Tumor | 3.987 | 1.935 | 371 | 1.675 (86,6%) |
| P32_Post-Tx | Tumor | 2.984 | 1.207 | 289 | 1.079 (89,4%) |
| P18_Pre-Tx | Tumor | 4.551 | 3.207 | 439 | 2.907 (90,6%) |
| P23_Pre-Tx | — | — | — | — | matrice senza TCR single-cell: escluso |
| P24_Pre-Tx | Tumor | 3.497 | 1.418 | 198 | 1.311 (92,5%) |
| P29_Pre-Tx | Tumor | 2.739 | 1.983 | 186 | 1.871 (94,4%) |
| P32_Pre-Tx | Tumor | 3.463 | 1.154 | 183 | 1.018 (88,2%) |
| P13_PBMC | PBMC | 6.746 | 3.119 | 221 | 2.983 (95,6%) |
| P14_PBMC | PBMC | 5.411 | 1.726 | 178 | 1.658 (96,1%) |
| P15_PBMC | PBMC | 6.171 | 3.037 | 640 | 2.893 (95,3%) |
| P16_PBMC | PBMC | 4.386 | 786 | 322 | 731 (93,0%) |
| P24_PBMC | PBMC | 4.591 | 1.199 | 120 | 1.130 (94,2%) |
| P26_PBMC | PBMC | 3.595 | 914 | 71 | 861 (94,2%) |
| P18_PBMC | PBMC | 1.906 | 869 | 163 | 788 (90,7%) |
| P29_PBMC | PBMC | 6.076 | 1.265 | 144 | 1.219 (96,4%) |
| P30_PBMC | PBMC | 3.846 | 779 | 60 | 750 (96,3%) |
| P19_PBMC | PBMC | 5.959 | 1.640 | 109 | 1.567 (95,5%) |
| P20_PBMC | PBMC | 6.634 | 2.078 | 136 | 1.968 (94,7%) |
| P21_PBMC | PBMC | 3.669 | 933 | 56 | 885 (94,9%) |
| P25_PBMC | PBMC | 5.232 | 1.887 | 165 | 1.825 (96,7%) |
| P28_PBMC | PBMC | 5.227 | 2.556 | 356 | 2.435 (95,3%) |
| P22_PBMC | PBMC | 5.000 | 1.653 | 156 | 1.576 (95,3%) |
| P23_PBMC | PBMC | 5.905 | 2.557 | 353 | 2.406 (94,1%) |
| P31_PBMC | PBMC | 3.916 | 1.400 | 333 | 1.293 (92,4%) |
| P32_PBMC | PBMC | 7.112 | 1.204 | 199 | 1.148 (95,3%) |

E3 sull'analisi principale (dopo la terapia): sangue 28.116/29.602
(95,0%), tumore
26.074/29.019
(89,9%); il campione più basso è sopra
il 50%. Con «almeno un contig produttivo TRA o TRB»: 96,7%
e 91,6%. **E3 soddisfatta.**

E1: 18 pazienti (P13, P14, P15, P16, P18, P19, P20, P21, P22, P23, P24, P25, P26, P28, P29, P30, P31, P32) hanno cellule di sangue assegnate a
un solo paziente e cellule di tumore dopo la terapia, con conteggi interi e TCR per cellula.
**E1 soddisfatta.**

| Condizione | Esito |
|---|---|
| E1 | soddisfatta: 18 pazienti (≥ 5) |
| E2 | soddisfatta: etichetta CD4/CD8 per cellula in entrambi i compartimenti, da trasferimento di etichette sull'espressione |
| E3 | soddisfatta: 95,0% (sangue) e 89,9% (tumore) con TRB produttiva |
| E4 | soddisfatta: matrici 10x intere, `filtered_contig_annotations` |
| E5 | soddisfatta: sangue CD45+ non selezionato per linea; corsie selezionate escluse |

**h5ad costruiti** (`build_gse200996.py`, colonne `patient`, `compartment`, `timepoint`, `label`,
`barcode`, `sample`, `lane`): principale 58.621 cellule × 36.601 geni, 18
pazienti ({'PBMC|CD4T': 21000, 'PBMC|CD8T': 8602, 'Tumor|CD4T': 16563, 'Tumor|CD8T': 12456}); descrittivo 12.299 cellule, pazienti
P18, P24, P29, P32 (tumore prima della terapia con matrice e TCR single-cell; P23 ha la matrice ma non il TCR ed è
escluso; P27 ha il tumore ma non il sangue).

## 2. Criteri × esito (esecuzione finale, dopo la correzione del margine)

| Criterio | Esito | Numeri misurati |
|---|---|---|
| B-E1 esecuzione | **superato** | `audit-sc tcr --h5ad GSE200996_principale.h5ad --vdj-manifest manifest_principale.csv --patient-col patient --compartment-col compartment --celltype-col label --barcode-col barcode --marker-map marcatori.json --reference-compartment PBMC`: exit 0, report HTML scritto; 0:46.73 reali, 1,14 GB di memoria massima (`/usr/bin/time -v`); barcode VDJ con TRB ritrovati: 84,3% |
| B-E2 implementazione indipendente | **superato** | conteggi identici; stime puntuali: differenza massima 1.7e-18 (tolleranza 1e-9); estremi: 3.5e-18 (tolleranza 1e-6); gradi di libertà: 0.0e+00 |
| B-E3 controllo negativo | **superato** | 20 permutazioni (semi 0-19): eccesso pesato medio +0,00004 (banda [−0,01; +0,01]; da −0,0040 a +0,0051); intervallo con lo zero in 20/20 (≥ 18) |
| B-E4 controllo positivo | **non valutabile** | unità clone-tumore con identità CD4 e maggioranza corretta: 3, meno delle 20 richieste; delle 114 unità nel tumore, 110 hanno identità CD8 e 4 CD4 |
| B-E5 riferimento proteico | **non valutabile su questo dataset** | nessuna proteina di superficie per cellula: nelle matrici del sangue le sole feature non geniche sono i 7 hashtag TotalSeq-C; nel tumore nessuna; rimandato alla voce 2 della roadmap |
| B-E6 sensibilità | **superato** | differenza massima dal default 0,0031 (soglia 0,05); nessun intervallo sotto 5 pazienti; tabella sotto |

**Controllo dello script delle repliche.** Sui dati originali lo script unico (tabella costruita una
volta, funzioni pubbliche del modulo) dà gli stessi numeri di `run_tcr_validation`, la funzione
chiamata da `cli.py`, con uguaglianza esatta (eccesso_pesato, eccesso_pesi_uguali, coppie, tasso_clone, pareggi),
e coincide con l'output della riga di comando di B-E1 arrotondato a 3 decimali (tasso per clone:
sì; eccesso: sì).
Tempo dello script (25 analisi): 596 s nella prima esecuzione, 1348 s
nella seconda, in parallelo con la suite completa dei test.

### B-E2 in dettaglio

`indipendente.py` è stato scritto da un sottoagente che ha letto solo `docs/METODI.md` §1-4 e non
importa nulla da `core/`; legge da sé h5ad, manifest e file VDJ (scelta del contig, clone,
positività ai marcatori, identità di riferimento, unità, eccesso, intervalli CR2 con gradi di
libertà di Bell-McCaffrey). Le assunzioni che ha dovuto fare dove METODI.md tace sono elencate nel
suo JSON (`results/BE2_indipendente.json`, campo `assunzioni`); nessuna ha prodotto una differenza
su questi dati, quindi nessuna correzione di METODI.md è richiesta dal criterio. Due erano
potenzialmente rilevanti: il filtro `full_length` dei contig (lo strumento lo applica, METODI non lo
dice; qui non cambia nessun numero) e il confronto con il margine 0,20 (sezione 5).

| Grandezza | Strumento | Indipendente |
|---|---|---|
| `clone_tumor.n_units` | 114 | 114 |
| `clone_tumor.n_ties` | 0 | 0 |
| `clone_tumor.n_errors` | 1 | 1 |
| `clone_tumor.n_patients` | 16 | 16 |
| `clone_tumor.rate.n_groups` | 16 | 16 |
| `clone_tumor.rate.n_obs` | 114 | 114 |
| `excess_PBMC_Tumor.n_pairs` | 198 | 198 |
| `excess_PBMC_Tumor.n_patients` | 18 | 18 |
| `excess_PBMC_Tumor.weighted.n_groups` | 18 | 18 |
| `excess_PBMC_Tumor.equal.n_groups` | 18 | 18 |
| `clone_tumor.rate.mean` | 0.008771929824561403 | 0.008771929824561403 |
| `excess_PBMC_Tumor.weighted.mean` | 0.014776250952296337 | 0.014776250952296338 |
| `excess_PBMC_Tumor.equal.mean` | 0.014982140799000371 | 0.014982140799000369 |
| `clone_tumor.rate.ci_low` | 0.00071906627052457776 | 0.00071906627052457776 |
| `clone_tumor.rate.ci_high` | 0.09815118259909579 | 0.09815118259909579 |
| `excess_PBMC_Tumor.weighted.ci_low` | -0.00036124911545377079 | -0.00036124911545376732 |
| `excess_PBMC_Tumor.weighted.ci_high` | 0.029913751020046446 | 0.029913751020046446 |
| `excess_PBMC_Tumor.equal.ci_low` | -0.00054586202673905654 | -0.00054586202673905654 |
| `excess_PBMC_Tumor.equal.ci_high` | 0.030510143624739798 | 0.030510143624739795 |

### B-E3 per replica (eccesso pesato sangue-tumore)

0: +0,0041; 1: +0,0021; 2: −0,0003; 3: −0,0039; 4: −0,0028; 5: +0,0029; 6: −0,0002; 7: +0,0018; 8: −0,0040; 9: −0,0035; 10: −0,0002; 11: +0,0050; 12: −0,0038; 13: −0,0021; 14: +0,0001; 15: +0,0049; 16: −0,0003; 17: −0,0029; 18: +0,0051; 19: −0,0011;

### B-E6

| Variante | Pazienti | Unità | Pareggi | Unità in errore | Tasso per clone nel tumore [IC95%] | Differenza dal default |
|---|---|---|---|---|---|---|
| TRB, pareggi esclusi (default) | 16 | 114 | 0 | 1 | 0,009 [0,001; 0,098] | +0,0000 |
| TRA+TRB, pareggi esclusi | 14 | 84 | 0 | 1 | 0,012 [0,001; 0,144] | +0,0031 |
| TRB, pareggi = errore | 16 | 114 | 0 | 1 | 0,009 [0,001; 0,098] | +0,0000 |
| TRB, pareggi = corretti | 16 | 114 | 0 | 1 | 0,009 [0,001; 0,098] | +0,0000 |

Nessuna unità in pareggio: le varianti «pareggi = errore» e «pareggi = corretti» coincidono per
costruzione con il default. Variante TRA+TRB: 48.963 cellule con entrambe le
catene su 54.190 con TRB.

## 3. Prima e dopo la correzione

La correzione del margine (sezione 5) è stata fatta dopo una prima esecuzione completa con il codice
0.9.0 (risultati in `results/prima_della_correzione/`, commit `be4d239`). Esito prima della
correzione: B-E1 superato (stesso output della CLI), B-E2
superato (differenza massima delle stime 1.7e-18),
B-E3 superato, B-E4 non valutabile
(3 candidate), B-E6 superato.
Dopo la correzione i numeri sono identici: output della CLI identico,
analisi descrittiva identica, script delle repliche
strumento_BE2 identico, BE3 identico, BE4 identico, BE6 identico, descrittivi identico. Nessun clone di
GSE200996 ha un margine esattamente pari a 0,20.

## 4. Risultati descrittivi (non sono criteri)

**Tumore dopo la terapia** (18 pazienti nell'h5ad principale):
- tasso d'errore per clone (convenzione del caso studio): 0,009 [0,001; 0,098], 114 unità da
  16 pazienti, 1 unità in errore, nessun pareggio;
- tasso d'errore per cellula: 0,006 [0,001; 0,033],
  1234 cellule da 18 pazienti;
- eccesso di discordanza sangue-tumore: pesato +0,015 [−0,000; +0,030], a pesi uguali
  +0,015 [−0,001; +0,031], su 198 coppie clone-compartimenti da 18 pazienti
  (l'intervallo include lo zero di poco);
- direzione degli errori (identità → etichetta di maggioranza nel tumore): CD4→CD8: 1; per paziente:
  P31 (1); nessun errore CD8→CD4;
- cloni con identità di riferimento: 394; cellule con TRB abbinata:
  54.190; flag per cellula valutabili: 2,1% delle cellule.

**Tumore prima della terapia** (h5ad descrittivo, `audit-sc tcr` eseguito una volta: exit 0;
barcode VDJ ritrovati 84,3%): 4 pazienti, quindi nessun intervallo. Tasso per
clone 0,000 su 21 unità (pareggi: 0); per cellula
0,000 su 175 cellule; eccesso non riportato dallo strumento con
meno di 5 pazienti.

**Come leggere questi numeri.** Il tasso d'errore è molto più basso di quello del caso studio su
GSE278694 (0,170 per clone nel tumore), ma il confronto non è alla pari, per tre ragioni (le prime due
scritte in `MAPPA_ETICHETTE.md` prima dell'esecuzione): (1) gli autori hanno tolto dalle etichette CD4 le cellule
con CD8A/CD8B rilevati e dalle CD8 quelle con CD4 rilevato, cioè hanno già ripulito le etichette con
gli stessi geni che il Modulo B usa per l'identità di riferimento; (2) le etichette vengono da un
trasferimento da un riferimento CITE-seq con soglia di predizione, non da cluster; (3) i metadati
contengono solo cellule T CD4/CD8: gli errori verso NK del caso studio qui non possono comparire, e il
confronto con il caso studio riguarda solo CD4↔CD8. Inoltre quasi tutte le unità giudicate nel tumore
appartengono a cloni CD8 (i cloni CD4 con almeno 3 cellule nel sangue e 3 nel tumore sono rari), e
per questo il controllo positivo B-E4 non è valutabile.

## 5. Difetti trovati e correzioni

**Difetto dello strumento: margine 0,20 confrontato in virgola mobile** (corretto nel commit
`8b8e761`). METODI §4 dice: identità = etichetta con la frazione più alta «se il margine sulla
seconda è almeno 0,20». `assign_reference_identity` confrontava `best - second < 0.20` in virgola
mobile: un margine esattamente 0,20 era rifiutato in alcuni casi (3/5 contro 2/5 = 0,19999999999999996;
7/10 contro 5/10; 3/10 contro 1/10; 6/10 contro 4/10) e accettato in altri (1/5 contro 0/5; 4/10
contro 2/10; 8/10 contro 6/10). L'ha fatto emergere `indipendente.py`, che confronta in aritmetica
esatta e aveva elencato il punto fra le sue assunzioni; verificato direttamente sulla funzione dello
strumento. Il caso studio (pdac-ml, `src/10_loco.py`, `abs(f8 - f4) < MIN_MARGIN`) ha lo stesso confronto,
quindi lo strumento la riproduceva anche in questo. Correzione: confronto con tolleranza 1e-9 (un
margine diverso da 0,20 ne dista almeno 1/(5n) per un clone di n cellule). Test di regressione
`tests/test_reference_margin.py`: fallisce sul codice precedente (5 casi su 8 con margine
esattamente 0,20), passa dopo. Effetto: nessuno su GSE200996 (sezione 3) e nessuno su GSE278694
(`tests/test_real_data_gse278694.py` passa con gli stessi numeri: 0,1697 su 277 unità).

**Difetti degli script di validazione** (non dello strumento), corretti prima di eseguire il modulo:
- nella prima versione di `build_gse200996.py` i percorsi del manifest erano relativi alla cartella
  di lavoro, mentre la guida (e `read_vdj_manifest`) li legge rispetto alla cartella del manifest:
  resi assoluti;
- `GSE200996_series_matrix.txt.gz` non esiste (404): la serie ha una series matrix per piattaforma.

**Lacune di METODI.md senza effetto sui numeri** (dalle assunzioni di `indipendente.py`): i filtri
dei contig (`is_cell`, `high_confidence`, `productive`, `full_length`), la regola di abbinamento dei
barcode e del suffisso `-1`, il trattamento delle cellule senza TRB e il quantile t con gradi di
libertà non interi non sono descritti. Il criterio B-E2 chiede di correggere METODI.md solo quando
una lacuna produce una differenza; qui non l'ha prodotta, e METODI.md non è stato modificato.

## 6. Verdetto con i criteri originali

**Validazione tecnica esterna superata nei criteri valutabili: B-E1, B-E2, B-E3 e B-E6 superati
(4 su 4); B-E4 non valutabile (3 unità candidate, ne servono 20); B-E5 non valutabile (nessuna
proteina di superficie).** Nessun criterio è stato modificato dopo aver visto i dati: non c'è una
sezione «Criteri modificati dopo aver visto i dati». Esito identico prima e dopo la correzione del
margine.

## 7. Che cosa questa validazione dimostra e che cosa no

Dimostra, su dati di un altro laboratorio (Dana-Farber, 10x 5′, sangue demultiplato con hashtag):
- che il Modulo B si esegue con i soli comandi e input della guida (B-E1), su 18 pazienti e circa
  59.000 cellule in meno di un minuto;
- che i suoi calcoli coincidono con un'implementazione scritta solo da METODI.md §1-4 (B-E2), e che
  quell'implementazione ha trovato un difetto vero, piccolo e ora corretto;
- che non segnala un eccesso di discordanza quando il legame fra clone e compartimento è rotto
  (controllo negativo, B-E3);
- che il tasso d'errore per clone è stabile alla definizione del clone (TRB o TRA+TRB) e alla regola
  dei pareggi (B-E6), qui senza pareggi.

Non dimostra:
- che il modulo rilevi errori iniettati su dati reali: B-E4 non è valutabile su GSE200996 (la
  sensibilità resta provata solo sui dati sintetici dei test);
- che le etichette segnalate come errate lo siano davvero: serve un riferimento proteico (B-E5,
  voce 2 della roadmap);
- che il tasso d'errore basso di GSE200996 descriva la qualità di un'annotazione qualsiasi: le
  etichette degli autori erano già state ripulite con gli stessi geni marcatori usati dal modulo.

## 8. Test

Codice cambiato (`core/tcr_validation.py`): suite completa su un clone pulito del branch
(`uv sync --extra dev`, dati di GSE278694 con `GSE278694_DIR`). Esito: 206 superati, 0 falliti, 7 fallimenti attesi (xfail strict, gli stessi della 0.9.0), nessun test saltato (inclusi i test su GSE278694), in 27 min 40 s con 6 processi (commit 8b8e761).
