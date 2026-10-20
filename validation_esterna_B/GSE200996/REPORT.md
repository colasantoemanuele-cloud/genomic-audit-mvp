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

**Eleggibilità: E1-E5 soddisfatte.** Questa versione del report è scritta prima di eseguire il Modulo B; i criteri B-E1…B-E6 seguono nelle sezioni successive dopo l'esecuzione.
