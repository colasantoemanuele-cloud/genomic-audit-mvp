# C2 — Prova della guida su un dataset mai usato

Dataset: **GSE144735** (Lee H.-O. et al., *Nature Genetics* 2020; coorte KUL3, carcinoma
colorettale: 6 pazienti, tumore, bordo e mucosa normale; 27.414 cellule, 33.694 geni). Mai usato
nello sviluppo: né i metadati né i conteggi. File scaricati da GEO
(`GSE144735_processed_KUL3_CRC_10X_raw_UMI_count_matrix.txt.gz`, 57 MB;
`GSE144735_processed_KUL3_CRC_10X_annotation.txt.gz`) in una cartella fuori dalla repo
(`../dati_esterni/GSE144735/`). Comandi eseguiti seguendo solo `docs/GUIDA_UTENTE.md`, con lo
strumento installato nell'ambiente del progetto. Nessun dato del dataset è versionato.

## Passi e lacune

| Passo | Che cosa diceva la guida | Che cosa è successo | Lacuna e correzione |
|---|---|---|---|
| 1. Preparare l'`.h5ad` | Solo conversione da Seurat (`.rds`) e cartelle 10x | I conteggi di GEO sono una tabella di testo geni × cellule: la guida non spiegava come leggerla | **Lacuna 1.** Aggiunta la sezione «Da una matrice di conteggi in testo»: lettura a blocchi (la matrice densa occuperebbe circa 7 GB). Seguita alla lettera: 6 min 47 s, 2,4 GB di picco, conteggi interi verificati con il controllo della guida |
| 2. Audit del disegno | `--h5ad` riduce le cellule alle combinazioni distinte | Eseguito senza errori; lo strumento avvisa che senza colonna libreria le librerie replicate collassano | **Lacuna 2.** La guida non diceva di dichiarare la colonna del campione (`Sample`). Aggiunto: dichiararla sempre con `--technical library=<colonna>`. Con la colonna dichiarata compare il fatto strutturale «unità tecnica» (ogni libreria = una coppia paziente-tessuto), in giallo |
| 2b. Livello scritto male (`Class:Tumor:Normale`) | Il confronto risulta «non valutato» in grigio, con i livelli presenti | Grigio come previsto, ma la riga della sintesi non elencava i livelli presenti (c'erano solo nel testo esteso) | **Lacuna 3.** Corretto il codice: la riga grigia riporta la frase con i livelli presenti |
| 3. Modulo A su `Cell_type` | Comando della sezione 3.2 | Eseguito in 34 s; macro-F1 0,999 in entrambi gli schemi, divario +0,001 entro la soglia di rumore, nessuna direzione indicata; training limitato a circa 20.000 cellule per fold e dichiarato | **Lacuna 4.** La guida non avvertiva che con tipi cellulari grossolani il controllo è poco informativo. Aggiunto |
| 3b. Modulo A su `Cell_subtype` | Comando della sezione 3.2 | Eseguito (40 sottotipi, 6 pazienti, 5 fold): split per paziente 0,732 ± 0,024, split casuale 0,845 ± 0,005, divario +0,114 oltre la soglia; classi assenti in 2 fold su 5, dichiarate; 895 s su una macchina occupata da altre esecuzioni | **Lacuna 5.** La guida non diceva come togliere o accorpare le classi rare se una classe compare in un solo paziente: aggiunto un esempio in Python (qui non è servito: ogni sottotipo compare in almeno 2 pazienti) |
| 4. Modulo B | Servono i file VDJ | GSE144735 non ha dati TCR: non applicabile, come la guida dice | nessuna |

## Esiti

Tutti i passi della guida applicabili a GSE144735 sono stati eseguiti dall'inizio alla fine senza
modifiche al motore, dopo aver corretto le lacune 1, 2, 4 e 5 nella guida e la lacuna 3 nel codice
(`core/verdict.py`). Output copiati in `revisione/log/c2/` (i report HTML restano fuori dalla repo, in `../dati_esterni/GSE144735/prova/risultati/`): `disegno_1.txt`,
`disegno_2.txt`, `moduloA.txt`, `moduloA_sottotipi.txt` e i report HTML in `risultati/`.

Nota: la prova è stata fatta da chi ha sviluppato lo strumento, seguendo il testo della guida; non
sostituisce una prova con un utente esterno.
