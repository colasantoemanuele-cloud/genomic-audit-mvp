# Dati della demo

Questa cartella contiene **dati pubblici di NCBI GEO, ridistribuiti in forma ridotta al solo
scopo della demo** («Demo immediata» della web app e test). Non contiene dati clinici
identificabili: solo conteggi di espressione e identificativi di campione e di paziente già
pubblicati dagli autori. I diritti sui dati restano degli autori originali: chi li usa deve
citare gli articoli di origine, non questa repo. Per qualunque analisi vanno scaricati i dataset
completi da GEO.

| File | Accessione | Articolo di origine | Contenuto e riduzione |
|---|---|---|---|
| `GSE125449_demo.h5ad` | [GSE125449](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE125449) | Ma L. et al., «Tumor Cell Biodiversity Drives Microenvironmental Reprogramming in Liver Cancer», *Cancer Cell* 36(4):418-430, 2019 | Tumori primitivi del fegato. Conteggi grezzi di **1.861 cellule su 9.946**: 10 pazienti estratti a caso fra i 19 (numpy `default_rng(0)`, senza reinserimento) e, per ciascuno, al massimo 200 cellule estratte a caso (stesso seme). 18.372 geni comuni ai due set del deposito; nessun filtro sui geni né sulle classi. Etichetta `Type` degli autori, invariata. |
| `GSE125449_samples.csv` | GSE125449 | come sopra | Metadati dei 19 campioni, dal *series matrix* (due piattaforme). Il paziente è ricavato dal titolo del campione (`LCP<id>`). |
| `GSE132465_samples.csv` | [GSE132465](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE132465) | Lee H.-O. et al., «Lineage-dependent gene expression programs influence the immune landscape of colorectal cancer», *Nature Genetics* 52:594-603, 2020 | Carcinoma colorettale. Solo metadati: 33 campioni, 23 pazienti, tumore e mucosa normale. Nessun dato di espressione. |
| `GSE131907_samples.csv` | [GSE131907](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE131907) | Kim N. et al., «Single-cell RNA sequencing demonstrates the molecular and cellular reprogramming of metastatic lung adenocarcinoma», *Nature Communications* 11:2285, 2020 | Adenocarcinoma polmonare. Solo metadati: 58 campioni, 44 pazienti, 7 tessuti. Nessun dato di espressione. |

## Come sono stati ridotti

1. I file sono stati scaricati da GEO in `validation_esterna/data/` (cartella non versionata).
2. `validation_esterna/parse_geo.py` estrae i metadati per campione dai *series matrix*.
3. `validation_esterna/build_gse125449.py` unisce i due set di GSE125449 sui geni comuni.
4. `validation_esterna/build_small.py` estrae la coorte ridotta (10 pazienti, al massimo 200
   cellule ciascuno, seme 0); `validation_esterna/build_demo.py` la copia qui in formato
   compresso.

Le regole della riduzione sono quelle della sezione 5 di `validation_esterna/CRITERI.md`. La
coorte ridotta coincide, cellula per cellula, con quella usata nella validazione esterna.

## Che cosa non c'è

Il Modulo B nella demo usa dati **sintetici** (`core/synthetic.py`): la repo non include dati TCR
reali. I dati di GSE278694 (Chen et al., *Cancer Cell* 2025), usati per validare il Modulo B,
non sono ridistribuiti qui.
