# Dati della demo

Dati pubblici di NCBI GEO, usati dalla modalità "Demo immediata" della web app. Nessun dato
clinico identificabile: solo identificativi di campione/paziente già pubblicati.

| File | Origine | Contenuto |
|---|---|---|
| `GSE125449_demo.h5ad` | [GSE125449](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE125449), Ma et al., *Cancer Cell* 2019 (tumori primitivi del fegato) | conteggi grezzi di 1.861 cellule (10 pazienti su 19, al massimo 200 cellule per paziente, seme 0), 18.372 geni comuni ai due set, etichetta `Type` degli autori |
| `GSE125449_samples.csv` | series matrix di GSE125449 (due piattaforme) | 19 campioni; paziente ricavato dal titolo (`LCP<id>`) |
| `GSE132465_samples.csv` | [GSE132465](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE132465), carcinoma colorettale | 33 campioni, 23 pazienti, tumore e mucosa normale |
| `GSE131907_samples.csv` | [GSE131907](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE131907), adenocarcinoma polmonare | 58 campioni, 44 pazienti, 7 tessuti |

Ricostruzione: `validation_esterna/parse_geo.py`, `build_gse125449.py`, `build_small.py`,
`build_demo.py` (dopo aver scaricato i file GEO in `validation_esterna/data/`).

Il Modulo B nella demo usa dati **sintetici** (`core/synthetic.py`): il progetto non include
dati TCR reali.
