# Guida per l'utente

Questa guida è per chi usa `genomic-audit` sui propri dati: un biologo, un clinico o un
bioinformatico che non ha sviluppato lo strumento. Spiega che cosa serve in ingresso, come
lanciare ciascuna parte, come leggere i risultati e che cosa i risultati **non** dicono.

> **Strumento per la ricerca. Non è un dispositivo medico e non va usato per decisioni cliniche
> su singoli pazienti.**

## 1. Che cosa fa, in una frase per parte

| Parte | Domanda a cui risponde | Che cosa serve |
|---|---|---|
| Audit del disegno | Quali confronti il disegno dello studio permette di stimare, e con quanti pazienti indipendenti? | Solo il foglio dei campioni (CSV/TSV) |
| Modulo A | Quanto cambia l'accuratezza di un classificatore cellulare se la valutazione mescola cellule dello stesso paziente fra training e test? | Un `.h5ad` con conteggi grezzi, paziente ed etichetta |
| Modulo B | Le cellule dello stesso clone T hanno la stessa etichetta CD4/CD8 in sangue, tumore e tessuto normale? | Il `.h5ad` più i file VDJ di Cell Ranger, un manifest e un compartimento di riferimento |

Tutto gira sul computer dell'utente: nessun dato esce dalla macchina (sezione 8).

## 2. Installazione

Le istruzioni sono nel [README](../README.md#avvio-rapido): con `uv` bastano `uv sync` e poi
`uv run audit-sc demo`. Per una prova senza dati propri: `audit-sc demo`, che scrive
`results/demo_report.html`.

In questa guida i comandi sono scritti come `audit-sc ...`. Se lo strumento è installato con
`uv`, si scrivono `uv run audit-sc ...`.

## 3. Che cosa serve in ingresso

### 3.1 Audit del disegno: il foglio dei campioni

Una tabella CSV o TSV con **una riga per campione o per libreria** (non per cellula). Le colonne
sono quelle del disegno; i nomi sono liberi, perché al comando si dice quale colonna ha quale
ruolo.

| patient | tissue | library | batch | protocol | response |
|---|---|---|---|---|---|
| P01 | Tumor | L01_T | run1 | scRNA | R |
| P01 | Normal | L01_N | run1 | scRNA | R |
| P02 | Tumor | L02_T | run2 | scRNA | NR |

Ruoli possibili: paziente (obbligatorio), tessuto o condizione, fattori tecnici (`library`,
`batch`, `chemistry`, `protocol`, `date`) ed esiti a livello di paziente (es. risposta alla
terapia). Le celle vuote diventano un livello esplicito `NA`, dichiarato nel report. Se una colonna
identifica il campione o la libreria (spesso si chiama `sample`, `Sample` o `library`), va sempre
dichiarata con `--technical library=<colonna>`: lo strumento segnala allora se ogni libreria
coincide con una coppia paziente-tessuto, un limite comune che altrimenti resta nascosto.

Se la tabella ha molte righe ripetute con gli stessi valori, lo strumento si ferma: è il segno di
una tabella per cellula, oppure di librerie replicate senza la colonna libreria. In quel caso si
indica la colonna libreria, oppure si usa `--h5ad`, che riduce da solo le cellule alle
combinazioni distinte (e avvisa che, senza colonna libreria, le librerie replicate collassano).

```bash
audit-sc design --meta campioni.csv --patient-col patient --tissue-col tissue \
    --technical library=library --technical batch=batch --technical protocol=protocol \
    --outcome-col response \
    --compare tissue:Tumor:Normal --compare response:R:NR \
    --out risultati/disegno.html
```

`--compare colonna:livello_a:livello_b` chiede di valutare un confronto; la colonna deve essere
fra quelle dichiarate. Un livello scritto male (es. `Tumour` invece di `Tumor`) non blocca
l'esecuzione: il confronto risulta «non valutato» in grigio, con l'elenco dei livelli presenti.

### 3.2 Modulo A: matrice con conteggi grezzi

Un file AnnData (`.h5ad`) con:
- in `X` i **conteggi grezzi** (interi, non normalizzati, non log-trasformati): lo strumento
  normalizza da solo dentro ogni fold e si ferma se trova valori non interi o negativi;
- in `obs` una colonna con il paziente e una con l'etichetta da classificare (es. il tipo
  cellulare).

In alternativa, una cartella 10x (`matrix.mtx`, `barcodes.tsv`, `features.tsv`, anche `.gz`) più
una tabella per cellula con la colonna dei barcode (solo dall'app web).

```bash
audit-sc leakage --h5ad dati.h5ad --target-col celltype --patient-col patient \
    --out risultati/moduloA.html
```

Con almeno 8 pazienti lo strumento esegue anche un confronto fra tre modelli lasciando fuori un
paziente alla volta; richiede più tempo (sezione 9). `--rapido` lo salta.

Ogni classe deve comparire in almeno 2 pazienti, altrimenti lo strumento si ferma («il training
contiene una sola classe»). Per togliere o accorpare le classi rare si crea una nuova colonna, per
esempio in Python:

```python
import anndata as ad
a = ad.read_h5ad("studio.h5ad")
pazienti = a.obs.groupby("celltype", observed=True)["patient"].nunique()
rare = pazienti[pazienti < 2].index
a.obs["celltype_ridotto"] = a.obs["celltype"].astype(str).where(~a.obs["celltype"].isin(rare), "altro")
a.write_h5ad("studio.h5ad")      # poi: --target-col celltype_ridotto
```

Con tipi cellulari grossolani (epiteliali, T, B, mieloidi) la macro-F1 è spesso vicina a 1 in
entrambi gli schemi e il divario è nullo: il controllo è più informativo con etichette più fini
(sottotipi) o con etichette di condizione.

### 3.3 Modulo B: matrice, file VDJ, manifest e riferimento

Servono:
1. il `.h5ad` con conteggi grezzi (servono per i geni marcatori) e in `obs` le colonne
   paziente, compartimento (es. `PBMC`, `Tumor`, `Normal`), tipo cellulare e barcode;
2. i file `filtered_contig_annotations.csv` di Cell Ranger VDJ, uno per libreria;
3. un **manifest**, cioè una tabella che dice a quale paziente e compartimento appartiene ogni
   file VDJ (i file di Cell Ranger non lo contengono);
4. il compartimento di riferimento, da cui si stima l'identità del clone: di solito il sangue;
5. una mappa dei marcatori in JSON, per esempio `{"CD4T": ["CD4"], "CD8T": ["CD8A", "CD8B"]}`:
   le chiavi devono essere esattamente le etichette usate nella colonna del tipo cellulare.

Esempio di manifest (`vdj_manifest.csv`; i percorsi relativi sono letti rispetto alla cartella
del manifest):

```
path,patient,compartment
vdj/P01_PBMC_filtered_contig_annotations.csv,P01,PBMC
vdj/P01_Tumor_filtered_contig_annotations.csv,P01,Tumor
```

I valori di `patient` e `compartment` devono coincidere con quelli dell'`.h5ad`. Se lo stesso
paziente e compartimento hanno **più librerie**, si aggiunge al manifest una colonna `library` e
si indica la colonna libreria dell'`.h5ad` con `--library-col`: senza, lo strumento si ferma,
perché lo stesso barcode potrebbe stare in due librerie. Il suffisso `-1` dei barcode di Cell
Ranger viene tolto in automatico se serve, e il report lo dichiara.

```bash
audit-sc tcr --h5ad dati.h5ad --vdj-manifest vdj_manifest.csv \
    --patient-col patient --compartment-col compartment --celltype-col celltype \
    --barcode-col barcode --marker-map marcatori.json --reference-compartment PBMC \
    --out risultati/moduloB.html
```

Opzioni utili: `--export-flags` scrive una copia dell'`.h5ad` con due colonne di flag per cellula
e un CSV; `--compartment-order Tumor,Normal` fissa il verso delle differenze fra compartimenti;
`--convention cell|clone|both` sceglie la definizione del tasso d'errore (default: entrambe).
Nella convenzione per clone le unità con l'etichetta di maggioranza in pareggio sono escluse e
contate a parte (il report dice quante sono); `--tie-rule error` o `--tie-rule correct` le conta
come errore o come corrette: vale la pena eseguire le due sensibilità quando i pareggi sono
più di qualche unità.

## 4. Da Seurat (`.rds`) a `.h5ad`

Molti gruppi lavorano in R. La strada consigliata è **zellkonverter** (Bioconductor), che permette
di scegliere quale matrice va in `X`:

```r
library(Seurat)
library(zellkonverter)              # BiocManager::install("zellkonverter")
obj <- readRDS("studio.rds")
obj <- JoinLayers(obj)              # solo Seurat v5, se i conteggi sono divisi in più layer
obj$barcode <- sub(".*_", "", colnames(obj))   # barcode 10x puro, se i nomi hanno un prefisso
sce <- as.SingleCellExperiment(obj, assay = "RNA")
writeH5AD(sce, "studio.h5ad", X_name = "counts")   # i conteggi GREZZI in X
```

Con **SeuratDisk** (`SaveH5Seurat` e poi `Convert(..., dest = "h5ad")`) in `X` non finiscono i
conteggi: secondo il contenuto dell'oggetto, `X` contiene `data` (normalizzata) o `scale.data`, e
`raw` contiene `counts` oppure `data`. Lo strumento allora si ferma con «valori non interi». Se
`raw` contiene i conteggi, si possono rimettere in `X` in Python, verificando poi che siano interi
(controlli qui sotto):

```python
import anndata as ad
a = ad.read_h5ad("studio.h5ad")
c = a.raw.to_adata()        # da usare solo se raw contiene i conteggi grezzi
c.obs = a.obs
c.write_h5ad("studio_conteggi.h5ad")
```

### Da una matrice di conteggi in testo (es. file supplementari di GEO)

Molti dataset pubblici danno i conteggi grezzi come tabella di testo, una riga per gene e una
colonna per cellula, più una tabella di annotazione per cellula. La matrice va letta a blocchi:
caricata tutta insieme in forma densa può richiedere decine di GB.

```python
import anndata as ad, numpy as np, pandas as pd, scipy.sparse as sp

parts, genes = [], []
for chunk in pd.read_csv("conteggi.txt.gz", sep="\t", index_col=0, chunksize=2000):  # 2000 geni per volta
    parts.append(sp.csr_matrix(chunk.to_numpy(dtype=np.float32)))
    genes.extend(chunk.index.astype(str))
cells = chunk.columns.astype(str)
X = sp.vstack(parts).T.tocsr()                                  # cellule x geni
obs = pd.read_csv("annotazione.txt.gz", sep="\t", index_col=0, dtype=str)
obs = obs.loc[cells]                                            # stesse cellule, stesso ordine
a = ad.AnnData(X=X, obs=obs, var=pd.DataFrame(index=genes))
a.var_names_make_unique()
a.write_h5ad("studio.h5ad")
```

`obs.loc[cells]` si ferma con un errore se una cellula della matrice non ha annotazione: è voluto.
Se i nomi delle colonne (separatore, colonna indice) sono diversi, vanno adattati gli argomenti di
`read_csv`.

Dopo la conversione, controllare:
- che `X` contenga interi: `python -c "import anndata as ad, numpy as np; a=ad.read_h5ad('studio.h5ad'); x=a.X[:1000]; x=x.toarray() if hasattr(x,'toarray') else x; print(np.allclose(x, np.round(x)), x.min())"`
  deve stampare `True` e un minimo `>= 0`;
- i nomi delle colonne di `obs` (`print(list(a.obs.columns))`): paziente, compartimento, tipo
  cellulare e barcode devono esserci, con gli stessi valori del manifest VDJ;
- che i barcode in `obs` siano la sequenza 10x (con o senza `-1`), senza il prefisso del campione;
- che le etichette del tipo cellulare coincidano con le chiavi della mappa dei marcatori (es.
  `CD4T` e non `CD4+ T cells`).

## 5. Come leggere la «Sintesi dei controlli»

All'inizio di ogni report c'è una tabella con una riga per controllo. Ogni riga ha uno di quattro
stati. Non c'è un colore complessivo: lo strumento non dà un voto allo studio.

| Stato | Testo | Esempio | Che cosa fare |
|---|---|---|---|
| verde | nessun limite rilevato dalle regole | «Confronto tissue: Tumor vs PBMC — stimabile con 10 unità indipendenti» | Leggere comunque il numero e la sua definizione |
| giallo | stima con limiti | «Confronto tissue: Tumor vs Normal — stimabile con bassa potenza: 5 unità indipendenti» | La stima c'è, ma va riportata insieme al limite scritto nella riga |
| rosso | stima non possibile con questi dati | «Confronto protocol: scRNA vs snRNA — non stimabile: nessun paziente ha entrambi i livelli» | Non trarre conclusioni su quel confronto da questi dati |
| grigio | controllo non eseguito o non valutabile | «Confronto tissue: Tumor vs Tumour — non valutato, errore d'input» | Correggere l'input, oppure leggere perché il controllo non era possibile |

**Che cosa il verde NON significa.** Il verde vuol dire solo che nessuna delle regole scritte
nel report ha segnalato un limite. Non dice che la stima sia corretta, che le etichette siano
giuste o che lo studio sia ben disegnato. Le regole non vedono, per esempio, le ipotesi biologiche
del Modulo B (sezione 6), la qualità dell'annotazione di partenza o un errore di laboratorio. In
particolare:
- un eccesso di discordanza in verde («nessun eccesso rilevato») non dimostra che le etichette
  siano corrette: con pochi pazienti un errore moderato può non essere rilevato;
- un tasso d'errore in verde vuol dire che l'intervallo è stato prodotto, non che il tasso sia
  basso: il numero va letto;
- il confronto fra modelli del Modulo A, quando eseguito, è giallo per costruzione: è solo
  descrittivo.

Le regole esatte sono stampate nel report, sotto la tabella.

## 6. Limiti, in breve

- **Modulo B, ipotesi biologica.** Si assume che le cellule dello stesso clone T abbiano la stessa
  linea CD4/CD8. Eccezioni note: CD4 citotossiche, cellule CD4 con CD8αα, MAIT, NKT e γδ, doppietti
  T-NK, RNA ambientale. Il metodo ha senso solo con etichette di linea (CD4/CD8), non di stato
  (effettrice, esausta).
- **Modulo B, una sola prova su dati reali.** È stato eseguito su GSE278694, il dataset su cui è
  nato: lì si è verificata la fedeltà del porting, non la correttezza del metodo.
- **Convenzione «per clone».** Di default conta come errore ogni etichetta fuori dalla mappa dei
  marcatori, anche mieloide o B; `--clone-error-labels` limita le etichette ammesse.
- **Pochi pazienti.** Sotto 5 pazienti nessun intervallo; con 5-10 pazienti gli intervalli sono
  larghi, ed è l'informazione principale che danno.
- **Modulo A.** Il divario fra split casuale e split per paziente include tutto ciò che distingue
  i pazienti, anche piattaforma o batch quando coincidono con il paziente. La deviazione standard
  fra 5 fold non è una misura valida dell'incertezza.
- **Audit del disegno.** I fatti strutturali (annidamenti, coincidenze) sono controlli
  deterministici sulla tabella, utili come lista di verifica; l'unico elemento statistico è il
  Cramér V, che con poche righe spesso non è valutabile.
- **Frazione di CD8.** È sperimentale e disponibile solo da riga di comando
  (`--sperimentale-cd8`): sui dati reali provati non è utilizzabile.

## 7. Web app

`audit-sc serve` apre l'interfaccia su `http://localhost:8501`. La modalità «Carica studio»
accetta gli stessi file della riga di comando. L'app è pensata per una workstation con **un solo
utente**: «localhost» non protegge da altri utenti della stessa macchina o di un server condiviso.

## 8. Dove finiscono i dati

- I dati di partenza non vengono mai modificati.
- La riga di comando scrive solo dove si indica: il report (`--out`) e, con `--export-flags`, una
  copia `.h5ad` e un CSV nella stessa cartella del report. Per cancellare tutto basta cancellare
  quella cartella.
- La demo scrive in `results/` e in `data/synthetic/` (dati sintetici).
- L'app web legge i file caricati in memoria, senza scriverne copie su disco; i report si
  scaricano dal browser, nella cartella scelta dall'utente.
- Nessun dato e nessuna statistica d'uso escono dalla macchina: niente account, niente servizi
  esterni, telemetria di Streamlit disattivata.

## 9. Tempi e memoria

Ordini di grandezza misurati su una workstation con 12 core: la conversione di una matrice di testo
da 27.000 cellule e 34.000 geni richiede alcuni minuti e circa 2,5 GB; il Modulo A su 10.000 cellule
con il confronto fra modelli circa 7 minuti; il Modulo B su 200.000 cellule circa 2 minuti e 6 GB
(il file `.h5ad` viene caricato tutto in memoria).

## 10. Messaggi d'errore frequenti

| Messaggio (inizio) | Causa | Che cosa fare |
|---|---|---|
| `la matrice non contiene conteggi interi` / `valori non interi` | `X` è normalizzata | Rimettere i conteggi grezzi in `X` (sezione 4) |
| `solo il ...% dei barcode VDJ è stato ritrovato` | paziente, compartimento o barcode diversi fra manifest e `.h5ad` | Confrontare gli esempi stampati nel messaggio |
| `due file VDJ ... per paziente ...` | più librerie per la stessa coppia | Colonna `library` nel manifest e `--library-col` |
| `righe ... ripetono la stessa combinazione` | foglio per cellula o librerie senza colonna | Usare `--h5ad` o dichiarare la colonna libreria |
| `la colonna '...' non è dichiarata` | `--compare` su una colonna senza ruolo | Aggiungere `--tissue-col`, `--technical` o `--outcome-col` |
| `nel fold ... il training contiene una sola classe` | una classe compare in un solo paziente | Togliere o accorpare quella classe |
