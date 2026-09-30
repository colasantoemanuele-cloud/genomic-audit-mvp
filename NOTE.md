# Note: idee scartate, rimandate o emerse durante il lavoro

Nessuna delle voci qui sotto e' implementata. Sono annotate per non perderle, non come
impegno.

## Scartate (dalla valutazione strategica per l'IRCCS De Bellis)

- **Modulo CNV / field effect sui margini chirurgici.** Doppione di inferCNV, CopyKAT e
  Numbat; il passaggio "tessuto adiacente -> margine chirurgico" non e' sostenuto dai dati
  (4 pazienti valutabili in GSE278694). Se serve, basta un controllo "separabilita' entro
  paziente contro aggregata" su etichette di malignita' gia' calcolate, che e' lo stesso
  principio del Modulo A.
- **Correzione automatica delle etichette** (`adata_cleaned.h5ad`). L'audit diventerebbe un
  annotatore e perderebbe l'indipendenza. Sostituita dai flag per cellula (Intervento 2),
  che segnalano senza riscrivere.
- **BCR / strutture linfoidi terziarie.** L'assunto "stesso clone, stessa identita'" non
  vale per le cellule B (ipermutazione somatica, switch isotipico, differenziazione in
  plasmacellula); le TLS sono strutture spaziali che la dissociazione distrugge.
- **Endpoint time-to-event (Cox, DFS).** Con 8-15 eventi produce proprio le sovrastime che lo
  strumento dovrebbe prevenire.
- **Rilevatore di doppietti per cellula.** Compete con scDblFinder e Scrublet. In `pdac-ml` la
  stima dei doppietti T-NK (72-79%) e' una stima di POPOLAZIONE sulle cellule NK-TCR+, non un
  classificatore per cellula: trasformarla in un flag richiederebbe una regola e una
  calibrazione nuove.

## Rimandate

- **Scomposizione della varianza nel report.** Implementata e testata in
  `core/design_audit.py` (`variance_decomposition`), ma NON esposta: la copertura del suo
  intervallo per la quota del tessuto e' 200/200 = 1.000, sopra la banda [0.90, 0.99] (test
  `test_variance_share_bootstrap_coverage_tissue`, xfail strict). Possibili strade, da
  calibrare prima di esporre: intervallo diverso per effetti fissi e casuali (es. bootstrap
  parametrico sui residui per il tessuto), oppure riportare solo la quota del paziente.
- **Coerenza con la trascrittomica spaziale (Visium).** Richiede Visium e scRNA-seq sugli
  stessi pazienti; in GSE278694 la sovrapposizione non e' nota.
- **Input da R/Seurat (`.rds`).** Da fare solo se il laboratorio partner lavora in R.
- **Esporre `gene_col`, `min_cells` (3) e `min_margin` (0.20) del Modulo B in CLI e app.**
  Oggi sono fissi ai default.
- **Matrice di confusione per paziente o per dimensione del clone (Intervento 3).** Oggi e'
  pooled fra pazienti; stratificarla per dimensione del clone permetterebbe di verificare,
  almeno fra i cloni condivisi, se l'errore cresce con l'espansione (rho = +0.177 in
  `pdac-ml`), invece di affidarsi solo agli scenari 0.5x/2x.
- **Controllo "unita'-tecnica" generalizzato.** Oggi si verifica solo la corrispondenza 1:1
  di un fattore tecnico con la coppia paziente-tessuto; si potrebbe estendere a qualunque
  combinazione di fattori biologici.
