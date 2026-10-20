# Mappa delle etichette CD4/CD8 per GSE200996 (scritta PRIMA di eseguire il Modulo B)

Condizione E2 di `CRITERI.md`. Scritta usando solo i nomi dei file e delle colonne dei metadati degli
autori e i Metodi dell'articolo (Luoma et al., *Cell* 2022, PMC9508682), senza aver eseguito il
Modulo B né scaricato le matrici. Conteggi in `results/metadati.json` (`verifica_metadati.py`).

## Da dove viene l'etichetta (Metodi dell'articolo)

Gli autori non danno un'etichetta CD4/CD8 in una colonna: la danno come **appartenenza a due file
separati per compartimento**, `GSE200996_CD4.<PBMC|tumor>.single.cell.meta.data.txt.gz` e
`GSE200996_CD8.<PBMC|tumor>.single.cell.meta.data.txt.gz` (colonne `Patient_ID`, `Stage`, `Cohort`,
`Path_response`, `CellType_ID`, `UMAP_1`, `UMAP_2`; `CellType_ID` è il numero del sottocluster).

Metodi, «Annotation of T cells» (testo citato): le cellule T «were mapped to previously annotated
CITE-seq reference datasets in Seurat» (Hao et al., 2021: 162.000 PBMC con 228 anticorpi CITE-seq);
«Only 'CD4 T', 'CD8 T' and 'Other T' cell populations in the reference dataset were selected for the
prediction» (`FindTransferAnchors`, `MapQuery`); «CD4 T and CD8 T cells were selected based on a
prediction score >0.6. A small number of CD8A and CD8B positive cells in the CD4 T cell group and CD4
positive cells in CD8 T cells were also removed.» Blood: 33.446 CD4 e 14.207 CD8; tumor: 24.032 CD4 e
18.126 CD8 (gli stessi numeri delle righe dei file).

Quindi l'etichetta è **un'annotazione computazionale sull'espressione** (trasferimento da un
riferimento, non clustering e marcatori come nel caso studio), **non una separazione fisica**
(nessun sorting per linea nei campioni usati: vedi E5). Conseguenza da riportare, non un criterio:
gli autori hanno tolto le cellule etichettate CD4 con CD8A o CD8B rilevati e le cellule etichettate
CD8 con CD4 rilevato. I geni marcatori della mappa del Modulo B (CD4; CD8A, CD8B) sono quindi già
stati usati per ripulire le etichette: per costruzione nessuna cellula CD4T è positiva a CD8A/CD8B e
nessuna cellula CD8T è positiva a CD4. L'identità di riferimento del clone e le etichette non sono
indipendenti, e gli «errori» che il modulo può trovare sono soprattutto disaccordi di etichetta fra
le cellule dello stesso clone, non contraddizioni fra etichetta e marcatori della stessa cellula.

## Mappa

| Origine della cellula (file degli autori) | Etichetta nell'h5ad (`label`) |
|---|---|
| `GSE200996_CD4.PBMC...` e `GSE200996_CD4.tumor...` | `CD4T` |
| `GSE200996_CD8.PBMC...` e `GSE200996_CD8.tumor...` | `CD8T` |
| cellula presente in entrambi i file (CD4 e CD8) | escluso (0 cellule nei metadati) |
| cellula del file CD45 assente da CD4 e da CD8 (T non assegnate dagli autori, NK, B, mieloidi, mastociti) | escluso |

Le cellule escluse non entrano nell'h5ad (non sono rietichettate «altro»). Cellule T del file CD45
senza etichetta CD4/CD8: 6.801 nel sangue e 6.309 nel tumore (tutti i pazienti e i tempi; il
conteggio per paziente dei campioni usati è stampato da `build_gse200996.py`).

**Annotazioni discordanti fra i due livelli degli autori.** Alcune cellule dei file CD4/CD8 hanno nel
file CD45 un tipo diverso da «T cell» o «Cycling»: nel sangue 135 «NK cell»; nel tumore 294 «NK
cell», 52 «B cell», 50 «Myeloid», 3 «Mast». Decisione presa ora, prima di eseguire il modulo:
**restano con l'etichetta CD4T/CD8T**, perché l'etichetta valutata dal Modulo B è quella di linea del
file CD4/CD8 (l'annotazione fine che gli autori usano per tutte le analisi delle cellule T), non il
tipo grossolano del file CD45. Il loro numero è riportato nel report.

Mappa dei marcatori passata con `--marker-map` (la CLI non ha una mappa di default; è l'esempio di
`docs/GUIDA_UTENTE.md` §3.3 e della demo): `{"CD4T": ["CD4"], "CD8T": ["CD8A", "CD8B"]}`.

## Altre scelte fissate qui

- Sangue: solo le corsie non selezionate (GSM6048235-GSM6048246, CD45+ vive); le corsie selezionate
  CD38+ HLA-DR+ CD4/CD8 (GSM6048338-GSM6048345) non sono usate (E5). Tutti i prelievi B1, B2, B3 di
  un paziente formano il suo PBMC.
- Ogni cellula di sangue è assegnata al paziente e al prelievo dall'ID degli autori
  (`<barcode>_<paziente>_<prelievo>`), che è l'esito della loro demultiplazione degli hashtag
  (GMM-demux, singoletti con confidenza ≥ 0,8: Metodi). La corsia è quella della chiave
  `GSE200996_hashtag_info_sc_GEX_PBMC.csv.gz` che contiene quella coppia paziente-prelievo (ogni
  paziente è in una sola corsia). I barcode della corsia che non sono nei metadati (doppietti,
  negativi, cellule tolte dal controllo di qualità degli autori) sono esclusi e contati.
- Tumore: `Stage` = `Post-Tx` per l'h5ad principale, `Pre-Tx` per l'h5ad descrittivo.
