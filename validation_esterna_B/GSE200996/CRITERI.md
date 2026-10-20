# Criteri della validazione tecnica esterna del Modulo B su GSE200996 (scritti PRIMA di qualsiasi download)

Secondo tentativo della voce 1 di `docs/ROADMAP.md`, dopo GSE139555 (non eleggibile per E1:
`validation_esterna_B/REPORT.md`). Dataset: GSE200996 (Luoma et al., *Cell* 2022: carcinoma
squamoso testa-collo, terapia neoadiuvante anti-PD-1 ± anti-CTLA-4; scRNA-seq + scTCR-seq di tumore
e sangue periferico). Versione dello strumento: 0.9.0, branch `validazione-b-gse200996` (creato da
`validazione-b`, perché master non ne contiene ancora il merge).

Informazioni preliminari, da una ricerca sulle pagine GEO e NON verificate: 18 pazienti con sangue e
tumore (P13, P14, P15, P16, P18, P19, P20, P21, P22, P23, P24, P25, P26, P28, P29, P30, P31, P32);
campioni di sangue raggruppati con hashtag (più pazienti nella stessa corsia, es. GSM6048242 «P18,
P29 and P30»); file di hashtag «sc_GEX_PBMC» e «sorted_sc_GEX_PBMC»; metadati per cellula degli
autori separati per CD4 e CD8, tumore e sangue. Da verificare: contenuto di GSE200996_RAW.tar
(1,5 GB), colonne dei metadati, demultiplazione degli hashtag, eventuale selezione (sorting) del
sangue per un tipo di cellule, chimica 10x.

Questo file è stato scritto e committato prima di scaricare qualunque file di GSE200996 e non si
modifica dopo il commit. Un criterio inapplicabile si dichiara «non valutabile» con il motivo. Se si
ritiene necessario cambiarne uno, `REPORT.md` contiene una sezione «Criteri modificati dopo aver
visto i dati», con testo prima e dopo, motivo ed esito con entrambe le versioni; il verdetto
principale resta quello con i criteri originali.

Che cosa può dimostrare questa validazione: che il Modulo B si esegue su dati di un altro
laboratorio, che i suoi calcoli sono corretti (implementazione indipendente), che non segnala
errori dove non ce ne sono (controllo negativo), che rileva errori iniettati (controllo positivo) e
che è stabile alle scelte dichiarate. Non può dimostrare che le etichette che segnala come errate lo
siano davvero: per quello serve un riferimento proteico (B-E5). Nel report l'esito va chiamato
«validazione tecnica esterna».

Regole: nessuna modifica a `core/`, `cli.py`, `app.py` durante la validazione (un difetto vero si
documenta e si corregge in un commit separato con un test di regressione, riportando l'esito prima e
dopo); nessuna soglia del Modulo B regolata su questo dataset (margine 0,20, minimo 3 cellule,
minimo 5 pazienti per gli intervalli, regola dei pareggi): si usano i default della 0.9.0. Errore
per compartimento e direzione degli errori sono risultati descrittivi, non criteri. Se GSE200996
non è eleggibile ci si ferma e si riporta, senza passare a un altro dataset.

## Definizioni fissate in anticipo

- Compartimento di riferimento: PBMC, unendo per ogni paziente le cellule di sangue di tutti i tempi
  disponibili (la linea CD4/CD8 di un clone non cambia nel tempo), limitate ai campioni ammessi da E5.
- Analisi principale: un h5ad che contiene SOLO PBMC e Tumor, dove Tumor = campioni di tumore dopo la
  terapia. Tutti i criteri B-E1…B-E6 si valutano su questo h5ad.
- Analisi descrittiva separata: un secondo h5ad con PBMC e i campioni di tumore prima della terapia,
  dove esistono, eseguito una volta con `audit-sc tcr`; nessun criterio si valuta su di esso.
- Clone: CDR3 nucleotidica della catena TRB (default dello strumento); variante TRA+TRB in B-E6.
- Cellule di sangue: solo quelle che i file di hashtag degli autori assegnano a un singolo paziente
  (doppietti e negativi degli hashtag esclusi, con il conteggio riportato).

## Condizioni di eleggibilità

- E1. Dopo la demultiplazione degli hashtag, almeno 5 pazienti hanno sia cellule di sangue sia
  cellule di tumore dopo la terapia, ciascuna con conteggi grezzi interi e TCR per cellula.
- E2. I metadati degli autori danno un'etichetta CD4 o CD8 per cellula in entrambi i compartimenti.
  Documenta da dove viene, citando i Metodi dell'articolo: clustering e marcatori sull'espressione
  (come nel caso studio), oppure separazione fisica (sorting) delle cellule. La mappa etichetta →
  {CD4, CD8, escluso} va scritta in MAPPA_ETICHETTE.md e committata PRIMA di eseguire il Modulo B,
  usando solo i nomi delle colonne e dei cluster e i Metodi, senza guardare nessun risultato del
  modulo. Le cellule escluse vengono tolte dall'h5ad prima di eseguire il modulo, con il conteggio
  per paziente.
- E3. La regola di abbinamento dei barcode fra metadati, matrici, file VDJ e hashtag è dichiarata
  nello script di conversione; almeno il 50% delle cellule T etichettate, per compartimento e
  sull'insieme dei campioni, ha un TCR produttivo abbinato (riporta la percentuale per campione).
- E4. L'archivio contiene, per i campioni usati, matrici di conteggi grezzi interi e annotazioni VDJ
  per cellula in un formato che il manifest del Modulo B accetta (filtered_contig_annotations o
  equivalente convertibile senza perdere informazione); la conversione è dichiarata.
- E5. Il sangue non è selezionato per una sola linea: dai Metodi e dai metadati risulta che le
  cellule di sangue usate contengono sia cellule T CD4 sia CD8. Se esistono sia campioni di sangue
  selezionati (sorted) sia non selezionati, si usano solo i non selezionati; se esistono solo
  selezionati per una linea, la condizione non è soddisfatta. Per decidere si usano solo i metadati
  e i Metodi, non i risultati del modulo.

## Criteri

- B-E1 (esecuzione). `audit-sc tcr` con i soli comandi e input documentati in docs/GUIDA_UTENTE.md
  (h5ad con conteggi grezzi, colonne paziente/compartimento/etichetta, manifest VDJ, riferimento =
  PBMC) termina con codice 0 e produce il report HTML, senza modifiche manuali al codice. Le etichette
  usano i nomi attesi dalla mappa dei marcatori di default, oppure la mappa è passata con l'opzione
  documentata: verificalo nella guida e nel codice prima di costruire l'h5ad.
- B-E2 (implementazione indipendente). indipendente.py, scritto leggendo SOLO docs/METODI.md §1-4 e
  senza importare nulla da core/, ricalcola: unità, pareggi e pazienti per compartimento; tasso
  d'errore per clone nel tumore; eccesso di discordanza pesato e a pesi uguali sangue-tumore; i loro
  intervalli t cluster-robusti. Superato se i conteggi coincidono esattamente, le stime puntuali entro
  1e-9 e gli estremi degli intervalli entro 1e-6. Se una differenza nasce da un punto che METODI.md
  non specifica, è una lacuna della documentazione: documentala, correggi METODI.md in un commit
  separato e riporta l'esito prima e dopo.
- B-E3 (controllo negativo). In 20 repliche (semi 0-19) si permuta, dentro ogni clone di ogni
  paziente, l'etichetta di compartimento fra le cellule del clone, lasciando invariati i numeri di
  cellule per compartimento e le etichette di linea. Superato se la media delle 20 stime dell'eccesso
  pesato sangue-tumore è in [−0,01; +0,01] e l'intervallo contiene lo zero in almeno 18 repliche su 20.
- B-E4 (controllo positivo con errori iniettati). Fra le unità clone-tumore con identità di riferimento
  CD4 e maggioranza corretta, se ne scelgono a caso il 10% e poi il 20% (seme dichiarato) e si
  rietichettano CD8 tutte le loro cellule nel tumore. Superato se, in entrambe le dosi, il numero di
  unità in errore nel tumore aumenta esattamente del numero di unità iniettate (pareggi riportati a
  parte) e l'eccesso a pesi uguali sangue-tumore aumenta; l'eccesso pesato è solo descrittivo. Se le
  unità con identità CD4 e maggioranza corretta sono meno di 20, B-E4 è «non valutabile».
- B-E5 (riferimento proteico). Non valutabile su questo dataset salvo che contenga proteine di
  superficie per cellula: verificalo e dichiaralo.
- B-E6 (sensibilità). Il tasso d'errore per clone nel tumore cambia al massimo di 0,05 in valore
  assoluto passando dal clone TRB al clone TRA+TRB, e dalla regola dei pareggi di default a «pareggi
  = errore» e «pareggi = corretti». Ogni variante riporta pazienti, unità e pareggi, e nessun
  intervallo sotto 5 pazienti.
