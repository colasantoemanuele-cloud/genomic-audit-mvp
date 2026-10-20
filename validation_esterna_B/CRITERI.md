# Criteri della validazione esterna del Modulo B su GSE139555 (scritti PRIMA di qualsiasi download)

Voce 1 di `docs/ROADMAP.md`. Dataset: GSE139555 (Wu et al., *Nature* 2020: scRNA-seq + scTCR-seq di
tumore, tessuto adiacente normale e sangue periferico in più tipi di tumore). Versione dello
strumento: 0.9.0 (master = cb9c1ec, tag v0.9.0), branch `validazione-b`.

Questo file è stato scritto e committato prima di scaricare qualunque file di GSE139555 e non si
modifica dopo il commit. Se un criterio si rivela inapplicabile, viene dichiarato «non valutabile»
con il motivo nel report. Se si ritiene necessario cambiarne uno, `REPORT.md` contiene una sezione
«Criteri modificati dopo aver visto i dati», con testo prima e dopo, motivo ed esito con entrambe le
versioni; il verdetto principale resta quello con i criteri originali.

Che cosa può dimostrare questa validazione: che il Modulo B si esegue su dati di un altro
laboratorio, che i suoi calcoli sono corretti (implementazione indipendente), che non segnala
errori dove non ce ne sono (controllo negativo), che rileva errori iniettati (controllo positivo) e
che è stabile alle scelte dichiarate. Non può dimostrare che le etichette che segnala come errate lo
siano davvero: per quello serve un riferimento proteico (B-E5), assente in GSE139555. Nel report
l'esito va chiamato «validazione tecnica esterna», non «validazione» senza aggettivi.

Regole: nessuna modifica a `core/`, `cli.py`, `app.py` durante la validazione (un difetto vero si
documenta e si corregge in un commit separato con un test di regressione, riportando l'esito prima e
dopo la correzione); nessuna soglia del Modulo B regolata su questo dataset (margine 0,20, minimo 3
cellule, minimo 5 pazienti per gli intervalli, regola dei pareggi): si usano i default della 0.9.0.
L'errore per compartimento e la direzione degli errori sono risultati descrittivi, non criteri.

## Condizioni di eleggibilità (se una non è soddisfatta: fermati e riporta, senza eseguire i criteri)

- E1. Per almeno 5 pazienti esistono, sia per il sangue sia per il tumore, la matrice dei conteggi
  grezzi interi e il file VDJ filtrato (filtered_contig_annotations). Il tessuto adiacente, dove
  presente, si usa ma non è richiesto.
- E2. I metadati degli autori danno per le cellule T un'etichetta di cluster da cui si può derivare la
  linea CD4 o CD8. La mappa cluster → {CD4 o CD8, escluso} va scritta in
  validation_esterna_B/MAPPA_ETICHETTE.md e committata PRIMA di eseguire il Modulo B, usando solo i
  nomi dei cluster e la legenda dell'articolo, senza guardare nessun risultato del modulo. Cluster
  ambigui (misti, ciclanti senza linea, MAIT, γδ, ad alta frazione mitocondriale senza linea) →
  escluso. Le cellule dei cluster esclusi vengono tolte dall'h5ad prima di eseguire il modulo (non
  rietichettate «altro»: la convenzione per clone conterebbe un'etichetta fuori mappa come errore);
  riporta quante cellule e quanti cloni si perdono per paziente.
- E3. La regola di abbinamento dei barcode fra metadati, matrici e file VDJ è dichiarata nello script
  di conversione; almeno il 50% delle cellule T etichettate dei metadati, sull'insieme dei campioni,
  ha un TCR produttivo abbinato (riporta la percentuale per campione).

## Criteri

- B-E1 (esecuzione). `audit-sc tcr` con i soli comandi e input documentati in docs/GUIDA_UTENTE.md
  (h5ad con conteggi grezzi, colonne paziente/compartimento/etichetta, manifest VDJ, riferimento =
  sangue) termina con codice 0 e produce il report HTML, senza modifiche manuali al codice. Le
  etichette CD4/CD8 usano i nomi che la mappa dei marcatori di default si aspetta, oppure la mappa è
  passata con l'opzione documentata: verificalo nella guida e nel codice prima di costruire l'h5ad.
- B-E2 (implementazione indipendente). validation_esterna_B/indipendente.py, scritto leggendo SOLO
  docs/METODI.md §1-4 e senza importare nulla da core/ (scritto da un sottoagente a cui è vietato
  leggere core/), ricalcola: unità, pareggi e pazienti per compartimento; tasso d'errore per clone
  nel tumore e nell'adiacente; eccesso di discordanza pesato e a pesi uguali per la coppia
  sangue-tumore; i loro intervalli t cluster-robusti (CR2, gradi di libertà di Bell-McCaffrey, scala
  logit per i tassi). Superato se i conteggi coincidono esattamente, le stime puntuali entro 1e-9 e
  gli estremi degli intervalli entro 1e-6. Se una differenza nasce da un punto che METODI.md non
  specifica, è una lacuna della documentazione: documentala, correggi METODI.md in un commit separato
  e riporta l'esito prima e dopo la correzione.
- B-E3 (controllo negativo). In 20 repliche (semi 0-19) si permuta, dentro ogni clone di ogni
  paziente, l'etichetta di compartimento fra le cellule del clone, lasciando invariati i numeri di
  cellule per compartimento e le etichette di linea. Sotto questa permutazione l'eccesso pesato ha
  valore atteso zero. Superato se la media delle 20 stime dell'eccesso pesato sangue-tumore è in
  [−0,01; +0,01] e l'intervallo contiene lo zero in almeno 18 repliche su 20.
- B-E4 (controllo positivo con errori iniettati). Fra le unità clone-tumore con identità di riferimento
  CD4 e maggioranza corretta, se ne scelgono a caso il 10% e poi il 20% (seme dichiarato) e si
  rietichettano CD8 tutte le loro cellule nel tumore. Superato se, in entrambe le dosi: il numero di
  unità in errore nel tumore aumenta esattamente del numero di unità iniettate (pareggi riportati a
  parte) e l'eccesso a pesi uguali sangue-tumore aumenta rispetto ai dati originali. L'eccesso pesato
  si riporta come descrittivo, senza soglia (METODI §3: può diminuire quando le proporzioni divergono).
- B-E5 (riferimento proteico, CITE-seq). Non applicabile a GSE139555, che non ha proteine di superficie:
  dichiarare «non valutabile su questo dataset», rimandato alla voce 2 della roadmap.
- B-E6 (sensibilità alle scelte). Il tasso d'errore per clone nel tumore cambia al massimo di 0,05 in
  valore assoluto passando dal clone definito dalla sola TRB al clone TRA+TRB, e passando dalla regola
  dei pareggi di default a «pareggi = errore» e «pareggi = corretti». Ogni variante riporta pazienti,
  unità e pareggi, e nessun intervallo sotto 5 pazienti.
