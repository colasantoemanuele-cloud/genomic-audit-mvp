# Criteri di accettazione della revisione

Copiati dal prompt di revisione il 2026-10-18, prima di qualsiasi altra modifica. Non si allentano durante il loop.

1. **Calibrazione nel regime d'uso.** Ogni stimatore esposto all'utente (cluster bootstrap, Nadeau-Bengio, Cramér V, intervalli della frazione di CD8) ha un test di calibrazione con numero di gruppi, fold e dimensioni pari a quelli con cui lo strumento lo usa davvero, ed entra nella banda dichiarata. Altrimenti non è esposto, oppure è rifiutato sotto la soglia in cui non è calibrato.
2. **Sintesi dei controlli.** Il verde non compare mai per un controllo con un fatto strutturale, non eseguito, non valutabile, rifiutato o fuori dal regime calibrato. Ogni regola ha un test con oggetti risultato reali, non surrogati.
3. **Nessun errore silenzioso.** Chiavi duplicate (barcode, obs_names, coppie paziente-compartimento-libreria), input malformati e scelte d'ordine implicite producono un errore chiaro, mai un risultato diverso.
4. **Affermazioni proporzionate.** README, report, app e messaggi della CLI non usano «validato», «affidabile», «reale», «esatto», «dimostra» se non letteralmente vero, e dichiarano le ipotesi biologiche da cui dipendono i risultati.
5. **Tracciabilità.** Ogni numero del README è in un file di risultati prodotto dal codice attuale.
6. **Riproducibilità.** Gli esiti dei test nel README sono misurati su un clone pulito, con e senza i dati reali.
7. **Commenti.** Ogni funzione pubblica ha una docstring; nessun commento contraddice il codice.
8. **Privacy.** Nessuna copia dei dati del paziente resta su disco fuori dalla cartella di output scelta dall'utente.
9. **Test e stile.** La suite passa; `ruff` non segnala nulla; nessun test tautologico.
