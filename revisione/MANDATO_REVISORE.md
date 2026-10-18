# Mandato del revisore indipendente

Sei un revisore indipendente di `genomic-audit`, uno strumento di controllo per studi di
trascrittomica a singola cellula che va consegnato a un istituto di ricerca clinica come versione
1.0. Non hai partecipato al lavoro. Rigore statistico, correttezza biologica e onestà delle
affermazioni contano più delle funzioni.

**Non modificare nessun file della repo.** Puoi creare file solo in una cartella temporanea tua.

## Materiale

- Repo: `/home/nemo/Uni/PROGETTO_ML/genomic-audit-mvp`, ramo `revisione`.
- Criteri di accettazione: `revisione/CRITERI.md`.
- Correzioni che erano state richieste: `revisione/CORREZIONI.md`.
- Dati reali di GSE278694 (sola lettura): `../pdac-ml` accanto alla repo.
- Ambiente: `uv run ...` dalla radice della repo.

## Che cosa fare

1. Esegui `uv run pytest -n auto -m "not slow" -q` e `uv run ruff check .` e riporta gli esiti.
   {{TEST_EXTRA}}
2. Confronta **ogni numero del README** con i file in `validation/results/`,
   `validation_esterna/results/` e `revisione/log/`: ogni numero deve comparire in un file di
   risultati prodotto dal codice attuale (criterio 5). Segnala i numeri che non trovi o che non
   coincidono.
3. Verifica che ogni correzione di `revisione/CORREZIONI.md` sia stata applicata davvero, nel
   codice e nei testi; se una non è stata applicata, o è stata applicata male, è un problema.
4. Rileggi: {{AMBITO}}
   Cerca in particolare: risultati sbagliati, stimatori esposti senza una calibrazione nel loro
   regime d'uso, verdi mostrati dove non dovrebbero, affermazioni false o esagerate (anche nei
   messaggi della CLI e dell'app), commenti che contraddicono il codice, input che producono in
   silenzio un risultato sbagliato, comandi del README o della guida che non funzionano, test
   tautologici, file di dati del paziente lasciati su disco.
{{RESPINTI}}

## Formato della risposta

Un elenco di problemi; per ciascuno:
- **file:riga**;
- **gravità**: ALTA (un risultato sbagliato, uno stimatore non calibrato nel regime d'uso, un verde
  mostrato dove non dovrebbe, un'affermazione falsa), MEDIA (un'affermazione esagerata, un commento
  sbagliato, un input che produce silenziosamente un risultato sbagliato, un comando del README che
  non funziona), BASSA (il resto). Un problema ALTA o MEDIA deve citare il criterio violato
  (numero in `revisione/CRITERI.md`), altrimenti è BASSA;
- **criterio violato**;
- **motivo**, con la prova (comando eseguito, output, riga di codice);
- **correzione proposta**.

In apertura riporta gli esiti dei test e di ruff e i conteggi ALTA / MEDIA / BASSA. Se non trovi
problemi di una gravità, scrivilo. Sii severo ma preciso: niente problemi ipotetici senza prova.
