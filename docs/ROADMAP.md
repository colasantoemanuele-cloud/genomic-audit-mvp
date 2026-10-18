# Roadmap verso la versione 1.0

La versione 0.9.0 è un **prototipo per studio pilota**. Le voci qui sotto erano previste per la
1.0 (Fase 4 della revisione) e sono **rimandate, non eseguite**. Le stime sono di lavoro effettivo
su una workstation come quella usata per la revisione (i tempi di calcolo sono indicati a parte
quando dominano).

| # | Voce | Che cosa serve | Stima |
|---|---|---|---|
| 1 | **Validazione esterna del Modulo B** su GSE139555 (Wu et al., *Nature* 2020: scRNA-seq + TCR di sangue, tumore e tessuto adiacente in più tipi di tumore) | Scrivere e committare `validation_esterna_B/CRITERI.md` (criteri B-E1…B-E6) **prima di qualsiasi download**, senza modificarli dopo; scaricare i dati; script `indipendente.py` che non importa nulla da `core/`; eseguire e riportare l'esito di ogni criterio | 6-10 h, di cui 1-2 h di download e conversione |
| 2 | **Criterio B-E5 su un dataset CITE-seq + TCR** (proteine di superficie CD4/CD8 come riferimento indipendente dai trascritti) | Individuare un dataset pubblico adatto; se non esiste, dichiarare il criterio non valutabile | 3-5 h |
| 3 | **Pacchetto `genomic_audit`** | Spostare `core/`, `cli.py` e `app.py` in un pacchetto importabile, mantenere il comando `audit-sc`, includere o incorporare `.streamlit/config.toml`, aggiornare test e documenti | 2-3 h |
| 4 | **Installazione senza rete** | Pacchetto di wheel per Linux, script di installazione, prova su una macchina senza connessione | 1-2 h |
| 5 | **Dockerfile** | Immagine che avvia l'app solo su localhost; prova di avvio; requisiti minimi misurati (RAM, disco, tempo) | 1-2 h |
| 6 | **Integrazione continua** | GitHub Actions su Ubuntu, Windows e macOS con Python 3.12 e 3.13; verifica al primo push | 1 h, più il tempo di correzione dei problemi specifici di Windows/macOS (non stimabile prima della prima esecuzione) |
| 7 | **Scalabilità** | Dataset sintetico di circa 300.000 cellule × 30.000 geni × 30 pazienti; tempi e memoria con `/usr/bin/time -v` per ogni modulo | 2-3 h, di cui 1-2 h di calcolo |
| 8 | **Documentazione 1.0** | `docs/VALIDAZIONE.md` (sintesi di tutte le verifiche), revisione di `GUIDA_UTENTE.md` e `METODI.md`, nota sulla protezione dei dati, report HTML di esempio su dati reali in `docs/esempi/` | 2-3 h |
| 9 | **Rilascio 1.0.0** | Versione 1.0.0, `CHANGELOG.md`, `CITATION.cff`, tag `v1.0.0` | 0,5 h |
| 10 | **Verdetto 1.0** | Un revisore indipendente scrive `revisione/VERDETTO.md` sulla lista di consegna (11 voci) e sui criteri B-E; «pronto» solo se tutte le voci sono soddisfatte | 1-2 h |
| 11 | **Giro di chiusura dopo le modifiche al codice delle voci 1-7** | Revisione completa, suite completa con i test lenti, riesecuzione delle validazioni | 2-3 h |

Totale stimato: **circa 22-35 h** di lavoro.

## Limiti della 0.9 che la 1.0 dovrebbe affrontare

- Il Modulo B è stato eseguito su dati reali solo su GSE278694, il dataset su cui è nato: è una
  verifica di fedeltà del porting, non una validazione (voci 1-2).
- Il confronto fra modelli del Modulo A è descrittivo: Nadeau-Bengio e Wilcoxon sono fuori banda
  in LeaveOneGroupOut e non sono esposti (`revisione/log/a2_calibrazione.txt`).
- Nessun utente esterno ha provato lo strumento; la prova della guida (C2, GSE144735) è stata
  fatta da chi lo ha sviluppato.
- Le prove di installazione coprono Linux con Python 3.12 e 3.13; Windows e macOS non sono provati
  (voce 6).
