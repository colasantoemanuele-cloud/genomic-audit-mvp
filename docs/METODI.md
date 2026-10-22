# Metodi

Documento per un revisore scientifico: per ogni stimatore esposto, la formula, le ipotesi, il
riferimento, la calibrazione (banda dichiarata prima dell'esecuzione e valore misurato) e i limiti.
I valori misurati sono nei file citati; le bande sono nei docstring dei test.

> Strumento per la ricerca. Non è un dispositivo medico e non va usato per decisioni cliniche su
> singoli pazienti.

## 1. Principio comune: il paziente è l'unità

Cellule, cloni e coppie clone-compartimento dello stesso paziente non sono osservazioni
indipendenti. Ogni intervallo è quindi calcolato sui pazienti, e ogni conteggio di «unità
indipendenti» conta pazienti. Sotto 5 pazienti lo strumento non produce intervalli
(`MIN_PATIENTS_FOR_CI_DEFAULT = 5`), perché la calibrazione parte da 5.

## 2. Intervallo t cluster-robusto (Modulo B)

**Dove.** Ogni intervallo del Modulo B: eccesso di discordanza (anche a pesi uguali e per coppia di
compartimenti), tasso d'errore per compartimento in entrambe le convenzioni, differenze appaiate
fra compartimenti (`core/stats.py::cluster_t_interval`).

**Stima puntuale.** Media di tutte le osservazioni, cioè media dei pazienti pesata per numero di
osservazioni: θ̂ = Σ_g S_g / N, con S_g la somma dei valori del paziente g, n_g il numero di
osservazioni, N = Σ n_g.

**Varianza.** Varianza cluster-robusta CR2 della regressione sulla sola intercetta (Bell e
McCaffrey 2002): V̂ = Σ_g e_g² / ((1 − n_g/N) N²), con e_g = S_g − θ̂ n_g la somma dei residui del
paziente g. Il fattore 1/(1 − n_g/N) è la correzione di leva del CR2 per questo modello.

**Gradi di libertà.** Bell-McCaffrey con modello di lavoro a errori indipendenti e
omoschedastici: ν = (tr M)² / tr(M²), M = A^½ K A^½, K = diag(n) − n nᵀ/N,
A = diag(1 / ((1 − n_g/N) N²)). Con pazienti di uguale dimensione ν = G − 1; con dimensioni diverse
ν < G − 1 (Imbens e Kolesár 2016).

**Intervallo.** θ̂ ± t_{ν, 0.975} √V̂. Per i tassi (valori in [0, 1]) l'intervallo è calcolato
per logit(θ̂) con il metodo delta (errore standard √V̂ / (θ̂ (1 − θ̂))) e riportato sulla scala
originale; con θ̂ = 0 o 1, o con variabilità nulla fra pazienti, nessun intervallo.

**Ipotesi.** Pazienti indipendenti; nessuna ipotesi sulla struttura entro paziente (cloni, cellule
dello stesso clone). L'approssimazione t peggiora con distribuzioni molto asimmetriche fra
pazienti: la scala logit serve a questo per i tassi.

**Calibrazione** (`tests/test_stats_calibration.py`, banda: frazione di intervalli al 95% che non
contengono il valore vero sotto H0 fra 2% e 8%; 1000 repliche per caso; 5, 8, 10 e 15 pazienti con
dimensioni diverse). Scenari: normale (3-14 osservazioni per paziente), con eterogeneità fra
pazienti, tasso per cellula modellato su GSE278694 (40-1000 cellule, tasso per paziente con media
0,09 e d.s. 0,06, cellule raggruppate in cloni), tasso per clone (5-60 unità, media 0,18, d.s.
0,125), eccesso di discordanza (13-285 coppie, metà zeri). Misurato: 3,6%-6,3% in tutti i 20 casi
(`revisione/log/a1_calibrazione.txt`). Potenza con effetto 0,5 d.s. e 10 pazienti: 0,97.

**Perché non il bootstrap percentile.** Il cluster bootstrap percentile (Davison e Hinkley 1997),
usato prima della revisione, è troppo stretto con pochi gruppi: 14,5% di intervalli che mancano il
valore vero con 5 pazienti nello scenario normale (test di calibrazione); nelle simulazioni
esplorative che hanno guidato la scelta del metodo, 18% con 5 pazienti e 9% con 10, e BCa 17% con 5
(`revisione/REGISTRO.md`, Fase 1). Il percentile resta nel codice
(`cluster_bootstrap`) solo per riprodurre i numeri del caso studio e per la frazione di CD8
sperimentale.

## 3. Eccesso di discordanza (Modulo B)

Per ogni paziente, clone T (CDR3 nucleotidica della TRB) e coppia di compartimenti A, B con almeno
2 cellule per lato:
- D_A = 1 − Σ_k n_k (n_k − 1) / (n (n − 1)): probabilità che due cellule DISTINTE del clone in A
  abbiano etichetta diversa (senza reinserimento, non distorto);
- D(A×B) = 1 − Σ_k p_{A,k} p_{B,k}: probabilità che una cellula di A e una di B abbiano etichetta
  diversa;
- eccesso pesato (stima principale): D(A×B) − (n_A D_A + n_B D_B)/(n_A + n_B);
- eccesso a pesi uguali (sensibilità): D(A×B) − (D_A + D_B)/2.

Valori attesi. L'eccesso a pesi uguali ha valore atteso ½‖p_A − p_B‖²: zero se le proporzioni di
etichette coincidono, positivo altrimenti. L'eccesso pesato ha valore atteso
w_A‖p_A‖² + w_B‖p_B‖² − p_A·p_B, zero sotto H0 ma negativo quando le proporzioni differiscono e
n_A ≠ n_B (esempio: 3 cellule nel sangue, tutte corrette, e 60 nel tumore con il 40% di errore:
circa −0,06 se le cellule sbagliate hanno tutte la stessa etichetta, circa −0,13 se si dividono in
parti uguali fra due etichette, come nel caso simulato di `tests/test_input_errors.py`, −0,131). Un intervallo pesato interamente sotto lo zero non
indica quindi assenza di errori, e la sintesi lo mostra in giallo.

**Ipotesi biologica.** Lo stesso clone ha la stessa linea CD4/CD8. Eccezioni: CD4 con CD8αα,
MAIT, NKT e γδ, doppietti T-NK, RNA ambientale specifico del compartimento. Le CD4 citotossiche, che esprimono programmi tipici delle CD8, NON sono un'eccezione: restano della linea CD4, e una loro etichetta CD8 è un errore di linea che il modulo conta correttamente; nel caso studio (progetto-ml) sono indicate come possibile causa di una parte delle chiamate CD8 sui cloni CD4 nel tumore di GSE278694. Solo
etichette di linea, non di stato.

## 4. Tasso d'errore via marcatori (Modulo B)

**Identità di riferimento del clone.** Solo dal compartimento di riferimento (di solito il sangue):
almeno 3 cellule del clone, frazione di cellule positive ai marcatori di ciascuna etichetta (conta
> 0 su conteggi grezzi interi), identità = etichetta con la frazione più alta se il margine sulla
seconda è almeno 0,20; altrimenti il clone non ha riferimento.

**Convenzione per cellula.** Unità = cellula di un clone con riferimento, fuori dal compartimento di
riferimento, con etichetta fra le chiavi della mappa dei marcatori; errore se l'etichetta è diversa
dall'identità.

**Convenzione per clone** (definizione del caso studio). Unità = coppia clone-compartimento con almeno 3
cellule; errore se l'etichetta di maggioranza è diversa dall'identità; le etichette fuori mappa
ammesse (default: tutte) contano come errore; con `--clone-marker-priority` la positività è
esclusiva in quell'ordine. Pareggi della maggioranza (due o più etichette con lo stesso conteggio
massimo): l'etichetta di maggioranza non è definita e l'unità è esclusa dal tasso, come nel caso studio
(progetto-ml, `10_loco.py`); il numero di pareggi è riportato accanto a ogni tasso. Come
sensibilità, `--tie-rule error` li conta come errore e `--tie-rule correct` come corretti: le due
sensibilità usano le stesse unità più i pareggi. La regola non dipende dall'ordine delle righe.

**Scelta del contig.** Se una cellula ha più contig della stessa catena, vince quello con più UMI,
poi quello con più reads, poi la CDR3 minore in ordine alfabetico. Il caso studio teneva il primo dopo un
ordinamento per UMI non stabile: su GSE278694 155 cellule hanno una CDR3 TRB diversa e 879 una CDR3
TRA diversa. Con la regola di esclusione dei pareggi il tasso per clone (TRB) coincide comunque con
quello del caso studio; la sola unità in più (PA12, tumore) è un pareggio
(`validation/results/a1_moduleB.json`, sezione `fedelta_S5`).

Intervalli: sezione 2, su scala logit. Calibrazione: scenari «cellula» e «clone» della sezione 2.

**Limiti delle sezioni 3-4 su dati reali.** Le formule sono calibrate su dati simulati (sezione 2),
e su dati reali il Modulo B è stato eseguito su GSE278694 (il dataset su cui è nato: verifica di
fedeltà) e su GSE200996 (validazione tecnica esterna: un'implementazione scritta solo da queste
sezioni dà gli stessi conteggi e le stesse stime, entro 1e-15; controllo negativo e sensibilità
superati; `validation_esterna_B/GSE200996/REPORT.md`). Il controllo positivo con errori iniettati
non è stato valutabile su dati reali (su GSE200996 le unità con identità CD4 erano 3), e nessun
riferimento proteico indipendente dai trascritti è stato ancora usato per verificare le etichette
segnalate come errate. GSE139555 non era eleggibile (`validation_esterna_B/REPORT.md`).

## 5. Modulo A

**Schemi.** Split per paziente: `StratifiedGroupKFold` a min(5, pazienti) fold, con shuffle e seme
0. Split casuale (controllo negativo): `StratifiedKFold` a min(5, cellule della classe più rara)
fold. Pipeline stimata sul solo training di ogni fold: CPM a 1e4 e log1p, 2000 geni a varianza più
alta, standardizzazione senza centratura; regressione logistica (C = 1, pesi bilanciati) sui geni;
random forest e gradient boosting su 50 componenti SVD. Training sottocampionato a circa 20.000
cellule per fold, con quote proporzionali per paziente (arrotondate per difetto, minimo 1).
Macro-F1 sulle sole classi presenti nel fold di test.

**Divario.** Media dello split casuale meno media dello split per paziente. Sotto ±0,05 è
considerato rumore e non se ne indica la direzione; sopra +0,05 la sintesi lo segnala in rosso,
sotto −0,05 in giallo (direzione inattesa: variabilità fra fold, classi assenti o pochi pazienti).
Con meno di 5 pazienti nessuna riga del Modulo A è verde: la soglia è stata provata solo con 19
pazienti (GSE125449). La soglia è una scelta dichiarata, non una stima.

**Nessun intervallo per la macro-F1.** I punteggi dei fold sono correlati (training sovrapposti): la
deviazione standard fra fold non è una stima valida dell'incertezza, e il rapporto fra le due
deviazioni standard è riportato come descrittivo.

**Confronto fra modelli.** Con almeno 8 pazienti, `LeaveOneGroupOut`; ogni modello contro la
regressione logistica, fissata come riferimento prima di vedere i dati. Esposto: differenza delle
medie e numero di pazienti in cui il modello fa meglio, peggio o uguale. **Nessun p-value esposto.**

**Calibrazione dei test nel regime d'uso** (`tests/test_stats_calibration.py`; banda di falsi
positivi sotto H0 fra 2% e 8%; 400 repliche; 8, 10, 15, 20 pazienti con 15-60 cellule; due
regressioni logistiche su metà diverse e scambiabili delle feature, quindi ugualmente buone in
aspettativa e addestrate sugli stessi training):

| Test | 8 | 10 | 15 | 20 | Esito |
|---|---|---|---|---|---|
| Nadeau e Bengio (2003), rapporto fra le medie n_test/n_train | 0,018 | 0,007 | 0,007 | 0,005 | fuori banda (conservativo): `xfail(strict=True)`, non esposto |
| Wilcoxon appaiato | 0,060 | 0,083 | 0,072 | 0,045 | fuori banda a 10: `xfail(strict=True)`, non esposto |
| t appaiato non corretto (solo contrasto) | 0,065 | 0,077 | 0,080 | 0,043 | in banda, non usato |

Valori in `revisione/log/a2_calibrazione.txt`. In LeaveOneGroupOut le differenze fra fold sono
quasi scorrelate, mentre la correzione di Nadeau-Bengio assume una correlazione n_test/n_train e
quasi raddoppia la varianza. Il metodo di Wilcoxon è scelto esplicitamente (esatto senza pareggi;
permutazione esatta dei segni con pareggi fino a 13 coppie; approssimazione normale oltre).

**Stabilità delle spiegazioni.** Jaccard fra i 50 geni con coefficiente più grande nei fold:
descrittivo, senza soglia.

## 6. Audit del disegno

**Fatti strutturali** (deterministici): costante, identificativo di riga, annidamento (ogni livello
di A in un solo livello di B), uno-a-uno, esito determinato da un fattore, unità tecnica (fattore
tecnico 1:1 con la coppia paziente-tessuto). Sono controlli sulla tabella, non stime.

**Cramér V** con la correzione del bias di Bergsma (2013):
φ̃² = max(0, χ²/n − (k−1)(r−1)/(n−1)), r̃ = r − (r−1)²/(n−1), k̃ = k − (k−1)²/(n−1),
V = √(φ̃² / min(k̃−1, r̃−1)). Valutabile solo con almeno 20 righe e in media almeno 2 righe per cella
(n/(r c) ≥ 2); allarme a V ≥ 0,5 se la coppia non ha un fatto strutturale. La soglia di 20 righe
è il foglio più piccolo calibrato: con 10-12 righe i falsi allarmi sotto fattori indipendenti
superano il 5% (0,096 a 10 righe, 2×2, simulazione del revisore indipendente, giro 1). Un
confronto «stimabile» è segnalato in giallo se un fattore tecnico è fortemente associato ai due
livelli o se quell'associazione non è valutabile.

**Calibrazione** (`tests/test_design_audit_calibration.py`, `revisione/log/design_calibrazione.txt`):
- fattori indipendenti, 80 righe: almeno un allarme in 1/200 repliche (banda 0-5%);
- fattori indipendenti, fogli da 20, 30 e 40 righe con 2 o 3 livelli: allarmi per coppia valutabile
  0,0%-0,7% (banda 0-5%, 400 repliche per caso);
- annidamento e confondimento paziente-esito iniettati: rilevati 200/200;
- potenza (3×3, 60 righe, una sola configurazione): 0,02 a V = 0,3; 0,46 a V = 0,5; 0,98 a V = 0,7.

**Confronti.** Appaiato se almeno un paziente ha entrambi i livelli (unità = pazienti con entrambi;
p-value minimo di Wilcoxon esatto 2/2^k); altrimenti fra pazienti solo per un esito a livello di
paziente (p-value minimo di Mann-Whitney esatto 2/C(n_a+n_b, n_a)), non stimabile per tessuto e
fattori tecnici. Bassa potenza sotto 5 unità o con p-value minimo sopra 0,05 (con 5 coppie la
potenza a 0,05 è nulla). Un fattore tecnico che separa perfettamente i due livelli rende il
confronto non stimabile.

**Scomposizione della varianza.** Implementata (metodo dei momenti su pseudobulk) ma non esposta:
la copertura dell'intervallo della quota del tessuto è 1,000, fuori dalla banda [0,90; 0,99]
(`xfail(strict=True)`).

## 7. Frazione di CD8 (sperimentale)

Inversione di una matrice di confusione 2×3 (vera CD4/CD8 → chiamata CD4/CD8/altro) stimata sulle
cellule con riferimento e aggregata fra pazienti; ricampionamento dei pazienti (schema del bootstrap
percentile) e Beta(n_CD8 + ½, n_CD4 + ½) per i conteggi; tre scenari di errore (0.5x, 1x, 2x).
Rifiuto se J < 0,2 o il suo intervallo include lo zero, e se la frazione osservata cade fuori
dall'intervallo invertibile [b4/(a4+b4), b8/(a8+b8)] (prima della revisione veniva troncata).
Calibrazione nella banda [0,90; 0,99] solo quando l'errore delle cellule con riferimento è
rappresentativo degli altri; con errore non rappresentativo nessuno scenario garantisce la
copertura (`tests/test_cd8_robustness.py`). Solo da riga di comando (`--sperimentale-cd8`).

## 8. Riferimenti

- Bell R. M., McCaffrey D. F. (2002). Bias reduction in standard errors for linear regression with
  multi-stage samples. *Survey Methodology* 28(2):169-181.
- Imbens G. W., Kolesár M. (2016). Robust standard errors in small samples: some practical advice.
  *Review of Economics and Statistics* 98(4):701-712.
- Cameron A. C., Miller D. L. (2015). A practitioner's guide to cluster-robust inference.
  *Journal of Human Resources* 50(2):317-372.
- Nadeau C., Bengio Y. (2003). Inference for the generalization error. *Machine Learning*
  52(3):239-281.
- Bergsma W. (2013). A bias-correction for Cramér's V and Tschuprow's T. *Journal of the Korean
  Statistical Society* 42(3):323-328.
- Davison A. C., Hinkley D. V. (1997). *Bootstrap Methods and their Application*. Cambridge
  University Press (bootstrap sui gruppi, intervalli percentile e BCa).
- Holm S. (1979). A simple sequentially rejective multiple test procedure. *Scandinavian Journal
  of Statistics* 6(2):65-70.
