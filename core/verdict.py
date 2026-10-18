"""Sintesi dei controlli: regole dichiarate, applicate ai risultati già calcolati.

La sintesi non è una nuova stima e non è un giudizio sullo studio né sul lavoro di chi lo ha
prodotto. Ogni riga descrive UN controllo con uno di quattro stati:

- verde,  «nessun limite rilevato dalle regole»: il controllo è stato eseguito e nessuna delle
  regole scritte qui sotto segnala un limite. NON significa che la stima sia corretta: le
  regole non vedono le ipotesi biologiche, la qualità delle etichette o i limiti dichiarati;
- giallo, «stima con limiti»: la stima esiste, ma va letta insieme al limite indicato
  (compresi i fatti strutturali del disegno);
- rosso,  «stima non possibile con questi dati»: i dati non contengono l'informazione che
  servirebbe per quella stima (compreso: meno pazienti del minimo per un intervallo);
- grigio, «controllo non eseguito o non valutabile»: non richiesto, non valutabile oppure
  bloccato da un errore d'input (es. un livello che non esiste).

Il verde non compare mai per un controllo con un fatto strutturale, non eseguito, non
valutabile, rifiutato o fuori dal regime calibrato. Non esiste un colore complessivo per
sezione o per studio: ogni sezione riporta solo quanti controlli sono in ciascuno stato. Le
regole sono scritte in RULES (mostrate nella web app e nel report) e non cambiano nessun
calcolo: leggono le classi e le soglie già prodotte dai moduli. La frazione di CD8, funzione
sperimentale disponibile solo da riga di comando, non entra nella sintesi.
"""

from __future__ import annotations

from dataclasses import dataclass

from core.design_audit import MIN_ROWS_FOR_V, MISSING_LEVEL_CLASS, DesignAuditResult
from core.leakage_audit import GAP_ALERT, MIN_PATIENTS_FOR_MODEL_COMPARISON_DEFAULT, LeakageAuditResult
from core.tcr_validation import BIOLOGICAL_ASSUMPTIONS, TcrValidationResult

VERDE, GIALLO, ROSSO, GRIGIO = "verde", "giallo", "rosso", "grigio"
# Modulo A: sotto questo numero di pazienti la soglia di rumore del divario (GAP_ALERT, provata con
# 19 pazienti su GSE125449) non è verificata: nessun verde.
MIN_PATIENTS_LEAKAGE = 5
STATES = (VERDE, GIALLO, ROSSO, GRIGIO)

STATE_LABEL = {
    VERDE: "nessun limite rilevato dalle regole",
    GIALLO: "stima con limiti",
    ROSSO: "stima non possibile con questi dati",
    GRIGIO: "controllo non eseguito o non valutabile",
}

TITLE = "Sintesi dei controlli"

# Fatti strutturali del disegno e stato con cui entrano nella sintesi
STRUCTURAL_STATE = {"uno-a-uno": ROSSO, "annidamento": GIALLO, "esito-determinato": GIALLO,
                    "unità-tecnica": GIALLO}
STRUCTURAL_NAME = {"uno-a-uno": "Fattori coincidenti", "annidamento": "Annidamento",
                   "esito-determinato": "Esito determinato da un fattore", "unità-tecnica": "Unità tecnica"}

RULES = {
    "Disegno": (
        "Un controllo per ogni confronto richiesto: verde se il confronto è «stimabile»; giallo se è "
        "«stimabile con bassa potenza» (meno di 5 unità indipendenti, oppure p-value minimo "
        "raggiungibile sopra 0.05); rosso se è «non stimabile» (nessun paziente con entrambi i "
        "livelli, oppure un fattore tecnico che separa i due livelli); grigio se un livello non "
        "compare nei metadati (errore d'input). Nessun confronto richiesto: grigio. Un controllo per "
        "ogni fatto strutturale (controlli deterministici sulla tabella): rosso per due fattori "
        "coincidenti, giallo per annidamento, esito determinato da un fattore o unità tecnica. "
        "Associazioni fra fattori (Cramér V, l'unico elemento inferenziale), escluse le coppie già "
        "coperte da un fatto strutturale: giallo per ogni coppia con V >= 0.5; grigio per le coppie "
        f"con V non valutabile (meno di {MIN_ROWS_FOR_V} righe oppure, in media, meno di 2 righe per "
        "cella); verde per le altre. Un confronto «stimabile» è giallo invece che verde se un fattore "
        "tecnico è fortemente associato ai due livelli (V >= 0.5) o se quell'associazione non è "
        "valutabile."),
    "Modulo A": (
        "Valutazione con split per paziente: giallo se almeno un fold di test ha classi assenti, se "
        "il training è stato sottocampionato o se il classificatore non è arrivato a convergenza in "
        f"almeno un fold, oppure con meno di {MIN_PATIENTS_LEAKAGE} pazienti; altrimenti verde. "
        "Valutazione con split casuale delle cellule: rosso se la sua macro-F1 supera quella dello "
        f"split per paziente di più di {GAP_ALERT:g}; giallo se è più bassa di più di {GAP_ALERT:g} "
        f"(direzione inattesa) o con meno di {MIN_PATIENTS_LEAKAGE} pazienti; altrimenti verde. "
        "Confronto fra modelli: grigio se non è stato eseguito (meno di "
        f"{MIN_PATIENTS_FOR_MODEL_COMPARISON_DEFAULT} pazienti oppure modalità rapida); giallo se è stato "
        f"eseguito, perché è solo descrittivo (nessun test calibrato con un paziente per fold)."),
    "Modulo B": (
        "Eccesso di discordanza fra compartimenti: rosso se i pazienti sono meno del minimo per un "
        "intervallo (5); giallo se l'intervallo non è stato prodotto per un altro motivo (variabilità "
        "nulla), se l'intervallo al 95% è interamente sopra o interamente sotto lo "
        "zero, oppure se lo è l'intervallo della sensibilità a pesi uguali; verde («nessun eccesso "
        "rilevato») se entrambi gli intervalli includono lo zero. Tasso d'errore per compartimento, per "
        "ogni convenzione: rosso se l'intervallo non è stato prodotto perché i pazienti sono meno di 5; "
        "giallo se la stima esiste ma l'intervallo non è stato prodotto per un altro motivo "
        "(variabilità nulla, tasso pari a 0 o 1); verde se è stato prodotto; grigio se il tasso non è "
        "stato calcolato (mappa dei marcatori o compartimento di riferimento non forniti) o se una "
        "convenzione non ha unità valutabili."),
}


@dataclass(frozen=True)
class Check:
    """Un controllo: che cosa è stato verificato, in quale stato e con quale descrizione."""
    name: str
    state: str
    text: str

    @property
    def label(self) -> str:
        """Testo dello stato (es. «stima con limiti»)."""
        return STATE_LABEL[self.state]


@dataclass(frozen=True)
class SectionSummary:
    """Controlli di una sezione, con la regola che li produce."""
    section: str
    checks: list[Check]
    rule: str

    @property
    def counts(self) -> dict[str, int]:
        """Numero di controlli per stato (solo gli stati presenti, nell'ordine di STATES)."""
        return {s: n for s in STATES if (n := sum(c.state == s for c in self.checks))}

    @property
    def counts_text(self) -> str:
        """Conteggi per stato in una riga di testo."""
        return "; ".join(f"{n} × {STATE_LABEL[s]}" for s, n in self.counts.items())


def _coppie(n: int) -> str:
    return "1 coppia" if n == 1 else f"{n} coppie"


def _is_structural(p) -> bool:
    return bool(p.a_nested_in_b or p.b_nested_in_a or p.one_to_one)


def design_summary(r: DesignAuditResult) -> SectionSummary:
    """Controlli dell'audit del disegno: confronti richiesti, fatti strutturali, associazioni."""
    checks = []
    for c in r.comparisons:
        name = f"Confronto {c.factor}: {c.level_a} vs {c.level_b}"
        limits = [f"'{f}' fortemente associato al confronto (Cramér V = {v:.2f})" for f, v in c.partial_confounders]
        if c.unassessed_technical:
            limits.append("associazione con " + ", ".join(repr(f) for f in c.unassessed_technical)
                          + " non valutabile (tabella troppo piccola)")
        lim = ("; " + "; ".join(limits)) if limits else ""
        if c.classification == "stimabile":
            checks.append(Check(name, GIALLO if limits else VERDE, f"stimabile con {c.n_units} unità indipendenti" + lim))
        elif c.classification == "stimabile con bassa potenza":
            checks.append(Check(name, GIALLO, f"stimabile con bassa potenza: {c.n_units} unità indipendenti" + lim))
        elif c.classification == MISSING_LEVEL_CLASS:
            checks.append(Check(name, GRIGIO, "errore d'input: " + c.sentence.split(": ", 1)[1]))
        else:
            checks.append(Check(name, ROSSO, "non stimabile: i due livelli non sono separabili dal paziente "
                                             "o da un fattore tecnico in questo disegno"))
    if not r.comparisons:
        checks.append(Check("Confronti", GRIGIO, "nessun confronto richiesto: controllo non eseguito"))
    for f in r.findings:
        if f.kind in STRUCTURAL_STATE:
            checks.append(Check(f"{STRUCTURAL_NAME[f.kind]}: {' / '.join(f.factors)}",
                                STRUCTURAL_STATE[f.kind], f.sentence))
    pairs = [p for p in r.pairs if not _is_structural(p)]
    strong = [p for p in pairs if p.v_alarm]
    not_evaluable = [p for p in pairs if p.cramer_v is None]
    evaluable = [p for p in pairs if p.cramer_v is not None and not p.v_alarm]
    for p in strong:
        checks.append(Check(f"Associazione {p.factor_a} × {p.factor_b}", GIALLO,
                            f"Cramér V = {p.cramer_v:.2f}, sopra la soglia di 0.5: l'effetto di un fattore "
                            f"è separabile dall'altro solo in parte"))
    if evaluable:
        checks.append(Check("Associazioni fra fattori", VERDE,
                            f"Cramér V sotto la soglia di 0.5 per {_coppie(len(evaluable))} "
                            f"fra quelle valutabili (non significa disegno bilanciato)"))
    if not_evaluable:
        checks.append(Check("Associazioni fra fattori", GRIGIO,
                            f"Cramér V non valutabile per {_coppie(len(not_evaluable))} (tabelle troppo "
                            f"piccole)"))
    return SectionSummary("Disegno", checks, RULES["Disegno"])


def leakage_summary(r: LeakageAuditResult) -> SectionSummary:
    """Controlli del Modulo A."""
    checks = []
    limits = []
    if r.grouped.n_folds_with_absent_classes:
        limits.append(f"classi assenti in {r.grouped.n_folds_with_absent_classes} fold di test")
    if r.grouped.capped or r.random.capped:
        limits.append("training sottocampionato per paziente")
    nc = r.grouped.n_not_converged + r.random.n_not_converged
    if nc:
        limits.append(f"classificatore non convergente in {nc} fold")
    few = r.n_patients < MIN_PATIENTS_LEAKAGE
    few_txt = (f"solo {r.n_patients} pazienti: con meno di {MIN_PATIENTS_LEAKAGE} la soglia di rumore "
               f"di {GAP_ALERT:g} non è verificata")
    if few:
        limits.append(few_txt)
    base = f"macro-F1 = {r.grouped.mean:.3f} ± {r.grouped.std:.3f} su {len(r.grouped.fold_scores)} fold"
    checks.append(Check("Valutazione con split per paziente", GIALLO if limits else VERDE,
                        base + ("; " + "; ".join(limits) if limits else "")))
    if r.gap > GAP_ALERT:
        checks.append(Check("Valutazione con split casuale delle cellule", ROSSO,
                            f"su questi dati lo split casuale dà una macro-F1 più alta di {r.gap:+.3f} "
                            f"rispetto allo split per paziente: non stima la prestazione su pazienti nuovi"))
    elif r.gap < -GAP_ALERT:
        checks.append(Check("Valutazione con split casuale delle cellule", GIALLO,
                            f"su questi dati lo split casuale dà una macro-F1 più bassa di {-r.gap:.3f} "
                            f"rispetto allo split per paziente, oltre la soglia di {GAP_ALERT:g}: direzione "
                            f"inattesa (variabilità fra fold, classi assenti o pochi pazienti); il divario "
                            f"non dice nulla sul leakage"))
    elif few:
        checks.append(Check("Valutazione con split casuale delle cellule", GIALLO,
                            f"divario rispetto allo split per paziente {r.gap:+.3f}; " + few_txt))
    else:
        checks.append(Check("Valutazione con split casuale delle cellule", VERDE,
                            f"divario rispetto allo split per paziente {r.gap:+.3f}, entro la soglia di "
                            f"±{GAP_ALERT:g}"))
    if r.model_comparison is None:
        checks.append(Check("Confronto fra modelli", GRIGIO,
                            f"non eseguito (meno di {MIN_PATIENTS_FOR_MODEL_COMPARISON_DEFAULT} pazienti "
                            f"oppure modalità rapida)"))
    else:
        mc = r.model_comparison
        checks.append(Check("Confronto fra modelli", GIALLO,
                            f"eseguito con un paziente alla volta fuori, ogni modello contro "
                            f"{mc.reference_model} (riferimento fissato prima); solo descrittivo: nessun "
                            f"test è calibrato in questo regime"))
    return SectionSummary("Modulo A", checks, RULES["Modulo A"])


def tcr_summary(r: TcrValidationResult) -> SectionSummary:
    """Controlli del Modulo B."""
    d = r.discordance
    e = d.equal_weights
    name = "Eccesso di discordanza fra compartimenti"
    n_min = r.min_patients_for_ci
    if not d.sufficient and d.n_patients < n_min:
        checks = [Check(name, ROSSO, f"intervallo non prodotto: {d.n_patients} pazienti, ne servono almeno {n_min}")]
    elif not d.sufficient:
        checks = [Check(name, GIALLO, f"stima {d.mean_excess:+.3f}, intervallo non prodotto: {d.reason}")]
    else:
        txt = f"{d.mean_excess:+.3f} [{d.ci_low:+.3f}, {d.ci_high:+.3f}]"
        eq_txt = (f"; pesi uguali {e.mean:+.3f} [{e.ci_low:+.3f}, {e.ci_high:+.3f}]"
                  if e is not None and e.sufficient else "")
        if d.ci_low > 0:
            checks = [Check(name, GIALLO, txt + eq_txt + ": entro lo stesso clone le etichette differiscono "
                            "fra compartimenti oltre il rumore entro compartimento; le quantità calcolate da "
                            "queste etichette vanno lette insieme al tasso d'errore stimato")]
        elif d.ci_high < 0:
            checks = [Check(name, GIALLO, txt + eq_txt + ": intervallo interamente sotto lo zero; con "
                            "compartimenti di dimensioni diverse lo stimatore pesato può essere negativo anche "
                            "in presenza di errori, quindi non indica assenza di errori")]
        elif e is not None and e.sufficient and (e.ci_low > 0 or e.ci_high < 0):
            checks = [Check(name, GIALLO, txt + eq_txt + ": l'intervallo pesato include lo zero, quello a "
                            "pesi uguali no")]
        else:
            checks = [Check(name, VERDE, txt + eq_txt + ": nessun eccesso rilevato (l'intervallo include lo "
                            "zero; non dimostra l'assenza di errori)")]
    if not r.conventions:
        checks.append(Check("Tasso d'errore per compartimento", GRIGIO,
                            "non calcolato: mappa dei marcatori o compartimento di riferimento non forniti"))
    else:
        for c in r.conventions.values():
            if not c.by_compartment:
                checks.append(Check(f"Tasso d'errore [{c.name}]", GRIGIO,
                                    "nessuna unità valutabile: nessun clone ha un'identità di riferimento "
                                    "oppure nessuna cellula giudicata fuori dal compartimento di riferimento"))
            for comp, b in c.by_compartment.items():
                n = f"Tasso d'errore [{c.name}] in {comp}"
                if b.sufficient:
                    checks.append(Check(n, VERDE, f"{b.mean:.3f} [{b.ci_low:.3f}, {b.ci_high:.3f}], "
                                                  f"{b.n_groups} pazienti"))
                elif b.n_groups < n_min:
                    checks.append(Check(n, ROSSO, f"stima {b.mean:.3f}, intervallo non prodotto: {b.reason}"))
                else:
                    checks.append(Check(n, GIALLO, f"stima {b.mean:.3f} su {b.n_groups} pazienti, intervallo "
                                                   f"non prodotto: {b.reason}"))
    return SectionSummary("Modulo B", checks, RULES["Modulo B"])


def summarize(design: DesignAuditResult | None = None, leakage: LeakageAuditResult | None = None,
              tcr: TcrValidationResult | None = None) -> list[SectionSummary]:
    """Sintesi delle sole sezioni eseguite, nell'ordine dell'interfaccia."""
    out = []
    if design is not None:
        out.append(design_summary(design))
    if leakage is not None:
        out.append(leakage_summary(leakage))
    if tcr is not None:
        out.append(tcr_summary(tcr))
    return out


USE_WARNING = ("Strumento per la ricerca. Non è un dispositivo medico e non va usato per decisioni "
               "cliniche su singoli pazienti.")


def standing_limits() -> list[str]:
    """Limiti dichiarati formalmente, sempre riportati insieme alla sintesi."""
    return [
        USE_WARNING,
        "Il verde significa che nessuna delle regole scritte segnala un limite, non che la stima sia "
        "corretta. I colori descrivono i singoli controlli: non valutano lo studio e non sostituiscono la "
        "lettura dei numeri e delle definizioni di ciascuna sezione.",
        "Il Modulo B è stato eseguito su dati reali solo in GSE278694 (PDAC), il dataset su cui è stato "
        "sviluppato: lì la riproduzione dei numeri della tesi verifica la fedeltà del porting, non la "
        "correttezza del metodo. La sua generalità ad altri tumori e laboratori non è verificata.",
        BIOLOGICAL_ASSUMPTIONS,
        "I fatti strutturali del disegno sono controlli deterministici sulla tabella dei metadati; "
        "l'unico elemento inferenziale dell'audit del disegno è il Cramér V.",
        "Con pochi pazienti la potenza statistica è bassa: «non significativo» non significa "
        "«equivalente».",
        "Nessuna sezione modifica le etichette originali: lo strumento segnala, non corregge.",
    ]
