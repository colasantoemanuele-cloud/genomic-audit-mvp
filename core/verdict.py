"""Sintesi dei controlli: regole dichiarate, applicate ai risultati già calcolati.

La sintesi non è una nuova stima e non è un giudizio sullo studio né sul lavoro di chi lo ha
prodotto. Ogni riga descrive UN controllo e dice che cosa i dati permettono di stimare, con uno
di quattro stati:

- verde,  «stima affidabile»: il controllo è stato eseguito e non segnala limiti;
- giallo, «stima con limiti»: la stima esiste, ma va letta insieme al limite indicato;
- rosso,  «stima non possibile con questi dati»: i dati non contengono l'informazione che
  servirebbe per quella stima;
- grigio, «controllo non eseguito o non valutabile»: non richiesto, rifiutato dallo strumento
  oppure non valutabile. Non è mai mostrato in verde.

Non esiste un colore complessivo per sezione o per studio: ogni sezione riporta solo quanti
controlli sono in ciascuno stato. Le regole sono scritte in RULES (mostrate nella web app e nel
report) e non cambiano nessun calcolo: leggono le classi e le soglie già prodotte dai moduli.
"""

from __future__ import annotations

from dataclasses import dataclass

from core.cd8_propagation import CD8PropagationResult
from core.design_audit import DesignAuditResult
from core.leakage_audit import GAP_ALERT, LeakageAuditResult
from core.tcr_validation import TcrValidationResult

VERDE, GIALLO, ROSSO, GRIGIO = "verde", "giallo", "rosso", "grigio"
STATES = (VERDE, GIALLO, ROSSO, GRIGIO)

STATE_LABEL = {
    VERDE: "stima affidabile",
    GIALLO: "stima con limiti",
    ROSSO: "stima non possibile con questi dati",
    GRIGIO: "controllo non eseguito o non valutabile",
}

TITLE = "Sintesi dei controlli"

RULES = {
    "Disegno": (
        "Un controllo per ogni confronto richiesto: verde se il confronto è «stimabile»; giallo se è "
        "«stimabile con bassa potenza» (meno di 5 unità indipendenti, oppure p-value minimo "
        "raggiungibile sopra 0.05); rosso se è «non stimabile» (nessun paziente con entrambi i "
        "livelli, oppure un fattore tecnico che separa i due livelli). Se non è stato richiesto "
        "nessun confronto: grigio. Un controllo per le associazioni fra fattori: giallo per ogni "
        "coppia con Cramér V >= 0.5; grigio per le coppie con V non valutabile (meno di 10 righe "
        "oppure meno di 2 righe attese per cella); verde per le coppie valutabili sotto la soglia."),
    "Modulo A": (
        "Valutazione con split per paziente: giallo se almeno un fold di test ha classi assenti, se "
        "il training è stato sottocampionato o se il classificatore non è arrivato a convergenza in "
        "almeno un fold; altrimenti verde. Valutazione con split casuale delle cellule: rosso se la "
        f"sua macro-F1 supera quella dello split per paziente di più di {GAP_ALERT:g}; altrimenti "
        "verde. Confronto fra modelli: grigio se non è stato eseguito (meno di 8 pazienti oppure "
        "modalità rapida); altrimenti verde."),
    "Modulo B": (
        "Eccesso di discordanza fra compartimenti: rosso se i pazienti sono meno del minimo per un "
        "intervallo (5); giallo se l'intervallo al 95% è interamente sopra lo zero; verde se "
        "l'intervallo include lo zero. Tasso d'errore per compartimento, per ogni convenzione: "
        "rosso se l'intervallo non è stato prodotto (meno di 5 pazienti); verde se è stato "
        "prodotto; grigio se il tasso non è stato calcolato (mappa dei marcatori o compartimento di "
        "riferimento non forniti)."),
    "Frazione di CD8": (
        "Funzione sperimentale: mai verde. Grigio se lo strumento ha rifiutato il calcolo (meno di "
        "5 pazienti con riferimento, cellule di riferimento mancanti) oppure se nello scenario 1x "
        "nessun intervallo è stato prodotto (indice J della matrice di confusione sotto 0.2); "
        "giallo se almeno un intervallo è stato prodotto, perché la copertura nominale non è "
        "garantita."),
}


@dataclass(frozen=True)
class Check:
    """Un controllo: che cosa è stato verificato, in quale stato e con quale descrizione."""
    name: str
    state: str
    text: str

    @property
    def label(self) -> str:
        return STATE_LABEL[self.state]


@dataclass(frozen=True)
class SectionSummary:
    section: str
    checks: list[Check]
    rule: str

    @property
    def counts(self) -> dict[str, int]:
        """Numero di controlli per stato (solo gli stati presenti, nell'ordine di STATES)."""
        return {s: n for s in STATES if (n := sum(c.state == s for c in self.checks))}

    @property
    def counts_text(self) -> str:
        return "; ".join(f"{n} × {STATE_LABEL[s]}" for s, n in self.counts.items())


def _coppie(n: int) -> str:
    return "1 coppia" if n == 1 else f"{n} coppie"


def design_summary(r: DesignAuditResult) -> SectionSummary:
    checks = []
    for c in r.comparisons:
        name = f"Confronto {c.factor}: {c.level_a} vs {c.level_b}"
        if c.classification == "stimabile":
            checks.append(Check(name, VERDE, f"stimabile con {c.n_units} unità indipendenti"))
        elif c.classification == "stimabile con bassa potenza":
            checks.append(Check(name, GIALLO, f"stimabile con bassa potenza: {c.n_units} unità indipendenti"))
        else:
            checks.append(Check(name, ROSSO, "non stimabile: i due livelli non sono separabili dal paziente "
                                             "o da un fattore tecnico in questo disegno"))
    if not r.comparisons:
        checks.append(Check("Confronti", GRIGIO, "nessun confronto richiesto: controllo non eseguito"))
    strong = [p for p in r.pairs if p.v_alarm]
    not_evaluable = [p for p in r.pairs if p.cramer_v is None]
    evaluable = [p for p in r.pairs if p.cramer_v is not None and not p.v_alarm]
    for p in strong:
        checks.append(Check(f"Associazione {p.factor_a} × {p.factor_b}", GIALLO,
                            f"Cramér V = {p.cramer_v:.2f}, sopra la soglia di 0.5: l'effetto di un fattore "
                            f"è separabile dall'altro solo in parte"))
    if evaluable:
        checks.append(Check("Associazioni fra fattori", VERDE,
                            f"Cramér V sotto la soglia di 0.5 per {_coppie(len(evaluable))} "
                            f"fra quelle valutabili"))
    if not_evaluable:
        checks.append(Check("Associazioni fra fattori", GRIGIO,
                            f"Cramér V non valutabile per {_coppie(len(not_evaluable))} (tabelle troppo "
                            f"piccole)"))
    return SectionSummary("Disegno", checks, RULES["Disegno"])


def leakage_summary(r: LeakageAuditResult) -> SectionSummary:
    checks = []
    limits = []
    if r.grouped.n_folds_with_absent_classes:
        limits.append(f"classi assenti in {r.grouped.n_folds_with_absent_classes} fold di test")
    if r.grouped.capped or r.random.capped:
        limits.append("training sottocampionato per paziente")
    nc = r.grouped.n_not_converged + r.random.n_not_converged
    if nc:
        limits.append(f"classificatore non convergente in {nc} fold")
    base = f"macro-F1 = {r.grouped.mean:.3f} ± {r.grouped.std:.3f} su {len(r.grouped.fold_scores)} fold"
    checks.append(Check("Valutazione con split per paziente", GIALLO if limits else VERDE,
                        base + ("; " + "; ".join(limits) if limits else "")))
    if r.gap > GAP_ALERT:
        checks.append(Check("Valutazione con split casuale delle cellule", ROSSO,
                            f"su questi dati lo split casuale dà una macro-F1 più alta di {r.gap:+.3f} "
                            f"rispetto allo split per paziente: non stima la prestazione su pazienti nuovi"))
    else:
        checks.append(Check("Valutazione con split casuale delle cellule", VERDE,
                            f"divario rispetto allo split per paziente {r.gap:+.3f}, entro la soglia di "
                            f"{GAP_ALERT:g}"))
    if r.model_comparison is None:
        checks.append(Check("Confronto fra modelli", GRIGIO,
                            "non eseguito (meno di 8 pazienti oppure modalità rapida)"))
    else:
        mc = r.model_comparison
        checks.append(Check("Confronto fra modelli", VERDE,
                            f"eseguito con un paziente alla volta fuori; macro-F1 più alta: {mc.best_model} "
                            f"({mc.scores[mc.best_model].mean:.3f})"))
    return SectionSummary("Modulo A", checks, RULES["Modulo A"])


def tcr_summary(r: TcrValidationResult) -> SectionSummary:
    d = r.discordance
    name = "Eccesso di discordanza fra compartimenti"
    if not d.sufficient:
        checks = [Check(name, ROSSO, f"intervallo non prodotto: {d.n_patients} pazienti, ne servono almeno 5")]
    elif d.ci_low > 0:
        checks = [Check(name, GIALLO,
                        f"{d.mean_excess:+.3f} [{d.ci_low:+.3f}, {d.ci_high:+.3f}]: entro lo stesso clone le "
                        f"etichette differiscono fra compartimenti oltre il rumore di base; le quantità "
                        f"calcolate da queste etichette vanno lette insieme al tasso d'errore stimato")]
    else:
        checks = [Check(name, VERDE,
                        f"{d.mean_excess:+.3f} [{d.ci_low:+.3f}, {d.ci_high:+.3f}]: l'intervallo include lo zero")]
    if not r.conventions:
        checks.append(Check("Tasso d'errore per compartimento", GRIGIO,
                            "non calcolato: mappa dei marcatori o compartimento di riferimento non forniti"))
    else:
        for c in r.conventions.values():
            for comp, b in c.by_compartment.items():
                n = f"Tasso d'errore [{c.name}] in {comp}"
                if b.sufficient:
                    checks.append(Check(n, VERDE, f"{b.mean:.3f} [{b.ci_low:.3f}, {b.ci_high:.3f}], "
                                                  f"{b.n_groups} pazienti"))
                else:
                    checks.append(Check(n, ROSSO, f"intervallo non prodotto: {b.n_groups} pazienti, ne "
                                                  f"servono almeno 5"))
    return SectionSummary("Modulo B", checks, RULES["Modulo B"])


def cd8_summary(r: CD8PropagationResult) -> SectionSummary:
    name = f"Intervalli della frazione di CD8 in {r.target_compartment}"
    if r.refused_reason:
        checks = [Check(name, GRIGIO, f"calcolo rifiutato dallo strumento: {r.refused_reason}")]
    else:
        n_ok = sum(1 for p in r.patients if p.scenarios[1.0].low is not None)
        if n_ok == 0:
            why = next((p.scenarios[1.0].refused_reason for p in r.patients
                        if p.scenarios[1.0].refused_reason), "nessun paziente valutabile")
            checks = [Check(name, GRIGIO, f"nessun intervallo prodotto nello scenario 1x: {why}")]
        else:
            checks = [Check(name, GIALLO, f"funzione sperimentale: intervalli prodotti per {n_ok} su "
                                          f"{len(r.patients)} pazienti nello scenario 1x; copertura nominale non "
                                          f"garantita")]
    return SectionSummary("Frazione di CD8", checks, RULES["Frazione di CD8"])


def summarize(design: DesignAuditResult | None = None, leakage: LeakageAuditResult | None = None,
              tcr: TcrValidationResult | None = None, cd8: CD8PropagationResult | None = None,
              ) -> list[SectionSummary]:
    """Sintesi delle sole sezioni eseguite, nell'ordine dell'interfaccia."""
    out = []
    if design is not None:
        out.append(design_summary(design))
    if leakage is not None:
        out.append(leakage_summary(leakage))
    if tcr is not None:
        out.append(tcr_summary(tcr))
    if cd8 is not None:
        out.append(cd8_summary(cd8))
    return out


def standing_limits(cd8: CD8PropagationResult | None = None) -> list[str]:
    """Limiti dichiarati formalmente, sempre riportati insieme alla sintesi."""
    out = [
        "Il Modulo B è validato su un solo dataset reale (GSE278694, PDAC): la sua generalità ad "
        "altri tumori non è verificata.",
        "Con pochi pazienti la potenza statistica è bassa: «non significativo» non significa "
        "«equivalente».",
        "I colori della sintesi descrivono i singoli controlli secondo regole fisse: non valutano lo "
        "studio e non sostituiscono la lettura dei numeri e delle definizioni di ciascuna sezione.",
        "Nessuna sezione modifica le etichette originali: lo strumento segnala, non corregge.",
    ]
    if cd8 is not None:
        out.append("La frazione di CD8 è una funzione sperimentale: la copertura nominale dei suoi "
                   "intervalli non è garantita.")
    return out
