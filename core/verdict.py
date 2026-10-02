"""Verdetto di audit per sezione: regole dichiarate, applicate ai risultati gia' calcolati.

Il verdetto non e' una nuova stima: riassume con un colore i risultati dei moduli, secondo
regole fisse e scritte qui (e mostrate all'utente nella web app e nel report). Non sostituisce
la lettura dei numeri e dei limiti di ciascuna sezione.
"""

from __future__ import annotations

from dataclasses import dataclass

from core.cd8_propagation import CD8PropagationResult
from core.design_audit import DesignAuditResult
from core.leakage_audit import GAP_ALERT, LeakageAuditResult
from core.tcr_validation import TcrValidationResult

VERDE, GIALLO, ROSSO = "verde", "giallo", "rosso"

RULES = {
    "Disegno": (
        "Rosso: almeno un confronto richiesto e' non stimabile (fattore confuso con il paziente "
        "o con un fattore tecnico). Giallo: almeno un confronto e' stimabile solo con bassa potenza, "
        "oppure due fattori hanno un'associazione forte (Cramér V >= 0.5). Verde: nessuna delle due."),
    "Modulo A": (
        f"Rosso: lo split casuale sovrastima la macro-F1 dello split per paziente di piu' di "
        f"{GAP_ALERT:g} (qualunque valutazione che non separa i pazienti e' inaffidabile su questi "
        f"dati). Giallo: classi assenti da almeno un fold di test, training sottocampionato, "
        f"classificatore non convergente in almeno un fold, oppure confronto fra modelli non "
        f"eseguito. Verde: nessuna delle condizioni precedenti."),
    "Modulo B": (
        "Rosso: l'eccesso di discordanza fra compartimenti ha un intervallo al 95% interamente sopra "
        "lo zero (errore sistematico di annotazione). Giallo: numerosita' insufficiente per un "
        "intervallo, oppure intervallo che include lo zero ma tasso d'errore per compartimento non "
        "stimabile. Verde: intervallo che include lo zero e tassi d'errore stimati."),
}


@dataclass(frozen=True)
class SectionVerdict:
    section: str
    color: str
    reasons: list[str]
    rule: str


def design_verdict(r: DesignAuditResult) -> SectionVerdict:
    reasons, color = [], VERDE
    bad = [c for c in r.comparisons if c.classification == "non stimabile"]
    low = [c for c in r.comparisons if c.classification == "stimabile con bassa potenza"]
    strong = [p for p in r.pairs if p.v_alarm]
    if bad:
        color = ROSSO
        reasons += [f"non stimabile: {c.factor} {c.level_a} vs {c.level_b}" for c in bad]
    if low or strong:
        color = color if color == ROSSO else GIALLO
        reasons += [f"bassa potenza: {c.factor} {c.level_a} vs {c.level_b} ({c.n_units} unita')" for c in low]
        reasons += [f"associazione forte: {p.factor_a} x {p.factor_b} (V = {p.cramer_v:.2f})" for p in strong]
    if not reasons:
        reasons.append("nessun confronto non stimabile o a bassa potenza, nessuna associazione forte")
    return SectionVerdict("Disegno", color, reasons, RULES["Disegno"])


def leakage_verdict(r: LeakageAuditResult) -> SectionVerdict:
    reasons, color = [], VERDE
    if r.gap > GAP_ALERT:
        color = ROSSO
        reasons.append(f"lo split casuale sovrastima di {r.gap:+.3f} la macro-F1 per paziente")
    yellow = []
    if r.grouped.n_folds_with_absent_classes:
        yellow.append(f"classi assenti in {r.grouped.n_folds_with_absent_classes} fold di test")
    if r.grouped.capped or r.random.capped:
        yellow.append("training sottocampionato per paziente")
    if r.model_comparison is None:
        yellow.append("confronto fra modelli non eseguito")
    nc = r.grouped.n_not_converged + r.random.n_not_converged
    if nc:
        yellow.append(f"classificatore non convergente in {nc} fold")
    if yellow:
        color = color if color == ROSSO else GIALLO
        reasons += yellow
    if not reasons:
        reasons.append(f"divario {r.gap:+.3f}, nessuna condizione di attenzione")
    return SectionVerdict("Modulo A", color, reasons, RULES["Modulo A"])


def tcr_verdict(r: TcrValidationResult) -> SectionVerdict:
    d = r.discordance
    if not d.sufficient:
        return SectionVerdict("Modulo B", GIALLO, [f"numerosita' insufficiente ({d.n_patients} pazienti)"],
                              RULES["Modulo B"])
    if d.ci_low > 0:
        return SectionVerdict("Modulo B", ROSSO,
                              [f"eccesso di discordanza {d.mean_excess:+.3f} [{d.ci_low:+.3f}, {d.ci_high:+.3f}]"],
                              RULES["Modulo B"])
    unresolved = []
    if r.conventions:
        for c in r.conventions.values():
            unresolved += [f"{c.name}/{comp}" for comp, b in c.by_compartment.items() if not b.sufficient]
    if unresolved:
        return SectionVerdict("Modulo B", GIALLO,
                              ["tasso d'errore senza intervallo: " + ", ".join(unresolved)], RULES["Modulo B"])
    return SectionVerdict("Modulo B", VERDE,
                          [f"eccesso di discordanza {d.mean_excess:+.3f}, intervallo che include lo zero"],
                          RULES["Modulo B"])


def standing_limits(cd8: CD8PropagationResult | None = None) -> list[str]:
    """Limiti dichiarati formalmente, sempre riportati nel verdetto."""
    out = [
        "Il Modulo B e' validato su un solo dataset reale (GSE278694, PDAC): la sua generalita' ad "
        "altri tumori non e' verificata.",
        "Con pochi pazienti la potenza statistica e' bassa: 'non significativo' non significa "
        "'equivalente'.",
        "Il verdetto a colori riassume regole fisse: non sostituisce la lettura dei numeri e delle "
        "definizioni di ciascuna sezione.",
        "Nessuna sezione modifica le etichette originali: lo strumento segnala, non corregge.",
    ]
    if cd8 is not None:
        out.append("La frazione di CD8 e' una funzione sperimentale: la copertura nominale dei suoi "
                   "intervalli non e' garantita.")
    return out
