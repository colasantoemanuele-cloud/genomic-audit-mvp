"""Propagazione dell'errore di annotazione sulla frazione di CD8 (Intervento 3).

Traduce l'errore misurato dal Modulo B in un intervallo sulla frazione di CD8 per
paziente: "frazione CD8 nel tumore del paziente X: riportata 0.22, plausibile fra A e B".

Metodo
------
1. Matrice di confusione CON DIREZIONE, stimata solo sulle cellule del compartimento
   bersaglio che hanno un'identita' di riferimento (``audit_reference_label`` del Modulo
   B, schema solo sangue): da identita' vera {CD4, CD8} a etichetta chiamata {CD4, CD8,
   altro}. "Altro" e' ogni etichetta diversa da CD4/CD8 (es. NK); serve solo qui,
   ``marker_error_rate`` resta invariata. La matrice e' POOLED fra pazienti.
2. Frazione riportata del paziente: r = n_CD8 / (n_CD4 + n_CD8), contando le cellule del
   compartimento bersaglio etichettate CD4 o CD8.
3. Inversione. Con x, y = cellule vere CD4 e CD8 del paziente e a4 = P(chiamata CD4 |
   vera CD4), b4 = P(chiamata CD8 | vera CD4), a8 = P(chiamata CD4 | vera CD8),
   b8 = P(chiamata CD8 | vera CD8):
       n_CD4 = a4 x + a8 y,   n_CD8 = b4 x + b8 y
   da cui, con (n_CD4, n_CD8) proporzionali a (1 - r, r):
       x ∝ b8 (1 - r) - a8 r,   y ∝ a4 r - b4 (1 - r),   frazione vera = y / (x + y)
   troncata a [0, 1]. Le chiamate "altro" non entrano nel denominatore: la perdita
   differenziale verso "altro" e' gia' nelle a/b, che per riga non sommano a 1.
   Condizionamento: J = b8/(a8+b8) - b4/(a4+b4), differenza fra le probabilita' di essere
   chiamata CD8 (fra le chiamate CD4/CD8) per una vera CD8 e per una vera CD4. Con J
   vicino a zero l'inversione amplifica il rumore senza limite: intervallo rifiutato se
   J < 0.2 oppure se l'intervallo bootstrap di J include lo zero.
4. Incertezza: ricampionamento dei PAZIENTI con lo stesso schema (stesse chiamate al
   generatore casuale) di ``core.stats.cluster_bootstrap`` -- che restituisce solo media
   e IC, non le repliche, e non va modificato -- per la matrice; per ogni replica una
   frazione riportata estratta da Beta(n_CD8 + 1/2, n_CD4 + 1/2) (Monte Carlo sui
   conteggi del paziente); intervallo = percentili 2.5-97.5 delle frazioni invertite.
5. Sensibilita': l'errore di ciascuna riga della matrice (1 - chiamata corretta) e'
   moltiplicato per 0.5, 1 e 2 mantenendone la composizione; tre intervalli sempre
   riportati, nessuno indicato come "il risultato".

Assunzioni
----------
- L'identita' di riferimento dal sangue e' considerata corretta.
- I cloni condivisi con il sangue (quelli con riferimento) sono rappresentativi dei non
  condivisi: e' l'assunzione che gli scenari 0.5x/2x mettono alla prova.
- Le cellule chiamate CD4/CD8 sono vere cellule T CD4 o CD8 (le non-T chiamate CD4/CD8,
  es. doppietti, non sono modellate).
"""

from __future__ import annotations

from collections.abc import Iterator
from dataclasses import dataclass

import numpy as np
import pandas as pd

from core.tcr_validation import MIN_PATIENTS_FOR_CI_DEFAULT

SCENARIOS = (0.5, 1.0, 2.0)
J_MIN = 0.2
MIN_CELLS_PER_PATIENT = 20
MAX_INVALID_FRACTION = 0.10

ASSUMPTIONS_TEXT = (
    "Assunzioni. (1) L'identita' di riferimento stimata dal sangue e' considerata corretta. "
    "(2) La matrice di confusione e' stimata solo sui cloni condivisi con il sangue, che per "
    "costruzione sono quelli espansi (almeno 3 cellule nel sangue e identita' netta), ed e' "
    "unica per tutti i pazienti (pooled), non per paziente. Nei dati originali (PDAC, "
    "GSE278694) l'errore di annotazione cresce con la dimensione del clone (rho = +0.177): "
    "e' quindi possibile che la matrice SOVRASTIMI l'errore dei cloni non condivisi, ma non "
    "e' verificabile, perche' quelle cellule non hanno riferimento. In simulazione nessuno "
    "scenario copre il valore vero in tutti i casi: se i cloni condivisi sbagliano 2-4 volte "
    "piu' degli altri lo scenario 1x scende fino al 28% di copertura e lo 0.5x resta intorno "
    "al 92-95%; se sbagliano quanto gli altri, lo 0.5x scende all'88%. Per questo sono "
    "riportati sempre tre scenari (errore 0.5x, 1x, 2x), nessuno dei quali e' \"il "
    "risultato\". (3) Le cellule etichettate CD4 o CD8 sono vere cellule T; doppietti "
    "e altre cellule non-T etichettate CD4/CD8 non sono modellati."
)


EXPERIMENTAL_NOTE = (
    "FUNZIONE SPERIMENTALE. Nei test di robustezza (cloni condivisi con il sangue che sbagliano "
    "piu' o meno degli altri) nessuno dei tre scenari garantisce la copertura nominale "
    "dell'intervallo al 95%. Sui dati reali su cui e' stata provata (GSE278694) gli intervalli "
    "sono risultati molto ampi, e lo scenario 2x non era calcolabile. Gli intervalli non vanno "
    "letti come una stima della frazione vera di CD8."
)

SCENARIO_WARNING = (
    "Per ogni paziente: frazione di cellule etichettate CD8 sul totale delle cellule etichettate "
    "CD4 o CD8 nel compartimento (denominatore dichiarato), intervallo dei soli conteggi (senza "
    "correzione) e intervalli plausibili al 95% dopo la correzione per l'errore di annotazione, in "
    "tre scenari di errore sui cloni non condivisi con il sangue (0.5x, 1x, 2x). Nessuno scenario "
    "e' \"il risultato\": la loro distanza mostra quanto la conclusione dipende dall'assunzione."
)


@dataclass(frozen=True)
class ScenarioInterval:
    factor: float
    low: float | None
    high: float | None
    refused_reason: str | None


@dataclass(frozen=True)
class PatientCD8:
    patient: str
    n_cd4_called: int
    n_cd8_called: int
    reported: float | None
    naive_low: float | None
    naive_high: float | None
    scenarios: dict[float, ScenarioInterval]
    sentence: str


@dataclass(frozen=True)
class CD8PropagationResult:
    target_compartment: str
    n_reference_cells: int
    n_reference_patients: int
    matrix: pd.DataFrame | None
    youden_j: float | None
    refused_reason: str | None
    patients: list[PatientCD8]
    assumptions: str
    narrative: str


def patient_resample_draws(groups: np.ndarray, n_boot: int, seed: int) -> Iterator[np.ndarray]:
    """Indici delle osservazioni per ogni replica del bootstrap sui pazienti. Riproduce
    esattamente lo schema di ``core.stats.cluster_bootstrap`` (stesso generatore, stesse
    chiamate, stesso ordine di concatenazione): la media dei valori su questi indici da'
    le stesse repliche, e quindi lo stesso IC, di quella funzione. Verificato in
    tests/test_cd8_propagation_calibration.py."""
    g = np.asarray(groups)
    unique = np.unique(g)
    by_group = {grp: np.flatnonzero(g == grp) for grp in unique}
    rng = np.random.default_rng(seed)
    for _ in range(n_boot):
        draw = rng.choice(unique, size=len(unique), replace=True)
        yield np.concatenate([by_group[grp] for grp in draw])


def _confusion(true_is_cd8: np.ndarray, called: np.ndarray) -> np.ndarray | None:
    """Matrice 2x3 di probabilita': righe vera CD4 / vera CD8, colonne chiamata
    CD4 / CD8 / altro. ``called``: 0 = CD4, 1 = CD8, 2 = altro. None se una riga e' vuota."""
    m = np.zeros((2, 3))
    np.add.at(m, (true_is_cd8.astype(int), called), 1.0)
    tot = m.sum(axis=1, keepdims=True)
    if (tot == 0).any():
        return None
    return m / tot


def _scale(P: np.ndarray, k: float) -> np.ndarray | None:
    """Moltiplica per k l'errore di ogni riga (1 - chiamata corretta), mantenendone la
    composizione fra 'chiamata sbagliata' e 'altro'. None se l'errore scalato supera 1."""
    out = P.copy()
    for t in (0, 1):
        err = 1.0 - P[t, t]
        new_err = k * err
        if new_err >= 1.0:
            return None
        if err > 0:
            out[t] = P[t] * (new_err / err)
        out[t, t] = 1.0 - new_err
    return out


def _youden(P: np.ndarray) -> float:
    a4, b4, a8, b8 = P[0, 0], P[0, 1], P[1, 0], P[1, 1]
    if a8 + b8 <= 0 or a4 + b4 <= 0:
        return float("nan")
    return float(b8 / (a8 + b8) - b4 / (a4 + b4))


def _invert(P: np.ndarray, r: np.ndarray) -> np.ndarray:
    a4, b4, a8, b8 = P[..., 0, 0], P[..., 0, 1], P[..., 1, 0], P[..., 1, 1]
    x = b8 * (1 - r) - a8 * r
    y = a4 * r - b4 * (1 - r)
    with np.errstate(divide="ignore", invalid="ignore"):
        frac = y / (x + y)
    frac = np.where((x + y) > 0, frac, np.nan)
    return np.clip(frac, 0.0, 1.0)


def _fmt_interval(iv: ScenarioInterval) -> str:
    label = f"{iv.factor:g}x"
    if iv.low is None:
        return f"{label}: non prodotto ({iv.refused_reason})"
    return f"{label}: fra {iv.low:.2f} e {iv.high:.2f}"


def cd8_fraction_intervals(
    obs: pd.DataFrame,
    patient_col: str,
    compartment_col: str,
    celltype_col: str,
    reference_col: str,
    target_compartment: str,
    cd4_label: str = "CD4T",
    cd8_label: str = "CD8T",
    n_boot: int = 1000,
    seed: int = 0,
    min_patients: int = MIN_PATIENTS_FOR_CI_DEFAULT,
    scenarios: tuple[float, ...] = SCENARIOS,
) -> CD8PropagationResult:
    """Intervalli sulla frazione di CD8 per paziente nel compartimento bersaglio.

    ``obs``: una riga per cellula, con l'etichetta assegnata (``celltype_col``) e
    l'identita' di riferimento del clone (``reference_col``, tipicamente la colonna
    ``audit_reference_label`` prodotta dal Modulo B; NA se assente)."""
    df = obs[[patient_col, compartment_col, celltype_col, reference_col]].copy()
    df.columns = ["patient", "compartment", "celltype", "ref"]
    df = df[df["compartment"].astype(str) == str(target_compartment)]
    df["patient"] = df["patient"].astype(str)
    lab = df["celltype"].astype(str)

    ref = df["ref"].astype("object")
    ref_ok = ref.notna() & ref.astype(str).isin([cd4_label, cd8_label])
    rdf = df[ref_ok.values]
    r_true_cd8 = (rdf["ref"].astype(str) == cd8_label).values
    r_lab = rdf["celltype"].astype(str).values
    r_called = np.where(r_lab == cd4_label, 0, np.where(r_lab == cd8_label, 1, 2))
    r_groups = rdf["patient"].values
    n_ref_pat = int(pd.unique(r_groups).size)

    refused = None
    P_hat, J_hat = None, None
    if n_ref_pat < min_patients:
        refused = (f"solo {n_ref_pat} pazienti hanno cellule con identita' di riferimento nel "
                   f"compartimento '{target_compartment}', ne servono almeno {min_patients}")
    else:
        P_hat = _confusion(r_true_cd8, r_called)
        if P_hat is None:
            refused = ("mancano cellule di riferimento per una delle due identita' (CD4 o CD8): "
                       "la matrice di confusione non e' stimabile")
        else:
            J_hat = _youden(P_hat)

    boot_P = None
    if refused is None:
        mats = []
        for idx in patient_resample_draws(r_groups, n_boot, seed):
            m = _confusion(r_true_cd8[idx], r_called[idx])
            if m is not None:
                mats.append(m)
        if len(mats) < (1 - MAX_INVALID_FRACTION) * n_boot:
            refused = (f"{n_boot - len(mats)} repliche bootstrap su {n_boot} senza una delle due "
                       f"identita': la matrice e' troppo instabile")
        else:
            boot_P = np.stack(mats)

    # condizionamento per scenario (globale, la matrice e' pooled)
    scen_reason: dict[float, str | None] = {}
    scen_P: dict[float, tuple[np.ndarray, np.ndarray] | None] = {}
    for k in scenarios:
        if refused is not None:
            scen_reason[k], scen_P[k] = refused, None
            continue
        P_k = _scale(P_hat, k)
        if P_k is None:
            scen_reason[k], scen_P[k] = (f"con errore {k:g}x la probabilita' di errore stimata "
                                         f"supera 1"), None
            continue
        j_k = _youden(P_k)
        if not (j_k >= J_MIN):
            scen_reason[k] = (f"matrice di confusione mal condizionata (J = {j_k:.2f}; servono "
                              f"J >= {J_MIN:g} e un intervallo di J che escluda lo zero)")
            scen_P[k] = None
            continue
        B_k = [b for b in (_scale(m, k) for m in boot_P) if b is not None]
        if len(B_k) < (1 - MAX_INVALID_FRACTION) * len(boot_P):
            scen_reason[k], scen_P[k] = (f"con errore {k:g}x la probabilita' di errore supera 1 in "
                                         f"{len(boot_P) - len(B_k)} repliche bootstrap su "
                                         f"{len(boot_P)}"), None
            continue
        B_k = np.stack(B_k)
        j_boot = np.array([_youden(b) for b in B_k])
        if np.nanpercentile(j_boot, 2.5) <= 0:
            scen_reason[k] = (f"matrice di confusione mal condizionata (J = {j_k:.2f}; servono "
                              f"J >= {J_MIN:g} e un intervallo di J che escluda lo zero)")
            scen_P[k] = None
        else:
            scen_reason[k], scen_P[k] = None, (P_k, B_k)

    rng = np.random.default_rng(seed + 1)
    patients: list[PatientCD8] = []
    for pid, g in df.groupby("patient", sort=True):
        glab = g["celltype"].astype(str)
        n4, n8 = int((glab == cd4_label).sum()), int((glab == cd8_label).sum())
        if n4 + n8 < MIN_CELLS_PER_PATIENT:
            reason = (f"troppo poche cellule etichettate {cd4_label}/{cd8_label} "
                      f"({n4 + n8}, servono almeno {MIN_CELLS_PER_PATIENT})")
            patients.append(PatientCD8(
                pid, n4, n8, (n8 / (n4 + n8)) if n4 + n8 else None, None, None,
                {k: ScenarioInterval(k, None, None, reason) for k in scenarios},
                f"Paziente {pid}: {reason}; nessuna frazione riportata o corretta."))
            continue
        reported = n8 / (n4 + n8)
        r_draws = rng.beta(n8 + 0.5, n4 + 0.5, size=boot_P.shape[0] if boot_P is not None else n_boot)
        naive_lo, naive_hi = (float(x) for x in np.percentile(r_draws, (2.5, 97.5)))
        ivs: dict[float, ScenarioInterval] = {}
        for k in scenarios:
            if scen_P[k] is None:
                ivs[k] = ScenarioInterval(k, None, None, scen_reason[k])
                continue
            _, B_k = scen_P[k]
            fr = _invert(B_k, r_draws[: len(B_k)])
            if np.isnan(fr).mean() > MAX_INVALID_FRACTION:
                ivs[k] = ScenarioInterval(k, None, None, "inversione non definita in troppe repliche")
                continue
            lo, hi = np.nanpercentile(fr, (2.5, 97.5))
            ivs[k] = ScenarioInterval(k, float(lo), float(hi), None)
        sentence = (
            f"Frazione CD8 nel compartimento '{target_compartment}' del paziente {pid} "
            f"({cd8_label} sul totale delle cellule etichettate {cd4_label} o {cd8_label}, "
            f"n = {n4 + n8}): riportata {reported:.2f} (IC95% dei soli conteggi "
            f"{naive_lo:.2f}-{naive_hi:.2f}, senza correzione). Intervallo plausibile al 95% "
            f"tenendo conto dell'errore di annotazione misurato, per scenario di errore sui cloni "
            f"non condivisi con il sangue -- " + "; ".join(_fmt_interval(ivs[k]) for k in scenarios)
            + ".")
        patients.append(PatientCD8(pid, n4, n8, reported, naive_lo, naive_hi, ivs, sentence))

    matrix = None
    if P_hat is not None:
        matrix = pd.DataFrame(P_hat, index=[f"vera {cd4_label}", f"vera {cd8_label}"],
                              columns=[f"chiamata {cd4_label}", f"chiamata {cd8_label}", "chiamata altro"])
    if refused:
        head = (f"Intervalli sulla frazione di CD8 non prodotti: {refused}. Sono riportate solo le "
                f"frazioni osservate, senza correzione.")
    else:
        head = (f"Matrice di confusione stimata su {len(rdf)} cellule con identita' di riferimento "
                f"da {n_ref_pat} pazienti (pooled fra pazienti), J = {J_hat:.2f}.")
    narrative = " ".join([head, SCENARIO_WARNING, ASSUMPTIONS_TEXT])
    return CD8PropagationResult(
        target_compartment=str(target_compartment), n_reference_cells=int(len(rdf)),
        n_reference_patients=n_ref_pat, matrix=matrix, youden_j=J_hat, refused_reason=refused,
        patients=patients, assumptions=ASSUMPTIONS_TEXT, narrative=narrative)


def refusal_notes(result: CD8PropagationResult) -> dict[str, int]:
    """Motivi di rifiuto distinti -> numero della nota, nell'ordine in cui compaiono."""
    notes: dict[str, int] = {}
    for p in result.patients:
        for k in sorted(p.scenarios):
            r = p.scenarios[k].refused_reason
            if r is not None and r not in notes:
                notes[r] = len(notes) + 1
    return notes


def format_cd8_text(result: CD8PropagationResult) -> str:
    """Testo per la CLI: intestazione e avvertenze UNA volta, poi una tabella con una riga per
    paziente; i motivi dei rifiuti come note numerate sotto la tabella."""
    notes = refusal_notes(result)

    def cell(iv: ScenarioInterval) -> str:
        return (f"{iv.low:.2f}-{iv.high:.2f}" if iv.low is not None
                else f"non prodotto [{notes[iv.refused_reason]}]")

    scen = sorted({k for p in result.patients for k in p.scenarios})
    header = ["paziente", "n CD4+CD8", "riportata", "IC95% conteggi"] + [f"{k:g}x" for k in scen]
    rows = []
    for p in result.patients:
        rows.append([p.patient, str(p.n_cd4_called + p.n_cd8_called),
                     "-" if p.reported is None else f"{p.reported:.2f}",
                     "-" if p.naive_low is None else f"{p.naive_low:.2f}-{p.naive_high:.2f}"]
                    + [cell(p.scenarios[k]) for k in scen])
    widths = [max(len(r[i]) for r in [header] + rows) for i in range(len(header))]
    line = lambda r: "  ".join(v.ljust(w) for v, w in zip(r, widths))  # noqa: E731
    head = result.narrative.split(SCENARIO_WARNING)[0].strip()
    out = [f"Frazione di CD8 nel compartimento '{result.target_compartment}'", EXPERIMENTAL_NOTE, "",
           head, SCENARIO_WARNING, "",
           line(header), line(["-" * w for w in widths])] + [line(r) for r in rows]
    if notes:
        out.append("")
        out += [f"[{n}] {reason}" for reason, n in notes.items()]
    out += ["", result.assumptions]
    return "\n".join(out)
