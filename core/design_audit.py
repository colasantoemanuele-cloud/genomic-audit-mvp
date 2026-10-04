"""Audit del disegno e del confondimento (Intervento 1).

Lavora sui soli METADATI (una riga per campione/libreria; se l'input è a livello di
cellula, le righe vengono prima ridotte alle combinazioni distinte dei fattori indicati).
Tre output:

1. Fattori confusi fra loro: per ogni coppia, tabella di contingenza, Cramér V con la
   correzione di Bergsma (2013) e i fatti STRUTTURALI -- annidamento, mappa uno-a-uno,
   esito determinato da un fattore. I fatti strutturali sono rilevati in modo
   deterministico: sono proprietà del disegno, non stime.
2. Confronti stimabili: per ogni confronto richiesto, unità indipendenti (i pazienti)
   e classe `stimabile` / `stimabile con bassa potenza` / `non stimabile`.
3. (Opzionale, richiede la matrice) scomposizione della varianza su pseudobulk
   paziente-tessuto, con intervallo bootstrap sui pazienti; se i fattori sono confusi
   la scomposizione NON viene prodotta ("non identificabile"). NON ESPOSTA in CLI, app o
   report: la copertura dell'intervallo della quota del tessuto è fuori banda (1.000 >
   0.99) dopo due correzioni del metodo -- test marcato xfail(strict=True).

Nessuno di questi output è un giudizio sul lavoro di chi ha disegnato lo studio: molte
strutture segnalate (es. una libreria per coppia paziente-tessuto) sono la norma nei
disegni a singola cellula. Lo strumento le rende esplicite perché limitano quali
conclusioni i dati possono sostenere.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from math import comb

import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse as sp
from scipy.stats import chi2_contingency

from core.leakage_audit import LogCPM
from core.stats import wilcoxon_min_pvalue
from core.tcr_validation import MIN_PATIENTS_FOR_CI_DEFAULT

# Stessa soglia del Modulo B: sotto 5 pazienti nessun intervallo di confidenza.
MIN_UNITS_DEFAULT = MIN_PATIENTS_FOR_CI_DEFAULT
CRAMER_V_THRESHOLD = 0.5
# Il V è "valutabile" solo con almeno 10 righe e almeno 2 righe attese per cella in
# media (n / (r*c) >= 2): sotto, il valore oscilla troppo per essere riportato.
MIN_ROWS_FOR_V = 10
MIN_EXPECTED_PER_CELL = 2.0
ALPHA = 0.05
MIN_PVALUE_NOTE = ("Il p-value minimo non misura la potenza: con questa numerosità solo effetti "
                   "grandi sono rilevabili.")
MIN_CELLS_PER_UNIT_DEFAULT = 20
N_HVG_DEFAULT = 1000
TECHNICAL_ROLES = ("library", "batch", "chemistry", "protocol", "date")

PROTOCOL_NOTE = (
    "Nota fissa sui protocolli: lo snRNA-seq (nuclei) sottorappresenta le cellule "
    "immunitarie e cambia la composizione misurata rispetto allo scRNA-seq; i protocolli "
    "per tessuto fissato (FFPE, 10x Flex) misurano l'RNA tramite sonde e danno profili non "
    "direttamente sovrapponibili a quelli da tessuto fresco. Un confronto fra protocolli è "
    "interpretabile solo se gli stessi pazienti sono misurati con entrambi."
)


# --------------------------------------------------------------------------- #
# Strutture dei risultati
# --------------------------------------------------------------------------- #
@dataclass(frozen=True)
class PairAssociation:
    factor_a: str
    factor_b: str
    n_rows: int
    n_levels_a: int
    n_levels_b: int
    cramer_v: float | None  # None = non valutabile
    v_alarm: bool
    a_nested_in_b: bool
    b_nested_in_a: bool
    one_to_one: bool
    table: pd.DataFrame
    sentence: str


@dataclass(frozen=True)
class StructuralFinding:
    kind: str  # annidamento | uno-a-uno | esito-determinato | unità-tecnica | costante | identificatore
    factors: tuple[str, ...]
    sentence: str


@dataclass(frozen=True)
class ComparisonAssessment:
    factor: str
    level_a: str
    level_b: str
    design: str  # "appaiato (entro paziente)" | "fra pazienti" | "-"
    n_units: int
    n_patients_a: int
    n_patients_b: int
    classification: str  # stimabile | stimabile con bassa potenza | non stimabile
    reasons: tuple[str, ...]
    sentence: str
    # p-value minimo raggiungibile da un test esatto con queste unità (None se non definito)
    min_pvalue: float | None = None
    min_pvalue_test: str | None = None


@dataclass(frozen=True)
class VarianceShare:
    factor: str
    share: float
    ci_low: float
    ci_high: float


@dataclass(frozen=True)
class VarianceDecomposition:
    identifiable: bool
    reason: str
    n_units: int
    n_patients: int
    n_genes: int
    shares: list[VarianceShare]
    residual_share: float | None
    sentence: str


@dataclass(frozen=True)
class DesignAuditResult:
    n_rows: int
    roles: dict[str, str]
    pairs: list[PairAssociation]
    findings: list[StructuralFinding]
    comparisons: list[ComparisonAssessment]
    notes: list[str] = field(default_factory=list)
    narrative: str = ""


# --------------------------------------------------------------------------- #
# Input
# --------------------------------------------------------------------------- #
MISSING_LEVEL = "NA"


def as_levels(df: pd.DataFrame) -> tuple[pd.DataFrame, dict[str, int]]:
    """Converte i fattori in stringhe con i valori mancanti come livello esplicito "NA".
    Senza questo passaggio (con pandas >= 3, ``astype(str)`` lascia i mancanti come
    mancanti) le righe con un valore mancante venivano scartate in silenzio da crosstab e
    groupby, ma contate nel numero di unità. Ritorna anche il conteggio dei mancanti per
    colonna, da dichiarare."""
    missing = {c: int(df[c].isna().sum()) for c in df.columns if df[c].isna().any()}
    out = df.astype(object).where(df.notna(), MISSING_LEVEL).astype(str)
    return out, missing


def sample_sheet_from_obs(obs: pd.DataFrame, cols: list[str]) -> pd.DataFrame:
    """Riduce una tabella a livello di cellula (es. ``adata.obs``) alle combinazioni
    distinte dei fattori indicati: le cellule non sono unità indipendenti del disegno.
    Se fra le colonne c'è un identificativo di libreria, ogni riga risultante è una
    libreria."""
    missing = [c for c in cols if c not in obs.columns]
    if missing:
        raise ValueError(f"colonne non trovate nei metadati: {missing}")
    return as_levels(obs[cols])[0].drop_duplicates().reset_index(drop=True)


# --------------------------------------------------------------------------- #
# 1. Associazioni e fatti strutturali
# --------------------------------------------------------------------------- #
def cramers_v_bergsma(table: np.ndarray) -> float | None:
    """Cramér V con correzione del bias di Bergsma (2013). None se non definito
    (meno di 2 livelli effettivi per lato dopo la correzione)."""
    t = np.asarray(table, dtype=float)
    t = t[t.sum(axis=1) > 0][:, t.sum(axis=0) > 0]
    r, k = t.shape
    n = t.sum()
    if r < 2 or k < 2 or n < 2:
        return None
    chi2 = chi2_contingency(t, correction=False)[0]
    phi2 = chi2 / n
    phi2c = max(0.0, phi2 - (k - 1) * (r - 1) / (n - 1))
    rc = r - (r - 1) ** 2 / (n - 1)
    kc = k - (k - 1) ** 2 / (n - 1)
    denom = min(kc - 1, rc - 1)
    if denom <= 0:
        return None
    return float(np.sqrt(phi2c / denom))


def _nested(sheet: pd.DataFrame, a: str, b: str) -> bool:
    """A annidato in B: ogni livello di A compare in un solo livello di B."""
    return bool(sheet.groupby(a, observed=True)[b].nunique().max() == 1)


def _pair_sentence(a: str, b: str, n: int, r: int, c: int, v: float | None, alarm: bool,
                   v_threshold: float, structural: bool = False) -> str:
    if v is not None and structural:
        return (f"'{a}' e '{b}': Cramér V corretto = {v:.2f} su {n} unità. La relazione fra i due "
                f"fattori è strutturale (annidamento o coincidenza, vedi i fatti strutturali): "
                f"l'associazione è riportata lì, non come allarme separato.")
    if v is None:
        return (f"'{a}' e '{b}': Cramér V non valutabile ({n} unità per {r}x{c} celle; servono "
                f"almeno {MIN_ROWS_FOR_V} unità e {MIN_EXPECTED_PER_CELL:g} unità attese per "
                f"cella). Un numero su una tabella così piccola sarebbe instabile.")
    if alarm:
        return (f"'{a}' e '{b}': associazione forte (Cramér V corretto = {v:.2f}, soglia "
                f"{v_threshold:g}, su {n} unità). I due fattori sono in parte confusi: un effetto "
                f"attribuito all'uno può dipendere in parte dall'altro. Il V descrive la "
                f"struttura del disegno, non un effetto biologico.")
    return (f"'{a}' e '{b}': Cramér V corretto = {v:.2f} su {n} unità, sotto la soglia di "
            f"{v_threshold:g}: nessuna segnalazione (non significa che il disegno sia "
            f"perfettamente bilanciato).")


def _pairs_and_findings(
    sheet: pd.DataFrame, roles: dict[str, str], outcome_like: list[str], v_threshold: float,
) -> tuple[list[PairAssociation], list[StructuralFinding]]:
    cols = list(roles)
    findings: list[StructuralFinding] = []
    n = len(sheet)

    constant = [c for c in cols if sheet[c].nunique() < 2]
    identifiers = [c for c in cols if c not in constant and sheet[c].nunique() == n]
    for c in constant:
        findings.append(StructuralFinding(
            "costante", (c,), f"Il fattore '{c}' ha un solo livello: non entra in nessun confronto."))
    for c in identifiers:
        if roles[c] == "patient":
            continue
        findings.append(StructuralFinding(
            "identificatore", (c,),
            f"Il fattore '{c}' ha un livello diverso per ogni riga (identifica il campione): "
            f"non ha senso cercarne l'associazione con gli altri fattori."))

    active = [c for c in cols if c not in constant and c not in identifiers]
    pairs: list[PairAssociation] = []
    seen_structural: set[tuple[str, str]] = set()
    for i in range(len(active)):
        for j in range(i + 1, len(active)):
            a, b = active[i], active[j]
            table = pd.crosstab(sheet[a], sheet[b])
            r, c = table.shape
            evaluable = n >= MIN_ROWS_FOR_V and n / (r * c) >= MIN_EXPECTED_PER_CELL
            v = cramers_v_bergsma(table.values) if evaluable else None
            a_in_b, b_in_a = _nested(sheet, a, b), _nested(sheet, b, a)
            one_to_one = a_in_b and b_in_a
            alarm = v is not None and v >= v_threshold and not (a_in_b or b_in_a)
            pairs.append(PairAssociation(
                factor_a=a, factor_b=b, n_rows=n, n_levels_a=int(sheet[a].nunique()),
                n_levels_b=int(sheet[b].nunique()), cramer_v=v, v_alarm=alarm,
                a_nested_in_b=a_in_b, b_nested_in_a=b_in_a, one_to_one=one_to_one, table=table,
                sentence=_pair_sentence(a, b, n, r, c, v, alarm, v_threshold,
                                        structural=a_in_b or b_in_a),
            ))
            if one_to_one:
                findings.append(StructuralFinding(
                    "uno-a-uno", (a, b),
                    f"Nel disegno attuale '{a}' e '{b}' coincidono: ogni livello dell'uno corrisponde "
                    f"a un solo livello dell'altro e viceversa. Statisticamente sono lo stesso "
                    f"fattore e i loro effetti non sono separabili."))
                seen_structural.add((a, b))
                continue
            for inner, outer in ((a, b), (b, a)):
                if not _nested(sheet, inner, outer):
                    continue
                if outer in outcome_like and roles[inner] != "outcome":
                    extra = (" Il confronto fra i valori di '" + outer + "' può essere fatto solo "
                             "fra pazienti diversi, mai entro lo stesso paziente."
                             if roles[inner] == "patient" else
                             " Un effetto attribuito a '" + outer + "' non è separabile da quello "
                             "di '" + inner + "'.")
                    findings.append(StructuralFinding(
                        "esito-determinato", (inner, outer),
                        f"Nel disegno attuale ogni livello di '{inner}' ha un solo valore di "
                        f"'{outer}': '{outer}' è interamente determinato da '{inner}'." + extra))
                else:
                    findings.append(StructuralFinding(
                        "annidamento", (inner, outer),
                        f"Nel disegno attuale il fattore '{inner}' è annidato in '{outer}': ogni "
                        f"livello di '{inner}' compare in un solo livello di '{outer}'. Un "
                        f"confronto fra livelli di '{outer}' non distingue l'effetto di '{outer}' "
                        f"da quello dei livelli di '{inner}' che contiene."))
    return pairs, findings


def _technical_unit_findings(sheet: pd.DataFrame, roles: dict[str, str], patient_col: str,
                              tissue_col: str | None) -> list[StructuralFinding]:
    """Fattori tecnici in corrispondenza uno-a-uno con la coppia paziente-tessuto (es.
    una libreria per coppia, come in GSE278694)."""
    if tissue_col is None:
        return []
    combo = sheet[patient_col].astype(str) + "|" + sheet[tissue_col].astype(str)
    out = []
    for c, role in roles.items():
        if role not in TECHNICAL_ROLES:
            continue
        tmp = pd.DataFrame({"t": sheet[c].astype(str), "u": combo})
        if tmp.t.nunique() < 2:
            continue
        if _nested(tmp, "t", "u") and _nested(tmp, "u", "t"):
            out.append(StructuralFinding(
                "unità-tecnica", (c,),
                f"Ogni livello di '{c}' corrisponde a una sola coppia paziente-tessuto e viceversa: "
                f"l'effetto tecnico di '{c}' non è separabile da quello della coppia "
                f"paziente-tessuto. È la norma senza multiplexing o librerie replicate, ma "
                f"significa che ogni differenza fra tessuti dello stesso paziente include anche "
                f"la differenza fra due librerie."))
    return out


# --------------------------------------------------------------------------- #
# 2. Confronti stimabili
# --------------------------------------------------------------------------- #
def _technical_confounders(sub: pd.DataFrame, factor: str, roles: dict[str, str]) -> list[str]:
    """Fattori tecnici (non identificativi di riga) che separano PERFETTAMENTE i due
    livelli confrontati nelle righe usate: ogni loro livello contiene un solo livello
    del confronto, con almeno 2 livelli presenti."""
    out = []
    for c, role in roles.items():
        if role not in TECHNICAL_ROLES or c == factor:
            continue
        k = sub[c].nunique()
        if k < 2 or k == len(sub):
            continue
        if _nested(sub, c, factor):
            out.append(c)
    return out


def _partial_confounders(sub: pd.DataFrame, factor: str, roles: dict[str, str],
                         v_threshold: float) -> list[tuple[str, float]]:
    out = []
    for c, role in roles.items():
        if role not in TECHNICAL_ROLES or c == factor:
            continue
        k = sub[c].nunique()
        if k < 2 or k == len(sub):
            continue
        table = pd.crosstab(sub[c], sub[factor])
        r, cc = table.shape
        if len(sub) >= MIN_ROWS_FOR_V and len(sub) / (r * cc) >= MIN_EXPECTED_PER_CELL:
            v = cramers_v_bergsma(table.values)
            if v is not None and v >= v_threshold:
                out.append((c, v))
    return out


def _mann_whitney_min_pvalue(n_a: int, n_b: int) -> float:
    """P-value minimo raggiungibile dal test di Mann-Whitney esatto (bilaterale) con
    gruppi di n_a e n_b pazienti: 2 / C(n_a+n_b, n_a)."""
    if n_a < 1 or n_b < 1:
        return 1.0
    return min(1.0, 2.0 / comb(n_a + n_b, n_a))


def assess_comparison(
    sheet: pd.DataFrame, factor: str, level_a: str, level_b: str, patient_col: str,
    roles: dict[str, str], min_units: int = MIN_UNITS_DEFAULT,
    v_threshold: float = CRAMER_V_THRESHOLD,
) -> ComparisonAssessment:
    """Classifica un confronto 'level_a vs level_b' del fattore ``factor``.

    - Se almeno un paziente ha entrambi i livelli: disegno appaiato, unità = pazienti
      con entrambi. Bassa potenza se meno di ``min_units`` oppure se anche un test di
      Wilcoxon esatto non potrebbe scendere sotto 0.05 (con 5 coppie il minimo è 0.0625).
    - Se nessun paziente ha entrambi i livelli: per tessuto e fattori tecnici il
      confronto coincide con un confronto fra gruppi di pazienti -> non stimabile; per un
      esito a livello di paziente (colonna dichiarata come esito) è un confronto fra
      pazienti, con unità = pazienti per gruppo.
    - In ogni caso, un fattore tecnico che separa perfettamente i due livelli rende il
      confronto non stimabile.
    """
    s = as_levels(sheet)[0]
    level_a, level_b = str(level_a), str(level_b)
    if factor not in s.columns:
        raise ValueError(f"fattore '{factor}' non trovato")
    if factor == patient_col:
        raise ValueError("il fattore del confronto non può essere la colonna paziente")
    rows = s[s[factor].isin([level_a, level_b])]
    pats_a = set(rows.loc[rows[factor] == level_a, patient_col])
    pats_b = set(rows.loc[rows[factor] == level_b, patient_col])
    label = f"'{factor}': {level_a} vs {level_b}"
    reasons: list[str] = []

    if not pats_a or not pats_b:
        missing = level_a if not pats_a else level_b
        return ComparisonAssessment(
            factor, level_a, level_b, "-", 0, len(pats_a), len(pats_b), "non stimabile",
            (f"livello '{missing}' assente",),
            f"Confronto {label}: non stimabile, il livello '{missing}' non compare nei metadati.")

    both = pats_a & pats_b
    role = roles.get(factor, "outcome")
    if both:
        design = "appaiato (entro paziente)"
        units = len(both)
        sub = rows[rows[patient_col].isin(both)]
        excluded = len((pats_a | pats_b) - both)
        if excluded:
            reasons.append(f"{excluded} pazienti con un solo livello esclusi dal confronto appaiato")
        min_p, min_p_test = wilcoxon_min_pvalue(units), "Wilcoxon appaiato"
        low_power = units < min_units or min_p > ALPHA
        power_txt = (f"con {units} pazienti appaiati anche un test di Wilcoxon esatto non può "
                     f"scendere sotto p = {min_p:.3f}" if low_power else "")
    else:
        sub = rows
        if role == "outcome":
            design = "fra pazienti"
            units = min(len(pats_a), len(pats_b))
            floor = _mann_whitney_min_pvalue(len(pats_a), len(pats_b))
            min_p, min_p_test = floor, "Mann-Whitney"
            low_power = units < min_units or floor > ALPHA
            power_txt = (f"con {len(pats_a)} contro {len(pats_b)} pazienti il p-value minimo "
                         f"raggiungibile da un test di Mann-Whitney esatto è {floor:.3f}"
                         if floor > ALPHA else
                         f"meno di {min_units} pazienti nel gruppo più piccolo")
        else:
            reason = (f"nessun paziente ha entrambi i livelli: la differenza fra '{level_a}' e "
                      f"'{level_b}' coincide con la differenza fra due gruppi di pazienti diversi "
                      f"(fattore confuso con il paziente)")
            return ComparisonAssessment(
                factor, level_a, level_b, "-", 0, len(pats_a), len(pats_b), "non stimabile",
                (reason,), f"Confronto {label}: non stimabile -- {reason}. Unità indipendenti: "
                           f"0; p-value minimo non definito (nessun test possibile).")

    tail = (f" Unità indipendenti: {units}; p-value minimo raggiungibile con un test esatto di "
            f"{min_p_test}: {min_p:.3f}. {MIN_PVALUE_NOTE}")
    confounders = _technical_confounders(sub, factor, roles)
    if confounders:
        reason = (f"il fattore tecnico {', '.join(repr(c) for c in confounders)} separa "
                  f"perfettamente i due livelli (fattore confuso con il batch)")
        return ComparisonAssessment(
            factor, level_a, level_b, design, units, len(pats_a), len(pats_b), "non stimabile",
            tuple(reasons + [reason]), f"Confronto {label}: non stimabile -- {reason}." + tail,
            min_p, min_p_test)

    for c, v in _partial_confounders(sub, factor, roles, v_threshold):
        reasons.append(f"'{c}' è fortemente associato al confronto (Cramér V = {v:.2f})")

    unit_word = "pazienti con entrambi i livelli" if design.startswith("appaiato") else \
        "pazienti nel gruppo più piccolo"
    if low_power:
        cls = "stimabile con bassa potenza"
        reasons.append(power_txt)
        sentence = (f"Confronto {label}: stimabile con bassa potenza, disegno {design}, {units} "
                    f"{unit_word} (soglia dichiarata: {min_units} pazienti, la stessa del "
                    f"Modulo B; inoltre {power_txt}). Un risultato non significativo qui non "
                    f"indica assenza di effetto.")
    else:
        cls = "stimabile"
        sentence = (f"Confronto {label}: stimabile, disegno {design}, {units} {unit_word}. "
                    f"L'unità indipendente è il paziente, non la cellula.")
    extra = [r for r in reasons if r != power_txt]
    if extra:
        sentence += " Attenzione: " + "; ".join(extra) + "."
    sentence += tail
    return ComparisonAssessment(factor, level_a, level_b, design, units, len(pats_a),
                                len(pats_b), cls, tuple(reasons), sentence, min_p, min_p_test)


# --------------------------------------------------------------------------- #
# 3. Scomposizione della varianza su pseudobulk
# --------------------------------------------------------------------------- #
def _design_matrix(factors: pd.DataFrame, cols: list[str]) -> np.ndarray:
    parts = [np.ones((len(factors), 1))]
    for c in cols:
        d = pd.get_dummies(factors[c].astype(str), drop_first=True).values.astype(float)
        if d.shape[1]:
            parts.append(d)
    return np.hstack(parts)


def _rss(X: np.ndarray, Y: np.ndarray) -> tuple[np.ndarray, int]:
    beta, *_ = np.linalg.lstsq(X, Y, rcond=None)
    resid = Y - X @ beta
    return (resid ** 2).sum(axis=0), int(np.linalg.matrix_rank(X))


def variance_shares(Y: np.ndarray, factors: pd.DataFrame,
                    fixed_cols: list[str] | None = None) -> dict[str, float] | None:
    """Quota di varianza attribuita a ciascun fattore, stimata con il metodo dei momenti
    sul contributo UNICO del fattore nel modello lineare additivo con tutti i fattori.

    Per ogni gene: SS_f = RSS(modello senza f) - RSS(modello completo), df_f i gradi di
    libertà propri di f, MS_res = RSS completo / df residui. Componente di f:
        C_f = (SS_f - df_f * MS_res) / D_f
    con D_f = N per un fattore FISSO (livelli fissati, es. il tessuto: stima la varianza
    fra i livelli, sum(alpha^2)/k) e D_f = N (L_f - 1) / L_f per un fattore CASUALE con L_f
    livelli campionati da una popolazione (paziente, batch: stima sigma^2_f). In un
    disegno bilanciato questi sono gli stimatori non distorti delle componenti; con
    disegni sbilanciati sono approssimati. Quota di f = media sui geni di C_f (troncata a
    0) divisa per la somma delle componenti medie più MS_res medio.

    Ritorna None se non identificabile: un fattore senza gradi di libertà propri
    (confuso con gli altri) o nessun grado di libertà residuo. ``Y``: unità x geni, su
    scala log."""
    fixed = set(fixed_cols or [])
    cols = list(factors.columns)
    Y = np.asarray(Y, dtype=float)
    n = Y.shape[0]
    X_full = _design_matrix(factors, cols)
    rss_full, rank_full = _rss(X_full, Y)
    df_res = n - rank_full
    if df_res <= 0:
        return None
    ms_res = rss_full / df_res
    comps = {}
    for c in cols:
        X_minus = _design_matrix(factors, [x for x in cols if x != c])
        rss_minus, rank_minus = _rss(X_minus, Y)
        df_f = rank_full - rank_minus
        if df_f <= 0:
            return None
        n_levels = factors[c].nunique()
        divisor = n if c in fixed else n * (n_levels - 1) / n_levels
        comps[c] = max(0.0, float(np.mean((rss_minus - rss_full - df_f * ms_res) / divisor)))
    total = sum(comps.values()) + float(np.mean(ms_res))
    return {c: v / total for c, v in comps.items()}


def _refused(reason: str, n_units: int, n_patients: int, n_genes: int) -> VarianceDecomposition:
    return VarianceDecomposition(
        identifiable=False, reason=reason, n_units=n_units, n_patients=n_patients,
        n_genes=n_genes, shares=[], residual_share=None,
        sentence=f"Scomposizione della varianza non identificabile: {reason}. Nessuna "
                 f"percentuale viene riportata, perché sarebbe arbitraria.")


def variance_decomposition_from_units(
    Y: np.ndarray, units: pd.DataFrame, patient_col: str, factor_cols: list[str],
    fixed_cols: list[str] | None = None,
    n_boot: int = 500, seed: int = 0, min_patients: int = MIN_UNITS_DEFAULT,
    z: float = 1.959964,
) -> VarianceDecomposition:
    """Quote di varianza con intervallo bootstrap ricampionando i PAZIENTI (con tutte le
    loro unità), non le unità: le unità dello stesso paziente non sono indipendenti.

    L'intervallo è normale: stima +- z * SE, con SE = deviazione standard delle stime
    bootstrap. Non si usano i quantili bootstrap perché la distribuzione bootstrap di una
    componente di varianza è spostata verso il basso (i pazienti estratti più volte
    riducono la varianza fra pazienti) di più della propria ampiezza: in calibrazione
    il percentile copriva il vero 0/200 volte e il "basic" lo ribaltava senza coprirlo.
    Lo stimatore puntuale è invece non distorto per costruzione (metodo dei momenti).
    ``fixed_cols``: fattori a livelli fissi (tipicamente il tessuto); gli altri sono
    trattati come casuali."""
    Y = np.asarray(Y, dtype=float)
    units = units.reset_index(drop=True).astype(str)
    n_units, n_genes = Y.shape
    patients = units[patient_col].unique()
    n_pat = len(patients)
    if n_pat < min_patients:
        return _refused(f"{n_pat} pazienti, ne servono almeno {min_patients} per un intervallo",
                        n_units, n_pat, n_genes)
    point = variance_shares(Y, units[factor_cols], fixed_cols)
    if point is None:
        return _refused("almeno un fattore coincide con altri fattori o con le unità stesse "
                        "(nessun grado di libertà proprio o residuo)", n_units, n_pat, n_genes)

    rng = np.random.default_rng(seed)
    idx_by_pat = {p: np.flatnonzero(units[patient_col].values == p) for p in patients}
    boot = {c: [] for c in factor_cols}
    n_fail = 0
    for _ in range(n_boot):
        draw = rng.choice(patients, size=n_pat, replace=True)
        idx, new_ids = [], []
        for i, p in enumerate(draw):
            idx.append(idx_by_pat[p])
            new_ids += [f"{p}#{i}"] * len(idx_by_pat[p])
        idx = np.concatenate(idx)
        f = units.iloc[idx].reset_index(drop=True).copy()
        f[patient_col] = new_ids
        s = variance_shares(Y[idx], f[factor_cols], fixed_cols)
        if s is None:
            n_fail += 1
            continue
        for c in factor_cols:
            boot[c].append(s[c])
    if n_fail > 0.1 * n_boot:
        return _refused(f"{n_fail} ricampionamenti bootstrap su {n_boot} non identificabili: "
                        f"il disegno è troppo fragile per un intervallo", n_units, n_pat, n_genes)

    shares = []
    for c in factor_cols:
        se = float(np.std(boot[c], ddof=1))
        lo = max(0.0, point[c] - z * se)
        hi = min(1.0, point[c] + z * se)
        shares.append(VarianceShare(c, point[c], float(lo), float(hi)))
    residual = max(0.0, 1.0 - sum(point.values()))
    parts = ", ".join(f"'{s.factor}' {s.share:.0%} (IC95% {s.ci_low:.0%}-{s.ci_high:.0%})"
                      for s in shares)
    sentence = (f"Quota della variabilità dell'espressione (pseudobulk, {n_units} unità, "
                f"{n_pat} pazienti, {n_genes} geni ad alta varianza) attribuibile in modo "
                f"univoco a ciascun fattore: {parts}; il resto ({residual:.0%}) è variabilità "
                f"residua o condivisa fra fattori. Gli intervalli ricampionano i pazienti; con "
                f"pochi pazienti sono larghi, ed è questa la loro informazione principale.")
    return VarianceDecomposition(True, "", n_units, n_pat, n_genes, shares, residual, sentence)


def pseudobulk(adata: ad.AnnData, unit_cols: list[str],
               min_cells: int = MIN_CELLS_PER_UNIT_DEFAULT) -> tuple[np.ndarray, pd.DataFrame, int]:
    """Somma dei conteggi per combinazione distinta di ``unit_cols``. Scarta le unità
    con meno di ``min_cells`` cellule e ne ritorna il numero."""
    obs = adata.obs[unit_cols].astype(str).reset_index(drop=True)
    key = obs.agg("|".join, axis=1)
    codes, uniques = pd.factorize(key)
    sizes = np.bincount(codes)
    keep = np.flatnonzero(sizes >= min_cells)
    ind = sp.csr_matrix((np.ones(len(codes)), (codes, np.arange(len(codes)))),
                        shape=(len(uniques), len(codes)))
    X = adata.X if sp.issparse(adata.X) else sp.csr_matrix(adata.X)
    summed = (ind @ X)[keep]
    first = obs.groupby(codes).first().iloc[keep].reset_index(drop=True)
    return summed, first, int(len(uniques) - len(keep))


def variance_decomposition(
    adata: ad.AnnData, patient_col: str, tissue_col: str | None = None,
    batch_cols: list[str] | None = None, n_hvg: int = N_HVG_DEFAULT,
    min_cells: int = MIN_CELLS_PER_UNIT_DEFAULT, n_boot: int = 500, seed: int = 0,
    min_patients: int = MIN_UNITS_DEFAULT,
) -> VarianceDecomposition:
    """Pseudobulk per combinazione paziente-tessuto(-batch), CPM + log1p (``LogCPM`` del
    Modulo A, con target 1e6), geni ad alta varianza, quote per fattore."""
    factor_cols = [patient_col] + ([tissue_col] if tissue_col else []) + list(batch_cols or [])
    summed, units, n_dropped = pseudobulk(adata, factor_cols, min_cells=min_cells)
    Y = LogCPM(target_sum=1e6).fit_transform(summed)
    Y = Y.toarray() if sp.issparse(Y) else np.asarray(Y)
    var = Y.var(axis=0)
    order = np.argsort(var)[::-1]
    order = order[var[order] > 0][:n_hvg]
    dec = variance_decomposition_from_units(
        Y[:, order], units, patient_col=patient_col, factor_cols=factor_cols,
        fixed_cols=[tissue_col] if tissue_col else [], n_boot=n_boot, seed=seed, min_patients=min_patients)
    if n_dropped and dec.identifiable:
        dec = VarianceDecomposition(
            dec.identifiable, dec.reason, dec.n_units, dec.n_patients, dec.n_genes, dec.shares,
            dec.residual_share,
            dec.sentence + f" {n_dropped} unità con meno di {min_cells} cellule sono state escluse.")
    return dec


# --------------------------------------------------------------------------- #
# Orchestratore
# --------------------------------------------------------------------------- #
def run_design_audit(
    sheet: pd.DataFrame,
    patient_col: str,
    tissue_col: str | None = None,
    technical_cols: dict[str, str] | None = None,
    outcome_cols: list[str] | None = None,
    comparisons: list[tuple[str, str, str]] | None = None,
    min_units: int = MIN_UNITS_DEFAULT,
    v_threshold: float = CRAMER_V_THRESHOLD,
) -> DesignAuditResult:
    """Esegue l'audit del disegno.

    ``sheet``: una riga per campione/libreria (usa ``sample_sheet_from_obs`` per partire
    da ``adata.obs``). ``technical_cols``: ruolo -> colonna, ruoli in
    ``TECHNICAL_ROLES``. ``comparisons``: lista di (colonna, livello_a, livello_b).

    La scomposizione della varianza (``variance_decomposition``) NON è inclusa qui: la
    calibrazione del suo intervallo è fuori banda per la quota del tessuto (vedi
    tests/test_design_audit_calibration.py, xfail) e per regola non viene esposta.
    """
    technical_cols = dict(technical_cols or {})
    outcome_cols = list(outcome_cols or [])
    bad = set(technical_cols) - set(TECHNICAL_ROLES)
    if bad:
        raise ValueError(f"ruoli tecnici non riconosciuti: {sorted(bad)}; ammessi: {TECHNICAL_ROLES}")

    roles: dict[str, str] = {patient_col: "patient"}
    if tissue_col:
        roles[tissue_col] = "tissue"
    for role, col in technical_cols.items():
        roles.setdefault(col, role)
    for col in outcome_cols:
        roles.setdefault(col, "outcome")
    missing = [c for c in roles if c not in sheet.columns]
    if missing:
        raise ValueError(f"colonne non trovate nei metadati: {missing}")

    s, missing = as_levels(sheet[list(roles)].reset_index(drop=True))
    outcome_like = ([tissue_col] if tissue_col else []) + outcome_cols
    pairs, findings = _pairs_and_findings(s, roles, outcome_like, v_threshold)
    findings += _technical_unit_findings(s, roles, patient_col, tissue_col)

    comps = [assess_comparison(s, f, a, b, patient_col, roles, min_units, v_threshold)
             for f, a, b in (comparisons or [])]

    notes = []
    if missing:
        notes.append("Valori mancanti trattati come livello esplicito \"NA\": " + ", ".join(
            f"'{c}' ({n} righe)" for c, n in missing.items()) + ". Un livello \"NA\" condiviso "
            "da più righe è trattato come un valore uguale per tutte.")
    if "protocol" in technical_cols:
        notes.append(PROTOCOL_NOTE)
    notes.append("Le unità di questo audit sono le righe della tabella dei metadati (campioni o "
                 "librerie), non le cellule: le cellule dello stesso campione non sono "
                 "osservazioni indipendenti del disegno.")

    structural = [f for f in findings if f.kind in ("annidamento", "uno-a-uno", "esito-determinato",
                                                     "unità-tecnica")]
    n_alarm = sum(p.v_alarm for p in pairs)
    parts = [f"Audit del disegno su {len(s)} unità (righe dei metadati) e {len(roles)} fattori."]
    parts.append(f"{len(structural)} fatti strutturali rilevati (annidamenti, coincidenze, esiti "
                 f"determinati da un fattore)" + (":" if structural else "."))
    parts += [f.sentence for f in structural]
    parts.append(f"{n_alarm} coppie di fattori con associazione forte (Cramér V >= {v_threshold:g}).")
    parts += [c.sentence for c in comps]
    return DesignAuditResult(n_rows=len(s), roles=roles, pairs=pairs, findings=findings,
                             comparisons=comps, notes=notes,
                             narrative=" ".join(parts))
