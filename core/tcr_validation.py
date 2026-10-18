"""Modulo B -- Verifica dell'annotazione di tipo cellulare via TCR.

Logica di discordanza a coppie portata da pdac-ml/src/07_robustness.py e pdac-ml/src/10_loco.py,
generalizzata: nomi di colonna (paziente, compartimento, tipo cellulare, barcode) configurabili,
nessuna dipendenza dal dataset PDAC originale. Ogni intervallo esposto è l'intervallo t
cluster-robusto sui pazienti di ``core.stats.cluster_t_interval`` (revisione, A1).

Ipotesi biologica da cui dipendono TUTTI i risultati del modulo: le cellule dello stesso clone T
(stessa sequenza CDR3 della catena TRB) hanno la stessa linea CD4/CD8. Eccezioni note, che il
modulo non distingue da un errore di annotazione:
- cellule CD4 che esprimono CD8αα (o CD8 su cellule attivate);
- cellule T CD4 citotossiche, che esprimono programmi tipici delle CD8 (probabile causa reale di
  parte delle chiamate CD8 sulle vere CD4 nel tumore di GSE278694);
- MAIT, NKT e cellule γδ: un clonotipo definito dalla sola TRB non le separa in modo affidabile
  (le γδ non hanno TRB; MAIT e NKT hanno catene α semi-invarianti che qui non sono usate);
- doppietti T-NK o T-mieloidi, che portano il TCR in una cellula con un'altra etichetta;
- RNA ambientale specifico del compartimento, che sposta i marcatori fra compartimenti.
L'eccesso di discordanza è interpretabile solo con etichette di LINEA (CD4/CD8), non di STATO
(es. effettrici contro esauste): lo stesso clone cambia legittimamente stato fra sangue e tumore.

Estimatore D_A/D_B (discordanza entro compartimento, termine di confronto dell'eccesso): versione
SENZA reinserimento (coppie di cellule distinte dello stesso gruppo), non distorta
(07_robustness.py, punto A); la versione con proporzioni (con reinserimento) sovrastima D sui
gruppi piccoli. D(A×B) fra compartimenti resta con proporzioni: confronta due gruppi diversi
(una cellula per lato), senza problemi di reinserimento.

Eccesso pesato (stima principale) e a pesi uguali (sensibilità). L'eccesso principale è
D(A×B) - [n_A D_A + n_B D_B] / (n_A + n_B). Sotto H0 (stesse proporzioni di etichette nei due
compartimenti) ha media zero, ma NON è una divergenza non negativa: con n_A molto diverso da n_B
il suo valore atteso può essere negativo anche quando le etichette differiscono davvero (es. 3
cellule CD8 nel sangue e 60 nel tumore con il 40% di errore: circa -0,06). Per questo è riportata
anche la versione a pesi uguali, D(A×B) - (D_A + D_B) / 2, il cui valore atteso è ½‖p_A − p_B‖²:
zero sotto H0 e positivo appena le proporzioni differiscono, a prescindere da n_A e n_B.
"""

from __future__ import annotations

import re
import warnings
from dataclasses import dataclass
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse as sp

from core.stats import ClusterInterval, cluster_t_interval

MIN_CELLS_DEFAULT = 2
MIN_PATIENTS_FOR_CI_DEFAULT = 5
MIN_CELLS_MARKER_DEFAULT = 3
MIN_MARGIN_DEFAULT = 0.20

BIOLOGICAL_ASSUMPTIONS = (
    "Ipotesi: le cellule dello stesso clone T (stessa CDR3 della catena TRB) hanno la stessa linea "
    "CD4/CD8. Eccezioni note che il modulo non distingue da un errore di annotazione: cellule CD4 con "
    "CD8αα, CD4 citotossiche, MAIT, NKT e cellule γδ (non separabili da un clonotipo basato solo sulla "
    "TRB), doppietti T-NK, RNA ambientale specifico del compartimento. L'eccesso di discordanza è "
    "interpretabile solo con etichette di linea (CD4/CD8), non di stato (es. effettrici contro esauste)."
)


# --------------------------------------------------------------------------- #
# Ingestione VDJ (Cell Ranger)
# --------------------------------------------------------------------------- #
def parse_vdj_contigs(files: list[tuple]) -> pd.DataFrame:
    """``files``: lista di (percorso_csv, patient_id, compartment) oppure (percorso_csv,
    patient_id, compartment, library). L'associazione file -> paziente -> compartimento (->
    libreria) è fornita dal chiamante (CLI/app), non dedotta dal nome del file: i CSV Cell Ranger
    non contengono questa informazione e non c'è una convenzione di nome universale.

    Più file sulla stessa coppia paziente-compartimento sono ammessi solo se ciascuno indica la
    propria libreria: senza, lo stesso barcode in due librerie sarebbe attribuito alla cellula
    sbagliata, quindi è un errore. Le librerie vanno indicate per tutti i file o per nessuno.

    Filtra is_cell/high_confidence/productive/full_length (se le colonne sono presenti) e le
    catene diverse da TRA/TRB. Per (paziente, compartimento, libreria, barcode, catena) tiene un
    solo contig, con una regola che non dipende dall'ordine delle righe: più UMI, poi più reads
    (se la colonna c'è), poi la sequenza CDR3 nucleotidica minore in ordine alfabetico. La
    colonna ``library`` vale "" se le librerie non sono indicate.
    """
    if not files:
        raise ValueError("nessun file VDJ fornito")
    with_lib = {len(f) >= 4 and f[3] not in (None, "") for f in files}
    if len(with_lib) > 1:
        raise ValueError("la libreria VDJ va indicata per tutti i file del manifest o per nessuno")
    seen: dict[tuple[str, str, str], str] = {}
    frames = []
    for f in files:
        path, patient, compartment = f[0], str(f[1]), str(f[2])
        library = str(f[3]) if len(f) >= 4 and f[3] not in (None, "") else ""
        key = (patient, compartment, library)
        if key in seen:
            raise ValueError(
                f"due file VDJ ({seen[key]} e {path}) per paziente '{patient}', compartimento "
                f"'{compartment}'" + (f", libreria '{library}'" if library else "") + ": lo stesso "
                "barcode potrebbe comparire in entrambi. Aggiungi al manifest una colonna 'library' "
                "e indica la colonna libreria dei metadati (--library-col).")
        seen[key] = str(path)
        df = pd.read_csv(path, low_memory=False)
        for col in ("is_cell", "high_confidence", "productive", "full_length"):
            if col in df.columns:
                df = df[df[col].astype(str).str.lower().isin(["true", "1"])]
        required = {"barcode", "chain", "cdr3_nt"}
        missing = required - set(df.columns)
        if missing:
            raise ValueError(f"{path}: colonne mancanti nel CSV VDJ: {sorted(missing)}")
        df = df[df.chain.isin(["TRA", "TRB"])]
        df = df[df.cdr3_nt.notna() & (df.cdr3_nt.astype(str) != "None")]
        df = df.assign(patient=patient, compartment=compartment, library=library,
                       barcode=df["barcode"].astype(str), cdr3_nt=df["cdr3_nt"].astype(str))
        if "umis" not in df.columns:
            df = df.assign(umis=1)
        df = df.assign(reads=df["reads"] if "reads" in df.columns else 0)
        frames.append(df[["patient", "compartment", "library", "barcode", "chain", "cdr3_nt",
                          "umis", "reads"]])
    c = pd.concat(frames, ignore_index=True)
    c = c.sort_values(["umis", "reads", "cdr3_nt"], ascending=[False, False, True], kind="mergesort")
    c = c.drop_duplicates(["patient", "compartment", "library", "barcode", "chain"])
    return c.drop(columns="reads").reset_index(drop=True)


def build_clonotypes(contigs: pd.DataFrame) -> pd.DataFrame:
    """Il clonotipo è definito dalla sequenza nucleotidica CDR3 della catena TRB. Una riga per
    cellula con TRB rilevato (``parse_vdj_contigs`` ne tiene già una sola per cellula)."""
    if "library" not in contigs.columns:
        contigs = contigs.assign(library="")
    trb = contigs[contigs.chain == "TRB"][["patient", "compartment", "library", "barcode", "cdr3_nt"]]
    trb = trb.rename(columns={"cdr3_nt": "clone_id"})
    return trb.drop_duplicates(["patient", "compartment", "library", "barcode"]).reset_index(drop=True)


# --------------------------------------------------------------------------- #
# Controllo del match dei barcode
# --------------------------------------------------------------------------- #
MIN_BARCODE_MATCH = 0.5
_SUFFIX_RE = re.compile(r"-\d+$")
_KEY = ["patient", "compartment", "library", "barcode"]


@dataclass(frozen=True)
class BarcodeMatch:
    """Esito della corrispondenza fra barcode VDJ e metadati."""
    n_vdj: int
    n_matched: int
    fraction: float
    normalized_suffixes: tuple[str, ...]  # suffissi rimossi; vuoto se nessuna normalizzazione
    sentence: str


def _match_fraction(obs: pd.DataFrame, clones: pd.DataFrame) -> tuple[int, float]:
    if clones.empty:
        return 0, 0.0
    n = len(clones.merge(obs[_KEY].drop_duplicates(), on=_KEY, how="inner"))
    return n, n / len(clones)


def _strip(b: pd.Series) -> pd.Series:
    return b.astype(str).str.replace(_SUFFIX_RE, "", regex=True)


def _duplicate_key_error(obs: pd.DataFrame) -> None:
    dup = obs.duplicated(_KEY, keep=False)
    if dup.any():
        ex = obs.loc[dup, _KEY].head(3).to_dict(orient="records")
        raise ValueError(
            f"{int(dup.sum())} righe dei metadati condividono paziente, compartimento"
            + (", libreria" if (obs["library"] != "").any() else "") + f" e barcode (es. {ex}): "
            "il TCR non può essere attribuito a una sola cellula. Succede con più librerie per la "
            "stessa coppia paziente-compartimento: indica la colonna libreria (--library-col) e la "
            "colonna 'library' nel manifest VDJ.")


def match_barcodes(obs: pd.DataFrame, clones: pd.DataFrame,
                   min_fraction: float = MIN_BARCODE_MATCH) -> tuple[pd.DataFrame, pd.DataFrame, BarcodeMatch]:
    """Frazione di barcode VDJ (cellule con TRB) ritrovati nei metadati, per (paziente,
    compartimento, libreria, barcode). Righe dei metadati con la stessa chiave sono un errore.
    Se rimuovere il suffisso "-N" (es. "-1" di Cell Ranger) da entrambi i lati aumenta il match
    senza creare collisioni (barcode diversi che diventano uguali, es. -1 e -2 dello stesso
    campione aggregato), la normalizzazione viene applicata e dichiarata. Errore esplicito se il
    match finale è sotto ``min_fraction``."""
    obs = obs if "library" in obs.columns else obs.assign(library="")
    clones = clones if "library" in clones.columns else clones.assign(library="")
    _duplicate_key_error(obs)
    n_raw, f_raw = _match_fraction(obs, clones)
    used_obs, used_clones, suffixes = obs, clones, ()
    obs_n = obs.assign(barcode=_strip(obs["barcode"]))
    clones_n = clones.assign(barcode=_strip(clones["barcode"]))
    n_norm, f_norm = _match_fraction(obs_n, clones_n)
    collisions = (clones_n.duplicated(_KEY).sum() > clones.duplicated(_KEY).sum()
                  or obs_n.duplicated(_KEY).sum() > obs.duplicated(_KEY).sum())
    if f_norm > f_raw and not collisions:
        found = pd.concat([obs["barcode"], clones["barcode"]]).astype(str).str.extract(r"(-\d+)$")[0]
        suffixes = tuple(sorted(found.dropna().unique()))
        used_obs, used_clones = obs_n, clones_n
        n_final, f_final = n_norm, f_norm
    else:
        n_final, f_final = n_raw, f_raw
    sentence = (f"Barcode VDJ (cellule con TRB) ritrovati nei metadati per paziente, compartimento e "
                f"barcode: {n_final:,}/{len(clones):,} ({f_final:.1%}).")
    if suffixes:
        n_v = int(clones["barcode"].astype(str).str.contains(_SUFFIX_RE).sum())
        n_o = int(obs["barcode"].astype(str).str.contains(_SUFFIX_RE).sum())
        sentence += (f" Normalizzato il suffisso {', '.join(suffixes)}: rimosso da {n_v:,} barcode VDJ e da "
                     f"{n_o:,} barcode dei metadati (prima della normalizzazione: {f_raw:.1%}).")
    elif f_norm > f_raw and collisions:
        sentence += (" La rimozione del suffisso -N avrebbe aumentato il match ma avrebbe fuso barcode "
                     "diversi: non applicata.")
    if f_final < min_fraction:
        ex_v = clones["barcode"].astype(str).head(3).tolist()
        ex_o = obs["barcode"].astype(str).head(3).tolist()
        raise ValueError(
            f"solo il {f_final:.1%} dei barcode VDJ è stato ritrovato nei metadati (soglia "
            f"{min_fraction:.0%}), anche dopo l'eventuale normalizzazione del suffisso -N. Esempi "
            f"VDJ: {ex_v}; esempi metadati: {ex_o}. Controlla che paziente e compartimento del "
            f"manifest coincidano con quelli dell'AnnData e che i barcode abbiano lo stesso formato.")
    return used_obs, used_clones, BarcodeMatch(len(clones), n_final, f_final, suffixes, sentence)


# --------------------------------------------------------------------------- #
# (A) Discordanza a coppie
# --------------------------------------------------------------------------- #
def within_compartment_discordance(labels: np.ndarray) -> float:
    """D_A = 1 - sum(n_k(n_k-1)) / (n(n-1)): probabilità che DUE cellule distinte dello
    stesso gruppo (senza reinserimento) abbiano etichetta diversa. Stimatore non distorto
    (07_robustness.py, punto A) -- la versione con proporzioni (1-sum(p_i^2)) sovrastima
    sistematicamente D su gruppi piccoli."""
    _, counts = np.unique(labels, return_counts=True)
    n = int(counts.sum())
    if n < 2:
        return float("nan")
    return 1.0 - float((counts * (counts - 1)).sum() / (n * (n - 1)))


def _label_proportions(labels: np.ndarray, categories: np.ndarray) -> np.ndarray:
    counts = pd.Series(labels).value_counts().reindex(categories, fill_value=0).values
    total = counts.sum()
    return counts / total if total else counts.astype(float)


def cross_compartment_discordance(p_a: np.ndarray, p_b: np.ndarray) -> float:
    """D(A×B) = 1 - sum(p_A_i * p_B_i): una cellula per lato, gruppi diversi -> nessun
    problema di reinserimento, la versione con proporzioni è corretta qui."""
    return 1.0 - float(np.dot(p_a, p_b))


def pairwise_excess_discordance(
    cells: pd.DataFrame,
    patient_col: str = "patient",
    compartment_col: str = "compartment",
    clone_col: str = "clone_id",
    label_col: str = "celltype",
    min_cells: int = MIN_CELLS_DEFAULT,
) -> pd.DataFrame:
    """Una riga per (paziente, clone, coppia di compartimenti) con >= min_cells cellule per lato.
    ``excess`` = D(A×B) - [n_A D_A + n_B D_B] / (n_A + n_B) (stima principale);
    ``excess_equal`` = D(A×B) - (D_A + D_B) / 2 (sensibilità a pesi uguali, valore atteso
    ½‖p_A − p_B‖² >= 0). Vedi il docstring del modulo per la differenza fra le due."""
    categories = np.sort(cells[label_col].astype(str).unique())
    rows = []
    for (pat, clone), g in cells.groupby([patient_col, clone_col], observed=True):
        by_comp = {
            comp: sub[label_col].astype(str).values
            for comp, sub in g.groupby(compartment_col, observed=True)
            if len(sub) >= min_cells
        }
        comps = sorted(by_comp)
        for i in range(len(comps)):
            for j in range(i + 1, len(comps)):
                a_labels, b_labels = by_comp[comps[i]], by_comp[comps[j]]
                n_a, n_b = len(a_labels), len(b_labels)
                p_a = _label_proportions(a_labels, categories)
                p_b = _label_proportions(b_labels, categories)
                d_across = cross_compartment_discordance(p_a, p_b)
                d_a = within_compartment_discordance(a_labels)
                d_b = within_compartment_discordance(b_labels)
                d_within = (n_a * d_a + n_b * d_b) / (n_a + n_b)
                rows.append({
                    "patient": pat, "clone": clone, "comp_a": comps[i], "comp_b": comps[j],
                    "n_a": n_a, "n_b": n_b, "d_across": d_across, "d_within": d_within,
                    "excess": d_across - d_within,
                    "excess_equal": d_across - (d_a + d_b) / 2,
                })
    return pd.DataFrame(rows, columns=["patient", "clone", "comp_a", "comp_b", "n_a", "n_b",
                                       "d_across", "d_within", "excess", "excess_equal"])


@dataclass(frozen=True)
class DiscordanceResult:
    """Eccesso di discordanza complessivo (stima principale pesata), sensibilità a pesi uguali
    e righe per coppia di compartimenti."""
    n_pairs: int
    n_patients: int
    mean_excess: float
    ci_low: float
    ci_high: float
    sufficient: bool
    by_compartment_pair: pd.DataFrame
    equal_weights: ClusterInterval | None = None


def _summarize_discordance(pw: pd.DataFrame, min_patients_for_ci: int) -> DiscordanceResult:
    if pw.empty:
        return DiscordanceResult(n_pairs=0, n_patients=0, mean_excess=float("nan"),
                                 ci_low=float("nan"), ci_high=float("nan"), sufficient=False,
                                 by_compartment_pair=pd.DataFrame(), equal_weights=None)
    overall = cluster_t_interval(pw["excess"].values, pw["patient"].values, min_patients_for_ci)
    equal = cluster_t_interval(pw["excess_equal"].values, pw["patient"].values, min_patients_for_ci)
    rows = []
    for (a, b), g in pw.groupby(["comp_a", "comp_b"]):
        res = cluster_t_interval(g["excess"].values, g["patient"].values, min_patients_for_ci)
        rows.append({"comp_a": a, "comp_b": b, "n_pairs": len(g), "n_patients": res.n_groups,
                     "excess": res.mean, "ci_low": res.ci_low, "ci_high": res.ci_high,
                     "sufficient": res.sufficient, "excess_equal": float(g["excess_equal"].mean())})
    return DiscordanceResult(n_pairs=len(pw), n_patients=int(pw.patient.nunique()),
                             mean_excess=overall.mean, ci_low=overall.ci_low,
                             ci_high=overall.ci_high, sufficient=overall.sufficient,
                             by_compartment_pair=pd.DataFrame(rows), equal_weights=equal)


# --------------------------------------------------------------------------- #
# (B, opzionale) Tasso d'errore per compartimento via marcatori canonici
# --------------------------------------------------------------------------- #
def _check_counts(values: np.ndarray, what: str) -> None:
    if values.size and (values < 0).any():
        raise ValueError(f"{what}: valori negativi; servono conteggi grezzi (interi >= 0)")
    if values.size and not np.allclose(values, np.round(values)):
        raise ValueError(f"{what}: valori non interi (dati normalizzati o log-trasformati?); la "
                         "positività ai marcatori va calcolata sui conteggi grezzi. Usa il layer dei "
                         "conteggi (es. adata.layers['counts']) come adata.X.")


def compute_marker_positivity(
    adata: ad.AnnData, marker_map: dict[str, list[str]], gene_col: str | None = None
) -> pd.DataFrame:
    """Una colonna booleana 'pos_<etichetta>' per etichetta in marker_map: True se almeno uno
    dei geni marcatori di quell'etichetta ha conta > 0 in quella cellula. Stesso ordine di riga
    di adata.obs. Le colonne dei geni marcatori devono contenere conteggi interi non negativi:
    altrimenti errore (su dati normalizzati "> 0" non significa "molecola rilevata")."""
    gene_names = (adata.var[gene_col] if gene_col else adata.var_names).astype(str).values
    gene_to_idx = {}
    for i, g in enumerate(gene_names):
        gene_to_idx.setdefault(g, i)
    X = adata.X

    out = {}
    for label, genes in marker_map.items():
        idx = [gene_to_idx[g] for g in genes if g in gene_to_idx]
        missing = [g for g in genes if g not in gene_to_idx]
        if missing:
            warnings.warn(f"geni marcatori non trovati per '{label}': {missing}", stacklevel=2)
        if not idx:
            out[f"pos_{label}"] = np.zeros(adata.n_obs, dtype=bool)
            continue
        sub = X[:, idx]
        sub = sub.toarray() if sp.issparse(sub) else np.asarray(sub)
        _check_counts(sub, f"geni marcatori di '{label}'")
        out[f"pos_{label}"] = (sub > 0).any(axis=1)
    return pd.DataFrame(out, index=adata.obs_names)


def assign_reference_identity(
    cells_with_markers: pd.DataFrame,
    marker_labels: list[str],
    patient_col: str = "patient",
    clone_col: str = "clone_id",
    compartment_col: str = "compartment",
    reference_compartment: str = "PBMC",
    min_cells: int = MIN_CELLS_MARKER_DEFAULT,
    min_margin: float = MIN_MARGIN_DEFAULT,
) -> pd.DataFrame:
    """Identità di riferimento del clone dai marcatori, stimata SOLO dal compartimento
    di riferimento (tipicamente il sangue): niente circolarità con i compartimenti poi
    giudicati. Cloni ambigui (margine fra l'etichetta vincente e la seconda < min_margin,
    pareggi compresi) sono scartati, non forzati."""
    ref = cells_with_markers[cells_with_markers[compartment_col] == reference_compartment]
    rows = []
    for (pat, clone), g in ref.groupby([patient_col, clone_col], observed=True):
        if len(g) < min_cells:
            continue
        fracs = {label: float(g[f"pos_{label}"].mean()) for label in marker_labels}
        ranked = sorted(fracs.items(), key=lambda kv: kv[1], reverse=True)
        best_label, best_frac = ranked[0]
        second_frac = ranked[1][1] if len(ranked) > 1 else 0.0
        if best_frac - second_frac < min_margin:
            continue
        rows.append({patient_col: pat, clone_col: clone, "reference_identity": best_label,
                     "n_ref_cells": len(g)})
    return pd.DataFrame(rows, columns=[patient_col, clone_col, "reference_identity", "n_ref_cells"])


@dataclass(frozen=True)
class MarkerErrorResult:
    """Tasso d'errore per cellula (convenzione "cell") per compartimento giudicato."""
    reference_compartment: str
    n_resolved_clones: int
    by_compartment: dict[str, ClusterInterval]


def _judged_cells(cells: pd.DataFrame, reference_identity: pd.DataFrame, marker_labels: list[str],
                  reference_compartment: str, patient_col: str = "patient", clone_col: str = "clone_id",
                  compartment_col: str = "compartment", label_col: str = "celltype") -> pd.DataFrame:
    """Cellule giudicate dalla convenzione per cellula: clone con identità di riferimento,
    compartimento diverso dal riferimento, etichetta fra quelle della marker_map. Il filtro è
    sulle etichette della mappa, non su quelle risolte: un'etichetta della mappa che nessun clone
    ha come identità resta giudicata (ed è un errore), invece di uscire dal calcolo."""
    j = cells.merge(reference_identity, on=[patient_col, clone_col], how="inner")
    j = j[j[compartment_col] != reference_compartment]
    j = j[j[label_col].astype(str).isin([str(x) for x in marker_labels])]
    return j.assign(error=(j[label_col].astype(str) != j["reference_identity"].astype(str)))


def marker_error_rate(
    cells: pd.DataFrame,
    reference_identity: pd.DataFrame,
    marker_labels: list[str],
    patient_col: str = "patient",
    clone_col: str = "clone_id",
    compartment_col: str = "compartment",
    label_col: str = "celltype",
    reference_compartment: str = "PBMC",
    min_patients_for_ci: int = MIN_PATIENTS_FOR_CI_DEFAULT,
) -> MarkerErrorResult:
    """Confronta l'etichetta assegnata con l'identità di riferimento del clone (dai marcatori,
    stimata solo nel compartimento di riferimento) su tutti gli ALTRI compartimenti. Tasso
    d'errore per compartimento con intervallo cluster-robusto sui pazienti (scala logit).
    Sono giudicate le cellule con etichetta in ``marker_labels`` (le chiavi della marker_map)."""
    if reference_identity.empty:
        return MarkerErrorResult(reference_compartment=reference_compartment,
                                 n_resolved_clones=0, by_compartment={})
    j = _judged_cells(cells, reference_identity, marker_labels, reference_compartment,
                      patient_col, clone_col, compartment_col, label_col)
    by_compartment: dict[str, ClusterInterval] = {}
    for comp, g in j.groupby(compartment_col, observed=True):
        by_compartment[str(comp)] = cluster_t_interval(
            g["error"].values.astype(float), g[patient_col].values, min_patients_for_ci, scale="logit")
    return MarkerErrorResult(reference_compartment=reference_compartment,
                             n_resolved_clones=len(reference_identity),
                             by_compartment=by_compartment)


# --------------------------------------------------------------------------- #
# (B2) Convenzioni del tasso d'errore: "per cellula" e "per clone" (tesi)
# --------------------------------------------------------------------------- #
CONVENTIONS = ("cell", "clone")


@dataclass(frozen=True)
class ConventionResult:
    """Tasso d'errore per compartimento secondo una convenzione dichiarata, più le
    differenze appaiate fra compartimenti (stessi cloni)."""
    name: str
    definition: str
    by_compartment: dict[str, ClusterInterval]
    paired_differences: dict[tuple[str, str], tuple[ClusterInterval, int]]


def _ordered_compartments(present: list[str], order: list[str] | None) -> list[str]:
    """Ordine dei compartimenti per le differenze appaiate: quello dichiarato dall'utente
    (``order``), poi gli altri in ordine alfabetico. Nessun compartimento ha un posto fisso."""
    order = [str(x) for x in (order or [])]
    first = [c for c in order if c in present]
    return first + sorted(c for c in present if c not in first)


def _paired_differences(units: pd.DataFrame, min_patients_for_ci: int,
                        compartment_order: list[str] | None = None,
                        ) -> dict[tuple[str, str], tuple[ClusterInterval, int]]:
    """``units``: colonne patient, clone, compartment, error (una riga per unità).
    Per ogni coppia di compartimenti giudicati (A prima di B secondo ``compartment_order``):
    errore medio del clone in A meno quello in B, sui cloni presenti in entrambi; intervallo
    cluster-robusto sui pazienti. Valore: (risultato, numero di cloni)."""
    out = {}
    if units.empty:
        return out
    per = units.groupby(["patient", "clone", "compartment"], observed=True).error.mean().unstack("compartment")
    comps = _ordered_compartments([str(c) for c in per.columns], compartment_order)
    for i in range(len(comps)):
        for j in range(i + 1, len(comps)):
            a, b = comps[i], comps[j]
            d = per.dropna(subset=[a, b]).reset_index()
            if d.empty:
                continue
            res = cluster_t_interval((d[a] - d[b]).values, d["patient"].values, min_patients_for_ci)
            out[(a, b)] = (res, len(d))
    return out


def _cell_units(cells: pd.DataFrame, reference_identity: pd.DataFrame, marker_labels: list[str],
                reference_compartment: str) -> pd.DataFrame:
    """Unità della convenzione "per cellula": stessi filtri di ``marker_error_rate``."""
    if reference_identity.empty:
        return pd.DataFrame(columns=["patient", "clone", "compartment", "error"])
    j = _judged_cells(cells, reference_identity, marker_labels, reference_compartment)
    return pd.DataFrame({"patient": j["patient"].values, "clone": j["clone_id"].values,
                         "compartment": j["compartment"].values,
                         "error": j["error"].values.astype(float)})


def _majority_error(labels: pd.Series, identity: str) -> float:
    """Errore dell'unità per clone: 0 o 1 se l'etichetta di maggioranza è unica; in caso di
    pareggio, la probabilità di errore scegliendo a caso fra le etichette a pari merito
    (1 - 1/k se l'identità è fra le k etichette, 1 altrimenti). Non dipende dall'ordine delle
    righe."""
    vc = labels.astype(str).value_counts()
    top = vc.index[vc.values == vc.values.max()]
    return 1.0 - (1.0 / len(top) if identity in set(top) else 0.0)


def clone_level_units(
    cells: pd.DataFrame,
    marker_labels: list[str],
    reference_compartment: str,
    error_labels: list[str] | None = None,
    marker_priority: list[str] | None = None,
    min_cells: int = MIN_CELLS_MARKER_DEFAULT,
    min_margin: float = MIN_MARGIN_DEFAULT,
) -> pd.DataFrame:
    """Unità della convenzione "per clone" (definizione della tesi, 10_loco.py schema blood).

    - Si considerano solo le cellule con etichetta in ``marker_labels`` + ``error_labels``
      (``error_labels=None``: tutte le etichette, anche quelle non T come mieloidi e B).
    - Positività ai marcatori: se ``marker_priority`` è dato, una cellula positiva ai
      marcatori di un'etichetta è considerata negativa per tutte le etichette che seguono
      nell'ordine (tesi: CD8T prima di CD4T, cioè CD4 positivo solo se CD8A e CD8B negativi).
    - Identità del clone dal solo compartimento di riferimento (>= ``min_cells`` cellule,
      margine >= ``min_margin``), come ``assign_reference_identity``.
    - Un'unità per coppia (clone, compartimento giudicato) con >= ``min_cells`` cellule;
      errore se l'etichetta di maggioranza è diversa dall'identità. Le etichette fuori mappa
      ammesse contano come errore. Pareggi: vedi ``_majority_error`` (la tesi usava l'ordine
      delle righe: su GSE278694 11 unità su 405 sono pareggi).
    """
    x = cells if error_labels is None else cells[cells["celltype"].isin(list(marker_labels) + list(error_labels))]
    x = x.copy()
    if marker_priority:
        taken = np.zeros(len(x), dtype=bool)
        for lab in marker_priority:
            col = f"pos_{lab}"
            x[col] = x[col].values & ~taken
            taken |= x[col].values
    ref = assign_reference_identity(x, marker_labels, reference_compartment=reference_compartment,
                                    min_cells=min_cells, min_margin=min_margin)
    if ref.empty:
        return pd.DataFrame(columns=["patient", "clone", "compartment", "error"])
    t = x[x["compartment"] != reference_compartment].merge(ref[["patient", "clone_id", "reference_identity"]],
                                                            on=["patient", "clone_id"])
    rows = []
    for (pat, clone, comp), g in t.groupby(["patient", "clone_id", "compartment"], observed=True):
        if len(g) < min_cells:
            continue
        rows.append({"patient": pat, "clone": clone, "compartment": comp,
                     "error": _majority_error(g["celltype"], str(g["reference_identity"].iloc[0]))})
    return pd.DataFrame(rows, columns=["patient", "clone", "compartment", "error"])


def _summarize_units(units: pd.DataFrame, min_patients_for_ci: int) -> dict[str, ClusterInterval]:
    out = {}
    for comp, g in units.groupby("compartment", observed=True):
        out[str(comp)] = cluster_t_interval(g["error"].values.astype(float), g["patient"].values,
                                            min_patients_for_ci, scale="logit")
    return out


def convention_definition(name: str, error_labels: list[str] | None = None,
                          marker_priority: list[str] | None = None,
                          min_cells: int = MIN_CELLS_MARKER_DEFAULT) -> str:
    """Definizione leggibile di una convenzione del tasso d'errore, con i parametri usati."""
    if name == "cell":
        return ("per cellula: ogni cellula giudicata è un'unità; etichette fuori dalla mappa "
                "dei marcatori escluse")
    labels = ("tutte, anche non T (es. mieloidi, B)" if error_labels is None
              else (", ".join(error_labels) or "nessuna"))
    prio = ("positività non esclusiva" if not marker_priority else
            "positività esclusiva in ordine " + " > ".join(marker_priority))
    return (f"per clone: un'unità per coppia clone-compartimento con almeno {min_cells} cellule, "
            f"etichetta di maggioranza (pareggi: probabilità di errore con scelta a caso fra le "
            f"etichette a pari merito); etichette fuori mappa contate come errore (etichette ammesse "
            f"fuori mappa: {labels}); {prio}")


# --------------------------------------------------------------------------- #
# (C) Flag per cellula
# --------------------------------------------------------------------------- #
FLAG_COLUMNS = ("audit_reference_label", "audit_label_vs_reference")
FLAG_DEFINITIONS = {
    "audit_reference_label": (
        f"Identità del clone T della cellula stimata dai marcatori canonici SOLO nel "
        f"compartimento di riferimento (>= {MIN_CELLS_MARKER_DEFAULT} cellule del clone nel "
        f"riferimento, margine >= {MIN_MARGIN_DEFAULT:.2f} fra la prima e la seconda identità). NA se "
        f"la cellula non ha TCR o il suo clone non ha un'identità di riferimento."),
    "audit_label_vs_reference": (
        "True se l'etichetta assegnata alla cellula discorda dall'identità di riferimento del "
        "clone, False se concorda. Stessa logica di marker_error_rate: NA se il clone non ha "
        "riferimento, se l'etichetta non è una chiave della marker_map, o se la cellula sta nel "
        "compartimento di riferimento (da cui il riferimento è stimato). È un segnale, non "
        "una correzione: l'etichetta originale non viene modificata."),
}


def cell_flags(
    cells: pd.DataFrame,
    reference_identity: pd.DataFrame,
    n_obs: int,
    marker_labels: list[str],
    reference_compartment: str = "PBMC",
    row_col: str = "_row",
) -> pd.DataFrame:
    """Flag per cellula, indicizzati per posizione di riga (0..n_obs-1) in adata.obs.

    ``audit_label_vs_reference`` ripete ESATTAMENTE i filtri di ``marker_error_rate``
    (``_judged_cells``), così la media dei flag valutabili di un compartimento coincide con la
    stima puntuale del tasso d'errore di quel compartimento."""
    ref_label = pd.Series(pd.NA, index=pd.RangeIndex(n_obs), dtype="object")
    vs_ref = pd.Series(pd.NA, index=pd.RangeIndex(n_obs), dtype="boolean")
    if not reference_identity.empty:
        j = cells.merge(reference_identity, on=["patient", "clone_id"], how="inner")
        ref_label.loc[j[row_col].values] = j["reference_identity"].astype(str).values
        judged = _judged_cells(cells, reference_identity, marker_labels, reference_compartment)
        vs_ref.loc[judged[row_col].values] = judged["error"].values
    return pd.DataFrame({"audit_reference_label": ref_label, "audit_label_vs_reference": vs_ref})


def export_audited(adata: ad.AnnData, result: "TcrValidationResult",
                   out_prefix: str | Path) -> tuple[Path, Path]:
    """Scrive ``<prefisso>_audited.h5ad`` (COPIA compressa gzip di adata con le due colonne
    di flag in obs e le definizioni in uns['genomic_audit_flags']) e ``<prefisso>_audit_flags.csv``
    (obs_names + flag). L'AnnData passato non viene modificato; le colonne esistenti non
    vengono sovrascritte (errore se una colonna di flag esiste già). Gli obs_names devono essere
    unici: il CSV identifica le cellule per nome."""
    if result.cell_flags is None:
        raise ValueError("nessun flag da esportare: servono marker_map e reference_compartment")
    if not adata.obs_names.is_unique:
        raise ValueError("adata.obs_names non sono unici: il CSV dei flag non identificherebbe le "
                         "cellule. Rendili unici (es. adata.obs_names_make_unique()) prima dell'export.")
    clash = [c for c in FLAG_COLUMNS if c in adata.obs.columns]
    if clash:
        raise ValueError(f"adata.obs contiene già le colonne {clash}: non vengono sovrascritte")
    if list(result.cell_flags.index) != list(adata.obs_names):
        raise ValueError("i flag non corrispondono alle righe di questo AnnData")
    out_prefix = Path(out_prefix)
    out_prefix.parent.mkdir(parents=True, exist_ok=True)
    h5ad_path = out_prefix.with_name(out_prefix.name + "_audited.h5ad")
    csv_path = out_prefix.with_name(out_prefix.name + "_audit_flags.csv")
    copy = adata.copy()
    for c in FLAG_COLUMNS:
        copy.obs[c] = result.cell_flags[c].values
    copy.obs["audit_reference_label"] = copy.obs["audit_reference_label"].astype("category")
    copy.uns["genomic_audit_flags"] = dict(FLAG_DEFINITIONS)
    # gzip: senza compressione la copia di GSE278694 occupava 2.9 GB (input 845 MB)
    copy.write_h5ad(h5ad_path, compression="gzip")
    result.cell_flags.to_csv(csv_path, index_label="obs_name")
    return h5ad_path, csv_path


# --------------------------------------------------------------------------- #
# Orchestratore
# --------------------------------------------------------------------------- #
@dataclass(frozen=True)
class TcrValidationResult:
    """Risultato completo del Modulo B."""
    n_cells_with_tcr: int
    n_clones_total: int
    discordance: DiscordanceResult
    marker_error: MarkerErrorResult | None
    narrative: str
    # Flag per cellula: una riga per riga di adata.obs, stesso indice (obs_names). None se
    # marker_map/reference_compartment non sono forniti.
    cell_flags: pd.DataFrame | None = None
    # Frazione delle cellule di adata con audit_label_vs_reference valutabile (non NA).
    flag_coverage: float | None = None
    # Tasso d'errore secondo le convenzioni richieste, con nome e definizione.
    conventions: dict[str, ConventionResult] | None = None
    # Esito del controllo del match dei barcode VDJ <-> metadati.
    barcode_match: BarcodeMatch | None = None
    min_patients_for_ci: int = MIN_PATIENTS_FOR_CI_DEFAULT


def _fmt_ci(res: ClusterInterval, signed: bool = False) -> str:
    f = "+.3f" if signed else ".3f"
    return f"IC95% [{res.ci_low:{f}}, {res.ci_high:{f}}]"


def _discordance_narrative(d: DiscordanceResult, min_patients_for_ci: int) -> str:
    if not d.sufficient:
        return (f"Eccesso di discordanza: intervallo non prodotto ({d.n_patients} pazienti "
                f"contribuiscono cloni comparabili, ne servono almeno {min_patients_for_ci}; "
                f"numerosità insufficiente).")
    if d.ci_low > 0:
        esito = ("l'intervallo è interamente sopra lo zero: entro lo stesso clone le etichette "
                 "differiscono fra compartimenti più che entro un compartimento. È compatibile con "
                 "un errore di annotazione, ma anche con le eccezioni biologiche elencate fra le "
                 "ipotesi")
    elif d.ci_high < 0:
        esito = ("l'intervallo è interamente sotto lo zero (direzione NEGATIVA). Con compartimenti di "
                 "dimensioni molto diverse lo stimatore pesato può avere valore atteso negativo anche "
                 "quando le etichette differiscono davvero: un valore negativo non indica assenza di "
                 "errori. Va letto con la sensibilità a pesi uguali")
    else:
        esito = "l'intervallo include lo zero: nessun eccesso rilevato rispetto al rumore entro compartimento"
    text = (f"Le cellule dello stesso clone T (stessa sequenza CDR3 della catena TRB) ricevono etichette "
            f"diverse fra compartimenti con un eccesso di discordanza di {d.mean_excess:+.3f} "
            f"(IC95% [{d.ci_low:+.3f}, {d.ci_high:+.3f}], "
            f"su {d.n_pairs:,} coppie clone-compartimenti da {d.n_patients} pazienti) rispetto al "
            f"rumore entro lo stesso compartimento: {esito}.")
    e = d.equal_weights
    if e is not None:
        ci = _fmt_ci(e, True) if e.sufficient else f"intervallo non prodotto ({e.reason})"
        text += f" Sensibilità a pesi uguali (valore atteso ½‖pA − pB‖² ≥ 0): {e.mean:+.3f} ({ci})."
    return text


def _conventions_narrative(conv: dict[str, ConventionResult] | None, reference_compartment: str,
                           min_patients_for_ci: int) -> str:
    if not conv:
        return ""
    parts = [f" Tasso d'errore dell'annotazione rispetto all'identità clonale dai marcatori "
             f"(riferimento: compartimento '{reference_compartment}'), secondo "
             f"{len(conv)} convenzion{'i' if len(conv) > 1 else 'e'} dichiarat{'e' if len(conv) > 1 else 'a'}:"]
    for c in conv.values():
        parts.append(f" [{c.name}] {c.definition} --")
        if not c.by_compartment:
            parts.append(" nessuna unità valutabile.")
        for comp, res in c.by_compartment.items():
            if not res.sufficient:
                parts.append(f" '{comp}': {res.mean:.3f}, IC non prodotto ({res.reason}).")
            else:
                parts.append(f" '{comp}': {res.mean:.3f} ({_fmt_ci(res)}, {res.n_groups} pazienti).")
        for (a, b), (res, n_cl) in c.paired_differences.items():
            ci = _fmt_ci(res, True) if res.sufficient else f"IC non prodotto ({res.reason})"
            parts.append(f" Differenza '{a}' - '{b}' sugli stessi cloni: {res.mean:+.3f}, {ci}, {n_cl} cloni.")
    return "".join(parts)


def run_tcr_validation(
    adata: ad.AnnData,
    contigs: pd.DataFrame,
    patient_col: str,
    compartment_col: str,
    celltype_col: str,
    barcode_col: str,
    min_cells_per_group: int = MIN_CELLS_DEFAULT,
    min_patients_for_ci: int = MIN_PATIENTS_FOR_CI_DEFAULT,
    marker_map: dict[str, list[str]] | None = None,
    reference_compartment: str | None = None,
    gene_col: str | None = None,
    convention: str = "both",
    clone_error_labels: list[str] | None = None,
    clone_marker_priority: list[str] | None = None,
    library_col: str | None = None,
    compartment_order: list[str] | None = None,
) -> TcrValidationResult:
    """Esegue il Modulo B. ``contigs`` è l'output di ``parse_vdj_contigs``.

    ``marker_map`` e ``reference_compartment`` sono opzionali: se entrambi forniti, viene
    calcolato anche il tasso d'errore per compartimento via marcatori canonici, secondo
    ``convention`` ("cell", "clone" o "both", default: entrambe affiancate). ``library_col``:
    colonna di adata.obs con la libreria, da indicare quando il manifest VDJ ha più librerie per
    la stessa coppia paziente-compartimento. ``compartment_order``: ordine dei compartimenti
    nelle differenze appaiate (A - B con A prima di B); senza, ordine alfabetico.
    """
    if convention not in ("cell", "clone", "both"):
        raise ValueError(f"convention deve essere cell, clone o both, non '{convention}'")
    if bool(marker_map) != bool(reference_compartment):
        raise ValueError("marker_map e reference_compartment vanno indicati insieme")
    cols = [patient_col, compartment_col, celltype_col, barcode_col] + ([library_col] if library_col else [])
    missing = [c for c in dict.fromkeys(cols) if c not in adata.obs.columns]
    if missing:
        raise ValueError(f"colonne non trovate in adata.obs: {missing}")
    if "library" not in contigs.columns:
        contigs = contigs.assign(library="")
    has_lib = (contigs["library"].astype(str) != "").any()
    if has_lib and not library_col:
        raise ValueError("il manifest VDJ indica le librerie: indica anche la colonna libreria dei "
                         "metadati (--library-col)")
    if library_col and not has_lib:
        raise ValueError("è indicata la colonna libreria dei metadati ma il manifest VDJ non ha la "
                         "colonna 'library'")
    clones = build_clonotypes(contigs)

    obs = pd.DataFrame({
        "patient": adata.obs[patient_col].astype(str).values,
        "compartment": adata.obs[compartment_col].astype(str).values,
        "celltype": adata.obs[celltype_col].astype(str).values,
        "barcode": adata.obs[barcode_col].astype(str).values,
        "library": adata.obs[library_col].astype(str).values if library_col else "",
    })
    # Posizione di riga originale: serve solo a riscrivere i flag per cellula sulle righe
    # giuste di adata.obs dopo la merge (che non preserva l'ordine). Non entra in nessun calcolo.
    obs["_row"] = np.arange(len(obs))

    # Le colonne pos_<etichetta> vengono allineate per POSIZIONE (stessa riga di adata.obs), non
    # per join sul barcode: il barcode 10x da solo NON è univoco a livello globale (si ripete fra
    # pazienti/compartimenti diversi), quindi un join per barcode mescolerebbe cellule diverse.
    if marker_map:
        pos_df = compute_marker_positivity(adata, marker_map, gene_col=gene_col).reset_index(drop=True)
        obs = pd.concat([obs, pos_df], axis=1)

    obs, clones, bmatch = match_barcodes(obs, clones)
    cells = obs.merge(clones, on=_KEY, how="inner")
    if cells.empty:
        raise ValueError(
            "nessuna cellula con TCR rilevato corrisponde ai metadati forniti: controlla che "
            "paziente/compartimento/barcode combacino fra AnnData e file VDJ"
        )

    pw = pairwise_excess_discordance(cells, min_cells=min_cells_per_group)
    discordance = _summarize_discordance(pw, min_patients_for_ci)

    marker_error = None
    flags, coverage, conventions = None, None, None
    if marker_map and reference_compartment:
        labels = list(marker_map)
        ref_identity = assign_reference_identity(cells, labels, reference_compartment=reference_compartment)
        marker_error = marker_error_rate(cells, ref_identity, labels,
                                         reference_compartment=reference_compartment,
                                         min_patients_for_ci=min_patients_for_ci)
        flags = cell_flags(cells, ref_identity, n_obs=adata.n_obs, marker_labels=labels,
                           reference_compartment=reference_compartment)
        flags.index = adata.obs_names
        coverage = float(flags["audit_label_vs_reference"].notna().mean()) if adata.n_obs else 0.0
        conventions = {}
        if convention in ("cell", "both"):
            units = _cell_units(cells, ref_identity, labels, reference_compartment)
            conventions["cell"] = ConventionResult(
                "cell", convention_definition("cell"), dict(marker_error.by_compartment),
                _paired_differences(units, min_patients_for_ci, compartment_order))
        if convention in ("clone", "both"):
            units = clone_level_units(cells, labels, reference_compartment,
                                      error_labels=clone_error_labels,
                                      marker_priority=clone_marker_priority)
            conventions["clone"] = ConventionResult(
                "clone", convention_definition("clone", clone_error_labels, clone_marker_priority),
                _summarize_units(units, min_patients_for_ci),
                _paired_differences(units, min_patients_for_ci, compartment_order))

    narrative = _discordance_narrative(discordance, min_patients_for_ci) + (
        _conventions_narrative(conventions, reference_compartment, min_patients_for_ci)
        if conventions else "")

    return TcrValidationResult(
        n_cells_with_tcr=len(cells), n_clones_total=int(cells["clone_id"].nunique()),
        discordance=discordance, marker_error=marker_error, narrative=narrative,
        cell_flags=flags, flag_coverage=coverage, conventions=conventions,
        barcode_match=bmatch, min_patients_for_ci=min_patients_for_ci,
    )


# --------------------------------------------------------------------------- #
# Descrittivo: cellule con TCR per etichetta (web app)
# --------------------------------------------------------------------------- #
T_LABEL_HINTS = ("CD4", "CD8", "T cell", "Treg", "MAIT", "gdT", "T_")


def tcr_by_celltype(adata: ad.AnnData, contigs: pd.DataFrame, patient_col: str,
                    compartment_col: str, celltype_col: str, barcode_col: str,
                    library_col: str | None = None) -> pd.DataFrame:
    """Per ogni etichetta di tipo cellulare: cellule totali, cellule con catena TRB rilevata e
    quante di queste appartengono a un clone espanso (>= 2 cellule nello stesso paziente).

    È una tabella DESCRITTIVA. Un TCR in un'etichetta non-T (NK, mieloidi, stromali) può
    essere un doppietto, RNA ambientale o un errore di annotazione: la tabella non distingue fra
    questi casi e non classifica le singole cellule. Usa la stessa corrispondenza dei barcode
    del Modulo B (``match_barcodes``)."""
    obs = pd.DataFrame({
        "patient": adata.obs[patient_col].astype(str).values,
        "compartment": adata.obs[compartment_col].astype(str).values,
        "celltype": adata.obs[celltype_col].astype(str).values,
        "barcode": adata.obs[barcode_col].astype(str).values,
        "library": adata.obs[library_col].astype(str).values if library_col else "",
    })
    clones = build_clonotypes(contigs)
    obs, clones, _ = match_barcodes(obs, clones)
    size = clones.groupby(["patient", "clone_id"]).barcode.transform("size")
    clones = clones.assign(expanded=size >= 2)
    j = obs.merge(clones, on=_KEY, how="left")
    out = j.groupby("celltype").agg(
        cellule=("barcode", "size"),
        con_TRB=("clone_id", lambda s: int(s.notna().sum())),
        in_cloni_espansi=("expanded", lambda s: int(s.fillna(False).astype(bool).sum())),
    )
    out["% con TRB"] = (100 * out["con_TRB"] / out["cellule"]).round(1)
    out["etichetta T"] = [any(h.lower() in str(c).lower() for h in T_LABEL_HINTS) for c in out.index]
    return out.sort_values(["etichetta T", "% con TRB"], ascending=[True, False])
