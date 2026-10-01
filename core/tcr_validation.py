"""Modulo B -- Validazione dell'annotazione di tipo cellulare via TCR.

Logica di discordanza a coppie e cluster bootstrap sui pazienti portata da
pdac-ml/src/07_robustness.py e pdac-ml/src/10_loco.py, generalizzata: nomi di colonna
(paziente, compartimento, tipo cellulare, barcode) configurabili, nessuna dipendenza dal
dataset PDAC originale.

Scelta sull'estimatore D_A/D_B (discordanza entro compartimento, baseline dell'eccesso):
usa la versione SENZA reinserimento (coppie di cellule distinte dallo stesso gruppo),
validata come non distorta in 07_robustness.py punto (A) -- non la versione con
proporzioni (con reinserimento), che sovrastima sistematicamente D su gruppi piccoli e
quindi sottostima l'eccesso proprio nei cloni piu' piccoli (il regime tipico qui, dato il
min_cells di default = 2). D(A×B) cross-compartimento resta con proporzioni: confronta
due gruppi DIVERSI (una cellula per lato), nessun problema di reinserimento.
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

from core.stats import BootstrapResult, cluster_bootstrap

MIN_CELLS_DEFAULT = 2
MIN_PATIENTS_FOR_CI_DEFAULT = 5
MIN_CELLS_MARKER_DEFAULT = 3
MIN_MARGIN_DEFAULT = 0.20


# --------------------------------------------------------------------------- #
# Ingestione VDJ (Cell Ranger)
# --------------------------------------------------------------------------- #
def parse_vdj_contigs(files: list[tuple[Path, str, str]]) -> pd.DataFrame:
    """``files``: lista di (percorso_csv, patient_id, compartment) -- l'associazione
    file->paziente->compartimento e' fornita dal chiamante (CLI/app), non dedotta dal
    nome del file: i CSV Cell Ranger non contengono questa informazione al loro interno,
    e non c'e' una convenzione di nome universale.

    Filtra is_cell/high_confidence/productive/full_length (se le colonne sono presenti),
    catene limitate a TRA/TRB, e per (paziente, compartimento, barcode, catena) mantiene
    il contig con piu' UMI.
    """
    frames = []
    for path, patient, compartment in files:
        df = pd.read_csv(path)
        for col in ("is_cell", "high_confidence", "productive", "full_length"):
            if col in df.columns:
                df = df[df[col].astype(str).str.lower().isin(["true", "1"])]
        required = {"barcode", "chain", "cdr3_nt"}
        missing = required - set(df.columns)
        if missing:
            raise ValueError(f"{path}: colonne mancanti nel CSV VDJ: {sorted(missing)}")
        df = df[df.chain.isin(["TRA", "TRB"])]
        df = df[df.cdr3_nt.notna() & (df.cdr3_nt.astype(str) != "None")]
        df = df.assign(patient=str(patient), compartment=str(compartment))
        if "umis" not in df.columns:
            df = df.assign(umis=1)
        frames.append(df[["patient", "compartment", "barcode", "chain", "cdr3_nt", "umis"]])
    if not frames:
        raise ValueError("nessun file VDJ fornito")
    c = pd.concat(frames, ignore_index=True)
    c = c.sort_values("umis", ascending=False).drop_duplicates(
        ["patient", "compartment", "barcode", "chain"]
    )
    return c


def build_clonotypes(contigs: pd.DataFrame) -> pd.DataFrame:
    """Il clonotipo e' definito dalla sequenza nucleotidica CDR3 della catena TRB. Una
    riga per cellula con TRB rilevato."""
    trb = contigs[contigs.chain == "TRB"][["patient", "compartment", "barcode", "cdr3_nt"]]
    trb = trb.rename(columns={"cdr3_nt": "clone_id"})
    return trb.drop_duplicates(["patient", "compartment", "barcode"]).reset_index(drop=True)


# --------------------------------------------------------------------------- #
# Controllo del match dei barcode (Pulizia 2)
# --------------------------------------------------------------------------- #
MIN_BARCODE_MATCH = 0.5
_SUFFIX_RE = re.compile(r"-\d+$")
_KEY = ["patient", "compartment", "barcode"]


@dataclass(frozen=True)
class BarcodeMatch:
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


def match_barcodes(obs: pd.DataFrame, clones: pd.DataFrame,
                   min_fraction: float = MIN_BARCODE_MATCH) -> tuple[pd.DataFrame, pd.DataFrame, BarcodeMatch]:
    """Frazione di barcode VDJ (cellule con TRB) ritrovati nei metadati, per (paziente,
    compartimento, barcode). Se rimuovere il suffisso "-N" (es. "-1" di Cell Ranger) da
    entrambi i lati aumenta il match senza creare collisioni (barcode diversi che diventano
    uguali, es. -1 e -2 dello stesso campione aggregato), la normalizzazione viene applicata
    e dichiarata. Errore esplicito se il match finale e' sotto ``min_fraction``."""
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
            f"solo il {f_final:.1%} dei barcode VDJ e' stato ritrovato nei metadati (soglia "
            f"{min_fraction:.0%}), anche dopo l'eventuale normalizzazione del suffisso -N. Esempi "
            f"VDJ: {ex_v}; esempi metadati: {ex_o}. Controlla che paziente e compartimento del "
            f"manifest coincidano con quelli dell'AnnData e che i barcode abbiano lo stesso formato.")
    return used_obs, used_clones, BarcodeMatch(len(clones), n_final, f_final, suffixes, sentence)


# --------------------------------------------------------------------------- #
# (A) Discordanza a coppie
# --------------------------------------------------------------------------- #
def within_compartment_discordance(labels: np.ndarray) -> float:
    """D_A = 1 - sum(n_k(n_k-1)) / (n(n-1)): probabilita' che DUE cellule distinte dello
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
    problema di reinserimento, la versione con proporzioni e' corretta qui."""
    return 1.0 - float(np.dot(p_a, p_b))


def pairwise_excess_discordance(
    cells: pd.DataFrame,
    patient_col: str = "patient",
    compartment_col: str = "compartment",
    clone_col: str = "clone_id",
    label_col: str = "celltype",
    min_cells: int = MIN_CELLS_DEFAULT,
) -> pd.DataFrame:
    """Una riga per (paziente, clone, coppia di compartimenti) con >=min_cells cellule
    per lato: eccesso = D(A×B) - [n_A D_A + n_B D_B] / (n_A+n_B)."""
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
                })
    return pd.DataFrame(rows)


@dataclass(frozen=True)
class DiscordanceResult:
    n_pairs: int
    n_patients: int
    mean_excess: float
    ci_low: float
    ci_high: float
    sufficient: bool
    by_compartment_pair: pd.DataFrame


def _summarize_discordance(pw: pd.DataFrame, n_boot: int, seed: int,
                            min_patients_for_ci: int) -> DiscordanceResult:
    if pw.empty:
        return DiscordanceResult(n_pairs=0, n_patients=0, mean_excess=float("nan"),
                                  ci_low=float("nan"), ci_high=float("nan"), sufficient=False,
                                  by_compartment_pair=pd.DataFrame())
    overall = cluster_bootstrap(pw["excess"].values, pw["patient"].values, n_boot=n_boot,
                                 seed=seed, min_groups=min_patients_for_ci)
    rows = []
    for (a, b), g in pw.groupby(["comp_a", "comp_b"]):
        res = cluster_bootstrap(g["excess"].values, g["patient"].values, n_boot=n_boot,
                                 seed=seed, min_groups=min_patients_for_ci)
        rows.append({"comp_a": a, "comp_b": b, "n_pairs": len(g), "n_patients": res.n_groups,
                     "excess": res.mean, "ci_low": res.ci_low, "ci_high": res.ci_high,
                     "sufficient": res.sufficient})
    return DiscordanceResult(n_pairs=len(pw), n_patients=int(pw.patient.nunique()),
                              mean_excess=overall.mean, ci_low=overall.ci_low,
                              ci_high=overall.ci_high, sufficient=overall.sufficient,
                              by_compartment_pair=pd.DataFrame(rows))


# --------------------------------------------------------------------------- #
# (B, opzionale) Tasso d'errore per compartimento via marcatori canonici
# --------------------------------------------------------------------------- #
def compute_marker_positivity(
    adata: ad.AnnData, marker_map: dict[str, list[str]], gene_col: str | None = None
) -> pd.DataFrame:
    """Una colonna booleana 'pos_<etichetta>' per etichetta in marker_map: True se
    almeno uno dei geni marcatori di quell'etichetta ha conta > 0 in quella cellula.
    Stesso ordine di riga di adata.obs."""
    gene_names = (adata.var[gene_col] if gene_col else adata.var_names).astype(str).values
    gene_to_idx = {g: i for i, g in enumerate(gene_names)}
    X = adata.X

    out = {}
    for label, genes in marker_map.items():
        idx = [gene_to_idx[g] for g in genes if g in gene_to_idx]
        missing = [g for g in genes if g not in gene_to_idx]
        if missing:
            warnings.warn(f"geni marcatori non trovati per '{label}': {missing}")
        if not idx:
            out[f"pos_{label}"] = np.zeros(adata.n_obs, dtype=bool)
            continue
        sub = X[:, idx]
        sub = sub.toarray() if sp.issparse(sub) else np.asarray(sub)
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
    """Identita' di riferimento del clone dai marcatori, stimata SOLO dal compartimento
    di riferimento (tipicamente il sangue): niente circolarita' con i compartimenti poi
    giudicati. Cloni ambigui (margine fra l'etichetta vincente e la seconda < min_margin)
    sono scartati, non forzati."""
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
    return pd.DataFrame(rows)


@dataclass(frozen=True)
class MarkerErrorResult:
    reference_compartment: str
    n_resolved_clones: int
    by_compartment: dict[str, BootstrapResult]


def marker_error_rate(
    cells: pd.DataFrame,
    reference_identity: pd.DataFrame,
    patient_col: str = "patient",
    clone_col: str = "clone_id",
    compartment_col: str = "compartment",
    label_col: str = "celltype",
    reference_compartment: str = "PBMC",
    n_boot: int = 2000,
    seed: int = 0,
    min_patients_for_ci: int = MIN_PATIENTS_FOR_CI_DEFAULT,
) -> MarkerErrorResult:
    """Confronta l'etichetta assegnata (celltype_col) con l'identita' di riferimento del
    clone (dai marcatori, stimata solo nel compartimento di riferimento) su tutti gli
    ALTRI compartimenti. Tasso d'errore per compartimento con cluster bootstrap sui
    pazienti."""
    if reference_identity.empty:
        return MarkerErrorResult(reference_compartment=reference_compartment,
                                  n_resolved_clones=0, by_compartment={})

    j = cells.merge(reference_identity, on=[patient_col, clone_col], how="inner")
    j = j[j[compartment_col] != reference_compartment]
    j = j[j[label_col].isin(reference_identity["reference_identity"].unique())]
    j["error"] = j[label_col].astype(str) != j["reference_identity"].astype(str)

    by_compartment: dict[str, BootstrapResult] = {}
    for comp, g in j.groupby(compartment_col, observed=True):
        by_compartment[str(comp)] = cluster_bootstrap(
            g["error"].values.astype(float), g[patient_col].values,
            n_boot=n_boot, seed=seed, min_groups=min_patients_for_ci,
        )
    return MarkerErrorResult(reference_compartment=reference_compartment,
                              n_resolved_clones=reference_identity[clone_col].nunique(),
                              by_compartment=by_compartment)


# --------------------------------------------------------------------------- #
# (B2) Convenzioni del tasso d'errore: "per cellula" (MVP) e "per clone" (tesi)
# --------------------------------------------------------------------------- #
CONVENTIONS = ("cell", "clone")


@dataclass(frozen=True)
class ConventionResult:
    """Tasso d'errore per compartimento secondo una convenzione dichiarata, piu' le
    differenze appaiate fra compartimenti (stessi cloni)."""
    name: str
    definition: str
    by_compartment: dict[str, BootstrapResult]
    paired_differences: dict[tuple[str, str], tuple[BootstrapResult, int]]


def _paired_differences(units: pd.DataFrame, n_boot: int, seed: int,
                        min_patients_for_ci: int) -> dict[tuple[str, str], tuple[BootstrapResult, int]]:
    """``units``: colonne patient, clone, compartment, error (una riga per unita').
    Per ogni coppia di compartimenti giudicati: errore medio del clone in A meno quello in B,
    sui cloni presenti in entrambi; cluster bootstrap sui pazienti. Valore: (risultato,
    numero di cloni)."""
    out = {}
    per = units.groupby(["patient", "clone", "compartment"], observed=True).error.mean().unstack("compartment")
    comps = sorted(per.columns, key=lambda c: (c != "Tumor", c))
    for i in range(len(comps)):
        for j in range(i + 1, len(comps)):
            a, b = comps[i], comps[j]
            d = per.dropna(subset=[a, b]).reset_index()
            if d.empty:
                continue
            res = cluster_bootstrap((d[a] - d[b]).values, d["patient"].values, n_boot=n_boot,
                                    seed=seed, min_groups=min_patients_for_ci)
            out[(a, b)] = (res, len(d))
    return out


def _cell_units(cells: pd.DataFrame, reference_identity: pd.DataFrame,
                reference_compartment: str) -> pd.DataFrame:
    """Unita' della convenzione "per cellula": stessi filtri di ``marker_error_rate``."""
    if reference_identity.empty:
        return pd.DataFrame(columns=["patient", "clone", "compartment", "error"])
    j = cells.merge(reference_identity, on=["patient", "clone_id"], how="inner")
    j = j[j["compartment"] != reference_compartment]
    j = j[j["celltype"].isin(reference_identity["reference_identity"].unique())]
    return pd.DataFrame({"patient": j["patient"].values, "clone": j["clone_id"].values,
                         "compartment": j["compartment"].values,
                         "error": (j["celltype"].astype(str) != j["reference_identity"].astype(str)).values.astype(float)})


def clone_level_units(
    cells: pd.DataFrame,
    marker_labels: list[str],
    reference_compartment: str,
    error_labels: list[str] | None = None,
    marker_priority: list[str] | None = None,
    min_cells: int = MIN_CELLS_MARKER_DEFAULT,
    min_margin: float = MIN_MARGIN_DEFAULT,
) -> pd.DataFrame:
    """Unita' della convenzione "per clone" (definizione della tesi, 10_loco.py schema blood).

    - Si considerano solo le cellule con etichetta in ``marker_labels`` + ``error_labels``
      (``error_labels=None``: tutte le etichette).
    - Positivita' ai marcatori: se ``marker_priority`` e' dato, una cellula positiva ai
      marcatori di un'etichetta e' considerata negativa per tutte le etichette che seguono
      nell'ordine (tesi: CD8T prima di CD4T, cioe' CD4 positivo solo se CD8A e CD8B negativi).
    - Identita' del clone dal solo compartimento di riferimento (>= ``min_cells`` cellule,
      margine >= ``min_margin``), come ``assign_reference_identity``.
    - Un'unita' per coppia (clone, compartimento giudicato) con >= ``min_cells`` cellule;
      etichetta = etichetta di maggioranza; errore se diversa dall'identita'. Le etichette
      fuori mappa contano come errore.
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
        lab = g["celltype"].value_counts().index[0]
        rows.append({"patient": pat, "clone": clone, "compartment": comp,
                     "error": float(lab != g["reference_identity"].iloc[0])})
    return pd.DataFrame(rows, columns=["patient", "clone", "compartment", "error"])


def _summarize_units(units: pd.DataFrame, n_boot: int, seed: int,
                     min_patients_for_ci: int) -> dict[str, BootstrapResult]:
    out = {}
    for comp, g in units.groupby("compartment", observed=True):
        out[str(comp)] = cluster_bootstrap(g["error"].values.astype(float), g["patient"].values,
                                           n_boot=n_boot, seed=seed, min_groups=min_patients_for_ci)
    return out


def convention_definition(name: str, error_labels: list[str] | None = None,
                          marker_priority: list[str] | None = None) -> str:
    if name == "cell":
        return ("per cellula: ogni cellula giudicata e' un'unita'; etichette fuori dalla mappa "
                "dei marcatori escluse")
    labels = "tutte" if error_labels is None else (", ".join(error_labels) or "nessuna")
    prio = ("positivita' non esclusiva" if not marker_priority else
            "positivita' esclusiva in ordine " + " > ".join(marker_priority))
    return ("per clone: un'unita' per coppia clone-compartimento con almeno 3 cellule, etichetta "
            f"di maggioranza; etichette fuori mappa contate come errore (etichette ammesse fuori "
            f"mappa: {labels}); {prio}")


# --------------------------------------------------------------------------- #
# (C) Flag per cellula (Intervento 2)
# --------------------------------------------------------------------------- #
FLAG_COLUMNS = ("audit_reference_label", "audit_label_vs_reference")
FLAG_DEFINITIONS = {
    "audit_reference_label": (
        "Identita' del clone T della cellula stimata dai marcatori canonici SOLO nel "
        "compartimento di riferimento (schema solo sangue: >= 3 cellule del clone nel "
        "riferimento, margine >= 0.20 fra la prima e la seconda identita'). NA se la cellula "
        "non ha TCR o il suo clone non ha un'identita' di riferimento."),
    "audit_label_vs_reference": (
        "True se l'etichetta assegnata alla cellula discorda dall'identita' di riferimento del "
        "clone, False se concorda. Stessa logica di marker_error_rate: NA se il clone non ha "
        "riferimento, se l'etichetta e' fuori dalla marker_map, o se la cellula sta nel "
        "compartimento di riferimento (da cui il riferimento e' stimato). E' un segnale, non "
        "una correzione: l'etichetta originale non viene modificata."),
}


def cell_flags(
    cells: pd.DataFrame,
    reference_identity: pd.DataFrame,
    n_obs: int,
    patient_col: str = "patient",
    clone_col: str = "clone_id",
    compartment_col: str = "compartment",
    label_col: str = "celltype",
    reference_compartment: str = "PBMC",
    row_col: str = "_row",
) -> pd.DataFrame:
    """Flag per cellula, indicizzati per posizione di riga (0..n_obs-1) in adata.obs.

    ``audit_label_vs_reference`` ripete ESATTAMENTE i filtri di ``marker_error_rate``
    (join interno sull'identita' di riferimento, esclusione del compartimento di
    riferimento, etichette limitate a quelle presenti fra le identita' di riferimento),
    cosi' la media dei flag valutabili di un compartimento coincide con la stima puntuale
    del tasso d'errore di quel compartimento."""
    ref_label = pd.Series(pd.NA, index=pd.RangeIndex(n_obs), dtype="object")
    vs_ref = pd.Series(pd.NA, index=pd.RangeIndex(n_obs), dtype="boolean")
    if not reference_identity.empty:
        j = cells.merge(reference_identity, on=[patient_col, clone_col], how="inner")
        ref_label.loc[j[row_col].values] = j["reference_identity"].astype(str).values
        j = j[j[compartment_col] != reference_compartment]
        j = j[j[label_col].isin(reference_identity["reference_identity"].unique())]
        err = j[label_col].astype(str) != j["reference_identity"].astype(str)
        vs_ref.loc[j[row_col].values] = err.values
    return pd.DataFrame({"audit_reference_label": ref_label, "audit_label_vs_reference": vs_ref})


def export_audited(adata: ad.AnnData, result: "TcrValidationResult",
                   out_prefix: str | Path) -> tuple[Path, Path]:
    """Scrive ``<prefisso>_audited.h5ad`` (COPIA compressa gzip di adata con le due colonne
    di flag in obs e le definizioni in uns['genomic_audit_flags']) e ``<prefisso>_audit_flags.csv``
    (obs_names + flag). L'AnnData passato non viene modificato; le colonne esistenti non
    vengono sovrascritte (errore se una colonna di flag esiste gia')."""
    if result.cell_flags is None:
        raise ValueError("nessun flag da esportare: servono marker_map e reference_compartment")
    clash = [c for c in FLAG_COLUMNS if c in adata.obs.columns]
    if clash:
        raise ValueError(f"adata.obs contiene gia' le colonne {clash}: non vengono sovrascritte")
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
    n_cells_with_tcr: int
    n_clones_total: int
    discordance: DiscordanceResult
    marker_error: MarkerErrorResult | None
    narrative: str
    # Flag per cellula (Intervento 2): una riga per riga di adata.obs, stesso indice
    # (obs_names). None se marker_map/reference_compartment non sono forniti.
    cell_flags: pd.DataFrame | None = None
    # Frazione delle cellule di adata con audit_label_vs_reference valutabile (non NA).
    flag_coverage: float | None = None
    # Tasso d'errore secondo le convenzioni richieste (Pulizia 1), con nome e definizione.
    conventions: dict[str, ConventionResult] | None = None
    # Esito del controllo del match dei barcode VDJ <-> metadati (Pulizia 2).
    barcode_match: BarcodeMatch | None = None


def _discordance_narrative(d: DiscordanceResult, min_patients_for_ci: int) -> str:
    if not d.sufficient:
        return (f"Numerosita' insufficiente per una stima affidabile dell'eccesso di discordanza "
                f"({d.n_patients} pazienti contribuiscono cloni comparabili, ne servono almeno "
                f"{min_patients_for_ci}): nessun intervallo di confidenza riportato.")
    if d.ci_low > 0:
        esito = "un effetto reale (l'intervallo esclude lo zero)"
    elif d.ci_high < 0:
        esito = ("l'intervallo esclude lo zero ma in direzione NEGATIVA -- l'opposto di un errore "
                 "di annotazione (che produrrebbe eccesso >=0). Non e' la firma cercata: probabile "
                 "rumore campionario su pochi pazienti/cloni, da non interpretare come un effetto")
    else:
        esito = "non distinguibile dal rumore di base entro compartimento (l'intervallo include lo zero)"
    return (f"Le cellule dello stesso clone T (stessa sequenza CDR3 della catena TRB) cambiano "
            f"etichetta di tipo cellulare fra compartimenti tissutali con un eccesso di discordanza "
            f"di {d.mean_excess:+.3f} (IC95% [{d.ci_low:+.3f}, {d.ci_high:+.3f}], su {d.n_pairs:,} "
            f"coppie clone-compartimenti da {d.n_patients} pazienti) rispetto al rumore di base "
            f"entro lo stesso compartimento: {esito}.")


def _marker_narrative(m: MarkerErrorResult | None, min_patients_for_ci: int) -> str:
    if m is None:
        return ""
    if not m.by_compartment:
        return " Nessun clone risolvibile per il tasso d'errore via marcatori (identita' ambigua o troppe poche cellule di riferimento)."
    parts = [f" Tasso d'errore dell'annotazione rispetto all'identita' clonale dai marcatori "
             f"(riferimento: compartimento '{m.reference_compartment}', {m.n_resolved_clones} cloni risolti):"]
    for comp, res in m.by_compartment.items():
        if not res.sufficient:
            parts.append(f" '{comp}': numerosita' insufficiente ({res.n_groups} pazienti, "
                         f"ne servono almeno {min_patients_for_ci}).")
        else:
            parts.append(f" '{comp}': {res.mean:.3f} (IC95% [{res.ci_low:.3f}, {res.ci_high:.3f}], "
                         f"{res.n_groups} pazienti).")
    return "".join(parts)


def _conventions_narrative(conv: dict[str, ConventionResult] | None, reference_compartment: str,
                           min_patients_for_ci: int) -> str:
    if not conv:
        return ""
    parts = [f" Tasso d'errore dell'annotazione rispetto all'identita' clonale dai marcatori "
             f"(riferimento: compartimento '{reference_compartment}'), secondo "
             f"{len(conv)} convenzion{'i' if len(conv) > 1 else 'e'} dichiarat{'e' if len(conv) > 1 else 'a'}:"]
    for c in conv.values():
        parts.append(f" [{c.name}] {c.definition} --")
        if not c.by_compartment:
            parts.append(" nessuna unita' valutabile.")
        for comp, res in c.by_compartment.items():
            if not res.sufficient:
                parts.append(f" '{comp}': {res.mean:.3f}, IC non prodotto ({res.n_groups} pazienti, "
                             f"ne servono almeno {min_patients_for_ci}).")
            else:
                parts.append(f" '{comp}': {res.mean:.3f} (IC95% [{res.ci_low:.3f}, {res.ci_high:.3f}], "
                             f"{res.n_groups} pazienti).")
        for (a, b), (res, n_cl) in c.paired_differences.items():
            ci = (f"IC95% [{res.ci_low:+.3f}, {res.ci_high:+.3f}]" if res.sufficient
                  else f"IC non prodotto ({res.n_groups} pazienti)")
            parts.append(f" Differenza '{a}' - '{b}' sugli stessi cloni: {res.mean:+.3f}, {ci}, {n_cl} cloni.")
    return "".join(parts)


def run_tcr_validation(
    adata: ad.AnnData,
    contigs: pd.DataFrame,
    patient_col: str,
    compartment_col: str,
    celltype_col: str,
    barcode_col: str,
    n_boot: int = 2000,
    seed: int = 0,
    min_cells_per_group: int = MIN_CELLS_DEFAULT,
    min_patients_for_ci: int = MIN_PATIENTS_FOR_CI_DEFAULT,
    marker_map: dict[str, list[str]] | None = None,
    reference_compartment: str | None = None,
    gene_col: str | None = None,
    convention: str = "both",
    clone_error_labels: list[str] | None = None,
    clone_marker_priority: list[str] | None = None,
) -> TcrValidationResult:
    """Esegue il Modulo B. ``contigs`` e' l'output gia' filtrato di ``parse_vdj_contigs``.

    ``marker_map`` e ``reference_compartment`` sono opzionali: se entrambi forniti, viene
    calcolato anche il tasso d'errore per compartimento via marcatori canonici, secondo
    ``convention`` ("cell", "clone" o "both", default: entrambe affiancate). ``marker_error``
    (convenzione per cellula) e' calcolato comunque, invariato.
    """
    if convention not in ("cell", "clone", "both"):
        raise ValueError(f"convention deve essere cell, clone o both, non '{convention}'")
    clones = build_clonotypes(contigs)

    obs = adata.obs[[patient_col, compartment_col, celltype_col, barcode_col]].copy()
    obs.columns = ["patient", "compartment", "celltype", "barcode"]
    for c in ("patient", "compartment", "celltype", "barcode"):
        obs[c] = obs[c].astype(str)
    obs = obs.reset_index(drop=True)
    # Posizione di riga originale: serve solo a riscrivere i flag per cellula sulle righe
    # giuste di adata.obs dopo la merge (che non preserva l'ordine). Non entra in nessun
    # calcolo.
    obs["_row"] = np.arange(len(obs))

    # Le colonne pos_<etichetta> vengono allineate per POSIZIONE (stessa riga di
    # adata.obs), non per join sul barcode: il barcode 10x da solo NON e' univoco a
    # livello globale (si ripete fra pazienti/compartimenti diversi), quindi un join
    # per barcode rischierebbe di mescolare cellule di pazienti diversi.
    if marker_map:
        pos_df = compute_marker_positivity(adata, marker_map, gene_col=gene_col).reset_index(drop=True)
        obs = pd.concat([obs, pos_df], axis=1)

    obs, clones, bmatch = match_barcodes(obs, clones)
    cells = obs.merge(clones, on=["patient", "compartment", "barcode"], how="inner")
    if cells.empty:
        raise ValueError(
            "nessuna cellula con TCR rilevato corrisponde ai metadati forniti: controlla che "
            "paziente/compartimento/barcode combacino fra AnnData e file VDJ"
        )

    pw = pairwise_excess_discordance(cells, min_cells=min_cells_per_group)
    discordance = _summarize_discordance(pw, n_boot, seed, min_patients_for_ci)

    marker_error = None
    flags, coverage = None, None
    if marker_map and reference_compartment:
        ref_identity = assign_reference_identity(
            cells, list(marker_map), reference_compartment=reference_compartment,
        )
        marker_error = marker_error_rate(
            cells, ref_identity, reference_compartment=reference_compartment,
            n_boot=n_boot, seed=seed, min_patients_for_ci=min_patients_for_ci,
        )
        flags = cell_flags(cells, ref_identity, n_obs=adata.n_obs,
                           reference_compartment=reference_compartment)
        flags.index = adata.obs_names
        coverage = float(flags["audit_label_vs_reference"].notna().mean()) if adata.n_obs else 0.0

    conventions = None
    if marker_map and reference_compartment:
        conventions = {}
        if convention in ("cell", "both"):
            units = _cell_units(cells, ref_identity, reference_compartment)
            conventions["cell"] = ConventionResult(
                "cell", convention_definition("cell"), dict(marker_error.by_compartment),
                _paired_differences(units, n_boot, seed, min_patients_for_ci))
        if convention in ("clone", "both"):
            units = clone_level_units(cells, list(marker_map), reference_compartment,
                                      error_labels=clone_error_labels,
                                      marker_priority=clone_marker_priority)
            conventions["clone"] = ConventionResult(
                "clone", convention_definition("clone", clone_error_labels, clone_marker_priority),
                _summarize_units(units, n_boot, seed, min_patients_for_ci),
                _paired_differences(units, n_boot, seed, min_patients_for_ci))

    narrative = _discordance_narrative(discordance, min_patients_for_ci) + (
        _conventions_narrative(conventions, reference_compartment, min_patients_for_ci)
        if conventions is not None else _marker_narrative(marker_error, min_patients_for_ci))

    return TcrValidationResult(
        n_cells_with_tcr=len(cells), n_clones_total=int(cells["clone_id"].nunique()),
        discordance=discordance, marker_error=marker_error, narrative=narrative,
        cell_flags=flags, flag_coverage=coverage, conventions=conventions,
        barcode_match=bmatch,
    )


# --------------------------------------------------------------------------- #
# Descrittivo: cellule con TCR per etichetta (web app)
# --------------------------------------------------------------------------- #
T_LABEL_HINTS = ("CD4", "CD8", "T cell", "Treg", "MAIT", "gdT", "T_")


def tcr_by_celltype(adata: ad.AnnData, contigs: pd.DataFrame, patient_col: str,
                    compartment_col: str, celltype_col: str, barcode_col: str) -> pd.DataFrame:
    """Per ogni etichetta di tipo cellulare: cellule totali, cellule con catena TRB rilevata e
    quante di queste appartengono a un clone espanso (>= 2 cellule nello stesso paziente).

    E' una tabella DESCRITTIVA. Un TCR in un'etichetta non-T (NK, mieloidi, stromali) puo'
    essere un doppietto, RNA ambientale o un errore di annotazione: la tabella non distingue fra
    questi casi e non classifica le singole cellule. Usa la stessa corrispondenza dei barcode
    del Modulo B (``match_barcodes``)."""
    obs = adata.obs[[patient_col, compartment_col, celltype_col, barcode_col]].astype(str).reset_index(drop=True)
    obs.columns = ["patient", "compartment", "celltype", "barcode"]
    clones = build_clonotypes(contigs)
    obs, clones, _ = match_barcodes(obs, clones)
    size = clones.groupby(["patient", "clone_id"]).barcode.transform("size")
    clones = clones.assign(expanded=size >= 2)
    j = obs.merge(clones, on=["patient", "compartment", "barcode"], how="left")
    out = j.groupby("celltype").agg(
        cellule=("barcode", "size"),
        con_TRB=("clone_id", lambda s: int(s.notna().sum())),
        in_cloni_espansi=("expanded", lambda s: int(s.fillna(False).astype(bool).sum())),
    )
    out["% con TRB"] = (100 * out["con_TRB"] / out["cellule"]).round(1)
    out["etichetta T"] = [any(h.lower() in str(c).lower() for h in T_LABEL_HINTS) for c in out.index]
    return out.sort_values(["etichetta T", "% con TRB"], ascending=[True, False])
