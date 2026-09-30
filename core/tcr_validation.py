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
    """Scrive ``<prefisso>_audited.h5ad`` (COPIA di adata con le due colonne di flag in
    obs e le definizioni in uns['genomic_audit_flags']) e ``<prefisso>_audit_flags.csv``
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
    copy.write_h5ad(h5ad_path)
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
) -> TcrValidationResult:
    """Esegue il Modulo B. ``contigs`` e' l'output gia' filtrato di ``parse_vdj_contigs``.

    ``marker_map`` e ``reference_compartment`` sono opzionali: se entrambi forniti, viene
    calcolato anche il tasso d'errore per compartimento via marcatori canonici.
    """
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

    narrative = _discordance_narrative(discordance, min_patients_for_ci) + \
        _marker_narrative(marker_error, min_patients_for_ci)

    return TcrValidationResult(
        n_cells_with_tcr=len(cells), n_clones_total=int(cells["clone_id"].nunique()),
        discordance=discordance, marker_error=marker_error, narrative=narrative,
        cell_flags=flags, flag_coverage=coverage,
    )
