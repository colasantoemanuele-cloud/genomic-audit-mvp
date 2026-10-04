"""Generatori di dati sintetici per la demo e per i test.

Modulo A: dataset con struttura per paziente e un effetto di leakage iniettato
deliberatamente -- un "fingerprint" genico casuale per paziente, non generalizzabile a
pazienti mai visti, che uno split casuale sulle cellule può sfruttare (le cellule dello
stesso paziente finiscono sia in train sia in test) ma uno split per paziente no.

Modulo B: cloni TCR con discordanza cross-compartimento iniettata deliberatamente in un
compartimento "rumoroso", ed espressione dei geni marcatori legata all'identità VERA del
clone (non all'etichetta assegnata, che può essere sbagliata) -- così il tasso d'errore
via marcatori ha qualcosa di reale da rilevare.
"""

from __future__ import annotations

from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse as sp


# --------------------------------------------------------------------------- #
# Modulo A
# --------------------------------------------------------------------------- #
def make_leakage_dataset(
    n_patients: int = 20,
    cells_per_patient: int = 40,
    n_genes: int = 200,
    n_fingerprint_genes: int = 30,
    fingerprint_strength: float = 2.5,
    n_signal_genes: int = 10,
    true_effect: float = 0.3,
    seed: int = 0,
    patient_col: str = "patient_id",
    label_col: str = "label",
) -> ad.AnnData:
    """Ogni paziente ha un'etichetta assegnata a caso (nessun legame biologico vero con
    l'identità del paziente) e un "fingerprint" genico casuale forte e specifico di
    quel paziente, condiviso da tutte le sue cellule. Un classificatore che vede cellule
    dello stesso paziente sia in train sia in test (split casuale) può riconoscere il
    fingerprint e "ricordare" l'etichetta di quel paziente -- senza aver imparato nulla
    di generalizzabile. Un piccolo segnale biologico vero e debole (n_signal_genes,
    effetto true_effect, condiviso fra pazienti della stessa classe) è anche presente,
    così lo split per paziente non è a livello del caso puro.
    """
    rng = np.random.default_rng(seed)
    classes = np.array(["ClassA", "ClassB"])
    patient_ids = [f"PT{i:03d}" for i in range(n_patients)]
    patient_label = classes[rng.integers(0, 2, n_patients)]

    signal_idx = rng.choice(n_genes, n_signal_genes, replace=False)
    fp_pool = np.setdiff1d(np.arange(n_genes), signal_idx)
    base_log_mean = rng.normal(1.0, 0.5, n_genes)

    rows_X, rows_patient, rows_label = [], [], []
    for pid, lab in zip(patient_ids, patient_label):
        fp_idx = rng.choice(fp_pool, n_fingerprint_genes, replace=False)
        fingerprint = np.zeros(n_genes)
        fingerprint[fp_idx] = rng.normal(0, fingerprint_strength, n_fingerprint_genes)

        signal = np.zeros(n_genes)
        signal[signal_idx] = true_effect if lab == "ClassA" else -true_effect

        for _ in range(cells_per_patient):
            cell_noise = rng.normal(0, 0.3, n_genes)
            log_mean = base_log_mean + fingerprint + signal + cell_noise
            mean = np.exp(np.clip(log_mean, -5, 8))
            counts = rng.poisson(mean)
            rows_X.append(counts)
            rows_patient.append(pid)
            rows_label.append(lab)

    X = sp.csr_matrix(np.vstack(rows_X).astype(np.float32))
    obs = pd.DataFrame({patient_col: rows_patient, label_col: rows_label})
    var = pd.DataFrame(index=[f"gene_{i}" for i in range(n_genes)])
    return ad.AnnData(X=X, obs=obs, var=var)


# --------------------------------------------------------------------------- #
# Modulo B
# --------------------------------------------------------------------------- #
def make_tcr_validation_dataset(
    n_patients: int = 12,
    n_clones_per_patient: int = 15,
    cells_per_clone_range: tuple[int, int] = (3, 10),
    compartments: tuple[str, str] = ("PBMC", "Tumor"),
    noisy_compartment: str = "Tumor",
    baseline_mislabel_rate: float = 0.05,
    injected_excess: float = 0.35,
    identities: tuple[str, str] = ("CD4T", "CD8T"),
    distractor_label: str = "NK",
    marker_genes: dict[str, list[str]] | None = None,
    seed: int = 0,
    patient_col: str = "patient_id",
    compartment_col: str = "tissue",
    celltype_col: str = "celltype",
    barcode_col: str = "barcode",
) -> tuple[ad.AnnData, pd.DataFrame]:
    """Ritorna (adata, contigs) già nel formato atteso da
    ``core.tcr_validation.run_tcr_validation`` (contigs ha lo stesso schema
    dell'output di ``parse_vdj_contigs``: patient, compartment, barcode, chain,
    cdr3_nt, umis).

    Ogni clone ha un'identità VERA (CD4T o CD8T) fissata alla nascita del clone.
    L'etichetta assegnata alla cellula (che simula l'output del clustering) sbaglia con
    probabilità ``baseline_mislabel_rate`` in tutti i compartimenti, PIÙ
    ``injected_excess`` nel compartimento rumoroso. L'espressione dei geni marcatori
    segue l'identità VERA del clone, non l'etichetta assegnata -- come nella realtà,
    dove il trascrittoma riflette la biologia vera indipendentemente da come il
    clustering ha etichettato la cellula.
    """
    if marker_genes is None:
        marker_genes = {"CD4T": ["CD4"], "CD8T": ["CD8A", "CD8B"]}
    rng = np.random.default_rng(seed)
    all_labels = list(identities) + [distractor_label]

    marker_gene_names = sorted({g for genes in marker_genes.values() for g in genes})
    n_background = 30
    gene_names = marker_gene_names + [f"bg_{i}" for i in range(n_background)]
    gene_idx = {g: i for i, g in enumerate(gene_names)}
    n_genes = len(gene_names)

    obs_patient, obs_comp, obs_celltype, obs_barcode = [], [], [], []
    contig_rows = []
    x_rows = []
    barcode_counter = 0

    for p_i in range(n_patients):
        pid = f"PT{p_i:03d}"
        for c_i in range(n_clones_per_patient):
            clone_id = f"CLONE_{p_i:03d}_{c_i:03d}"
            true_identity = identities[rng.integers(0, len(identities))]
            for comp in compartments:
                n_cells = rng.integers(cells_per_clone_range[0], cells_per_clone_range[1] + 1)
                mislabel_p = baseline_mislabel_rate + (injected_excess if comp == noisy_compartment else 0.0)
                for _ in range(n_cells):
                    barcode = f"BC{barcode_counter:08d}"
                    barcode_counter += 1
                    if rng.random() < mislabel_p:
                        wrong_options = [lab for lab in all_labels if lab != true_identity]
                        assigned = wrong_options[rng.integers(0, len(wrong_options))]
                    else:
                        assigned = true_identity

                    obs_patient.append(pid)
                    obs_comp.append(comp)
                    obs_celltype.append(assigned)
                    obs_barcode.append(barcode)

                    for chain, cdr3_suffix, umi_mean in (("TRA", "A", 4), ("TRB", "B", 6)):
                        contig_rows.append({
                            "patient": pid, "compartment": comp, "barcode": barcode,
                            "chain": chain, "cdr3_nt": f"{clone_id}_{cdr3_suffix}",
                            "umis": max(1, int(rng.poisson(umi_mean))),
                        })

                    expr = np.zeros(n_genes)
                    for gene in marker_gene_names:
                        is_true_marker = gene in marker_genes.get(true_identity, [])
                        mean = 6.0 if is_true_marker else 0.1
                        expr[gene_idx[gene]] = rng.poisson(mean)
                    expr[len(marker_gene_names):] = rng.poisson(1.0, n_background)
                    x_rows.append(expr)

    X = sp.csr_matrix(np.vstack(x_rows).astype(np.float32))
    obs = pd.DataFrame({
        patient_col: obs_patient, compartment_col: obs_comp,
        celltype_col: obs_celltype, barcode_col: obs_barcode,
    })
    var = pd.DataFrame(index=gene_names)
    adata = ad.AnnData(X=X, obs=obs, var=var)

    contigs = pd.DataFrame(contig_rows)
    return adata, contigs


def write_vdj_csvs(contigs: pd.DataFrame, out_dir: Path) -> list[tuple[Path, str, str]]:
    """Scrive un CSV in formato Cell Ranger per ogni (paziente, compartimento), per
    dimostrare/testare il percorso di ingestione reale (``parse_vdj_contigs``). Ritorna
    la lista (percorso, paziente, compartimento) da passare a ``parse_vdj_contigs``."""
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    files = []
    for (pat, comp), g in contigs.groupby(["patient", "compartment"]):
        df = g.copy()
        df["is_cell"] = "True"
        df["high_confidence"] = "True"
        df["productive"] = "True"
        df["full_length"] = "True"
        path = out_dir / f"{pat}_{comp}_filtered_contig_annotations.csv"
        df[["barcode", "chain", "cdr3_nt", "umis", "is_cell", "high_confidence",
            "productive", "full_length"]].to_csv(path, index=False)
        files.append((path, pat, comp))
    return files


# --------------------------------------------------------------------------- #
# Audit del disegno (Intervento 1)
# --------------------------------------------------------------------------- #
def make_design_sheet_independent(
    n_patients: int = 10,
    samples_per_patient: int = 8,
    factor_levels: dict[str, int] | None = None,
    seed: int = 0,
) -> pd.DataFrame:
    """Tabella di campioni con fattori INDIPENDENTI fra loro e dal paziente: ogni
    campione riceve un livello estratto a caso, in modo uniforme e indipendente, per
    ciascun fattore. Con ``samples_per_patient`` alto un annidamento casuale (TUTTI i
    campioni di TUTTI i pazienti nello stesso livello) è di fatto impossibile. Una
    colonna ``library`` identifica ogni riga (un campione = una libreria)."""
    if factor_levels is None:
        factor_levels = {"tissue": 3, "batch": 4, "chemistry": 2, "protocol": 2}
    rng = np.random.default_rng(seed)
    rows = []
    for p in range(n_patients):
        for s in range(samples_per_patient):
            row = {"patient": f"P{p:02d}", "library": f"L{p:02d}_{s:02d}"}
            for name, k in factor_levels.items():
                row[name] = f"{name}{int(rng.integers(0, k))}"
            rows.append(row)
    return pd.DataFrame(rows)


def make_design_sheet_nested(
    n_patients: int = 10,
    samples_per_patient: int = 4,
    patients_per_batch: int = 2,
    seed: int = 0,
) -> pd.DataFrame:
    """Annidamento iniettato: ogni paziente è processato interamente in un solo batch
    (paziente annidato nel batch), ogni batch contiene ``patients_per_batch`` pazienti.
    Il tessuto resta estratto a caso."""
    rng = np.random.default_rng(seed)
    order = rng.permutation(n_patients)
    batch_of = {int(p): f"B{i // patients_per_batch}" for i, p in enumerate(order)}
    rows = []
    for p in range(n_patients):
        for s in range(samples_per_patient):
            rows.append({"patient": f"P{p:02d}", "library": f"L{p:02d}_{s:02d}",
                         "batch": batch_of[p],
                         "tissue": ("Tumor", "Normal")[int(rng.integers(0, 2))]})
    return pd.DataFrame(rows)


def make_design_sheet_outcome_confounded(
    n_patients: int = 10,
    samples_per_patient: int = 3,
    seed: int = 0,
) -> pd.DataFrame:
    """Confondimento completo paziente-esito: ogni paziente ha un solo valore dell'esito
    (es. pazienti che hanno solo tumore, altri solo normale). Entrambi i valori sono
    presenti nella coorte (almeno un paziente per valore)."""
    rng = np.random.default_rng(seed)
    labels = np.array(["Tumor"] * (n_patients // 2) + ["Normal"] * (n_patients - n_patients // 2))
    labels = rng.permutation(labels)
    rows = []
    for p in range(n_patients):
        for s in range(samples_per_patient):
            rows.append({"patient": f"P{p:02d}", "library": f"L{p:02d}_{s:02d}",
                         "batch": f"B{int(rng.integers(0, 3))}", "tissue": labels[p]})
    return pd.DataFrame(rows)


def make_design_sheet_association(
    target_v: float,
    n_samples: int = 60,
    n_levels: int = 3,
    seed: int = 0,
) -> pd.DataFrame:
    """Due fattori con ``n_levels`` livelli ciascuno e Cramér V di POPOLAZIONE pari a
    ``target_v``: A uniforme; B = A con probabilità ``target_v``, altrimenti uniforme e
    indipendente. Per la mistura p_ij = s*delta_ij/k + (1-s)/k^2 vale phi^2 = s^2 (k-1),
    quindi V = s esattamente."""
    rng = np.random.default_rng(seed)
    a = rng.integers(0, n_levels, n_samples)
    copy = rng.random(n_samples) < target_v
    b = np.where(copy, a, rng.integers(0, n_levels, n_samples))
    return pd.DataFrame({"library": [f"L{i:03d}" for i in range(n_samples)],
                         "factor_a": [f"a{x}" for x in a], "factor_b": [f"b{x}" for x in b]})


# Struttura delle librerie di GSE278694 (dalla tabella qc_by_library_sc.csv di pdac-ml):
# scRNA-seq, 14 pazienti, 29 librerie, ogni libreria = una coppia paziente-tessuto;
# 5 pazienti (04, 06, 08, 09, 10) con tumore e tessuto adiacente.
_GSE278694_SC_LIBRARIES = {
    "PDAC01": ("Tumor",), "PDAC02": ("Tumor",),
    "PDAC04": ("Tumor", "Adjacent_normal"),
    "PDAC05": ("Tumor", "PBMC"),
    "PDAC06": ("Tumor", "Adjacent_normal", "PBMC"),
    "PDAC07": ("Tumor", "PBMC"),
    "PDAC08": ("Tumor", "Adjacent_normal", "PBMC"),
    "PDAC09": ("Tumor", "Adjacent_normal", "PBMC"),
    "PDAC10": ("Tumor", "Adjacent_normal", "PBMC"),
    "PDAC11": ("Tumor", "PBMC"), "PDAC12": ("Tumor", "PBMC"),
    "PDAC13": ("Tumor", "PBMC"), "PDAC14": ("Tumor", "PBMC"),
    "PDAC15": ("PBMC",),
}


def make_gse278694_like_sheet(n_sn_patients: int = 8) -> pd.DataFrame:
    """Riproduce la STRUTTURA del disegno di GSE278694 (nessun dato di espressione):
    coorte scRNA-seq come sopra, più una coorte snRNA-seq di ``n_sn_patients`` pazienti
    DISGIUNTI (nessun paziente ha entrambe le modalità), un campione tumorale ciascuno.
    Gli identificativi dei pazienti snRNA sono fittizi."""
    rows = []
    for pat, tissues in _GSE278694_SC_LIBRARIES.items():
        for t in tissues:
            rows.append({"patient": pat, "tissue": t, "protocol": "scRNA",
                         "library": f"{pat}-{t}-GEX"})
    for i in range(n_sn_patients):
        pat = f"SN{i + 1:02d}"
        rows.append({"patient": pat, "tissue": "Tumor", "protocol": "snRNA",
                     "library": f"{pat}-Tumor-snRNA"})
    return pd.DataFrame(rows)


def make_variance_units(
    n_patients: int = 12,
    tissues: tuple[str, ...] = ("Tumor", "Normal"),
    shares: tuple[float, float, float] = (0.5, 0.3, 0.2),
    n_genes: int = 300,
    seed: int = 0,
) -> tuple[np.ndarray, pd.DataFrame]:
    """Matrice unità x geni (già su scala log) con quote di varianza NOTE per gene:
    ``shares`` = (paziente, tessuto, residuo). Effetto paziente ~ N(0, s_p) per gene;
    effetto tessuto a media nulla fra i livelli con varianza fra livelli pari a s_t (con
    2 tessuti: +-sqrt(s_t), segno casuale per gene); residuo ~ N(0, s_e). Disegno
    bilanciato: ogni paziente ha un'unità per tessuto."""
    s_p, s_t, s_e = shares
    rng = np.random.default_rng(seed)
    k = len(tissues)
    # effetti di tessuto centrati con varianza (1/k) sum(alpha^2) = s_t
    raw = rng.normal(0, 1, (k, n_genes))
    raw -= raw.mean(axis=0, keepdims=True)
    raw /= np.sqrt((raw ** 2).mean(axis=0, keepdims=True))
    tissue_eff = raw * np.sqrt(s_t)
    patient_eff = rng.normal(0, np.sqrt(s_p), (n_patients, n_genes))
    rows, Y = [], []
    for p in range(n_patients):
        for ti, t in enumerate(tissues):
            rows.append({"patient": f"P{p:02d}", "tissue": t})
            Y.append(patient_eff[p] + tissue_eff[ti] + rng.normal(0, np.sqrt(s_e), n_genes))
    return np.vstack(Y), pd.DataFrame(rows)


def make_pseudobulk_adata(
    n_patients: int = 10,
    tissues: tuple[str, ...] = ("Tumor", "Normal"),
    cells_per_unit: int = 60,
    n_genes: int = 400,
    shares: tuple[float, float, float] = (0.5, 0.3, 0.2),
    seed: int = 0,
) -> ad.AnnData:
    """AnnData a livello di cellula le cui medie di espressione per unità
    paziente-tessuto seguono le quote di ``make_variance_units`` (scala log naturale),
    per il test end-to-end della pseudobulk. Una colonna ``library`` 1:1 con la coppia
    paziente-tessuto, come in GSE278694."""
    rng = np.random.default_rng(seed)
    Y, units = make_variance_units(n_patients, tissues, shares, n_genes, seed=seed)
    base = rng.normal(2.0, 0.5, n_genes)
    X_rows, obs = [], []
    for u, row in units.iterrows():
        mean = np.exp(base + 0.6 * Y[u])
        counts = rng.poisson(np.tile(mean, (cells_per_unit, 1)))
        X_rows.append(counts)
        obs += [{"patient_id": row.patient, "tissue": row.tissue,
                 "library": f"{row.patient}-{row.tissue}"}] * cells_per_unit
    X = sp.csr_matrix(np.vstack(X_rows).astype(np.float32))
    return ad.AnnData(X=X, obs=pd.DataFrame(obs),
                      var=pd.DataFrame(index=[f"g{i}" for i in range(n_genes)]))


# --------------------------------------------------------------------------- #
# Propagazione dell'errore sulla frazione di CD8 (Intervento 3)
# --------------------------------------------------------------------------- #
def make_cd8_fraction_dataset(
    n_patients: int = 10,
    p_cd4_to_cd8: float = 0.15,
    p_cd8_to_cd4: float = 0.15,
    p_to_other: float = 0.05,
    cells_range: tuple[int, int] = (150, 300),
    fraction_range: tuple[float, float] = (0.2, 0.7),
    reference_fraction: float = 0.3,
    seed: int = 0,
    reference_error_factor: float = 1.0,
) -> tuple[pd.DataFrame, dict[str, float]]:
    """Cellule T del compartimento 'Tumor' con frazione VERA di CD8 nota per paziente
    (``truth``: paziente -> frazione latente) ed errore di annotazione iniettato con
    direzione: una cellula vera CD4 è chiamata CD8T con prob. ``p_cd4_to_cd8``, "NK"
    (altro) con prob. ``p_to_other``; simmetricamente per le vere CD8. Una frazione
    ``reference_fraction`` delle cellule, scelta a caso e indipendentemente da identità ed
    errore, ha un'identità di riferimento (= identità vera: simula i cloni condivisi con
    il sangue, con l'assunzione che siano rappresentativi). Colonne come l'output dei flag
    del Modulo B: patient, compartment, celltype, audit_reference_label.

    ``reference_error_factor`` (B3) rompe di proposito quell'assunzione: l'errore delle
    cellule CON riferimento è moltiplicato per questo fattore, quello delle altre no."""
    rng = np.random.default_rng(seed)
    rows, truth = [], {}
    for p in range(n_patients):
        pid = f"P{p:02d}"
        pi = float(rng.uniform(*fraction_range))
        truth[pid] = pi
        n = int(rng.integers(cells_range[0], cells_range[1] + 1))
        is_cd8 = rng.random(n) < pi
        u = rng.random(n)
        has_ref_draw = rng.random(n)
        has_ref = has_ref_draw < reference_fraction
        # errore delle cellule con riferimento moltiplicato per reference_error_factor (B3):
        # nessuna estrazione casuale in più, quindi con fattore 1 i dati sono identici a prima
        f = np.where(has_ref, reference_error_factor, 1.0)
        e84, e48, eo = p_cd8_to_cd4 * f, p_cd4_to_cd8 * f, p_to_other * f
        called = np.where(
            is_cd8,
            np.where(u < e84, "CD4T", np.where(u < e84 + eo, "NK", "CD8T")),
            np.where(u < e48, "CD8T", np.where(u < e48 + eo, "NK", "CD4T")),
        )
        ref = np.where(is_cd8, "CD8T", "CD4T").astype(object)
        ref[~has_ref] = pd.NA
        for c, r in zip(called, ref):
            rows.append({"patient": pid, "compartment": "Tumor", "celltype": c,
                         "audit_reference_label": r})
    return pd.DataFrame(rows), truth
