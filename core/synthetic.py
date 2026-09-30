"""Generatori di dati sintetici per la demo e per i test.

Modulo A: dataset con struttura per paziente e un effetto di leakage iniettato
deliberatamente -- un "fingerprint" genico casuale per paziente, non generalizzabile a
pazienti mai visti, che uno split casuale sulle cellule puo' sfruttare (le cellule dello
stesso paziente finiscono sia in train sia in test) ma uno split per paziente no.

Modulo B: cloni TCR con discordanza cross-compartimento iniettata deliberatamente in un
compartimento "rumoroso", ed espressione dei geni marcatori legata all'identita' VERA del
clone (non all'etichetta assegnata, che puo' essere sbagliata) -- cosi' il tasso d'errore
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
    l'identita' del paziente) e un "fingerprint" genico casuale forte e specifico di
    quel paziente, condiviso da tutte le sue cellule. Un classificatore che vede cellule
    dello stesso paziente sia in train sia in test (split casuale) puo' riconoscere il
    fingerprint e "ricordare" l'etichetta di quel paziente -- senza aver imparato nulla
    di generalizzabile. Un piccolo segnale biologico vero e debole (n_signal_genes,
    effetto true_effect, condiviso fra pazienti della stessa classe) e' anche presente,
    cosi' lo split per paziente non è a livello del caso puro.
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
    """Ritorna (adata, contigs) gia' nel formato atteso da
    ``core.tcr_validation.run_tcr_validation`` (contigs ha lo stesso schema
    dell'output di ``parse_vdj_contigs``: patient, compartment, barcode, chain,
    cdr3_nt, umis).

    Ogni clone ha un'identita' VERA (CD4T o CD8T) fissata alla nascita del clone.
    L'etichetta assegnata alla cellula (che simula l'output del clustering) sbaglia con
    probabilita' ``baseline_mislabel_rate`` in tutti i compartimenti, PIU'
    ``injected_excess`` nel compartimento rumoroso. L'espressione dei geni marcatori
    segue l'identita' VERA del clone, non l'etichetta assegnata -- come nella realta',
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
