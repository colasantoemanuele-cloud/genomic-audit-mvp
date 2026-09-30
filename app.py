"""Interfaccia Streamlit -- file sottile: tutta la logica sta in core/, qui solo upload,
scelte dell'utente e visualizzazione. Eseguibile con: streamlit run app.py
"""

from __future__ import annotations

import json
import tempfile
from pathlib import Path

import anndata as ad
import pandas as pd
import streamlit as st

from core.leakage_audit import run_leakage_audit
from core.report import render_report
from core.synthetic import make_leakage_dataset, make_tcr_validation_dataset
from core.tcr_validation import parse_vdj_contigs, run_tcr_validation

st.set_page_config(page_title="Audit genomico — coorti piccole", layout="wide")
st.title("Audit genomico per coorti cliniche piccole")
st.caption(
    "Prototipo per un singolo studio pilota. Tutto gira in locale: nessun dato lascia questa macchina."
)

for key in ("leakage_result", "tcr_result"):
    if key not in st.session_state:
        st.session_state[key] = None


@st.cache_resource(show_spinner=False)
def _read_h5ad(raw_bytes: bytes) -> ad.AnnData:
    with tempfile.NamedTemporaryFile(suffix=".h5ad", delete=False) as tmp:
        tmp.write(raw_bytes)
        tmp_path = tmp.name
    return ad.read_h5ad(tmp_path)


use_demo = st.checkbox(
    "Usa dati sintetici di esempio (nessun upload richiesto)",
    help="Genera al volo due dataset sintetici, uno per modulo, con un effetto noto "
         "gia' iniettato -- utile per vedere subito cosa fa lo strumento.",
)

adata = None
if use_demo:
    st.info(
        "Modalita' demo: il Modulo A e il Modulo B usano ciascuno un dataset sintetico "
        "dedicato, costruito apposta per mostrare chiaramente l'effetto che quel modulo misura."
    )
else:
    uploaded = st.file_uploader("File AnnData (.h5ad)", type=["h5ad"])
    if uploaded is not None:
        with st.spinner("Carico il file..."):
            adata = _read_h5ad(uploaded.getvalue())
        st.success(f"{adata.n_obs:,} cellule, {adata.n_vars:,} geni.")

tab_a, tab_b = st.tabs(["Modulo A — Leakage per paziente", "Modulo B — Validazione via TCR"])

# --------------------------------------------------------------------------- #
with tab_a:
    st.markdown(
        "Confronta una valutazione onesta (split per paziente) con un controllo negativo "
        "(split casuale sulle cellule, che ignora il paziente) per un task di classificazione "
        "a tua scelta."
    )
    if use_demo:
        leakage_adata = make_leakage_dataset(seed=0)
        target_col, patient_col = "label", "patient_id"
        st.caption(f"Dataset sintetico: {leakage_adata.n_obs:,} cellule, {len(leakage_adata.obs.patient_id.unique())} pazienti.")
    elif adata is not None:
        leakage_adata = adata
        cols = list(adata.obs.columns)
        target_col = st.selectbox("Colonna target (etichetta da classificare)", cols, key="target_col")
        patient_col = st.selectbox("Colonna identificativo paziente", cols, key="patient_col")
    else:
        leakage_adata = None
        st.info("Carica un file .h5ad o attiva i dati sintetici di esempio per procedere.")

    if leakage_adata is not None:
        n_folds = st.slider("Numero di fold", min_value=3, max_value=10, value=5, key="n_folds")
        if st.button("Esegui Modulo A", type="primary"):
            with st.spinner("Eseguo l'audit del leakage (puo' richiedere qualche minuto se ci sono molti pazienti)..."):
                try:
                    st.session_state.leakage_result = run_leakage_audit(
                        leakage_adata, target_col=target_col, patient_col=patient_col, n_folds=n_folds,
                    )
                except ValueError as e:
                    st.error(str(e))

    result = st.session_state.leakage_result
    if result is not None:
        st.markdown(f"**{result.narrative}**")
        c1, c2, c3 = st.columns(3)
        c1.metric("macro-F1 — split per paziente (onesto)", f"{result.grouped.mean:.3f}", f"± {result.grouped.std:.3f}")
        c2.metric("macro-F1 — split casuale (controllo negativo)", f"{result.random.mean:.3f}", f"± {result.random.std:.3f}")
        c3.metric("Divario sulla media", f"{result.gap:+.3f}")

        if result.model_comparison is not None:
            st.subheader("Confronto fra modelli (LeaveOneGroupOut)")
            mc = result.model_comparison
            rows = [{"modello": name, "macro-F1 medio": s.mean, "dev. standard": s.std}
                    for name, s in mc.scores.items()]
            st.dataframe(pd.DataFrame(rows).sort_values("macro-F1 medio", ascending=False), hide_index=True)
            st.caption(f"Modello migliore: **{mc.best_model}**")
            comp_rows = [{"modello": c.model, "Δ vs migliore": c.delta_vs_best,
                          "p (Nadeau-Bengio)": c.p_nadeau_bengio, "p (Wilcoxon)": c.p_wilcoxon,
                          "nota": c.note} for c in mc.comparisons]
            st.dataframe(pd.DataFrame(comp_rows), hide_index=True)
        else:
            st.caption(
                f"Confronto multi-modello non eseguito (servono almeno 8 pazienti; questo "
                f"dataset ne ha {result.n_patients})."
            )

# --------------------------------------------------------------------------- #
with tab_b:
    st.markdown(
        "Misura, tramite il repertorio T-cell receptor, quanto le etichette di tipo cellulare "
        "assegnate dal clustering sono coerenti per uno stesso clone attraverso i compartimenti "
        "tissutali (es. sangue vs tumore)."
    )
    if use_demo:
        tcr_adata, tcr_contigs = make_tcr_validation_dataset(seed=0)
        patient_col_b, compartment_col_b, celltype_col_b, barcode_col_b = (
            "patient_id", "tissue", "celltype", "barcode")
        marker_map = {"CD4T": ["CD4"], "CD8T": ["CD8A", "CD8B"]}
        reference_compartment = "PBMC"
        st.caption(
            f"Dataset sintetico: {tcr_adata.n_obs:,} cellule, discordanza cross-compartimento "
            f"iniettata nel compartimento 'Tumor'."
        )
        run_ready = True
    elif adata is not None:
        cols = list(adata.obs.columns)
        patient_col_b = st.selectbox("Colonna paziente", cols, key="patient_col_b")
        compartment_col_b = st.selectbox("Colonna compartimento tissutale", cols, key="compartment_col_b")
        celltype_col_b = st.selectbox("Colonna etichetta di tipo cellulare", cols, key="celltype_col_b")
        barcode_col_b = st.selectbox("Colonna barcode cellula", cols, key="barcode_col_b")

        st.markdown("**File VDJ Cell Ranger** (uno o piu' CSV)")
        uploaded_vdj = st.file_uploader("CSV VDJ", type=["csv"], accept_multiple_files=True, key="vdj_files")
        run_ready = False
        tcr_adata, tcr_contigs = adata, None
        if uploaded_vdj:
            manifest = pd.DataFrame({
                "file": [f.name for f in uploaded_vdj],
                "patient": ["" for _ in uploaded_vdj],
                "compartment": ["" for _ in uploaded_vdj],
            })
            st.caption("Indica paziente e compartimento per ciascun file (i CSV Cell Ranger non li contengono).")
            edited = st.data_editor(manifest, hide_index=True, key="vdj_manifest_editor")
            patient_filled = edited["patient"].fillna("").astype(str).str.strip() != ""
            compartment_filled = edited["compartment"].fillna("").astype(str).str.strip() != ""
            if patient_filled.all() and compartment_filled.all():
                with tempfile.TemporaryDirectory() as tmp_dir:
                    files = []
                    for f, row in zip(uploaded_vdj, edited.itertuples()):
                        p = Path(tmp_dir) / f.name
                        p.write_bytes(f.getvalue())
                        files.append((p, str(row.patient), str(row.compartment)))
                    tcr_contigs = parse_vdj_contigs(files)
                run_ready = True

        st.markdown("**Marcatori canonici (opzionale)** -- per il tasso d'errore per compartimento")
        marker_json = st.text_area(
            "Mappa etichetta -> geni marcatori, in JSON",
            value='{"CD8T": ["CD8A", "CD8B"], "CD4T": ["CD4"]}',
            help="Lascia vuoto per saltare questa parte opzionale.",
        )
        marker_map = None
        if marker_json.strip():
            try:
                marker_map = json.loads(marker_json)
            except json.JSONDecodeError as e:
                st.error(f"JSON non valido nella mappa marcatori: {e}")
        compartment_options = [""] + sorted(adata.obs[compartment_col_b].astype(str).unique())
        reference_compartment = st.selectbox(
            "Compartimento di riferimento (tipicamente il sangue)", compartment_options,
        ) or None
    else:
        run_ready = False
        tcr_adata = tcr_contigs = None
        st.info("Carica un file .h5ad o attiva i dati sintetici di esempio per procedere.")

    if run_ready and st.button("Esegui Modulo B", type="primary"):
        with st.spinner("Costruisco i clonotipi e calcolo la discordanza (cluster bootstrap: puo' richiedere qualche secondo)..."):
            try:
                st.session_state.tcr_result = run_tcr_validation(
                    tcr_adata, tcr_contigs, patient_col=patient_col_b, compartment_col=compartment_col_b,
                    celltype_col=celltype_col_b, barcode_col=barcode_col_b,
                    marker_map=marker_map, reference_compartment=reference_compartment,
                )
            except ValueError as e:
                st.error(str(e))

    result_b = st.session_state.tcr_result
    if result_b is not None:
        st.markdown(f"**{result_b.narrative}**")
        d = result_b.discordance
        if d.sufficient:
            c1, c2 = st.columns(2)
            c1.metric("Eccesso di discordanza", f"{d.mean_excess:+.3f}")
            c2.metric("IC95%", f"[{d.ci_low:+.3f}, {d.ci_high:+.3f}]")
            st.dataframe(d.by_compartment_pair, hide_index=True)
        else:
            st.warning("Numerosita' insufficiente per una stima affidabile (servono almeno 5 pazienti).")

        if result_b.marker_error is not None and result_b.marker_error.by_compartment:
            st.subheader("Tasso d'errore per compartimento (marcatori canonici)")
            rows = [{"compartimento": comp, "tasso d'errore": r.mean, "IC95% basso": r.ci_low,
                     "IC95% alto": r.ci_high, "pazienti": r.n_groups, "sufficiente": r.sufficient}
                    for comp, r in result_b.marker_error.by_compartment.items()]
            st.dataframe(pd.DataFrame(rows), hide_index=True)

# --------------------------------------------------------------------------- #
st.divider()
if st.session_state.leakage_result is not None or st.session_state.tcr_result is not None:
    report_html = render_report(
        leakage_result=st.session_state.leakage_result, tcr_result=st.session_state.tcr_result,
        dataset_name="demo sintetico" if use_demo else "dataset caricato",
    )
    st.download_button("Scarica report HTML completo", report_html, file_name="report_audit.html",
                        mime="text/html")
