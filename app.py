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

from core.cd8_propagation import EXPERIMENTAL_NOTE, cd8_fraction_intervals
from core.design_audit import TECHNICAL_ROLES, run_design_audit, sample_sheet_from_obs
from core.leakage_audit import run_leakage_audit
from core.report import render_report
from core.synthetic import (
    make_gse278694_like_sheet,
    make_leakage_dataset,
    make_tcr_validation_dataset,
)
from core.tcr_validation import export_audited, parse_vdj_contigs, run_tcr_validation

st.set_page_config(page_title="Audit genomico — coorti piccole", layout="wide")
st.title("Audit genomico per coorti cliniche piccole")
st.caption(
    "Prototipo per un singolo studio pilota. Tutto gira in locale: nessun dato lascia questa macchina."
)

for key in ("design_result", "leakage_result", "tcr_result", "cd8_result"):
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

tab_d, tab_a, tab_b = st.tabs(["Audit del disegno", "Modulo A — Leakage per paziente",
                               "Modulo B — Validazione via TCR"])

# --------------------------------------------------------------------------- #
with tab_d:
    st.markdown(
        "Legge **solo i metadati** (una riga per campione o libreria) e segnala quali fattori "
        "sono confusi fra loro e quali confronti il disegno permette davvero, con quante unita' "
        "indipendenti. Non usa l'espressione genica."
    )
    sheet = None
    if use_demo:
        sheet = make_gse278694_like_sheet()
        st.caption("Metadati sintetici con la struttura di GSE278694: 14 pazienti scRNA-seq "
                   "(una libreria per coppia paziente-tessuto) e 8 pazienti snRNA-seq disgiunti.")
        d_patient, d_tissue = "patient", "tissue"
        d_technical = {"protocol": "protocol", "library": "library"}
        d_outcomes: list[str] = []
        d_comparisons = [("tissue", "Tumor", "Adjacent_normal"), ("protocol", "scRNA", "snRNA")]
    else:
        meta_file = st.file_uploader("CSV dei metadati (opzionale se hai caricato un .h5ad)",
                                     type=["csv"], key="meta_csv")
        if meta_file is not None:
            sheet = pd.read_csv(meta_file, dtype=str)
        elif adata is not None:
            sheet = adata.obs.astype(str)
        if sheet is not None:
            cols = list(sheet.columns)
            d_patient = st.selectbox("Colonna paziente", cols, key="d_patient")
            d_tissue = st.selectbox("Colonna tessuto/compartimento (opzionale)", [""] + cols,
                                    key="d_tissue") or None
            d_technical = {}
            for role in TECHNICAL_ROLES:
                c = st.selectbox(f"Colonna '{role}' (opzionale)", [""] + cols, key=f"d_{role}")
                if c:
                    d_technical[role] = c
            d_outcomes = st.multiselect("Colonne di esito (opzionali)", cols, key="d_outcomes")
            comp_text = st.text_area(
                "Confronti da valutare, uno per riga: colonna:livello_a:livello_b",
                value="", key="d_comparisons")
            d_comparisons = []
            for line in comp_text.splitlines():
                parts = [x.strip() for x in line.split(":")]
                if len(parts) == 3 and all(parts):
                    d_comparisons.append(tuple(parts))
                elif line.strip():
                    st.error(f"Riga non valida (serve colonna:livello_a:livello_b): {line}")
            used = [d_patient] + ([d_tissue] if d_tissue else []) + list(d_technical.values()) + d_outcomes
            if meta_file is None:
                sheet = sample_sheet_from_obs(sheet, list(dict.fromkeys(used)))
                st.caption(f"Metadati ridotti da cellule a {len(sheet)} combinazioni distinte dei fattori scelti.")
        else:
            st.info("Carica un CSV di metadati o un file .h5ad, oppure attiva i dati sintetici.")

    if sheet is not None and st.button("Esegui audit del disegno", type="primary"):
        try:
            st.session_state.design_result = run_design_audit(
                sheet, patient_col=d_patient, tissue_col=d_tissue, technical_cols=d_technical,
                outcome_cols=d_outcomes, comparisons=d_comparisons,
            )
        except ValueError as e:
            st.error(str(e))

    result_d = st.session_state.design_result
    if result_d is not None:
        if result_d.comparisons:
            st.subheader("Confronti richiesti")
            st.dataframe(pd.DataFrame([
                {"confronto": f"{c.factor}: {c.level_a} vs {c.level_b}", "classe": c.classification,
                 "disegno": c.design, "unita' indipendenti": c.n_units, "spiegazione": c.sentence}
                for c in result_d.comparisons]), hide_index=True)
        st.subheader("Fatti strutturali")
        structural = [f for f in result_d.findings
                      if f.kind in ("annidamento", "uno-a-uno", "esito-determinato", "unita'-tecnica")]
        for f in structural:
            st.markdown(f"- {f.sentence}")
        if not structural:
            st.markdown("Nessun annidamento, coincidenza o esito determinato da un singolo fattore.")
        st.subheader("Associazione fra coppie di fattori")
        st.dataframe(pd.DataFrame([
            {"coppia": f"{p.factor_a} × {p.factor_b}",
             "Cramér V": "non valutabile" if p.cramer_v is None else f"{p.cramer_v:.3f}",
             "spiegazione": p.sentence} for p in result_d.pairs]), hide_index=True)
        for note in result_d.notes:
            st.caption(note)

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
        if result_b.barcode_match is not None:
            st.caption(result_b.barcode_match.sentence)
        st.markdown(f"**{result_b.narrative}**")
        d = result_b.discordance
        if d.sufficient:
            c1, c2 = st.columns(2)
            c1.metric("Eccesso di discordanza", f"{d.mean_excess:+.3f}")
            c2.metric("IC95%", f"[{d.ci_low:+.3f}, {d.ci_high:+.3f}]")
            st.dataframe(d.by_compartment_pair, hide_index=True)
        else:
            st.warning("Numerosita' insufficiente per una stima affidabile (servono almeno 5 pazienti).")

        if result_b.conventions:
            st.subheader("Tasso d'errore per compartimento (marcatori canonici)")
            st.markdown("Ogni numero e' riportato con la sua convenzione; le convenzioni sono affiancate.")
            for c in result_b.conventions.values():
                st.markdown(f"- **{c.name}**: {c.definition}")
            rows = []
            for c in result_b.conventions.values():
                for comp, r in c.by_compartment.items():
                    rows.append({"convenzione": c.name, "compartimento": comp, "tasso d'errore": f"{r.mean:.3f}",
                                 "IC95%": f"[{r.ci_low:.3f}, {r.ci_high:.3f}]" if r.sufficient
                                 else f"non prodotto ({r.n_groups} pazienti)", "pazienti": r.n_groups})
                for (a, b), (r, n_cl) in c.paired_differences.items():
                    rows.append({"convenzione": c.name, "compartimento": f"{a} − {b} (stessi cloni, {n_cl})",
                                 "tasso d'errore": f"{r.mean:+.3f}",
                                 "IC95%": f"[{r.ci_low:+.3f}, {r.ci_high:+.3f}]" if r.sufficient
                                 else f"non prodotto ({r.n_groups} pazienti)", "pazienti": r.n_groups})
            st.dataframe(pd.DataFrame(rows), hide_index=True)

        if result_b.cell_flags is not None:
            st.subheader("Flag per cellula")
            st.markdown(
                f"Cellule valutabili: **{result_b.flag_coverage:.1%}**. Le altre sono NA: non "
                "verificabili (nessun TCR, clone senza riferimento, etichetta fuori dalla mappa "
                "dei marcatori, o compartimento di riferimento). NA non significa \"corretta\"; "
                "le etichette originali non vengono modificate."
            )
            st.download_button(
                "Scarica i flag per cellula (CSV)",
                result_b.cell_flags.to_csv(index_label="obs_name"),
                file_name="audit_flags.csv", mime="text/csv",
            )
            if st.button("Prepara la copia .h5ad con i flag"):
                with tempfile.TemporaryDirectory() as tmp_dir:
                    try:
                        h5ad_path, _ = export_audited(tcr_adata, result_b, Path(tmp_dir) / "dati")
                        st.session_state.audited_h5ad = h5ad_path.read_bytes()
                    except ValueError as e:
                        st.error(str(e))
            if st.session_state.get("audited_h5ad"):
                st.download_button("Scarica dati_audited.h5ad", st.session_state.audited_h5ad,
                                   file_name="dati_audited.h5ad")

            st.subheader("Frazione di CD8: effetto dell'errore di annotazione")
            st.warning(EXPERIMENTAL_NOTE)
            comp_options = sorted(tcr_adata.obs[compartment_col_b].astype(str).unique())
            cd8_comp = st.selectbox("Compartimento", comp_options,
                                    index=comp_options.index("Tumor") if "Tumor" in comp_options else 0,
                                    key="cd8_comp")
            labels = list(marker_map) if marker_map else ["CD4T", "CD8T"]
            c1, c2 = st.columns(2)
            cd4_label = c1.selectbox("Etichetta CD4", labels, index=0, key="cd4_label")
            cd8_label = c2.selectbox("Etichetta CD8", labels, index=min(1, len(labels) - 1), key="cd8_label")
            if st.button("Calcola gli intervalli sulla frazione di CD8"):
                obs_flags = tcr_adata.obs.join(result_b.cell_flags)
                st.session_state.cd8_result = cd8_fraction_intervals(
                    obs_flags, patient_col=patient_col_b, compartment_col=compartment_col_b,
                    celltype_col=celltype_col_b, reference_col="audit_reference_label",
                    target_compartment=cd8_comp, cd4_label=cd4_label, cd8_label=cd8_label)
            cd8_res = st.session_state.cd8_result
            if cd8_res is not None:
                if cd8_res.refused_reason:
                    st.warning(f"Intervalli non prodotti: {cd8_res.refused_reason}.")
                elif cd8_res.matrix is not None:
                    st.caption(f"Matrice di confusione (pooled fra {cd8_res.n_reference_patients} pazienti, "
                               f"J = {cd8_res.youden_j:.2f})")
                    st.dataframe(cd8_res.matrix.round(3))

                def _iv(iv):
                    return (f"{iv.low:.2f}–{iv.high:.2f}" if iv.low is not None
                            else f"non prodotto: {iv.refused_reason}")
                st.dataframe(pd.DataFrame([
                    {"paziente": p.patient, "n (CD4+CD8)": p.n_cd4_called + p.n_cd8_called,
                     "riportata": "—" if p.reported is None else f"{p.reported:.2f}",
                     "errore 0.5x": _iv(p.scenarios[0.5]), "errore 1x": _iv(p.scenarios[1.0]),
                     "errore 2x": _iv(p.scenarios[2.0])} for p in cd8_res.patients]), hide_index=True)
                st.caption(cd8_res.assumptions)

# --------------------------------------------------------------------------- #
st.divider()
if any(st.session_state[k] is not None for k in ("design_result", "leakage_result", "tcr_result", "cd8_result")):
    report_html = render_report(
        leakage_result=st.session_state.leakage_result, tcr_result=st.session_state.tcr_result,
        dataset_name="demo sintetico" if use_demo else "dataset caricato",
        design_result=st.session_state.design_result, cd8_result=st.session_state.cd8_result,
    )
    st.download_button("Scarica report HTML completo", report_html, file_name="report_audit.html",
                        mime="text/html")
