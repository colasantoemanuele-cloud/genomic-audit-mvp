"""Web app locale di genomic-audit. Avvio: `audit-sc serve` (oppure `python cli.py serve`).

Tutta la logica sta in core/: qui solo caricamento dei dati, scelte dell'utente e
visualizzazione. L'app gira su localhost con la telemetria di Streamlit disattivata
(.streamlit/config.toml e opzioni di `serve`): nessun dato e nessuna statistica d'uso
lasciano la macchina.
"""

from __future__ import annotations

import html
import json
import tempfile
from pathlib import Path

import altair as alt
import anndata as ad
import numpy as np
import pandas as pd
import streamlit as st

from core.cd8_propagation import EXPERIMENTAL_NOTE, cd8_fraction_intervals
from core.design_audit import TECHNICAL_ROLES, run_design_audit
from core.io import load_matrix_market, read_table
from core.leakage_audit import GAP_ALERT, run_leakage_audit
from core.report import render_markdown_report, render_report
from core.synthetic import make_tcr_validation_dataset
from core.tcr_validation import parse_vdj_contigs, run_tcr_validation, tcr_by_celltype
from core.verdict import (
    STATE_LABEL,
    TITLE,
    design_summary,
    leakage_summary,
    standing_limits,
    summarize,
    tcr_summary,
)

ROOT = Path(__file__).resolve().parent
DEMO = ROOT / "data" / "demo"

DEMO_DESIGN = {
    "GSE132465 — carcinoma colorettale (23 pazienti, tumore e mucosa normale)": dict(
        file="GSE132465_samples.csv", patient="patient_id", tissue="tissue type",
        technical={"library": "gsm", "batch": "platform"}, outcomes=["tumor stage", "region"],
        comparisons=[("tissue type", "Colorectal cancer", "Normal mucosa"), ("tumor stage", "2", "3")]),
    "GSE131907 — adenocarcinoma polmonare (44 pazienti, 7 tessuti)": dict(
        file="GSE131907_samples.csv", patient="patient id", tissue="tissue origin abbrevation",
        technical={"library": "gsm", "batch": "platform"}, outcomes=["tumor stage"],
        comparisons=[("tissue origin abbrevation", "tLung", "nLung"),
                     ("tissue origin abbrevation", "mLN", "nLN"),
                     ("tissue origin abbrevation", "tLung", "mBrain")]),
    "GSE125449 — tumori primitivi del fegato (19 pazienti, 2 piattaforme)": dict(
        file="GSE125449_samples.csv", patient="patient_from_title", tissue=None,
        technical={"library": "gsm", "batch": "platform"}, outcomes=["cancer type"],
        comparisons=[("cancer type", "Hepatocellular carcinoma", "Intrahepatic cholangiocarcinoma"),
                     ("platform", "GPL18573", "GPL20301")]),
}

st.set_page_config(page_title="Audit genomico — coorti piccole", layout="wide")
st.markdown("""
<style>
.badge{display:inline-block;padding:.15rem .6rem;border-radius:999px;font-weight:600;font-size:.85rem;color:#fff}
.b-verde{background:#2e7d4f}.b-giallo{background:#b7860b}.b-rosso{background:#b3261e}.b-grigio{background:#6b7785}
.card{background:#fff;border:1px solid #dfe4ea;border-radius:10px;padding:.9rem 1.1rem;margin:.4rem 0}
.small{color:#5b6573;font-size:.88rem}
table.sem{border-collapse:collapse;width:100%}
table.sem td,table.sem th{border-bottom:1px solid #e3e7ec;padding:.35rem .5rem;text-align:left;vertical-align:top}
</style>""", unsafe_allow_html=True)

COLOR_OF = {"stimabile": "verde", "stimabile con bassa potenza": "giallo", "non stimabile": "rosso"}


def badge(color: str, text: str) -> str:
    return f'<span class="badge b-{color}">{html.escape(text)}</span>'


for k in ("design", "leakage", "tcr", "cd8", "tcr_table", "dataset_name", "b_ctx", "a_ctx", "d_ctx",
          "autorun", "run_design", "run_a", "run_b"):
    if k not in st.session_state:
        st.session_state[k] = None


# --------------------------------------------------------------------------- #
# Barra laterale: modalità
# --------------------------------------------------------------------------- #
with st.sidebar:
    st.markdown("### Audit genomico")
    st.caption("Strumento locale per coorti cliniche piccole. Nessun dato lascia questa macchina.")
    mode = st.radio("Modalità", ["Demo immediata", "Carica studio"], key="mode")
    if mode == "Demo immediata":
        demo_design = st.selectbox("Disegno dimostrativo (metadati reali GEO)", list(DEMO_DESIGN))
        rapido_demo = st.checkbox("Modulo A rapido (senza confronto fra modelli)", value=False)
        if st.button("Carica la demo", type="primary", width="stretch"):
            st.session_state.update(design=None, leakage=None, tcr=None, cd8=None, tcr_table=None)
            cfg = DEMO_DESIGN[demo_design]
            st.session_state.d_ctx = dict(sheet=read_table(DEMO / cfg["file"]), **cfg)
            st.session_state.a_ctx = dict(adata=ad.read_h5ad(DEMO / "GSE125449_demo.h5ad"),
                                          target="Type", patient="patient", rapido=rapido_demo,
                                          label="GSE125449, sottoinsieme reale: 10 pazienti, 1.861 cellule")
            adata_b, contigs_b = make_tcr_validation_dataset(seed=0)
            st.session_state.b_ctx = dict(adata=adata_b, contigs=contigs_b, patient="patient_id",
                                          compartment="tissue", celltype="celltype", barcode="barcode",
                                          markers={"CD4T": ["CD4"], "CD8T": ["CD8A", "CD8B"]},
                                          reference="PBMC", label="dati SINTETICI (nessun dato TCR reale incluso)")
            st.session_state.dataset_name = "demo"
            st.session_state.autorun = True
        st.caption("Disegno e Modulo A usano dati reali pubblici (GEO). Il Modulo B usa dati sintetici "
                   "dichiarati: il progetto non include dati TCR reali.")
    else:
        st.caption("Indica i percorsi locali (consigliato per file grandi) oppure carica i file.")
        st.session_state.dataset_name = st.text_input("Nome dello studio", value="studio")

st.title("Audit genomico per coorti cliniche piccole")
st.markdown('<p class="small">Unità statistica indipendente: il paziente. Ogni numero è '
            'accompagnato dalla sua definizione e dai suoi limiti; lo strumento segnala, non corregge '
            'le etichette.</p>', unsafe_allow_html=True)

tab_d, tab_a, tab_b, tab_v = st.tabs(["1 · Disegno sperimentale", "2 · Modulo A — Leakage",
                                      "3 · Modulo B — Verifica TCR", f"{TITLE} e report"])


# --------------------------------------------------------------------------- #
# Caricamento dati in modalità studio
# --------------------------------------------------------------------------- #
def load_matrix_ui(key: str):
    kind = st.radio("Formato della matrice", ["AnnData (.h5ad)", "Cartella 10x (Matrix Market)"],
                    horizontal=True, key=f"{key}_kind")
    if kind.startswith("AnnData"):
        path = st.text_input("Percorso del file .h5ad", key=f"{key}_path")
        up = st.file_uploader("...oppure carica il file .h5ad", type=["h5ad"], key=f"{key}_up")
        if path:
            return ad.read_h5ad(path)
        if up is not None:
            with tempfile.NamedTemporaryFile(suffix=".h5ad", delete=False) as tmp:
                tmp.write(up.getvalue())
            return ad.read_h5ad(tmp.name)
        return None
    folder = st.text_input("Cartella con matrix.mtx, barcodes.tsv, features.tsv (anche .gz)", key=f"{key}_mtx")
    meta = st.text_input("Tabella dei metadati per cellula (CSV/TSV)", key=f"{key}_meta")
    if folder and meta:
        m = read_table(meta)
        bcol = st.selectbox("Colonna dei barcode nei metadati", list(m.columns), key=f"{key}_bcol")
        return load_matrix_market(folder, m, bcol)
    return None


def summary_line(s) -> None:
    """Conteggio dei controlli della sezione per stato: nessun colore complessivo."""
    badges = " ".join(badge(state, f"{n} × {STATE_LABEL[state]}") for state, n in s.counts.items())
    st.markdown(f'{badges}<br><span class="small">{TITLE}: un colore per controllo, dettaglio nell\'ultima '
                f'scheda. Regola: {html.escape(s.rule)}</span>', unsafe_allow_html=True)


# --------------------------------------------------------------------------- #
# Tab 1: disegno
# --------------------------------------------------------------------------- #
def factor_matrix_chart(res) -> alt.Chart:
    rows = []
    facts = {(f.factors[0], f.factors[1]): f.kind for f in res.findings if len(f.factors) == 2}
    for p in res.pairs:
        if p.one_to_one or p.a_nested_in_b or p.b_nested_in_a:
            stato = "struttura (annidamento/coincidenza/esito)"
        elif p.cramer_v is None:
            stato = "non valutabile"
        elif p.v_alarm:
            stato = "associazione forte (V >= 0.5)"
        else:
            stato = "associazione debole"
        v = "" if p.cramer_v is None else f"{p.cramer_v:.2f}"
        for a, b in ((p.factor_a, p.factor_b), (p.factor_b, p.factor_a)):
            rows.append({"fattore 1": a, "fattore 2": b, "stato": stato, "V": v,
                         "dettaglio": facts.get((a, b), facts.get((b, a), "")), "frase": p.sentence})
    df = pd.DataFrame(rows)
    scale = alt.Scale(domain=["struttura (annidamento/coincidenza/esito)", "associazione forte (V >= 0.5)",
                              "associazione debole", "non valutabile"],
                      range=["#b3261e", "#d08a1e", "#2e7d4f", "#aab4bf"])
    base = alt.Chart(df).encode(x=alt.X("fattore 1:N", title=None), y=alt.Y("fattore 2:N", title=None))
    return (base.mark_rect(stroke="white").encode(
        color=alt.Color("stato:N", scale=scale, legend=alt.Legend(orient="bottom", title=None)),
        tooltip=["fattore 1", "fattore 2", "stato", "V", "dettaglio", "frase"])
        + base.mark_text(color="white", fontWeight="bold").encode(text="V:N")).properties(height=320)


with tab_d:
    st.markdown("Legge **solo i metadati** (una riga per campione): quali fattori sono confusi e quali "
                "confronti il disegno permette, con quante unità indipendenti. Utilizzabile anche "
                "prima di sequenziare.")
    ctx = st.session_state.d_ctx if mode == "Demo immediata" else None
    if mode == "Carica studio":
        src = st.file_uploader("Metadati per campione (CSV/TSV)", type=["csv", "tsv", "txt"], key="d_up")
        path = st.text_input("...oppure percorso locale", key="d_path")
        sheet = read_table(src) if src is not None else (read_table(path) if path else None)
        if sheet is not None:
            cols = list(sheet.columns)
            c1, c2 = st.columns(2)
            patient = c1.selectbox("Colonna paziente", cols, key="d_pat")
            tissue = c2.selectbox("Colonna tessuto/condizione (opzionale)", [""] + cols, key="d_tis") or None
            technical = {}
            for col_ui, role in zip(st.columns(len(TECHNICAL_ROLES)), TECHNICAL_ROLES):
                val = col_ui.selectbox(role, [""] + cols, key=f"d_{role}")
                if val:
                    technical[role] = val
            outcomes = st.multiselect("Colonne di esito", cols, key="d_out")
            comp_txt = st.text_area("Confronti, uno per riga: colonna:livello_a:livello_b", key="d_cmp")
            comps = [tuple(x.strip() for x in line.split(":")) for line in comp_txt.splitlines()
                     if len(line.split(":")) == 3]
            if st.button("Esegui audit del disegno", type="primary"):
                st.session_state.d_ctx = dict(sheet=sheet, patient=patient, tissue=tissue, technical=technical,
                                              outcomes=outcomes, comparisons=comps)
                st.session_state.design = None
                st.session_state.run_design = True
                ctx = st.session_state.d_ctx
    if ctx is not None and st.session_state.design is None and (
            st.session_state.autorun or st.session_state.run_design):
        try:
            st.session_state.design = run_design_audit(
                ctx["sheet"], patient_col=ctx["patient"], tissue_col=ctx["tissue"],
                technical_cols=ctx["technical"], outcome_cols=ctx["outcomes"], comparisons=ctx["comparisons"])
        except ValueError as e:
            st.error(str(e))
        st.session_state.run_design = False
    res = st.session_state.design
    if res is None:
        st.info("Carica la demo dalla barra laterale oppure scegli 'Carica studio'.")
    else:
        summary_line(design_summary(res))
        st.subheader("Confronti: che cosa il disegno permette di stimare")
        rows = "".join(
            f"<tr><td>{html.escape(c.factor)}: {html.escape(c.level_a)} vs {html.escape(c.level_b)}</td>"
            f"<td>{badge(COLOR_OF[c.classification], c.classification)}</td><td>{c.n_units}</td>"
            f"<td>{'—' if c.min_pvalue is None else f'{c.min_pvalue:.3f}'}</td>"
            f"<td class='small'>{html.escape(c.sentence)}</td></tr>" for c in res.comparisons)
        st.markdown("<table class='sem'><tr><th>Confronto</th><th>Classe</th><th>Unità indipendenti</th>"
                    f"<th>p-value minimo</th><th>Spiegazione</th></tr>{rows}</table>", unsafe_allow_html=True)
        c1, c2 = st.columns([3, 2])
        with c1:
            st.subheader("Matrice dei fattori")
            if res.pairs:
                st.altair_chart(factor_matrix_chart(res), width="stretch")
            st.caption("Numero = Cramér V corretto (Bergsma); vuoto = non valutabile su una tabella troppo "
                       "piccola. Passa il mouse sulle celle per la frase completa.")
        with c2:
            st.subheader("Fatti strutturali")
            for f in res.findings:
                if f.kind in ("annidamento", "uno-a-uno", "esito-determinato", "unità-tecnica"):
                    st.markdown(f"- {f.sentence}")
            for n in res.notes:
                st.caption(n)


# --------------------------------------------------------------------------- #
# Tab 2: Modulo A
# --------------------------------------------------------------------------- #
def folds_chart(r) -> alt.Chart:
    df = pd.concat([
        pd.DataFrame({"schema": "per paziente (onesto)", "fold": range(1, len(r.grouped.fold_scores) + 1),
                      "macro-F1": r.grouped.fold_scores}),
        pd.DataFrame({"schema": "casuale (leakage)", "fold": range(1, len(r.random.fold_scores) + 1),
                      "macro-F1": r.random.fold_scores})])
    summ = df.groupby("schema")["macro-F1"].agg(["mean", "std"]).reset_index()
    summ["lo"], summ["hi"] = summ["mean"] - summ["std"], summ["mean"] + summ["std"]
    color = alt.Color("schema:N", scale=alt.Scale(domain=["per paziente (onesto)", "casuale (leakage)"],
                                                  range=["#1f3a5f", "#c0392b"]), legend=None)
    y = alt.Y("schema:N", title=None)
    pts = alt.Chart(df).mark_circle(size=90, opacity=.75).encode(
        x=alt.X("macro-F1:Q", scale=alt.Scale(zero=False)), y=y, color=color,
        tooltip=["schema", "fold", alt.Tooltip("macro-F1:Q", format=".3f")])
    bars = alt.Chart(summ).mark_rule(strokeWidth=3).encode(x="lo:Q", x2="hi:Q", y=y, color=color)
    mean = alt.Chart(summ).mark_tick(thickness=4, size=28).encode(
        x="mean:Q", y=y, color=color,
        tooltip=[alt.Tooltip("mean:Q", format=".3f", title="media"), alt.Tooltip("std:Q", format=".3f", title="dev. std")])
    return (bars + mean + pts).properties(height=180)


def jaccard_chart(m: np.ndarray, title: str) -> alt.Chart:
    n = m.shape[0]
    df = pd.DataFrame([{"fold i": f"F{i + 1}", "fold j": f"F{j + 1}", "Jaccard": float(m[i, j])}
                       for i in range(n) for j in range(n)])
    return alt.Chart(df).mark_rect().encode(
        x=alt.X("fold i:N", title=None), y=alt.Y("fold j:N", title=None),
        color=alt.Color("Jaccard:Q", scale=alt.Scale(domain=[0, 1], scheme="blues")),
        tooltip=["fold i", "fold j", alt.Tooltip("Jaccard:Q", format=".2f")]).properties(title=title, height=220)


with tab_a:
    st.markdown("Quanto una valutazione che non separa i pazienti fra training e test **sovrastima** "
                "l'accuratezza di un classificatore cellulare, e quanto ne **sottostima** l'incertezza.")
    a_ctx = st.session_state.a_ctx if mode == "Demo immediata" else None
    if mode == "Carica studio":
        adata = load_matrix_ui("a")
        if adata is not None:
            cols = list(adata.obs.columns)
            c1, c2, c3 = st.columns(3)
            target = c1.selectbox("Colonna target (etichetta da classificare)", cols, key="a_t")
            patient = c2.selectbox("Colonna paziente", cols, key="a_p")
            rapido = c3.checkbox("Rapido (senza confronto fra modelli)", key="a_r")
            if st.button("Esegui Modulo A", type="primary"):
                st.session_state.a_ctx = dict(adata=adata, target=target, patient=patient, rapido=rapido,
                                              label=st.session_state.dataset_name)
                st.session_state.leakage = None
                st.session_state.run_a = True
                a_ctx = st.session_state.a_ctx
    if a_ctx is not None and st.session_state.leakage is None and (
            st.session_state.autorun or st.session_state.run_a):
        bar = st.progress(0.0, text="Modulo A in esecuzione...")
        state = {"total": None}

        def cb(msg: str) -> None:
            if msg.startswith("Modulo A:"):
                state["total"] = int(msg.split("Addestramenti previsti: ")[1].split()[0].rstrip("."))
                bar.progress(0.0, text=msg)
            elif msg.startswith("[") and state["total"]:
                done = int(msg[1:msg.index("/")])
                bar.progress(min(done / state["total"], 1.0), text=msg)
        try:
            st.session_state.leakage = run_leakage_audit(
                a_ctx["adata"], target_col=a_ctx["target"], patient_col=a_ctx["patient"],
                benchmark=not a_ctx["rapido"], progress=cb)
        except ValueError as e:
            st.error(str(e))
        bar.empty()
        st.session_state.run_a = False
    r = st.session_state.leakage
    if r is None:
        st.info("Carica la demo dalla barra laterale oppure scegli 'Carica studio'.")
    else:
        summary_line(leakage_summary(r))
        if a_ctx:
            st.caption(f"Dataset: {a_ctx.get('label', '')}. {r.n_cells:,} cellule, {r.n_patients} pazienti, "
                       f"{r.n_classes} classi; tempo di calcolo {r.elapsed_seconds:.0f} s.")
        m1, m2, m3, m4 = st.columns(4)
        m1.metric("macro-F1 per paziente (onesta)", f"{r.grouped.mean:.3f}", f"± {r.grouped.std:.3f}", delta_color="off")
        m2.metric("macro-F1 casuale (leakage)", f"{r.random.mean:.3f}", f"± {r.random.std:.3f}", delta_color="off")
        m3.metric("Divario", f"{r.gap:+.3f}", f"soglia di allarme {GAP_ALERT:g}", delta_color="off")
        m4.metric("Incertezza sottostimata", f"{r.std_ratio:.1f}x" if np.isfinite(r.std_ratio) else "—",
                  "dev. std onesta / casuale", delta_color="off")
        st.altair_chart(folds_chart(r), width="stretch")
        st.caption("Punti = singoli fold; tacca = media; barra = media ± 1 deviazione standard fra fold. "
                   "Non è un intervallo di confidenza: i punteggi dei fold sono correlati (i training si "
                   "sovrappongono) e un intervallo non sarebbe calibrato. Definizione: macro-F1 sulle "
                   "sole classi presenti nel fold di test.")
        st.markdown(r.narrative)
        if r.grouped.n_folds_with_absent_classes:
            st.subheader("Classi assenti dai fold di test (split per paziente)")
            st.dataframe(pd.DataFrame({"fold": range(1, len(r.grouped.absent_classes) + 1),
                                       "classi assenti": [", ".join(a) or "—" for a in r.grouped.absent_classes]}),
                         hide_index=True)
        if r.xai is not None and r.xai.grouped_matrix.shape[0] > 1:
            st.subheader("Stabilità delle spiegazioni (descrittiva)")
            c1, c2 = st.columns(2)
            c1.altair_chart(jaccard_chart(r.xai.grouped_matrix, f"split per paziente — media {r.xai.grouped_mean:.2f}"),
                            width="stretch")
            c2.altair_chart(jaccard_chart(r.xai.random_matrix, f"split casuale — media {r.xai.random_mean:.2f}"),
                            width="stretch")
            st.caption(f"Jaccard fra i {r.xai.k} geni con coefficiente più grande della regressione logistica "
                       "nei diversi fold. Se togliere pazienti cambia i geni scelti più di quanto lo cambi "
                       "togliere cellule a caso, l'informazione è organizzata per paziente. Misura descrittiva, "
                       "senza test statistico.")
        if r.model_comparison is not None:
            mc = r.model_comparison
            st.subheader("Confronto fra modelli (un paziente alla volta fuori)")
            rows = [{"modello": mc.best_model + " (migliore)", "macro-F1": round(mc.scores[mc.best_model].mean, 3),
                     "p Nadeau-Bengio": "—", "p Wilcoxon": "—", "nota": ""}]
            rows += [{"modello": c.model, "macro-F1": round(c.mean_macro_f1, 3),
                      "p Nadeau-Bengio": f"{c.p_nadeau_bengio:.3f}", "p Wilcoxon": f"{c.p_wilcoxon:.3f}",
                      "nota": c.note} for c in mc.comparisons]
            st.dataframe(pd.DataFrame(rows), hide_index=True)
            st.caption("Regressione logistica sui 2.000 geni più variabili; random forest e gradient boosting "
                       "su 50 componenti SVD. Geni e componenti sono stimati sul solo training di ogni fold.")
            ref = mc.scores["logreg"]
            if ref.n_folds_with_absent_classes:
                with st.expander(f"Classi assenti dal paziente lasciato fuori: {ref.n_folds_with_absent_classes} "
                                 f"fold su {len(ref.fold_scores)} (macro-F1 calcolata sulle classi presenti)"):
                    st.dataframe(pd.DataFrame({"fold (paziente)": range(1, len(ref.absent_classes) + 1),
                                               "classi assenti": [", ".join(a) or "—" for a in ref.absent_classes]}),
                                 hide_index=True)


# --------------------------------------------------------------------------- #
# Tab 3: Modulo B
# --------------------------------------------------------------------------- #
with tab_b:
    st.markdown("Usa il repertorio **TCR** come identità indipendente dal trascrittoma: le cellule dello "
                "stesso clone T dovrebbero avere la stessa etichetta in ogni compartimento.")
    b_ctx = st.session_state.b_ctx if mode == "Demo immediata" else None
    if mode == "Carica studio":
        adata_b = load_matrix_ui("b")
        if adata_b is not None:
            cols = list(adata_b.obs.columns)
            c = st.columns(4)
            sel = dict(patient=c[0].selectbox("Paziente", cols, key="b_p"),
                       compartment=c[1].selectbox("Compartimento", cols, key="b_c"),
                       celltype=c[2].selectbox("Tipo cellulare", cols, key="b_t"),
                       barcode=c[3].selectbox("Barcode", cols, key="b_b"))
            man = st.text_input("Manifest VDJ (CSV con colonne path,patient,compartment)", key="b_man")
            markers = st.text_area("Mappa marcatori (JSON)", value='{"CD4T": ["CD4"], "CD8T": ["CD8A", "CD8B"]}', key="b_mk")
            comps = sorted(adata_b.obs[sel["compartment"]].astype(str).unique())
            ref = st.selectbox("Compartimento di riferimento (sangue)", comps, key="b_ref")
            if man and st.button("Esegui Modulo B", type="primary"):
                m = read_table(man)
                contigs = parse_vdj_contigs([(Path(x.path), str(x.patient), str(x.compartment)) for x in m.itertuples()])
                st.session_state.b_ctx = dict(adata=adata_b, contigs=contigs, markers=json.loads(markers),
                                              reference=ref, label=st.session_state.dataset_name, **sel)
                st.session_state.tcr = None
                st.session_state.run_b = True
                b_ctx = st.session_state.b_ctx
    if b_ctx is not None and st.session_state.tcr is None and (
            st.session_state.autorun or st.session_state.run_b):
        with st.spinner("Modulo B in esecuzione..."):
            try:
                st.session_state.tcr = run_tcr_validation(
                    b_ctx["adata"], b_ctx["contigs"], patient_col=b_ctx["patient"],
                    compartment_col=b_ctx["compartment"], celltype_col=b_ctx["celltype"],
                    barcode_col=b_ctx["barcode"], marker_map=b_ctx["markers"],
                    reference_compartment=b_ctx["reference"])
                st.session_state.tcr_table = tcr_by_celltype(
                    b_ctx["adata"], b_ctx["contigs"], b_ctx["patient"], b_ctx["compartment"],
                    b_ctx["celltype"], b_ctx["barcode"])
            except ValueError as e:
                st.error(str(e))
        st.session_state.run_b = False
    t = st.session_state.tcr
    if t is None:
        st.info("Carica la demo dalla barra laterale oppure scegli 'Carica studio'.")
    else:
        summary_line(tcr_summary(t))
        if b_ctx:
            st.caption(f"Dataset: {b_ctx.get('label', '')}.")
        if t.barcode_match is not None:
            st.caption(t.barcode_match.sentence)
        d = t.discordance
        m1, m2 = st.columns(2)
        m1.metric("Eccesso di discordanza fra compartimenti", f"{d.mean_excess:+.3f}")
        m2.metric("IC 95% (bootstrap sui pazienti)",
                  f"[{d.ci_low:+.3f}, {d.ci_high:+.3f}]" if d.sufficient else "non prodotto")
        st.markdown(t.narrative)
        if t.conventions:
            st.subheader("Tasso d'errore dell'annotazione per compartimento")
            for c in t.conventions.values():
                st.markdown(f"- **{c.name}**: {c.definition}")
            rows = []
            for c in t.conventions.values():
                for comp, b in c.by_compartment.items():
                    rows.append({"convenzione": c.name, "compartimento": comp, "tasso d'errore": f"{b.mean:.3f}",
                                 "IC 95%": f"[{b.ci_low:.3f}, {b.ci_high:.3f}]" if b.sufficient
                                 else f"non prodotto ({b.n_groups} pazienti)", "pazienti": b.n_groups})
            st.dataframe(pd.DataFrame(rows), hide_index=True)
        if st.session_state.tcr_table is not None:
            st.subheader("Cellule con TCR per etichetta (contaminazioni cross-lineage)")
            st.dataframe(st.session_state.tcr_table)
            st.caption("Tabella descrittiva. Un TCR in un'etichetta non-T (NK, mieloidi, stromali) può essere un "
                       "doppietto, RNA ambientale o un errore di annotazione: la tabella non li distingue e non "
                       "classifica le singole cellule.")
        if t.cell_flags is not None and b_ctx is not None:
            st.download_button("Scarica i flag per cellula (CSV)", t.cell_flags.to_csv(index_label="obs_name"),
                               file_name="audit_flags.csv", mime="text/csv")
            st.caption(f"Cellule con flag valutabile: {t.flag_coverage:.1%}. NA = non verificabile, non 'corretta'.")
            st.subheader("Frazione di CD8 validata dal repertorio")
            st.warning(EXPERIMENTAL_NOTE)
            comps = sorted(b_ctx["adata"].obs[b_ctx["compartment"]].astype(str).unique())
            comps = [c for c in comps if c != b_ctx["reference"]]
            cd8_comp = st.selectbox("Compartimento", comps, key="cd8_comp")
            if st.button("Calcola gli intervalli"):
                obs_flags = b_ctx["adata"].obs.join(t.cell_flags)
                labels = list(b_ctx["markers"])
                st.session_state.cd8 = cd8_fraction_intervals(
                    obs_flags, b_ctx["patient"], b_ctx["compartment"], b_ctx["celltype"],
                    "audit_reference_label", cd8_comp, cd4_label=labels[0], cd8_label=labels[-1])
            cd8 = st.session_state.cd8
            if cd8 is not None:
                if cd8.refused_reason:
                    st.error(f"Intervalli non prodotti: {cd8.refused_reason}.")

                def iv(x):
                    return f"{x.low:.2f}–{x.high:.2f}" if x.low is not None else f"non prodotto: {x.refused_reason}"
                st.dataframe(pd.DataFrame([{
                    "paziente": p.patient, "n CD4+CD8": p.n_cd4_called + p.n_cd8_called,
                    "riportata": "—" if p.reported is None else f"{p.reported:.2f}",
                    "errore 0.5x": iv(p.scenarios[0.5]), "errore 1x": iv(p.scenarios[1.0]),
                    "errore 2x": iv(p.scenarios[2.0])} for p in cd8.patients]), hide_index=True)
                st.caption(cd8.assumptions)


# --------------------------------------------------------------------------- #
# Sintesi dei controlli e report
# --------------------------------------------------------------------------- #
with tab_v:
    verdicts = summarize(st.session_state.design, st.session_state.leakage, st.session_state.tcr,
                         st.session_state.cd8)
    if not verdicts:
        st.info("Esegui almeno una sezione per ottenere la sintesi dei controlli.")
    else:
        st.markdown('<p class="small">Ogni riga descrive un controllo e dice che cosa i dati permettono di '
                    'stimare. Non è un giudizio sullo studio né sul lavoro di chi lo ha prodotto, e non '
                    'esiste un colore complessivo.</p>', unsafe_allow_html=True)
    for v in verdicts:
        rows = "".join(f"<tr><td>{html.escape(c.name)}</td><td>{badge(c.state, c.label)}</td>"
                       f"<td class='small'>{html.escape(c.text)}</td></tr>" for c in v.checks)
        st.markdown(f'<div class="card"><b>{html.escape(v.section)}</b>'
                    f"<table class='sem'><tr><th>Controllo</th><th>Stato</th><th>Descrizione</th></tr>{rows}</table>"
                    f'<span class="small">Regola: {html.escape(v.rule)}</span></div>', unsafe_allow_html=True)
    st.subheader("Limiti dichiarati")
    for x in standing_limits(st.session_state.cd8):
        st.markdown(f"- {x}")
    if verdicts:
        name = st.session_state.dataset_name or "studio"
        md = render_markdown_report(st.session_state.design, st.session_state.leakage, st.session_state.tcr,
                                    st.session_state.cd8, dataset_name=name)
        page = render_report(leakage_result=st.session_state.leakage, tcr_result=st.session_state.tcr,
                             dataset_name=name, design_result=st.session_state.design,
                             cd8_result=st.session_state.cd8)
        c1, c2 = st.columns(2)
        c1.download_button("Scarica il report (Markdown)", md, file_name=f"report_{name}.md", mime="text/markdown")
        c2.download_button("Scarica il report (HTML, stampabile in PDF dal browser)", page,
                           file_name=f"report_{name}.html", mime="text/html")

st.session_state.autorun = False
