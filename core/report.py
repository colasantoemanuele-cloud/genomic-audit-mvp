"""Report HTML autocontenuto (grafici inline come data URI, nessuna dipendenza esterna
al momento dell'apertura del file) e report Markdown, per un ricercatore clinico non tecnico:
ogni numero è accompagnato da una frase che ne spiega il significato pratico. La frazione di CD8
sperimentale non compare nei report (revisione, A7).
"""

from __future__ import annotations

import base64
import html
import io
from datetime import UTC, datetime
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from core.design_audit import MISSING_LEVEL_CLASS, DesignAuditResult
from core.leakage_audit import GAP_ALERT, LeakageAuditResult, ModelComparisonResult
from core.tcr_validation import BIOLOGICAL_ASSUMPTIONS, MarkerErrorResult, TcrValidationResult
from core.verdict import (
    GIALLO,
    GRIGIO,
    ROSSO,
    TITLE,
    USE_WARNING,
    SectionSummary,
    design_summary,
    leakage_summary,
    standing_limits,
    summarize,
    tcr_summary,
)

# Colore della classe di un confronto: lo stesso della riga corrispondente nella sintesi
COMPARISON_STATE = {"stimabile": "verde", "stimabile con bassa potenza": GIALLO, "non stimabile": ROSSO,
                    MISSING_LEVEL_CLASS: GRIGIO}

_CSS = """
body { font-family: -apple-system, "Segoe UI", Roboto, Helvetica, Arial, sans-serif;
       max-width: 900px; margin: 2rem auto; padding: 0 1.5rem; color: #1a1a1a; line-height: 1.5; }
h1 { font-size: 1.6rem; border-bottom: 3px solid #1b9e77; padding-bottom: .4rem; }
h2 { font-size: 1.25rem; margin-top: 2.5rem; color: #1b6e5c; }
h3 { font-size: 1.05rem; margin-top: 1.5rem; }
.card { background: #f7f9f9; border: 1px solid #dde5e3; border-radius: 8px; padding: 1rem 1.3rem; margin: 1rem 0; }
.verdict-yes { border-left: 5px solid #1b9e77; }
.verdict-no { border-left: 5px solid #999; }
.verdict-warn { border-left: 5px solid #d95f02; }
table { border-collapse: collapse; width: 100%; margin: .8rem 0; font-size: .92rem; }
th, td { text-align: left; padding: .35rem .6rem; border-bottom: 1px solid #e2e8e6; }
th { background: #eef3f2; }
.narrative { font-size: .97rem; color: #333; }
img { max-width: 100%; height: auto; display: block; margin: .8rem 0; }
.tag { display: inline-block; padding: .1rem .5rem; border-radius: 4px; font-size: .8rem;
       font-weight: 600; }
.tag-ok { background: #dff3ec; color: #0f6848; }
.tag-warn { background: #fdecdc; color: #a34c05; }
.tag-verde { background: #dff3ec; color: #0f6848; }
.tag-giallo { background: #fbf0cf; color: #7a5a00; }
.tag-rosso { background: #fbe0dd; color: #96261c; }
.tag-grigio { background: #e8ebee; color: #48525e; }
.footer { margin-top: 3rem; padding-top: 1rem; border-top: 1px solid #ddd; font-size: .85rem; color: #666; }
.warning { background: #fff4e5; border: 1px solid #f0c27b; border-radius: 8px; padding: .7rem 1rem; font-weight: 600; }
"""


def _fig_to_data_uri(fig) -> str:
    buf = io.BytesIO()
    fig.savefig(buf, format="png", dpi=140, bbox_inches="tight")
    plt.close(fig)
    return "data:image/png;base64," + base64.b64encode(buf.getvalue()).decode("ascii")


# --------------------------------------------------------------------------- #
# Audit del disegno
# --------------------------------------------------------------------------- #
def _card_class(summary: SectionSummary) -> str:
    """Colore del bordo della scheda: verde solo se TUTTI i controlli della sezione sono verdi;
    grigio se qualcuno non è stato eseguito o non è valutabile; arancione negli altri casi."""
    states = {c.state for c in summary.checks}
    if states & {ROSSO, GIALLO}:
        return "verdict-warn"
    return "verdict-no" if GRIGIO in states else "verdict-yes"


def _limits_section() -> str:
    items = "".join(f"<li>{html.escape(x)}</li>" for x in standing_limits())
    return f'<h2>Limiti dichiarati</h2><div class="card verdict-no"><ul class="narrative">{items}</ul></div>'


def _summary_section(sections: list[SectionSummary]) -> str:
    rows = "".join(
        f"<tr><td>{html.escape(s.section)}</td><td>{html.escape(c.name)}</td>"
        f'<td><span class="tag tag-{c.state}">{html.escape(c.label)}</span></td>'
        f"<td>{html.escape(c.text)}</td></tr>" for s in sections for c in s.checks)
    rules = "".join(f"<li><b>{html.escape(s.section)}</b>: {html.escape(s.rule)}</li>" for s in sections)
    return (f"<h2>{TITLE}</h2>"
            '<div class="card verdict-no"><p class="narrative">Ogni riga descrive un controllo e dice che cosa '
            "i dati permettono di stimare. Non è un giudizio sullo studio né sul lavoro di chi lo ha "
            "prodotto, e non esiste un colore complessivo.</p>"
            f"<table><tr><th>Sezione</th><th>Controllo</th><th>Stato</th><th>Descrizione</th></tr>{rows}</table>"
            f'<p class="narrative">Regole:</p><ul class="narrative">{rules}</ul></div>')


_KIND_LABEL = {
    "annidamento": "annidamento", "uno-a-uno": "fattori coincidenti",
    "esito-determinato": "esito determinato da un fattore", "unità-tecnica": "unità tecnica",
    "costante": "fattore costante", "identificatore": "identificativo di campione",
}


def _design_section(result: DesignAuditResult) -> str:
    structural = [f for f in result.findings
                  if f.kind in ("annidamento", "uno-a-uno", "esito-determinato", "unità-tecnica")]
    verdict_cls = "verdict-warn" if structural else _card_class(design_summary(result))
    roles = ", ".join(f"{html.escape(c)} ({html.escape(r)})" for c, r in result.roles.items())
    parts = [
        "<h2>Audit del disegno e del confondimento</h2>",
        f'<div class="card {verdict_cls}">',
        f'<p class="narrative">Unità analizzate: {result.n_rows} righe dei metadati. Fattori: {roles}. '
        "Questo controllo guarda solo la struttura del disegno (chi è stato misurato come, "
        "quando, in quale tessuto): non usa l'espressione genica e non giudica il lavoro di "
        "chi ha disegnato lo studio. Rende espliciti i limiti che il disegno pone alle "
        "conclusioni.</p>",
    ]
    if result.comparisons:
        parts.append("<h3>Confronti richiesti</h3>")
        parts.append("<table><tr><th>Confronto</th><th>Classe</th><th>Disegno</th>"
                     "<th>Unità indipendenti</th><th>Cosa significa</th></tr>")
        for c in result.comparisons:
            tag = f"tag-{COMPARISON_STATE.get(c.classification, GRIGIO)}"
            parts.append(
                f"<tr><td>{html.escape(c.factor)}: {html.escape(c.level_a)} vs {html.escape(c.level_b)}</td>"
                f'<td><span class="tag {tag}">{html.escape(c.classification)}</span></td>'
                f"<td>{html.escape(c.design)}</td><td>{c.n_units}</td>"
                f"<td>{html.escape(c.sentence)}</td></tr>")
        parts.append("</table>")
    parts.append("<h3>Fatti strutturali del disegno</h3>")
    parts.append('<p class="narrative">Controlli deterministici sulla tabella dei metadati, utili come '
                 "lista di verifica: non sono stime. L'unico elemento inferenziale di questa sezione è "
                 "il Cramér V.</p>")
    if structural:
        parts.append("<ul>" + "".join(
            f"<li><b>{html.escape(_KIND_LABEL[f.kind])}</b>: {html.escape(f.sentence)}</li>"
            for f in structural) + "</ul>")
    else:
        parts.append('<p class="narrative">Nessun annidamento, coincidenza o esito determinato '
                     "da un singolo fattore.</p>")
    other = [f for f in result.findings if f.kind in ("costante", "identificatore")]
    if other:
        parts.append('<p class="narrative">' + " ".join(html.escape(f.sentence) for f in other) + "</p>")
    if result.pairs:
        parts.append("<h3>Associazione fra coppie di fattori (Cramér V corretto di Bergsma)</h3>")
        parts.append("<table><tr><th>Coppia</th><th>V</th><th>Cosa significa</th></tr>")
        for pr in result.pairs:
            v = "non valutabile" if pr.cramer_v is None else f"{pr.cramer_v:.2f}"
            parts.append(f"<tr><td>{html.escape(pr.factor_a)} × {html.escape(pr.factor_b)}</td>"
                         f"<td>{v}</td><td>{html.escape(pr.sentence)}</td></tr>")
        parts.append("</table>")
    for note in result.notes:
        parts.append(f'<p class="narrative"><i>{html.escape(note)}</i></p>')
    parts.append("</div>")
    return "".join(parts)


# --------------------------------------------------------------------------- #
# Modulo A
# --------------------------------------------------------------------------- #
def _leakage_chart(result: LeakageAuditResult) -> str:
    fig, ax = plt.subplots(figsize=(5, 3.5))
    labels = ["Split per paziente", "Split casuale\n(controllo negativo)"]
    means = [result.grouped.mean, result.random.mean]
    stds = [result.grouped.std, result.random.std]
    ax.bar(labels, means, yerr=stds, capsize=6, color=["#1b9e77", "#d95f02"], alpha=0.85)
    ax.set_ylabel("macro-F1")
    ax.set_ylim(0, min(1.05, max(means) + max(stds) + 0.15))
    ax.set_title(f"Modulo A — {html.escape(result.task)}")
    fig.tight_layout()
    return _fig_to_data_uri(fig)


def _model_comparison_chart(mc: ModelComparisonResult) -> str:
    names = list(mc.scores)
    means = np.array([mc.scores[n].mean for n in names])
    stds = np.array([mc.scores[n].std for n in names])
    order = np.argsort(means)
    fig, ax = plt.subplots(figsize=(5.5, 0.6 * len(names) + 1.3))
    ax.barh(np.array(names)[order], means[order], xerr=stds[order], capsize=4,
            color="#7570b3", alpha=0.85)
    ax.set_xlabel("macro-F1 (LeaveOneGroupOut, un paziente alla volta)")
    fig.tight_layout()
    return _fig_to_data_uri(fig)


def _leakage_section(result: LeakageAuditResult) -> str:
    verdict_cls = _card_class(leakage_summary(result))
    tag = (f'<span class="tag tag-warn">divario oltre la soglia di ±{GAP_ALERT:g}</span>' if abs(result.gap) > GAP_ALERT
           else f'<span class="tag tag-ok">divario entro la soglia di ±{GAP_ALERT:g}</span>')
    html_parts = [
        "<h2>Modulo A — Audit del leakage per paziente</h2>",
        f'<div class="card {verdict_cls}">',
        f"<p>{tag}</p>",
        f'<p class="narrative">{html.escape(result.narrative)}</p>',
        "<table><tr><th>Quantità</th><th>Valore</th><th>Cosa significa</th></tr>",
        f"<tr><td>Cellule totali</td><td>{result.n_cells:,}</td>"
        f"<td>Numero di cellule usate per questo task, dopo aver escluso quelle senza etichetta.</td></tr>",
        f"<tr><td>Pazienti</td><td>{result.n_patients}</td>"
        f"<td>Unità statistica indipendente: con più pazienti la stima è meno variabile.</td></tr>",
        f"<tr><td>Classi</td><td>{result.n_classes}</td><td>Categorie del task di classificazione.</td></tr>",
        f"<tr><td>Divario sulla media</td><td>{result.gap:+.3f}</td>"
        f"<td>macro-F1 dello split casuale (non valido) meno quella dello split per paziente. Entro "
        f"±{GAP_ALERT:g} è considerato rumore e non se ne indica la direzione.</td></tr>",
        f"<tr><td>Rapporto delle deviazioni standard fra fold</td><td>{result.std_ratio:.2f}</td>"
        f"<td>Descrittivo: con pochi fold correlati la deviazione standard fra fold non è una stima "
        f"valida dell'incertezza sulla prestazione.</td></tr>",
        "</table>",
        f'<img src="{_leakage_chart(result)}" alt="Confronto macro-F1 fra split">',
        "</div>",
    ]
    if result.model_comparison is not None:
        mc = result.model_comparison
        html_parts.append("<h3>Confronto fra modelli (LeaveOneGroupOut)</h3>")
        threshold = (getattr(result, "settings", {}) or {}).get("min_patients_for_model_comparison", 8)
        html_parts.append(
            f'<p class="narrative">Con almeno {threshold} pazienti (questo dataset ne ha '
            f"{result.n_patients}), il confronto fra modelli lascia fuori un paziente alla volta. Ogni "
            f"modello è confrontato con <b>{html.escape(mc.reference_model)}</b>, il riferimento fissato "
            f"prima di vedere i dati. Il confronto è <b>descrittivo</b>: nessun p-value, perché i test "
            f"disponibili (Nadeau-Bengio, Wilcoxon) non sono calibrati con un paziente per fold e 8-20 "
            f"pazienti (vedi docs/METODI.md). Macro-F1 media più alta: "
            f"<b>{html.escape(mc.best_model)}</b> ({mc.scores[mc.best_model].mean:.3f}).</p>"
        )
        html_parts.append(f'<img src="{_model_comparison_chart(mc)}" alt="Confronto fra modelli">')
        ref = mc.scores[mc.reference_model]
        html_parts.append(
            f"<table><tr><th>Modello</th><th>macro-F1</th><th>Δ vs {html.escape(mc.reference_model)}</th>"
            "<th>Pazienti in cui fa meglio / peggio / uguale</th></tr>"
            f"<tr><td>{html.escape(mc.reference_model)} (riferimento)</td><td>{ref.mean:.3f}</td>"
            "<td>—</td><td>—</td></tr>"
        )
        for comp in mc.comparisons:
            html_parts.append(
                f"<tr><td>{html.escape(comp.model)}</td><td>{comp.mean_macro_f1:.3f}</td>"
                f"<td>{comp.delta_vs_reference:+.3f}</td><td>{comp.n_better} / {comp.n_worse} / "
                f"{comp.n_ties}</td></tr>"
            )
        html_parts.append("</table>")
    else:
        settings = getattr(result, "settings", {}) or {}
        if settings.get("min_patients_for_model_comparison", 8) <= result.n_patients:
            why = "escluso su richiesta (opzione --rapido)."
        else:
            why = (f"servono almeno {settings.get('min_patients_for_model_comparison', 8)} pazienti, "
                   f"questo dataset ne ha {result.n_patients}.")
        html_parts.append(f'<p class="narrative">Confronto fra modelli non eseguito: {why}</p>')
    return "".join(html_parts)


# --------------------------------------------------------------------------- #
# Modulo B
# --------------------------------------------------------------------------- #
def _discordance_chart(result: TcrValidationResult) -> str | None:
    d = result.discordance
    if d.by_compartment_pair.empty:
        return None
    rows = d.by_compartment_pair
    labels = [f"{a} ↔ {b}" for a, b in zip(rows.comp_a, rows.comp_b)]
    excess = rows.excess.values
    lo, hi = rows.ci_low.values, rows.ci_high.values
    y = np.arange(len(labels))
    fig, ax = plt.subplots(figsize=(6.2, 0.55 * len(labels) + 1.3))
    xerr = np.vstack([np.nan_to_num(excess - lo, nan=0), np.nan_to_num(hi - excess, nan=0)])
    ax.errorbar(excess, y, xerr=xerr, fmt="o", capsize=4, color="#d95f02")
    ax.axvline(0, ls="--", c="k", lw=0.8)
    ax.set_yticks(y, labels)
    ax.set_xlabel("eccesso di discordanza (IC95%, intervallo t cluster-robusto sui pazienti)")
    fig.tight_layout()
    return _fig_to_data_uri(fig)


def _marker_error_chart(m: MarkerErrorResult) -> str | None:
    comps = [c for c, r in m.by_compartment.items() if r.sufficient]
    if not comps:
        return None
    means = [m.by_compartment[c].mean for c in comps]
    lo = [m.by_compartment[c].ci_low for c in comps]
    hi = [m.by_compartment[c].ci_high for c in comps]
    yerr = [[max(0, mm - ll) for mm, ll in zip(means, lo)], [max(0, hh - mm) for mm, hh in zip(means, hi)]]
    fig, ax = plt.subplots(figsize=(5, 3.5))
    ax.bar(comps, means, yerr=yerr, capsize=5, color="#d95f02", alpha=0.85)
    ax.set_ylabel("tasso d'errore dell'etichetta")
    ax.set_title(f"riferimento: {html.escape(m.reference_compartment)}")
    fig.tight_layout()
    return _fig_to_data_uri(fig)


def _tcr_section(result: TcrValidationResult) -> str:
    d = result.discordance
    if not d.sufficient:
        tag = '<span class="tag tag-warn">numerosità insufficiente</span>'
    elif d.ci_low > 0:
        tag = '<span class="tag tag-warn">eccesso di discordanza sopra lo zero</span>'
    elif d.ci_high < 0:
        tag = '<span class="tag tag-warn">intervallo interamente sotto lo zero</span>'
    else:
        tag = '<span class="tag tag-ok">nessun eccesso rilevato</span>'
    # bordo verde solo se ogni controllo della sezione è verde
    verdict_cls = _card_class(tcr_summary(result))

    html_parts = [
        "<h2>Modulo B — Verifica dell'annotazione via TCR</h2>",
        f'<div class="card {verdict_cls}">',
        f"<p>{tag}</p>",
        f'<p class="narrative"><i>{html.escape(BIOLOGICAL_ASSUMPTIONS)}</i></p>',
        "<table><tr><th>Quantità</th><th>Valore</th><th>Cosa significa</th></tr>",
        (f"<tr><td>Match dei barcode VDJ</td><td>{result.barcode_match.fraction:.1%}</td>"
         f"<td>{html.escape(result.barcode_match.sentence)}</td></tr>" if result.barcode_match else ""),
        f"<tr><td>Cellule con TCR</td><td>{result.n_cells_with_tcr:,}</td>"
        f"<td>Cellule per cui è stato ricostruito un clonotipo (catena TRB rilevata).</td></tr>",
        f"<tr><td>Coppie clone-compartimenti confrontabili</td><td>{d.n_pairs:,}</td>"
        f"<td>Ogni riga confronta lo stesso clone T in due compartimenti diversi.</td></tr>",
        f"<tr><td>Pazienti che contribuiscono</td><td>{d.n_patients}</td>"
        f"<td>L'unità statistica indipendente è il paziente, non il clone: i cloni di uno stesso "
        f"paziente non sono osservazioni indipendenti.</td></tr>",
    ]
    if d.sufficient:
        html_parts.append(
            f"<tr><td>Eccesso di discordanza</td><td>{d.mean_excess:+.3f} "
            f"[{d.ci_low:+.3f}, {d.ci_high:+.3f}]</td>"
            f"<td>Quanto più spesso le cellule dello stesso clone T ricevono etichette diverse fra "
            f"compartimenti rispetto al rumore entro un compartimento. Zero nell'intervallo = nessun "
            f"eccesso rilevato. Con compartimenti di dimensioni molto diverse può essere negativo anche "
            f"in presenza di errori: vedi la riga sotto.</td></tr>"
        )
        e = d.equal_weights
        if e is not None:
            val = (f"{e.mean:+.3f} [{e.ci_low:+.3f}, {e.ci_high:+.3f}]" if e.sufficient
                   else f"{e.mean:+.3f} (intervallo non prodotto)")
            html_parts.append(
                f"<tr><td>Sensibilità a pesi uguali</td><td>{val}</td><td>Stessa misura con i due "
                f"compartimenti pesati allo stesso modo: il valore atteso è ½‖pA − pB‖², zero solo se le "
                f"proporzioni di etichette coincidono e mai negativo.</td></tr>")
    html_parts.append("</table>")
    html_parts.append(f'<p class="narrative">{html.escape(result.narrative)}</p>')

    chart = _discordance_chart(result)
    if chart:
        html_parts.append(f'<img src="{chart}" alt="Eccesso di discordanza per coppia di compartimenti">')
        html_parts.append(
            "<table><tr><th>Coppia di compartimenti</th><th>Coppie clone-comp.</th><th>Pazienti</th>"
            "<th>Eccesso [IC95%]</th></tr>"
        )
        for r in d.by_compartment_pair.itertuples():
            if r.sufficient:
                val = f"{r.excess:+.3f} [{r.ci_low:+.3f}, {r.ci_high:+.3f}]"
            else:
                val = "numerosità insufficiente"
            html_parts.append(
                f"<tr><td>{html.escape(r.comp_a)} ↔ {html.escape(r.comp_b)}</td>"
                f"<td>{r.n_pairs}</td><td>{r.n_patients}</td><td>{val}</td></tr>"
            )
        html_parts.append("</table>")

    if result.conventions:
        html_parts.append("<h3>Tasso d'errore per compartimento (marcatori canonici)</h3>")
        ref = result.marker_error.reference_compartment if result.marker_error else ""
        html_parts.append(
            f'<p class="narrative">Identità di riferimento del clone stimata SOLO dal compartimento '
            f"'{html.escape(ref)}', confrontata con l'etichetta assegnata negli altri compartimenti. "
            f"Lo stesso confronto può essere contato in modi diversi: ogni numero è riportato con la "
            f"sua convenzione, e le convenzioni sono affiancate.</p><ul>"
            + "".join(f"<li><b>{html.escape(c.name)}</b>: {html.escape(c.definition)}</li>"
                      for c in result.conventions.values()) + "</ul>")
        if result.marker_error is not None:
            chart2 = _marker_error_chart(result.marker_error)
            if chart2 and "cell" in result.conventions:
                html_parts.append(f'<img src="{chart2}" alt="Tasso d\'errore per compartimento (per cellula)">')
        names = list(result.conventions)
        comps = sorted({c for conv in result.conventions.values() for c in conv.by_compartment})

        def fmt(res, ties=None):
            if res is None:
                return "—"
            units = f"{res.n_obs} unità" + (f"; unità in pareggio: {ties}" if ties is not None else "")
            if not res.sufficient:
                return f"{res.mean:.3f} (IC non prodotto: {html.escape(res.reason)}; {units})"
            return f"{res.mean:.3f} [{res.ci_low:.3f}, {res.ci_high:.3f}], {res.n_groups} pazienti, {units}"
        html_parts.append("<table><tr><th>Compartimento</th>" + "".join(
            f"<th>convenzione {html.escape(n)}</th>" for n in names) + "</tr>")
        for comp in comps:
            html_parts.append(f"<tr><td>{html.escape(comp)}</td>" + "".join(
                f"<td>{fmt(result.conventions[n].by_compartment.get(comp), (result.conventions[n].n_ties or {}).get(comp))}</td>"
                for n in names) + "</tr>")
        pairs = sorted({p for conv in result.conventions.values() for p in conv.paired_differences})
        for a, b in pairs:
            cells_txt = []
            for n in names:
                v = result.conventions[n].paired_differences.get((a, b))
                cells_txt.append("—" if v is None else
                                 fmt(v[0]).replace(f"{v[0].mean:.3f}", f"{v[0].mean:+.3f}", 1)
                                 .replace(f"{v[0].n_obs} unità", f"{v[1]} cloni"))
            html_parts.append(f"<tr><td>differenza {html.escape(a)} − {html.escape(b)} (stessi cloni)</td>"
                              + "".join(f"<td>{t}</td>" for t in cells_txt) + "</tr>")
        html_parts.append("</table>")

    if result.cell_flags is not None:
        f = result.cell_flags
        n = len(f)
        n_ref = int(f["audit_reference_label"].notna().sum())
        n_eval = int(f["audit_label_vs_reference"].notna().sum())
        n_disc = int(f["audit_label_vs_reference"].fillna(False).astype(bool).sum())
        html_parts.append("<h3>Flag per cellula</h3>")
        html_parts.append(
            f'<p class="narrative">Il controllo produce due colonne per cellula, esportabili in una '
            f"copia dell'AnnData (<code>*_audited.h5ad</code>) e in un CSV: "
            f"<code>audit_reference_label</code> (identità del clone stimata nel sangue) e "
            f"<code>audit_label_vs_reference</code> (l'etichetta assegnata discorda da quell'identità). "
            f"Le etichette originali non vengono modificate. Copertura: {n_ref:,} cellule su {n:,} "
            f"({n_ref / n:.1%}) hanno un'identità di riferimento; {n_eval:,} ({result.flag_coverage:.1%}) "
            f"sono valutabili, e di queste {n_disc:,} sono discordanti. Tutte le altre sono NA: "
            f"lo strumento non poteva verificarle (nessun TCR, clone senza cellule sufficienti nel "
            f"riferimento, etichetta fuori dalla marker_map, o cellula del compartimento di "
            f"riferimento). NA non "
            f"significa \"corretta\".</p>")

    html_parts.append("</div>")
    return "".join(html_parts)


# --------------------------------------------------------------------------- #
def render_report(
    leakage_result: LeakageAuditResult | None = None,
    tcr_result: TcrValidationResult | None = None,
    dataset_name: str = "dataset caricato",
    design_result: DesignAuditResult | None = None,
) -> str:
    """Ritorna una stringa HTML autocontenuta: nessun asset esterno, i grafici sono
    immagini PNG incorporate come data URI. Contiene la sintesi dei controlli, le sezioni
    eseguite e i limiti dichiarati."""
    now = datetime.now(UTC).strftime("%Y-%m-%d %H:%M UTC")
    body = [
        f"<h1>Report di audit — {html.escape(dataset_name)}</h1>",
        f'<p class="warning">{html.escape(USE_WARNING)}</p>',
        f'<p class="footer" style="margin-top:0;border-top:none;padding-top:0;">'
        f"Generato il {now}. Strumento diagnostico: non è una certificazione, è un supporto per chi "
        f"analizza i dati.</p>",
    ]
    sections = summarize(design_result, leakage_result, tcr_result)
    if sections:
        body.append(_summary_section(sections))
    if design_result is not None:
        body.append(_design_section(design_result))
    if leakage_result is not None:
        body.append(_leakage_section(leakage_result))
    if tcr_result is not None:
        body.append(_tcr_section(tcr_result))
    if leakage_result is None and tcr_result is None and design_result is None:
        body.append("<p>Nessun risultato da mostrare: esegui almeno una sezione.</p>")
    body.append(_limits_section())
    body.append(
        '<div class="footer">Report generato da genomic-audit. Le sezioni possibili sono quattro: la '
        "sintesi dei controlli, l'audit del disegno (legge solo i metadati e segnala quali confronti la "
        "struttura dello studio permette), il Modulo A (quanto una valutazione che non raggruppa per "
        "paziente cambia la macro-F1 rispetto a una che lo fa) e il Modulo B (quanto l'annotazione di "
        "tipo cellulare è coerente per uno stesso clone T attraverso i compartimenti, con le ipotesi "
        "biologiche dichiarate). Nessuna sezione sostituisce una verifica biologica indipendente.</div>"
    )
    return f"<!doctype html><html lang=\"it\"><head><meta charset=\"utf-8\">" \
           f"<title>Report di audit — {html.escape(dataset_name)}</title>" \
           f"<style>{_CSS}</style></head><body>{''.join(body)}</body></html>"


def save_report(
    out_path: str | Path,
    leakage_result: LeakageAuditResult | None = None,
    tcr_result: TcrValidationResult | None = None,
    dataset_name: str = "dataset caricato",
    design_result: DesignAuditResult | None = None,
) -> Path:
    """Scrive il report HTML in ``out_path`` (crea la cartella se serve) e ne ritorna il percorso."""
    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(render_report(leakage_result, tcr_result, dataset_name, design_result),
                        encoding="utf-8")
    return out_path


# --------------------------------------------------------------------------- #
# Report Markdown (export della web app e della CLI)
# --------------------------------------------------------------------------- #
def render_markdown_report(
    design_result: DesignAuditResult | None = None,
    leakage_result: LeakageAuditResult | None = None,
    tcr_result: TcrValidationResult | None = None,
    dataset_name: str = "dataset",
) -> str:
    """Report in Markdown: sintesi dei controlli con le regole che la producono, numeri con la
    loro definizione e limiti dichiarati. Stesse frasi del report HTML."""

    now = datetime.now(UTC).strftime("%Y-%m-%d %H:%M UTC")
    md = [f"# Report di audit — {dataset_name}", "", f"**{USE_WARNING}**", "",
          f"Generato il {now}. Strumento diagnostico locale: non è una certificazione.", ""]
    sections = summarize(design_result, leakage_result, tcr_result)
    if sections:
        md += [f"## {TITLE}", "",
               "Ogni riga descrive un controllo e dice che cosa i dati permettono di stimare. Non è un "
               "giudizio sullo studio, e non esiste un colore complessivo. Il verde significa che nessuna "
               "regola segnala un limite, non che la stima sia corretta.", "",
               "| Sezione | Controllo | Stato | Descrizione |", "|---|---|---|---|"]
        md += [f"| {s.section} | {c.name} | {c.state}: {c.label} | {c.text} |" for s in sections for c in s.checks]
        md += [""] + [f"- *Regola {s.section}*: {s.rule}" for s in sections] + [""]
    if design_result is not None:
        md += ["## Audit del disegno", "", design_result.narrative, ""]
        if design_result.comparisons:
            md += ["| Confronto | Classe | Unità | p-value minimo |", "|---|---|---|---|"]
            for c in design_result.comparisons:
                mp = "—" if c.min_pvalue is None else f"{c.min_pvalue:.3f}"
                md.append(f"| {c.factor}: {c.level_a} vs {c.level_b} | {c.classification} | {c.n_units} | {mp} |")
            md.append("")
        md += [f"> {n}" for n in design_result.notes] + [""]
    if leakage_result is not None:
        r = leakage_result
        md += ["## Modulo A — Leakage per paziente", "", r.narrative, "",
               "| Schema | macro-F1 media | dev. std | fold |", "|---|---|---|---|",
               f"| split per paziente | {r.grouped.mean:.3f} | {r.grouped.std:.3f} | {len(r.grouped.fold_scores)} |",
               f"| split casuale (controllo negativo) | {r.random.mean:.3f} | {r.random.std:.3f} | {len(r.random.fold_scores)} |",
               "", "Definizione: macro-F1 calcolata sulle sole classi presenti nel fold di test; "
               "le classi assenti sono dichiarate. Nessun intervallo di confidenza: i punteggi dei "
               "fold sono correlati e un intervallo non sarebbe calibrato.", ""]
        if r.xai is not None:
            md += [f"Stabilità delle spiegazioni (descrittiva): Jaccard medio fra i {r.xai.k} geni "
                   f"principali dei fold = {r.xai.grouped_mean:.3f} con split per paziente, "
                   f"{r.xai.random_mean:.3f} con split casuale.", ""]
        if r.model_comparison is not None:
            mc = r.model_comparison
            md += [f"Ogni modello contro {mc.reference_model} (riferimento fissato prima di vedere i dati), "
                   "LeaveOneGroupOut. Confronto descrittivo: nessun p-value, perché Nadeau-Bengio e Wilcoxon "
                   "non sono calibrati in questo regime.", "",
                   f"| Modello | macro-F1 | Δ vs {mc.reference_model} | pazienti meglio / peggio / uguale |",
                   "|---|---|---|---|"]
            md += [f"| {mc.reference_model} (riferimento) | {mc.scores[mc.reference_model].mean:.3f} | — | — |"]
            md += [f"| {c.model} | {c.mean_macro_f1:.3f} | {c.delta_vs_reference:+.3f} | "
                   f"{c.n_better} / {c.n_worse} / {c.n_ties} |" for c in mc.comparisons]
            md.append("")
            ref = mc.scores[mc.reference_model]
            if ref.n_folds_with_absent_classes:
                md += [f"Classi assenti dal paziente lasciato fuori in {ref.n_folds_with_absent_classes} fold su "
                       f"{len(ref.fold_scores)} (macro-F1 calcolata sulle classi presenti): "
                       + "; ".join(f"fold {i + 1}: {', '.join(a)}" for i, a in enumerate(ref.absent_classes) if a),
                       ""]
    if tcr_result is not None:
        md += ["## Modulo B — Verifica dell'annotazione via TCR", "", f"*{BIOLOGICAL_ASSUMPTIONS}*", ""]
        if tcr_result.barcode_match is not None:
            md += [tcr_result.barcode_match.sentence, ""]
        md += [tcr_result.narrative, ""]
    md += ["## Limiti dichiarati", ""] + [f"- {x}" for x in standing_limits()] + [""]
    return "\n".join(md)
