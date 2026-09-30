"""Report HTML autocontenuto (grafici inline come data URI, nessuna dipendenza esterna
al momento dell'apertura del file) per un ricercatore clinico non tecnico: ogni numero e'
accompagnato da una frase che ne spiega il significato pratico.
"""

from __future__ import annotations

import base64
import html
import io
from datetime import datetime, timezone
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from core.design_audit import DesignAuditResult
from core.leakage_audit import LeakageAuditResult, ModelComparisonResult
from core.tcr_validation import MarkerErrorResult, TcrValidationResult

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
.footer { margin-top: 3rem; padding-top: 1rem; border-top: 1px solid #ddd; font-size: .85rem; color: #666; }
"""


def _fig_to_data_uri(fig) -> str:
    buf = io.BytesIO()
    fig.savefig(buf, format="png", dpi=140, bbox_inches="tight")
    plt.close(fig)
    return "data:image/png;base64," + base64.b64encode(buf.getvalue()).decode("ascii")


# --------------------------------------------------------------------------- #
# Audit del disegno
# --------------------------------------------------------------------------- #
_KIND_LABEL = {
    "annidamento": "annidamento", "uno-a-uno": "fattori coincidenti",
    "esito-determinato": "esito determinato da un fattore", "unita'-tecnica": "unita' tecnica",
    "costante": "fattore costante", "identificatore": "identificativo di campione",
}


def _design_section(result: DesignAuditResult) -> str:
    structural = [f for f in result.findings
                  if f.kind in ("annidamento", "uno-a-uno", "esito-determinato", "unita'-tecnica")]
    not_estimable = any(c.classification == "non stimabile" for c in result.comparisons)
    verdict_cls = "verdict-warn" if (structural or not_estimable) else "verdict-yes"
    roles = ", ".join(f"{html.escape(c)} ({html.escape(r)})" for c, r in result.roles.items())
    parts = [
        "<h2>Audit del disegno e del confondimento</h2>",
        f'<div class="card {verdict_cls}">',
        f'<p class="narrative">Unita\' analizzate: {result.n_rows} righe dei metadati. Fattori: {roles}. '
        "Questo controllo guarda solo la struttura del disegno (chi e' stato misurato come, "
        "quando, in quale tessuto): non usa l'espressione genica e non giudica il lavoro di "
        "chi ha disegnato lo studio. Rende espliciti i limiti che il disegno pone alle "
        "conclusioni.</p>",
    ]
    if result.comparisons:
        parts.append("<h3>Confronti richiesti</h3>")
        parts.append("<table><tr><th>Confronto</th><th>Classe</th><th>Disegno</th>"
                     "<th>Unita' indipendenti</th><th>Cosa significa</th></tr>")
        for c in result.comparisons:
            tag = ("tag-ok" if c.classification == "stimabile" else "tag-warn")
            parts.append(
                f"<tr><td>{html.escape(c.factor)}: {html.escape(c.level_a)} vs {html.escape(c.level_b)}</td>"
                f'<td><span class="tag {tag}">{html.escape(c.classification)}</span></td>'
                f"<td>{html.escape(c.design)}</td><td>{c.n_units}</td>"
                f"<td>{html.escape(c.sentence)}</td></tr>")
        parts.append("</table>")
    parts.append("<h3>Fatti strutturali del disegno</h3>")
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
    labels = ["Split per paziente\n(onesto)", "Split casuale\n(controllo negativo)"]
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
    verdict_cls = "verdict-warn" if result.gap > 0.05 else "verdict-yes"
    tag = ('<span class="tag tag-warn">leakage rilevabile</span>' if result.gap > 0.05
           else '<span class="tag tag-ok">nessun leakage evidente</span>')
    html_parts = [
        "<h2>Modulo A — Audit del leakage per paziente</h2>",
        f'<div class="card {verdict_cls}">',
        f"<p>{tag}</p>",
        f'<p class="narrative">{html.escape(result.narrative)}</p>',
        "<table><tr><th>Quantita'</th><th>Valore</th><th>Cosa significa</th></tr>",
        f"<tr><td>Cellule totali</td><td>{result.n_cells:,}</td>"
        f"<td>Numero di cellule usate per questo task, dopo aver escluso quelle senza etichetta.</td></tr>",
        f"<tr><td>Pazienti</td><td>{result.n_patients}</td>"
        f"<td>Unita' statistica indipendente: piu' pazienti = stima piu' affidabile.</td></tr>",
        f"<tr><td>Classi</td><td>{result.n_classes}</td><td>Categorie del task di classificazione.</td></tr>",
        f"<tr><td>Divario sulla media</td><td>{result.gap:+.3f}</td>"
        f"<td>Quanto lo split casuale (non valido) sovrastima l'accuratezza rispetto allo split "
        f"per paziente (onesto). Piu' vicino a zero, meglio e'.</td></tr>",
        f"<tr><td>Rapporto delle deviazioni standard</td><td>{result.std_ratio:.2f}</td>"
        f"<td>Quanto lo split casuale sottostima l'incertezza sulla performance "
        f"(valori &gt;1 indicano sottostima).</td></tr>",
        "</table>",
        f'<img src="{_leakage_chart(result)}" alt="Confronto macro-F1 fra split">',
        "</div>",
    ]
    if result.model_comparison is not None:
        mc = result.model_comparison
        html_parts.append("<h3>Confronto fra modelli (LeaveOneGroupOut)</h3>")
        html_parts.append(
            f'<p class="narrative">Con almeno {result.n_patients} pazienti, il confronto fra modelli '
            f"e' stato ripetuto lasciando fuori un paziente alla volta (la valutazione piu' onesta "
            f"possibile su una coorte piccola). Il modello con macro-F1 media piu' alta e' "
            f"<b>{html.escape(mc.best_model)}</b> ({mc.scores[mc.best_model].mean:.3f}).</p>"
        )
        html_parts.append(f'<img src="{_model_comparison_chart(mc)}" alt="Confronto fra modelli">')
        html_parts.append(
            "<table><tr><th>Modello</th><th>macro-F1</th><th>Δ vs migliore</th>"
            "<th>p (Nadeau-Bengio)</th><th>p (Wilcoxon)</th><th>Nota</th></tr>"
        )
        for comp in mc.comparisons:
            note = html.escape(comp.note) if comp.note else "—"
            html_parts.append(
                f"<tr><td>{html.escape(comp.model)}</td><td>{comp.mean_macro_f1:.3f}</td>"
                f"<td>{comp.delta_vs_best:+.3f}</td><td>{comp.p_nadeau_bengio:.3f}</td>"
                f"<td>{comp.p_wilcoxon:.3f}</td><td>{note}</td></tr>"
            )
        html_parts.append("</table>")
        html_parts.append(
            '<p class="narrative">Con pochi fold la potenza statistica e\' bassa: "non significativo" '
            "non vuol dire \"equivalente\". Se un modello piu' semplice non viene battuto in modo "
            "significativo, e' un'informazione utile, non un errore dello strumento.</p>"
        )
    else:
        html_parts.append(
            f'<p class="narrative">Confronto multi-modello non eseguito: servono almeno '
            f"{result.n_patients} pazienti in piu' per una valutazione LeaveOneGroupOut robusta "
            f"(soglia di default: 8 pazienti).</p>"
        )
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
    ax.set_xlabel("eccesso di discordanza (IC95%, cluster bootstrap sui pazienti)")
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
        verdict_cls, tag = "verdict-no", '<span class="tag tag-warn">numerosita\' insufficiente</span>'
    elif d.ci_low > 0:
        verdict_cls, tag = "verdict-warn", '<span class="tag tag-warn">discordanza reale rilevata</span>'
    elif d.ci_high < 0:
        verdict_cls, tag = "verdict-warn", '<span class="tag tag-warn">direzione inattesa (probabile rumore)</span>'
    else:
        verdict_cls, tag = "verdict-yes", '<span class="tag tag-ok">nessuna discordanza oltre il rumore</span>'

    html_parts = [
        "<h2>Modulo B — Validazione dell'annotazione via TCR</h2>",
        f'<div class="card {verdict_cls}">',
        f"<p>{tag}</p>",
        "<table><tr><th>Quantita'</th><th>Valore</th><th>Cosa significa</th></tr>",
        f"<tr><td>Cellule con TCR</td><td>{result.n_cells_with_tcr:,}</td>"
        f"<td>Cellule per cui e' stato ricostruito un clonotipo (catena TRB rilevata).</td></tr>",
        f"<tr><td>Coppie clone-compartimenti confrontabili</td><td>{d.n_pairs:,}</td>"
        f"<td>Ogni riga confronta lo stesso clone T in due compartimenti diversi.</td></tr>",
        f"<tr><td>Pazienti che contribuiscono</td><td>{d.n_patients}</td>"
        f"<td>L'unita' statistica indipendente e' il paziente, non il clone: i cloni di uno stesso "
        f"paziente non sono osservazioni indipendenti.</td></tr>",
    ]
    if d.sufficient:
        html_parts.append(
            f"<tr><td>Eccesso di discordanza</td><td>{d.mean_excess:+.3f} "
            f"[{d.ci_low:+.3f}, {d.ci_high:+.3f}]</td>"
            f"<td>Quanto piu' spesso le cellule dello stesso clone T ricevono etichette diverse fra "
            f"compartimenti rispetto a quanto ci si aspetterebbe dal solo rumore entro un "
            f"compartimento. Zero nell'intervallo = non distinguibile dal rumore.</td></tr>"
        )
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
                val = "numerosita' insufficiente"
            html_parts.append(
                f"<tr><td>{html.escape(r.comp_a)} ↔ {html.escape(r.comp_b)}</td>"
                f"<td>{r.n_pairs}</td><td>{r.n_patients}</td><td>{val}</td></tr>"
            )
        html_parts.append("</table>")

    if result.marker_error is not None and result.marker_error.by_compartment:
        html_parts.append("<h3>Tasso d'errore per compartimento (marcatori canonici)</h3>")
        html_parts.append(
            f'<p class="narrative">Identita\' di riferimento del clone stimata SOLO dal compartimento '
            f"'{html.escape(result.marker_error.reference_compartment)}' "
            f"({result.marker_error.n_resolved_clones} cloni risolti), confrontata con l'etichetta "
            f"assegnata negli altri compartimenti.</p>"
        )
        chart2 = _marker_error_chart(result.marker_error)
        if chart2:
            html_parts.append(f'<img src="{chart2}" alt="Tasso d\'errore per compartimento">')
        html_parts.append("<table><tr><th>Compartimento</th><th>Tasso d'errore [IC95%]</th><th>Pazienti</th></tr>")
        for comp, res in result.marker_error.by_compartment.items():
            val = f"{res.mean:.3f} [{res.ci_low:.3f}, {res.ci_high:.3f}]" if res.sufficient else "numerosita' insufficiente"
            html_parts.append(f"<tr><td>{html.escape(comp)}</td><td>{val}</td><td>{res.n_groups}</td></tr>")
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
            f"<code>audit_reference_label</code> (identita' del clone stimata nel sangue) e "
            f"<code>audit_label_vs_reference</code> (l'etichetta assegnata discorda da quell'identita'). "
            f"Le etichette originali non vengono modificate. Copertura: {n_ref:,} cellule su {n:,} "
            f"({n_ref / n:.1%}) hanno un'identita' di riferimento; {n_eval:,} ({result.flag_coverage:.1%}) "
            f"sono valutabili, e di queste {n_disc:,} sono discordanti. Tutte le altre sono NA: "
            f"lo strumento non poteva verificarle (nessun TCR, clone senza cellule sufficienti nel "
            f"riferimento, etichetta non CD4/CD8, o cellula del compartimento di riferimento). NA non "
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
    immagini PNG incorporate come data URI."""
    now = datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M UTC")
    body = [
        f"<h1>Report di audit — {html.escape(dataset_name)}</h1>",
        f'<p class="footer" style="margin-top:0;border-top:none;padding-top:0;">'
        f"Generato il {now}. Strumento diagnostico per un singolo studio pilota: "
        f"non e' una certificazione, e' un supporto alla decisione per chi analizza i dati.</p>",
    ]
    if design_result is not None:
        body.append(_design_section(design_result))
    if leakage_result is not None:
        body.append(_leakage_section(leakage_result))
    if tcr_result is not None:
        body.append(_tcr_section(tcr_result))
    if leakage_result is None and tcr_result is None and design_result is None:
        body.append("<p>Nessun risultato da mostrare: esegui almeno un modulo.</p>")
    body.append(
        '<div class="footer">Report generato da genomic-audit. L\'audit del disegno legge solo i '
        "metadati e segnala quali confronti la struttura dello studio permette. Il Modulo A misura quanto una "
        "valutazione che non raggruppa per paziente sovrastimerebbe l'accuratezza e sottostimerebbe "
        "l'incertezza. Il Modulo B misura, tramite il repertorio T-cell receptor, quanto l'annotazione "
        "di tipo cellulare da clustering e' internamente coerente per uno stesso clone attraverso i "
        "compartimenti tissutali. Nessuno dei due modulo sostituisce una validazione biologica "
        "indipendente.</div>"
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
    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(render_report(leakage_result, tcr_result, dataset_name, design_result),
                        encoding="utf-8")
    return out_path
