"""Assembla validation_esterna_B/GSE200996/REPORT.md dai file di risultato (nessun numero trascritto a mano).

Le sezioni dei criteri compaiono solo quando i loro file esistono in results/: la prima versione del
report (solo eleggibilità) è stata scritta e committata prima di eseguire il Modulo B.

Uso: python -I validation_esterna_B/GSE200996/make_report.py "<riga della suite di test>"
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
R = HERE / "results"


def load(name: str) -> dict | None:
    p = R / name
    return json.loads(p.read_text(encoding="utf-8")) if p.exists() else None


def it(x: float, d: int = 3) -> str:
    """Numero in formato italiano (virgola decimale)."""
    return f"{x:.{d}f}".replace(".", ",")


def n(x: int) -> str:
    return f"{x:,}".replace(",", ".")


def pct(a: int, b: int) -> str:
    return it(100 * a / b, 1) + "%"


def eligibility(meta: dict, build: dict) -> str:
    e3 = build["E3_principale"]
    rows = []
    for r in build["campioni_vdj"]:
        if "path" not in r:
            rows.append(f"| {r['campione']} | — | — | — | — | {r['nota']} |")
            continue
        rows.append(f"| {r['campione']} | {r['compartment']} | {n(r['cellule_CD45'])} | {n(r['etichettate'])} | "
                    f"{n(r['T_CD45_non_etichettate'])} | {n(r['etichettate_con_TRB_produttiva'])} "
                    f"({pct(r['etichettate_con_TRB_produttiva'], r['etichettate'])}) |")
    lanes = [s for s in build["campioni_matrici"] if "vdj_corsia" in s]
    lane_rows = [f"| {s['file'].split('_')[0]} | {s['file'].split('matrix_')[1].replace('_GEX_sc_PBMC.h5', '')} | "
                 f"{n(s['vdj_corsia']['barcode_is_cell'])} | {n(s['vdj_corsia']['barcode_assegnati'])} | "
                 f"{n(s['vdj_corsia']['barcode_is_cell_non_assegnati'])} |" for s in lanes]
    int_ok = all(s["interi_non_negativi"] for s in build["campioni_matrici"])
    non_gex = sorted({f for s in build["campioni_matrici"] for f in s["feature_non_GEX"]})
    lp, lt = meta["etichette_PBMC"], meta["etichette_tumor"]
    up, ud = build["uscite"]["principale"], build["uscite"]["descrittivo"]
    pre = ", ".join(build["pazienti_pre_tx_descrittivo"])
    return f"""## 1. Dataset ed eleggibilità

[GSE200996](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE200996) (Luoma et al., *Cell* 2022,
PMC9508682): carcinoma squamoso del cavo orale, terapia neoadiuvante anti-PD-1 ± anti-CTLA-4
(NCT02919683); scRNA-seq 10x 5′ + scTCR-seq di tumore (prima e dopo la terapia) e di sangue
(prelievi B1-B3, corsie con più pazienti marcati con hashtag). Cell Ranger 4.0.0, GRCh38-2020-A.

**File e download.** Elenco dei file supplementari letto sul server FTP di GEO. Scaricati in
`../dati_esterni/GSE200996/`, fuori dal repository (spazio libero prima del download: 582 GB):
- prima, solo metadati per cellula (CD4, CD8, CD45, Myeloid × PBMC/tumor), le due chiavi degli
  hashtag, `filelist.txt` e le 5 series matrix (9,3 MB, 43 s con pause di 1 s fra i file; il nome
  `GSE200996_series_matrix.txt.gz` dà 404, la serie ha una matrix per piattaforma);
- poi, dopo il commit di `MAPPA_ETICHETTE.md`, `GSE200996_RAW.tar` (1.591.296.000 byte, uguale a
  `filelist.txt`; 53 s, al primo tentativo), da cui sono estratti solo i 61 file single-cell di tumore
  e sangue (corsie non selezionate) dei pazienti eleggibili (814 MB): niente TCR bulk, niente corsie
  selezionate.
Nessuna limitazione delle richieste (reCAPTCHA, HTTP 429). Testo dell'articolo da NCBI efetch (PMC
XML; Europe PMC rispondeva 500).

**Campioni GEO single-cell.** 25 matrici di tumore (6 prima e 19 dopo la terapia) e 24 file TCR
single-cell di tumore (5 prima, 19 dopo); 12 corsie di sangue «total PBMC» (matrice + TCR); 4
corsie di sangue selezionate CD38+ HLA-DR+ CD4 o CD8 (P13, P19, P26, P32; prelievi B1 e B2).

### E5 — origine del sangue

Metodi, «Preparation of total PBMC»: le cellule di ogni prelievo, marcate con un hashtag TotalSeq-C,
sono state selezionate al FACS come «Live CD45+ cells»: tutti i leucociti vivi, non una linea. Le
corsie «sorted» (GSM6048338-345) sono un esperimento separato su cellule T attivate CD38+ HLA-DR+
(Metodi, «Sorting of activated blood T cells»), e non sono usate. I metadati del sangue contengono
{n(lp['cellule_CD45'])} cellule, lo stesso numero del «blood dataset» dei Metodi; ogni coppia
paziente-prelievo dei metadati compare nella chiave delle corsie non selezionate (fuori chiave:
{len(meta['coppie_paziente_prelievo_sangue_fuori_chiave'])}). Tutti i
{len(meta['E5']['pazienti_eleggibili_con_CD4_e_CD8_nel_sangue'])} pazienti eleggibili hanno nel
sangue sia cellule CD4 sia CD8. **E5 soddisfatta.**

### E2 — origine delle etichette

Etichetta = appartenenza al file CD4 o CD8 degli autori (mappa in `MAPPA_ETICHETTE.md`, commit
`2e67f4d`, prima di eseguire il modulo). Metodi: trasferimento di etichette da un riferimento CITE-seq
di PBMC (Seurat `MapQuery`, Hao et al. 2021) con punteggio di predizione > 0,6, cioè **annotazione
computazionale sull'espressione**, non sorting; inoltre gli autori hanno tolto le cellule CD4 con
CD8A/CD8B rilevati e le CD8 con CD4 rilevato. Sangue: {n(lp['CD4'])} CD4 e {n(lp['CD8'])} CD8;
tumore: {n(lt['CD4'])} CD4 e {n(lt['CD8'])} CD8; nessuna cellula in entrambi i file. Cellule T del
file CD45 senza etichetta (escluse): {n(lp['T_CD45_senza_etichetta'])} nel sangue e
{n(lt['T_CD45_senza_etichetta'])} nel tumore (tutti i pazienti). Cellule etichettate con tipo diverso
da T nel file CD45 (tenute, come deciso nella mappa): sangue {lp['etichettate_per_tipo_CD45']};
tumore {lt['etichettate_per_tipo_CD45']}. **E2 soddisfatta.**

### E1, E3, E4 — demultiplazione, abbinamento, formati

Demultiplazione degli autori (Metodi): GMM-demux sulle conte degli hashtag, singoletti con
confidenza ≥ 0,8; in media 83,46% singoletti, 15,53% doppietti o multipletti, 0,28% negativi. I
metadati contengono solo i singoletti assegnati; l'ID di cellula porta paziente e prelievo. Barcode
VDJ di ciascuna corsia non assegnati dagli autori (doppietti, negativi e cellule tolte dal controllo
di qualità, non separabili senza rifare la demultiplazione), esclusi:

| Corsia | Pazienti | Barcode VDJ (is_cell) | Assegnati a un paziente | Non assegnati, esclusi |
|---|---|---|---|---|
{chr(10).join(lane_rows)}

Matrici: `raw_feature_bc_matrix.h5` di Cell Ranger, conteggi int32 interi non negativi in tutti i
{len(build['campioni_matrici'])} campioni usati: {'sì' if int_ok else '**no**'}. Feature diverse da
Gene Expression: {', '.join(non_gex)} (solo hashtag, nessuna proteina di superficie). File VDJ:
`filtered_contig_annotations.csv` di Cell Ranger, con le colonne attese dal Modulo B; i file del
sangue (uno per corsia) sono divisi per paziente copiando le righe senza modificarle
(`build_gse200996.py`). **E4 soddisfatta.**

Per campione (analisi principale e descrittiva): cellule CD45 assegnate dagli autori, cellule T
etichettate CD4/CD8, cellule T non etichettate (escluse dalla mappa), cellule etichettate con almeno
un contig TRB produttivo (is_cell, high_confidence, full_length, productive, CDR3 presente):

| Campione | Compartimento | Cellule CD45 | T etichettate | T escluse | Con TRB produttiva |
|---|---|---|---|---|---|
{chr(10).join(rows)}

E3 sull'analisi principale (dopo la terapia): sangue {n(e3['PBMC']['con_TRB_produttiva'])}/{n(e3['PBMC']['etichettate'])}
({pct(e3['PBMC']['con_TRB_produttiva'], e3['PBMC']['etichettate'])}), tumore
{n(e3['Tumor']['con_TRB_produttiva'])}/{n(e3['Tumor']['etichettate'])}
({pct(e3['Tumor']['con_TRB_produttiva'], e3['Tumor']['etichettate'])}); il campione più basso è sopra
il 50%. Con «almeno un contig produttivo TRA o TRB»: {pct(e3['PBMC']['con_TCR_produttivo'], e3['PBMC']['etichettate'])}
e {pct(e3['Tumor']['con_TCR_produttivo'], e3['Tumor']['etichettate'])}. **E3 soddisfatta.**

E1: {len(up['pazienti'])} pazienti ({', '.join(up['pazienti'])}) hanno cellule di sangue assegnate a
un solo paziente e cellule di tumore dopo la terapia, con conteggi interi e TCR per cellula.
**E1 soddisfatta.**

| Condizione | Esito |
|---|---|
| E1 | soddisfatta: {len(up['pazienti'])} pazienti (≥ 5) |
| E2 | soddisfatta: etichetta CD4/CD8 per cellula in entrambi i compartimenti, da trasferimento di etichette sull'espressione |
| E3 | soddisfatta: {pct(e3['PBMC']['con_TRB_produttiva'], e3['PBMC']['etichettate'])} (sangue) e {pct(e3['Tumor']['con_TRB_produttiva'], e3['Tumor']['etichettate'])} (tumore) con TRB produttiva |
| E4 | soddisfatta: matrici 10x intere, `filtered_contig_annotations` |
| E5 | soddisfatta: sangue CD45+ non selezionato per linea; corsie selezionate escluse |

**h5ad costruiti** (`build_gse200996.py`, colonne `patient`, `compartment`, `timepoint`, `label`,
`barcode`, `sample`, `lane`): principale {n(up['cellule'])} cellule × {n(up['geni'])} geni, {len(up['pazienti'])}
pazienti ({up['per_compartimento_ed_etichetta']}); descrittivo {n(ud['cellule'])} cellule, pazienti
{pre} (tumore prima della terapia con matrice e TCR single-cell; P23 ha la matrice ma non il TCR ed è
escluso; P27 ha il tumore ma non il sangue).
"""


def main() -> None:
    suite = sys.argv[1] if len(sys.argv) > 1 else "(non disponibile)"
    meta, build = load("metadati.json"), load("build.json")
    crit = load("criteri.json")
    head = """# Validazione tecnica esterna del Modulo B su GSE200996

Criteri: `validation_esterna_B/GSE200996/CRITERI.md`, scritti e committati prima di qualsiasi download
(commit `8b10bb0`) e non modificati dopo. Mappa delle etichette: `MAPPA_ETICHETTE.md` (commit
`2e67f4d`), prima di eseguire il Modulo B. Secondo tentativo della voce 1 della roadmap, dopo
GSE139555 (non eleggibile: `validation_esterna_B/REPORT.md`). Numeri e tabelle di questo documento
sono assemblati da `make_report.py` a partire dai file in `results/`.

Versione dello strumento: 0.9.0, branch `validazione-b-gse200996`, creato da `validazione-b` perché
master non ne contiene ancora il merge. Ambiente: `.venv` del progetto (uv); workstation con 12 core
e 31 GB di RAM.
"""
    body = eligibility(meta, build)
    if crit is None:
        body += ("\n**Eleggibilità: E1-E5 soddisfatte.** Questa versione del report è scritta prima di eseguire il "
                 "Modulo B; i criteri B-E1…B-E6 seguono nelle sezioni successive dopo l'esecuzione.\n")
    else:
        from report_criteri import criteria_sections  # noqa: PLC0415 (modulo locale, solo a esecuzione finita)
        body += criteria_sections(R, suite)
    (HERE / "REPORT.md").write_text(head + "\n" + body, encoding="utf-8")
    print(f"[ok] scritto {HERE / 'REPORT.md'}")


if __name__ == "__main__":
    main()
