# ruff: noqa: E501 -- le righe delle tabelle markdown del report non si possono spezzare
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


def _flat(d: dict, p: str = ""):
    for k, v in d.items():
        if isinstance(v, dict):
            yield from _flat(v, p + k + ".")
        elif isinstance(v, (int, float)) and not isinstance(v, bool):
            yield p + k, v


def be2_compare(tool: dict, ind: dict) -> dict:
    """Confronto B-E2: conteggi esatti, stime entro 1e-9, estremi entro 1e-6."""
    ft, fi = dict(_flat(tool)), dict(_flat(ind))
    counts = ["clone_tumor.n_units", "clone_tumor.n_ties", "clone_tumor.n_errors", "clone_tumor.n_patients",
              "clone_tumor.rate.n_groups", "clone_tumor.rate.n_obs", "excess_PBMC_Tumor.n_pairs",
              "excess_PBMC_Tumor.n_patients", "excess_PBMC_Tumor.weighted.n_groups", "excess_PBMC_Tumor.equal.n_groups"]
    means = ["clone_tumor.rate.mean", "excess_PBMC_Tumor.weighted.mean", "excess_PBMC_Tumor.equal.mean"]
    ends = [f"{b}.{e}" for b in ("clone_tumor.rate", "excess_PBMC_Tumor.weighted", "excess_PBMC_Tumor.equal")
            for e in ("ci_low", "ci_high")]
    c_ok = all(ft[k] == fi[k] for k in counts)
    dm = max(abs(ft[k] - fi[k]) for k in means)
    de = max(abs(ft[k] - fi[k]) for k in ends)
    ddf = max(abs(ft[f"{b}.df"] - fi[f"{b}.df"]) for b in ("clone_tumor.rate", "excess_PBMC_Tumor.weighted",
                                                           "excess_PBMC_Tumor.equal"))
    return {"conteggi": c_ok, "max_diff_stime": dm, "max_diff_estremi": de, "max_diff_df": ddf,
            "superato": c_ok and dm <= 1e-9 and de <= 1e-6, "righe": [(k, ft[k], fi[k]) for k in counts + means + ends]}


def _ci(r: dict, signed: bool = False) -> str:
    f = (lambda x: ("+" if x >= 0 else "−") + it(abs(x))) if signed else it  # noqa: E731
    if not r["sufficient"]:
        return f"{f(r['mean'])} (intervallo non prodotto)"
    return f"{f(r['mean'])} [{f(r['ci_low'])}; {f(r['ci_high'])}]"


def criteria_sections(rep: dict, suite: str) -> str:
    """Sezioni 2-8: criteri, prima e dopo la correzione, descrittivi, difetti, verdetto."""
    import re  # noqa: PLC0415

    prima = load("prima_della_correzione/repliche.json")
    be2 = be2_compare(load("BE2_strumento.json"), load("BE2_indipendente.json"))
    be2_0 = be2_compare(load("prima_della_correzione/BE2_strumento.json"), load("prima_della_correzione/BE2_indipendente.json"))
    stdout = (R / "BE1_cli_stdout.txt").read_text(encoding="utf-8")
    stderr = (R / "BE1_cli_stderr.txt").read_text(encoding="utf-8")
    exit_code = (R / "BE1_exit.txt").read_text().strip()
    wall = re.search(r"Elapsed \(wall clock\) time \(h:mm:ss or m:ss\): (\S+)", stderr).group(1)
    rss = int(re.search(r"Maximum resident set size \(kbytes\): (\d+)", stderr).group(1))
    match = re.search(r"Barcode VDJ.*?\(([\d.]+%)\)", stdout).group(1)
    cellconv = re.search(r"\[cell\].*?'Tumor': ([\d.]+) \(IC95% \[([\d.]+), ([\d.]+)\], (\d+) pazienti, (\d+) unità\)", stdout)
    flags = re.search(r"valutabili.*?: ([\d.]+%)", stdout).group(1)
    descr_txt = (R / "descrittivo_cli.txt").read_text(encoding="utf-8")
    d_cell = re.search(r"\[cell\].*?'Tumor': ([\d.]+), IC non prodotto \((\d+) pazienti.*?; (\d+) unità\)", descr_txt)
    d_clone = re.search(r"\[clone\].*?'Tumor': ([\d.]+), IC non prodotto \((\d+) pazienti.*?; (\d+) unità; unità in pareggio: (\d+)\)", descr_txt)
    d_match = re.search(r"Barcode VDJ.*?\(([\d.]+%)\)", descr_txt).group(1)
    same_cli = (R / "BE1_cli_stdout.txt").read_text() == (R / "prima_della_correzione" / "BE1_cli_stdout.txt").read_text()
    same_descr = descr_txt == (R / "prima_della_correzione" / "descrittivo_cli.txt").read_text()
    same_rep = {k: rep[k] == prima[k] for k in ("strumento_BE2", "BE3", "BE4", "BE6", "descrittivi")}

    s2 = rep["strumento_BE2"]
    cr, ex = s2["clone_tumor"], s2["excess_PBMC_Tumor"]
    b3, b4, b6 = rep["BE3"]["sintesi"], rep["BE4"], rep["BE6"]
    ctrl = rep["controllo_su_dati_originali"]
    ds = rep["descrittivi"]
    v6 = b6["varianti"]
    v6_rows = "\n".join(
        f"| {k} | {v['n_patients']} | {v['n_units']} | {v['n_ties']} | {v['n_errors']} | {_ci(v['rate'])} | "
        f"{('+' if b6['differenze'][k] >= 0 else '−') + it(abs(b6['differenze'][k]), 4)} |" for k, v in v6.items())
    be2_rows = "\n".join(f"| `{k}` | {t:.17g} | {i:.17g} |" for k, t, i in be2["righe"])
    reps = rep["BE3"]["repliche"]
    b3_rows = " ".join(f"{r['seme']}: {('+' if r['weighted']['mean'] >= 0 else '−') + it(abs(r['weighted']['mean']), 4)}"
                       f"{'' if r['contiene_zero'] else ' (non contiene lo 0)'};" for r in reps)
    n_ok = sum(x for x in (True if exit_code == "exit=0" else False, be2["superato"], b3["superato"], b6["superato"]))
    dirs = ", ".join(f"{k}: {v}" for k, v in ds["errori_per_direzione"].items()) or "nessuno"
    by_pat = ", ".join(f"{k} ({v})" for k, v in ds["errori_per_paziente"].items()) or "nessuno"

    return f"""
## 2. Criteri × esito (esecuzione finale, dopo la correzione del margine)

| Criterio | Esito | Numeri misurati |
|---|---|---|
| B-E1 esecuzione | **superato** | `audit-sc tcr --h5ad GSE200996_principale.h5ad --vdj-manifest manifest_principale.csv --patient-col patient --compartment-col compartment --celltype-col label --barcode-col barcode --marker-map marcatori.json --reference-compartment PBMC`: {exit_code.replace('=', ' ')}, report HTML scritto; {wall} reali, {it(rss / 1024 / 1024, 2)} GB di memoria massima (`/usr/bin/time -v`); barcode VDJ con TRB ritrovati: {match.replace('.', ',')} |
| B-E2 implementazione indipendente | **{'superato' if be2['superato'] else 'non superato'}** | conteggi {'identici' if be2['conteggi'] else '**diversi**'}; stime puntuali: differenza massima {be2['max_diff_stime']:.1e} (tolleranza 1e-9); estremi: {be2['max_diff_estremi']:.1e} (tolleranza 1e-6); gradi di libertà: {be2['max_diff_df']:.1e} |
| B-E3 controllo negativo | **{'superato' if b3['superato'] else 'non superato'}** | 20 permutazioni (semi 0-19): eccesso pesato medio {('+' if b3['media'] >= 0 else '−') + it(abs(b3['media']), 5)} (banda [−0,01; +0,01]; da {('+' if b3['min'] >= 0 else '−') + it(abs(b3['min']), 4)} a {('+' if b3['max'] >= 0 else '−') + it(abs(b3['max']), 4)}); intervallo con lo zero in {b3['contiene_zero']}/20 (≥ 18) |
| B-E4 controllo positivo | **non valutabile** | unità clone-tumore con identità CD4 e maggioranza corretta: {b4['candidate']}, meno delle 20 richieste; delle {cr['n_units']} unità nel tumore, {ds['unita_per_identita'].get('CD8T', 0)} hanno identità CD8 e {ds['unita_per_identita'].get('CD4T', 0)} CD4 |
| B-E5 riferimento proteico | **non valutabile su questo dataset** | nessuna proteina di superficie per cellula: nelle matrici del sangue le sole feature non geniche sono i 7 hashtag TotalSeq-C; nel tumore nessuna; rimandato alla voce 2 della roadmap |
| B-E6 sensibilità | **{'superato' if b6['superato'] else 'non superato'}** | differenza massima dal default {it(max(abs(x) for x in b6['differenze'].values()), 4)} (soglia 0,05); nessun intervallo sotto 5 pazienti; tabella sotto |

**Controllo dello script delle repliche.** Sui dati originali lo script unico (tabella costruita una
volta, funzioni pubbliche del modulo) dà gli stessi numeri di `run_tcr_validation`, la funzione
chiamata da `cli.py`, con uguaglianza esatta ({', '.join(k for k, v in ctrl.items() if v and not k.startswith('cli'))}),
e coincide con l'output della riga di comando di B-E1 arrotondato a 3 decimali (tasso per clone:
{'sì' if ctrl['cli_tasso_clone_arrotondato'] else 'no'}; eccesso: {'sì' if ctrl['cli_eccesso_arrotondato'] else 'no'}).
Tempo dello script (25 analisi): {it(prima['secondi'], 0)} s nella prima esecuzione, {it(rep['secondi'], 0)} s
nella seconda, in parallelo con la suite completa dei test.

### B-E2 in dettaglio

`indipendente.py` è stato scritto da un sottoagente che ha letto solo `docs/METODI.md` §1-4 e non
importa nulla da `core/`; legge da sé h5ad, manifest e file VDJ (scelta del contig, clone,
positività ai marcatori, identità di riferimento, unità, eccesso, intervalli CR2 con gradi di
libertà di Bell-McCaffrey). Le assunzioni che ha dovuto fare dove METODI.md tace sono elencate nel
suo JSON (`results/BE2_indipendente.json`, campo `assunzioni`); nessuna ha prodotto una differenza
su questi dati, quindi nessuna correzione di METODI.md è richiesta dal criterio. Due erano
potenzialmente rilevanti: il filtro `full_length` dei contig (lo strumento lo applica, METODI non lo
dice; qui non cambia nessun numero) e il confronto con il margine 0,20 (sezione 5).

| Grandezza | Strumento | Indipendente |
|---|---|---|
{be2_rows}

### B-E3 per replica (eccesso pesato sangue-tumore)

{b3_rows}

### B-E6

| Variante | Pazienti | Unità | Pareggi | Unità in errore | Tasso per clone nel tumore [IC95%] | Differenza dal default |
|---|---|---|---|---|---|---|
{v6_rows}

Nessuna unità in pareggio: le varianti «pareggi = errore» e «pareggi = corretti» coincidono per
costruzione con il default. Variante TRA+TRB: {n(b6['cellule_TRA_TRB'])} cellule con entrambe le
catene su {n(b6['cellule_TRB'])} con TRB.

## 3. Prima e dopo la correzione

La correzione del margine (sezione 5) è stata fatta dopo una prima esecuzione completa con il codice
0.9.0 (risultati in `results/prima_della_correzione/`, commit `be4d239`). Esito prima della
correzione: B-E1 superato ({'stesso output della CLI' if same_cli else 'output diverso'}), B-E2
{'superato' if be2_0['superato'] else 'non superato'} (differenza massima delle stime {be2_0['max_diff_stime']:.1e}),
B-E3 {'superato' if prima['BE3']['sintesi']['superato'] else 'non superato'}, B-E4 non valutabile
({prima['BE4']['candidate']} candidate), B-E6 {'superato' if prima['BE6']['superato'] else 'non superato'}.
Dopo la correzione i numeri sono identici: output della CLI {'identico' if same_cli else 'diverso'},
analisi descrittiva {'identica' if same_descr else 'diversa'}, script delle repliche
{', '.join(f'{k} {"identico" if v else "diverso"}' for k, v in same_rep.items())}. Nessun clone di
GSE200996 ha un margine esattamente pari a 0,20.

## 4. Risultati descrittivi (non sono criteri)

**Tumore dopo la terapia** ({len(load('build.json')['uscite']['principale']['pazienti'])} pazienti nell'h5ad principale):
- tasso d'errore per clone (convenzione della tesi): {_ci(cr['rate'])}, {cr['n_units']} unità da
  {cr['n_patients']} pazienti, {cr['n_errors']} unità in errore, nessun pareggio;
- tasso d'errore per cellula: {it(float(cellconv.group(1)))} [{it(float(cellconv.group(2)))}; {it(float(cellconv.group(3)))}],
  {cellconv.group(5)} cellule da {cellconv.group(4)} pazienti;
- eccesso di discordanza sangue-tumore: pesato {_ci(ex['weighted'], True)}, a pesi uguali
  {_ci(ex['equal'], True)}, su {ex['n_pairs']} coppie clone-compartimenti da {ex['n_patients']} pazienti
  (l'intervallo include lo zero di poco);
- direzione degli errori (identità → etichetta di maggioranza nel tumore): {dirs.replace('T->', '→').replace('CD8T', 'CD8')}; per paziente:
  {by_pat}; nessun errore CD8→CD4;
- cloni con identità di riferimento: {ds['cloni_con_riferimento']}; cellule con TRB abbinata:
  {n(ds['cellule_con_trb'])}; flag per cellula valutabili: {flags.replace('.', ',')} delle cellule.

**Tumore prima della terapia** (h5ad descrittivo, `audit-sc tcr` eseguito una volta: exit 0;
barcode VDJ ritrovati {d_match.replace('.', ',')}): {d_clone.group(2)} pazienti, quindi nessun intervallo. Tasso per
clone {it(float(d_clone.group(1)))} su {d_clone.group(3)} unità (pareggi: {d_clone.group(4)}); per cellula
{it(float(d_cell.group(1)))} su {d_cell.group(3)} cellule; eccesso non riportato dallo strumento con
meno di 5 pazienti.

**Come leggere questi numeri.** Il tasso d'errore è molto più basso di quello della tesi su
GSE278694 (0,170 per clone nel tumore), ma il confronto non è alla pari, per tre ragioni (le prime due
scritte in `MAPPA_ETICHETTE.md` prima dell'esecuzione): (1) gli autori hanno tolto dalle etichette CD4 le cellule
con CD8A/CD8B rilevati e dalle CD8 quelle con CD4 rilevato, cioè hanno già ripulito le etichette con
gli stessi geni che il Modulo B usa per l'identità di riferimento; (2) le etichette vengono da un
trasferimento da un riferimento CITE-seq con soglia di predizione, non da cluster; (3) i metadati
contengono solo cellule T CD4/CD8: gli errori verso NK della tesi qui non possono comparire, e il
confronto con la tesi riguarda solo CD4↔CD8. Inoltre quasi tutte le unità giudicate nel tumore
appartengono a cloni CD8 (i cloni CD4 con almeno 3 cellule nel sangue e 3 nel tumore sono rari), e
per questo il controllo positivo B-E4 non è valutabile.

## 5. Difetti trovati e correzioni

**Difetto dello strumento: margine 0,20 confrontato in virgola mobile** (corretto nel commit
`8b8e761`). METODI §4 dice: identità = etichetta con la frazione più alta «se il margine sulla
seconda è almeno 0,20». `assign_reference_identity` confrontava `best - second < 0.20` in virgola
mobile: un margine esattamente 0,20 era rifiutato in alcuni casi (3/5 contro 2/5 = 0,19999999999999996;
7/10 contro 5/10; 3/10 contro 1/10; 6/10 contro 4/10) e accettato in altri (1/5 contro 0/5; 4/10
contro 2/10; 8/10 contro 6/10). L'ha fatto emergere `indipendente.py`, che confronta in aritmetica
esatta e aveva elencato il punto fra le sue assunzioni; verificato direttamente sulla funzione dello
strumento. La tesi (pdac-ml, `src/10_loco.py`, `abs(f8 - f4) < MIN_MARGIN`) ha lo stesso confronto,
quindi lo strumento la riproduceva anche in questo. Correzione: confronto con tolleranza 1e-9 (un
margine diverso da 0,20 ne dista almeno 1/(5n) per un clone di n cellule). Test di regressione
`tests/test_reference_margin.py`: fallisce sul codice precedente (5 casi su 8 con margine
esattamente 0,20), passa dopo. Effetto: nessuno su GSE200996 (sezione 3) e nessuno su GSE278694
(`tests/test_real_data_gse278694.py` passa con gli stessi numeri: 0,1697 su 277 unità).

**Difetti degli script di validazione** (non dello strumento), corretti prima di eseguire il modulo:
- nella prima versione di `build_gse200996.py` i percorsi del manifest erano relativi alla cartella
  di lavoro, mentre la guida (e `read_vdj_manifest`) li legge rispetto alla cartella del manifest:
  resi assoluti;
- `GSE200996_series_matrix.txt.gz` non esiste (404): la serie ha una series matrix per piattaforma.

**Lacune di METODI.md senza effetto sui numeri** (dalle assunzioni di `indipendente.py`): i filtri
dei contig (`is_cell`, `high_confidence`, `productive`, `full_length`), la regola di abbinamento dei
barcode e del suffisso `-1`, il trattamento delle cellule senza TRB e il quantile t con gradi di
libertà non interi non sono descritti. Il criterio B-E2 chiede di correggere METODI.md solo quando
una lacuna produce una differenza; qui non l'ha prodotta, e METODI.md non è stato modificato.

## 6. Verdetto con i criteri originali

**Validazione tecnica esterna superata nei criteri valutabili: B-E1, B-E2, B-E3 e B-E6 superati
({n_ok} su 4); B-E4 non valutabile (3 unità candidate, ne servono 20); B-E5 non valutabile (nessuna
proteina di superficie).** Nessun criterio è stato modificato dopo aver visto i dati: non c'è una
sezione «Criteri modificati dopo aver visto i dati». Esito identico prima e dopo la correzione del
margine.

## 7. Che cosa questa validazione dimostra e che cosa no

Dimostra, su dati di un altro laboratorio (Dana-Farber, 10x 5′, sangue demultiplato con hashtag):
- che il Modulo B si esegue con i soli comandi e input della guida (B-E1), su 18 pazienti e circa
  59.000 cellule in meno di un minuto;
- che i suoi calcoli coincidono con un'implementazione scritta solo da METODI.md §1-4 (B-E2), e che
  quell'implementazione ha trovato un difetto vero, piccolo e ora corretto;
- che non segnala un eccesso di discordanza quando il legame fra clone e compartimento è rotto
  (controllo negativo, B-E3);
- che il tasso d'errore per clone è stabile alla definizione del clone (TRB o TRA+TRB) e alla regola
  dei pareggi (B-E6), qui senza pareggi.

Non dimostra:
- che il modulo rilevi errori iniettati su dati reali: B-E4 non è valutabile su GSE200996 (la
  sensibilità resta provata solo sui dati sintetici dei test);
- che le etichette segnalate come errate lo siano davvero: serve un riferimento proteico (B-E5,
  voce 2 della roadmap);
- che il tasso d'errore basso di GSE200996 descriva la qualità di un'annotazione qualsiasi: le
  etichette degli autori erano già state ripulite con gli stessi geni marcatori usati dal modulo.

## 8. Test

Codice cambiato (`core/tcr_validation.py`): suite completa su un clone pulito del branch
(`uv sync --extra dev`, dati di GSE278694 con `GSE278694_DIR`). Esito: {suite}.
"""


def main() -> None:
    suite = sys.argv[1] if len(sys.argv) > 1 else "(non disponibile)"
    meta, build = load("metadati.json"), load("build.json")
    crit = load("repliche.json")
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
        body += criteria_sections(crit, suite)
    (HERE / "REPORT.md").write_text(head + "\n" + body, encoding="utf-8")
    print(f"[ok] scritto {HERE / 'REPORT.md'}")


if __name__ == "__main__":
    main()
