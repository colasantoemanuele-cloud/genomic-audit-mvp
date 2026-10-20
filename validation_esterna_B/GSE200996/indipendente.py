# ruff: noqa: I001, B007 -- file dell'implementazione indipendente, lasciato com'era scritto (solo stile)
"""Implementazione indipendente di alcune statistiche del Modulo B.

Scritta partendo solo da docs/METODI.md, sezioni 1-4, senza consultare il codice dello strumento.
Ricalcola, per il compartimento giudicato "Tumor" con riferimento "PBMC":
- il tasso d'errore con la convenzione per clone (METODI §4), con intervallo t cluster-robusto su
  scala logit (METODI §2);
- l'eccesso di discordanza pesato e a pesi uguali per la coppia PBMC-Tumor (METODI §3), con
  intervallo t cluster-robusto su scala originale.

Ogni punto in cui METODI.md non basta è dichiarato nella lista ASSUNZIONI, copiata nell'output.

Uso:
    python -I indipendente.py --h5ad X.h5ad --vdj-manifest manifest.csv \
        --marker-map mappa.json --reference-compartment PBMC --tie-rule exclude --out out.json
    python -I indipendente.py --selftest
"""

from __future__ import annotations

import argparse
import json
import math
import os
import re
import sys
import tempfile
from collections import Counter

import numpy as np
import pandas as pd
from scipy import sparse
from scipy import stats

MIN_PAZIENTI = 5          # METODI §1
MIN_CELLULE_RIF = 3       # METODI §4, identità di riferimento
MARGINE_NUM, MARGINE_DEN = 1, 5   # margine 0,20 come frazione esatta 1/5
MIN_CELLULE_CLONE = 3     # METODI §4, convenzione per clone
MIN_CELLULE_LATO = 2      # METODI §3
COMPARTIMENTO_GIUDICATO = "Tumor"
CATENA = "TRB"

ASSUNZIONI = [
    # --- contig e definizione del clone ---
    "Filtri sui contig: METODI non dice quali contig usare. Tengo solo le righe con chain == 'TRB', "
    "is_cell vero (se la colonna esiste), high_confidence vero, productive vero e cdr3_nt non vuota "
    "e diversa da 'None'. Non filtro su full_length. I valori booleani sono letti come veri solo se "
    "il testo è 'true' (senza distinzione di maiuscole); 'None' o vuoto valgono falso.",
    "Ordine filtri/scelta: la regola di scelta del contig (più UMI, poi più reads, poi CDR3 minore) "
    "è applicata DOPO i filtri, solo fra i contig TRB rimasti; METODI non dice se i filtri "
    "precedano la scelta.",
    "Nella regola di scelta del contig, 'la CDR3 minore in ordine alfabetico' è intesa come la "
    "cdr3_nt (nucleotidica), coerente con la definizione del clone; METODI non dice se sia la CDR3 "
    "aminoacidica o nucleotidica. UMI e reads sono le colonne 'umis' e 'reads', convertite in numero "
    "(valori non numerici = 0). L'ordine alfabetico è quello dei codici carattere (Python).",
    "Il clone è definito dalla sola cdr3_nt della TRB scelta; la TRA è ignorata e cellule senza TRA "
    "sono incluse. METODI non dice se cellule con TRB ma senza TRA, o con più TRB, siano escluse.",
    "Il clone è definito entro paziente: la stessa cdr3_nt in due pazienti dà due cloni distinti "
    "(METODI §3 dice 'per ogni paziente, clone T', non lo dice esplicitamente per §4).",
    "Le cellule dell'AnnData senza TRB valida (nessun contig che passi i filtri, o barcode assente "
    "dai file VDJ, o coppia paziente-compartimento assente dal manifest) sono escluse da tutto il "
    "Modulo B. I contig VDJ senza cellula corrispondente nell'AnnData sono ignorati.",
    # --- abbinamento dei barcode ---
    "Abbinamento dei barcode: METODI non lo descrive. Una cellula è la terna (patient, compartment, "
    "barcode); da entrambi i lati tolgo spazi e un eventuale suffisso finale '-<cifre>' (es. '-1'), "
    "e confronto le stringhe esattamente (maiuscole incluse). patient e compartment devono "
    "coincidere come stringhe fra AnnData e manifest.",
    "Terne (patient, compartment, barcode normalizzato) duplicate nell'AnnData o coppie "
    "paziente-compartimento ripetute nel manifest: errore bloccante, non una scelta silenziosa.",
    # --- etichette ---
    "Cellule con etichetta mancante (NaN o stringa vuota) nell'AnnData: escluse dalle convenzioni per "
    "clone e dall'eccesso di discordanza, ma incluse nel calcolo dell'identità di riferimento, che "
    "usa solo l'espressione. METODI non tratta le etichette mancanti.",
    "Identità di riferimento: le cellule del clone nel compartimento di riferimento sono contate "
    "tutte, qualunque sia la loro etichetta (anche fuori mappa); METODI non dice se le cellule con "
    "etichetta fuori mappa entrino nel denominatore.",
    # --- positività ai marcatori ---
    "Positività per un'etichetta con più marcatori (es. CD8T: CD8A, CD8B): una cellula è positiva se "
    "almeno uno dei geni ha conta > 0 (regola OR). METODI dice solo 'conta > 0' e non dice come "
    "combinare più geni.",
    "La positività non è esclusiva (senza --clone-marker-priority): una cellula positiva a CD4 e a "
    "CD8 conta per entrambe le etichette; una cellula negativa a tutto resta nel denominatore. La "
    "frazione di positivi è (cellule positive) / (tutte le cellule del clone nel riferimento).",
    "Geni marcatori assenti da var_names: ignorati (contano come zero); se nessun gene di "
    "un'etichetta è presente il calcolo si interrompe. Se un simbolo compare più volte in var_names "
    "vale l'OR su tutte le colonne con quel nome. Nessun controllo che X sia davvero intera: si usa "
    "X così com'è con il test > 0.",
    "Il margine 0,20 fra prima e seconda frazione è confrontato in aritmetica esatta: "
    "5·(pos_1 − pos_2) ≥ n. In virgola mobile, per esempio 0,6 − 0,4 = 0,19999999999999996 < 0,20, "
    "cioè un clone 3/5 contro 2/5 perderebbe il riferimento; METODI non specifica la precisione "
    "del confronto. Con una sola etichetta in mappa la seconda frazione vale 0. Se le prime due "
    "frazioni sono uguali il margine è 0 e il clone non ha riferimento.",
    "Etichette candidate per l'identità: solo le chiavi della mappa dei marcatori.",
    # --- convenzione per clone ---
    "Convenzione per clone: la soglia di 3 cellule della coppia clone-compartimento è contata sulle "
    "cellule con TRB valida ed etichetta non mancante; tutte le etichette osservate entrano nel "
    "conteggio della maggioranza (etichette fuori mappa ammesse: tutte, come il default di METODI); "
    "una maggioranza fuori mappa conta come errore. Unità solo per i cloni con riferimento.",
    "clone_tumor.n_patients = pazienti con almeno un'unità che entra nel tasso (uguale a "
    "rate.n_groups); i pazienti con soli pareggi esclusi non sono contati. METODI non definisce il "
    "conteggio dei pazienti accanto al tasso.",
    "Con tie-rule error/correct i pareggi entrano nel tasso (1 o 0) e sono comunque riportati in "
    "n_ties; n_errors include i pareggi contati come errore.",
    # --- eccesso di discordanza ---
    "Eccesso di discordanza: calcolato per tutti i cloni con almeno 2 cellule per lato, con o senza "
    "identità di riferimento (METODI §3 non richiede il riferimento). Le etichette k sono tutte le "
    "etichette osservate, anche fuori mappa. Il conteggio per lato usa le cellule con TRB valida ed "
    "etichetta non mancante.",
    "excess.n_patients = pazienti con almeno una coppia clone-compartimenti; n_pairs = numero di "
    "coppie (= n_obs di entrambi gli intervalli).",
    # --- intervallo ---
    "Quantile t con gradi di libertà non interi: scipy.stats.t.ppf(0.975, ν) con ν reale, senza "
    "arrotondare; METODI non dice se ν vada arrotondato.",
    "G = numero di pazienti con almeno un'osservazione; con G < 5 l'intervallo non è prodotto e df, "
    "se, ci_low, ci_high sono null, mentre mean, n_groups e n_obs sono riportati. Senza osservazioni "
    "anche mean è null.",
    "Quando l'intervallo non è prodotto (G < 5, θ̂ = 0 o 1, variabilità nulla) riporto null anche "
    "per df e se, benché calcolabili; METODI non dice cosa riportare.",
    "'Variabilità nulla fra pazienti' (solo per i tassi) è intesa come: tutte le medie per paziente "
    "uguali a θ̂ entro 1e-12, cioè V̂ = 0. Per l'eccesso di discordanza METODI non prevede questo "
    "caso: con V̂ = 0 l'intervallo prodotto è degenere [θ̂, θ̂] con se = 0.",
    "Per i tassi il campo 'se' è l'errore standard sulla scala originale, √V̂; l'errore standard "
    "sulla scala logit, √V̂ / (θ̂(1−θ̂)), è usato solo per costruire l'intervallo, poi riportato con "
    "la funzione logistica. METODI non dice quale dei due esporre.",
    "tr(M²) è calcolata come somma dei quadrati degli elementi di M (M è simmetrica).",
    "n_cells_with_trb = cellule dell'AnnData (tutti i compartimenti, tutte le etichette, anche "
    "mancanti) a cui è stata assegnata una TRB valida. n_reference_clones = cloni (paziente, "
    "cdr3_nt) con identità di riferimento definita (non solo con almeno 3 cellule).",
]


# ---------------------------------------------------------------------------
# Lettura dei dati
# ---------------------------------------------------------------------------

def normalizza_barcode(b) -> str:
    """Toglie spazi e un suffisso finale '-<cifre>'."""
    return re.sub(r"-\d+$", "", str(b).strip())


def _vero(serie: pd.Series) -> pd.Series:
    return serie.astype(str).str.strip().str.lower() == "true"


def scegli_contig_trb(contig: pd.DataFrame) -> pd.DataFrame:
    """Filtra i contig e ne tiene uno TRB per cellula (METODI §4, scelta del contig).

    `contig` ha le colonne di Cell Ranger più patient, compartment. Restituisce una riga per
    cellula con patient, compartment, barcode (normalizzato), cdr3_nt.
    """
    d = contig.copy()
    m = d["chain"].astype(str).str.strip() == CATENA
    if "is_cell" in d.columns:
        m &= _vero(d["is_cell"])
    m &= _vero(d["high_confidence"])
    m &= _vero(d["productive"])
    cdr3 = d["cdr3_nt"].astype(str).str.strip()
    m &= (cdr3 != "") & (cdr3.str.lower() != "none") & (cdr3.str.lower() != "nan")
    d = d.loc[m].copy()
    d["cdr3_nt"] = d["cdr3_nt"].astype(str).str.strip()
    d["barcode"] = d["barcode"].map(normalizza_barcode)
    d["umis"] = pd.to_numeric(d["umis"], errors="coerce").fillna(0)
    d["reads"] = pd.to_numeric(d["reads"], errors="coerce").fillna(0)
    d = d.sort_values(
        ["patient", "compartment", "barcode", "umis", "reads", "cdr3_nt"],
        ascending=[True, True, True, False, False, True],
        kind="mergesort",
    )
    d = d.drop_duplicates(["patient", "compartment", "barcode"], keep="first")
    return d[["patient", "compartment", "barcode", "cdr3_nt"]].reset_index(drop=True)


def leggi_vdj(manifest_path: str) -> pd.DataFrame:
    man = pd.read_csv(manifest_path, dtype=str, keep_default_na=False)
    for c in ("path", "patient", "compartment"):
        if c not in man.columns:
            raise ValueError(f"manifest senza colonna {c!r}")
    man["patient"] = man["patient"].str.strip()
    man["compartment"] = man["compartment"].str.strip()
    dup = man.duplicated(["patient", "compartment"])
    if dup.any():
        raise ValueError(f"coppie paziente-compartimento ripetute nel manifest: "
                         f"{man.loc[dup, ['patient', 'compartment']].values.tolist()}")
    base = os.path.dirname(os.path.abspath(manifest_path))
    parti = []
    for r in man.itertuples(index=False):
        f = pd.read_csv(os.path.join(base, r.path), dtype=str, keep_default_na=False)
        f["patient"] = r.patient
        f["compartment"] = r.compartment
        parti.append(f)
    if not parti:
        raise ValueError("manifest vuoto")
    return scegli_contig_trb(pd.concat(parti, ignore_index=True))


def positivita(X, var_names, marker_map: dict) -> dict:
    """Per ogni etichetta, vettore booleano: almeno un marcatore con conta > 0 (regola OR)."""
    nomi = np.asarray([str(v) for v in var_names])
    out = {}
    for lab, geni in marker_map.items():
        idx = np.flatnonzero(np.isin(nomi, list(geni)))
        if idx.size == 0:
            raise ValueError(f"nessun marcatore di {lab!r} presente in var_names: {geni}")
        sub = X[:, idx]
        if sparse.issparse(sub):
            pos = np.asarray((sub > 0).sum(axis=1)).ravel() > 0
        else:
            pos = (np.asarray(sub) > 0).any(axis=1)
        out[lab] = pos
    return out


def costruisci_celle(obs: pd.DataFrame, X, var_names, vdj: pd.DataFrame,
                     marker_map: dict) -> pd.DataFrame:
    """Tabella delle cellule con TRB valida: patient, compartment, barcode, label, cdr3_nt, pos_*."""
    o = pd.DataFrame({
        "riga": np.arange(obs.shape[0]),
        "patient": obs["patient"].astype(str).str.strip().values,
        "compartment": obs["compartment"].astype(str).str.strip().values,
        "barcode": obs["barcode"].map(normalizza_barcode).values,
    })
    lab = obs["label"]
    lab = lab.astype(object).where(lab.notna(), None)
    lab = lab.map(lambda v: None if v is None or str(v).strip() == "" else str(v).strip())
    o["label"] = lab.values
    dup = o.duplicated(["patient", "compartment", "barcode"], keep=False)
    if dup.any():
        raise ValueError(f"{int(dup.sum())} cellule con (patient, compartment, barcode) duplicato")
    pos = positivita(X, var_names, marker_map)
    for k, v in pos.items():
        o["pos_" + k] = v
    celle = o.merge(vdj, on=["patient", "compartment", "barcode"], how="inner")
    return celle


# ---------------------------------------------------------------------------
# Statistiche
# ---------------------------------------------------------------------------

def d_entro(conteggi) -> float:
    """D_A = 1 − Σ n_k(n_k − 1) / (n(n − 1)), senza reinserimento (METODI §3)."""
    c = np.asarray(list(conteggi), dtype=float)
    n = c.sum()
    return 1.0 - float(np.sum(c * (c - 1)) / (n * (n - 1)))


def d_tra(cnt_a: Counter, cnt_b: Counter) -> float:
    """D(A×B) = 1 − Σ_k p_{A,k} p_{B,k}."""
    na, nb = sum(cnt_a.values()), sum(cnt_b.values())
    return 1.0 - sum(cnt_a[k] / na * cnt_b.get(k, 0) / nb for k in cnt_a)


def eccessi(cnt_a: Counter, cnt_b: Counter) -> tuple[float, float]:
    """Restituisce (eccesso pesato, eccesso a pesi uguali) per una coppia clone-compartimenti."""
    na, nb = sum(cnt_a.values()), sum(cnt_b.values())
    da, db = d_entro(cnt_a.values()), d_entro(cnt_b.values())
    dab = d_tra(cnt_a, cnt_b)
    pesato = dab - (na * da + nb * db) / (na + nb)
    uguale = dab - (da + db) / 2.0
    return pesato, uguale


def _vuoto(n_groups=0, n_obs=0, mean=None) -> dict:
    return {"mean": mean, "df": None, "se": None, "ci_low": None, "ci_high": None,
            "n_groups": int(n_groups), "n_obs": int(n_obs)}


def bm_df(n: np.ndarray) -> float:
    """Gradi di libertà Bell-McCaffrey: ν = (tr M)² / tr(M²), M = A^½ K A^½ (METODI §2)."""
    n = np.asarray(n, dtype=float)
    N = n.sum()
    a = 1.0 / ((1.0 - n / N) * N ** 2)
    K = np.diag(n) - np.outer(n, n) / N
    sa = np.sqrt(a)
    M = sa[:, None] * K * sa[None, :]
    return float(np.trace(M) ** 2 / np.sum(M * M))


def cluster_t(valori, gruppi, logit: bool = False, min_pazienti: int = MIN_PAZIENTI) -> dict:
    """Intervallo t cluster-robusto CR2 con gradi di libertà Bell-McCaffrey (METODI §2)."""
    v = np.asarray(valori, dtype=float)
    g = np.asarray(gruppi, dtype=object)
    N = v.size
    if N == 0:
        return _vuoto()
    _, inv = np.unique(g.astype(str), return_inverse=True)
    S = np.bincount(inv, weights=v)
    n = np.bincount(inv).astype(float)
    G = S.size
    theta = float(S.sum() / N)
    if G < min_pazienti:
        return _vuoto(G, N, theta)
    e = S - theta * n
    V = float(np.sum(e ** 2 / ((1.0 - n / N) * N ** 2)))
    nu = bm_df(n)
    q = float(stats.t.ppf(0.975, nu))
    se = math.sqrt(V)
    out = {"mean": theta, "df": nu, "se": se, "ci_low": None, "ci_high": None,
           "n_groups": int(G), "n_obs": int(N)}
    if logit:
        nulla = bool(np.all(np.abs(S / n - theta) <= 1e-12))
        if theta <= 0.0 or theta >= 1.0 or nulla:
            return _vuoto(G, N, theta)
        se_l = se / (theta * (1.0 - theta))
        l0 = math.log(theta / (1.0 - theta))
        out["ci_low"] = 1.0 / (1.0 + math.exp(-(l0 - q * se_l)))
        out["ci_high"] = 1.0 / (1.0 + math.exp(-(l0 + q * se_l)))
    else:
        out["ci_low"] = theta - q * se
        out["ci_high"] = theta + q * se
    return out


# ---------------------------------------------------------------------------
# Modulo B
# ---------------------------------------------------------------------------

def identita_riferimento(celle: pd.DataFrame, rif: str, etichette: list) -> dict:
    """Identità di riferimento per clone (paziente, cdr3_nt), solo dal compartimento di riferimento."""
    ident = {}
    r = celle[celle["compartment"] == rif]
    for chiave, grp in r.groupby(["patient", "cdr3_nt"], sort=True):
        n = len(grp)
        if n < MIN_CELLULE_RIF:
            continue
        pos = sorted(((int(grp["pos_" + k].sum()), k) for k in etichette), reverse=True)
        p1, k1 = pos[0]
        p2 = pos[1][0] if len(pos) > 1 else 0
        # margine (p1 − p2)/n ≥ 1/5, in aritmetica intera
        if MARGINE_DEN * (p1 - p2) >= MARGINE_NUM * n and p1 > p2:
            ident[chiave] = k1
    return ident


def convenzione_clone(celle: pd.DataFrame, comp: str, ident: dict, tie_rule: str) -> dict:
    c = celle[(celle["compartment"] == comp) & celle["label"].notna()]
    valori, gruppi = [], []
    n_ties = 0
    for (paz, cdr3), grp in c.groupby(["patient", "cdr3_nt"], sort=True):
        if (paz, cdr3) not in ident or len(grp) < MIN_CELLULE_CLONE:
            continue
        cnt = Counter(grp["label"])
        massimo = max(cnt.values())
        testa = [k for k, v in cnt.items() if v == massimo]
        if len(testa) > 1:
            n_ties += 1
            if tie_rule == "exclude":
                continue
            err = 1 if tie_rule == "error" else 0
        else:
            err = int(testa[0] != ident[(paz, cdr3)])
        valori.append(err)
        gruppi.append(paz)
    rate = cluster_t(valori, gruppi, logit=True)
    return {"n_units": len(valori), "n_ties": n_ties, "n_errors": int(sum(valori)),
            "n_patients": len(set(gruppi)), "rate": rate}


def eccesso_coppia(celle: pd.DataFrame, comp_a: str, comp_b: str) -> dict:
    c = celle[celle["label"].notna()]
    pes, ugu, gruppi = [], [], []
    for (paz, cdr3), grp in c.groupby(["patient", "cdr3_nt"], sort=True):
        ca = Counter(grp.loc[grp["compartment"] == comp_a, "label"])
        cb = Counter(grp.loc[grp["compartment"] == comp_b, "label"])
        if sum(ca.values()) < MIN_CELLULE_LATO or sum(cb.values()) < MIN_CELLULE_LATO:
            continue
        w, u = eccessi(ca, cb)
        pes.append(w)
        ugu.append(u)
        gruppi.append(paz)
    return {"n_pairs": len(pes), "n_patients": len(set(gruppi)),
            "weighted": cluster_t(pes, gruppi), "equal": cluster_t(ugu, gruppi)}


def calcola(obs, X, var_names, vdj, marker_map, rif, tie_rule) -> dict:
    celle = costruisci_celle(obs, X, var_names, vdj, marker_map)
    etichette = list(marker_map.keys())
    ident = identita_riferimento(celle, rif, etichette)
    return {
        "clone_tumor": convenzione_clone(celle, COMPARTIMENTO_GIUDICATO, ident, tie_rule),
        "excess_PBMC_Tumor": eccesso_coppia(celle, rif, COMPARTIMENTO_GIUDICATO),
        "n_cells_with_trb": int(len(celle)),
        "n_reference_clones": int(len(ident)),
        "assunzioni": list(ASSUNZIONI),
    }


def esegui(args) -> dict:
    import anndata as ad
    a = ad.read_h5ad(args.h5ad)
    for c in ("patient", "compartment", "label", "barcode"):
        if c not in a.obs.columns:
            raise ValueError(f"obs senza colonna {c!r}")
    X = a.X
    if sparse.issparse(X):
        X = sparse.csr_matrix(X)
    with open(args.marker_map) as fh:
        marker_map = json.load(fh)
    vdj = leggi_vdj(args.vdj_manifest)
    return calcola(a.obs, X, a.var_names, vdj, marker_map,
                   args.reference_compartment, args.tie_rule)


def _json_default(o):
    if isinstance(o, (np.integer,)):
        return int(o)
    if isinstance(o, (np.floating,)):
        return float(o)
    raise TypeError(type(o))


# ---------------------------------------------------------------------------
# Autoverifica su dati sintetici
# ---------------------------------------------------------------------------

def _vicino(a, b, tol=1e-12):
    return a is not None and abs(a - b) <= tol


def selftest() -> None:
    ok = []

    def verifica(cond, nome):
        if not cond:
            raise AssertionError(nome)
        ok.append(nome)

    # 1. D_A senza reinserimento: etichette [a, a, b] -> 1 − 2/6 = 2/3
    verifica(_vicino(d_entro(Counter("aab").values()), 2 / 3), "D_A [a,a,b] = 2/3")
    verifica(_vicino(d_entro(Counter("aaaa").values()), 0.0), "D_A omogeneo = 0")
    verifica(_vicino(d_entro(Counter("ab").values()), 1.0), "D_A [a,b] = 1 (con reinserimento 0,5)")

    # 2. Eccesso: A = [a,a,b], B = [a,b,b,b]
    #    D_A = 2/3, D_B = 1 − 6/12 = 1/2, D(A×B) = 1 − (2/3·1/4 + 1/3·3/4) = 7/12
    #    pesato = 7/12 − (3·2/3 + 4·1/2)/7 = 7/12 − 4/7 = 1/84; uguale = 7/12 − 7/12 = 0
    A, B = Counter("aab"), Counter("abbb")
    verifica(_vicino(d_tra(A, B), 7 / 12), "D(A×B) = 7/12")
    w, u = eccessi(A, B)
    verifica(_vicino(w, 1 / 84), "eccesso pesato = 1/84")
    verifica(_vicino(u, 0.0), "eccesso a pesi uguali = 0")

    # 3. Stima puntuale = media di tutte le osservazioni (non media delle medie)
    val = [1, 0, 0, 0, 1, 1, 0, 1, 0, 0, 0, 1]
    grp = ["P1", "P1", "P1", "P1", "P2", "P3", "P3", "P4", "P5", "P5", "P5", "P5"]
    r = cluster_t(val, grp)
    verifica(_vicino(r["mean"], 5 / 12), "stima puntuale = 5/12 (media di tutte le osservazioni)")
    verifica(r["df"] < 4 - 1e-9, "dimensioni diverse: ν < G − 1")

    # 4. Pazienti di uguale dimensione: ν = G − 1; varianza a mano
    #    gruppi da 2: [0,1] [1,1] [0,0] [1,0] [1,1]; θ = 0,6; e = S − 1,2 = (−0,2, 0,8, −1,2, −0,2, 0,8)
    #    Σe² = 2,8; V = 2,8 / ((1 − 0,2)·100) = 0,035
    val = [0, 1, 1, 1, 0, 0, 1, 0, 1, 1]
    grp = ["A", "A", "B", "B", "C", "C", "D", "D", "E", "E"]
    r = cluster_t(val, grp)
    verifica(_vicino(r["df"], 4.0, 1e-9), "uguale dimensione: ν = G − 1 = 4")
    verifica(_vicino(r["se"], math.sqrt(0.035)), "V̂ CR2 = 0,035")
    q = stats.t.ppf(0.975, 4)
    verifica(_vicino(r["ci_low"], 0.6 - q * math.sqrt(0.035)), "intervallo scala originale")
    for G, m in ((7, 3), (12, 1), (5, 40)):
        verifica(_vicino(bm_df(np.full(G, m)), G - 1, 1e-9), f"ν = G − 1 per G={G}, n_g={m}")
    rl = cluster_t(val, grp, logit=True)
    sl = math.sqrt(0.035) / 0.24
    lo = 1 / (1 + math.exp(-(math.log(1.5) - q * sl)))
    verifica(_vicino(rl["ci_low"], lo), "intervallo su scala logit")
    verifica(cluster_t([0] * 10, grp, logit=True)["ci_low"] is None, "θ̂ = 0: nessun intervallo")
    verifica(cluster_t([0, 1] * 5, grp, logit=True)["ci_low"] is None,
             "variabilità nulla fra pazienti: nessun intervallo")
    verifica(cluster_t(val[:8], grp[:8])["ci_low"] is None, "4 pazienti: nessun intervallo")

    # 5. Prova completa su file sintetici (un paziente)
    _selftest_completo(verifica)

    for nome in ok:
        print("ok  ", nome)
    print(f"selftest superato: {len(ok)} controlli")


def _selftest_completo(verifica) -> None:
    import anndata as ad

    geni = ["CD4", "CD8A", "CD8B", "ACTB"]
    righe = []   # (compartment, barcode, label, cd4, cd8a, cd8b)
    # Clone X (AAA), PBMC 5 cellule: CD4+ in 3, CD8A+ in 1, CD8B+ in 1 -> 0,6 contro 0,4, margine
    # esattamente 0,20 -> identità CD4T (con il confronto in virgola mobile la perderebbe).
    for i, (c4, c8a, c8b) in enumerate([(2, 0, 0), (1, 0, 0), (3, 0, 0), (0, 1, 0), (0, 0, 1)]):
        righe.append(("PBMC", f"XP{i:014d}", "CD4T", c4, c8a, c8b, "AAA"))
    # Clone Y (CCC), PBMC 3 cellule CD8+ -> CD8T
    for i in range(3):
        righe.append(("PBMC", f"YP{i:014d}", "CD8T", 0, 2, 1, "CCC"))
    # Clone Z (GGG), PBMC 2 cellule -> nessun riferimento
    for i in range(2):
        righe.append(("PBMC", f"ZP{i:014d}", "CD4T", 1, 0, 0, "GGG"))
    # Tumore: X 4 cellule 2 CD4T + 2 CD8T -> pareggio
    for i, lab in enumerate(["CD4T", "CD4T", "CD8T", "CD8T"]):
        righe.append(("Tumor", f"XT{i:014d}", lab, 0, 0, 0, "AAA"))
    # Y 4 cellule: CD4T, CD4T, CD8T con CCC e una CD8T che per la scelta del contig diventa TTT
    for i, lab in enumerate(["CD4T", "CD4T", "CD8T", "CD8T"]):
        righe.append(("Tumor", f"YT{i:014d}", lab, 0, 0, 0, "CCC" if i < 3 else "TTT"))
    # Z 3 cellule CD4T
    for i in range(3):
        righe.append(("Tumor", f"ZT{i:014d}", "CD4T", 0, 0, 0, "GGG"))
    # cellula senza VDJ e cellula con sola TRA
    righe.append(("Tumor", "NOVDJ00000000000", "CD8T", 0, 0, 0, None))
    righe.append(("Tumor", "SOLOTRA000000000", "CD8T", 0, 0, 0, "TRA"))

    X = sparse.csr_matrix(np.array([[r[3], r[4], r[5], 5] for r in righe], dtype=np.int64))
    obs = pd.DataFrame({
        "patient": ["P1"] * len(righe),
        "compartment": [r[0] for r in righe],
        "label": [r[2] for r in righe],
        "barcode": [r[1] for r in righe],   # senza suffisso
    }, index=[f"c{i}" for i in range(len(righe))])
    adata = ad.AnnData(X=X, obs=obs, var=pd.DataFrame(index=geni))

    def contig(bc, chain, cdr3, umis, reads, prod="True", hc="True"):
        return {"barcode": bc + "-1", "is_cell": "True", "contig_id": bc + "_c",
                "high_confidence": hc, "length": "500", "chain": chain, "v_gene": "V",
                "d_gene": "", "j_gene": "J", "c_gene": "C", "full_length": "True",
                "productive": prod, "cdr3": "CASS", "cdr3_nt": cdr3, "reads": str(reads),
                "umis": str(umis), "raw_clonotype_id": "clonotype1", "raw_consensus_id": ""}

    vdj = {"PBMC": [], "Tumor": []}
    for comp, bc, lab, _, _, _, cl in righe:
        if cl is None:
            continue
        if cl == "TRA":
            vdj[comp].append(contig(bc, "TRA", "TTTTTT", 9, 99))
            continue
        if cl == "TTT":
            # due TRB a pari UMI: vince quella con più reads (TTT)
            vdj[comp].append(contig(bc, "TRB", "CCC", 5, 10))
            vdj[comp].append(contig(bc, "TRB", "TTT", 5, 20))
            continue
        vdj[comp].append(contig(bc, "TRB", cl, 3, 30))
        vdj[comp].append(contig(bc, "TRA", "GATGAT", 4, 40))
    # contig da ignorare: non produttivo e non high_confidence, con più UMI
    vdj["PBMC"].append(contig("XP00000000000000", "TRB", "AAAAAA", 50, 500, prod="False"))
    vdj["PBMC"].append(contig("XP00000000000001", "TRB", "AAAAAA", 50, 500, hc="False"))

    with tempfile.TemporaryDirectory() as tmp:
        h5 = os.path.join(tmp, "x.h5ad")
        adata.write_h5ad(h5)
        man = []
        for comp, rr in vdj.items():
            os.makedirs(os.path.join(tmp, comp))
            pd.DataFrame(rr).to_csv(os.path.join(tmp, comp, "filtered_contig_annotations.csv"),
                                    index=False)
            man.append({"path": f"{comp}/filtered_contig_annotations.csv",
                        "patient": "P1", "compartment": comp})
        pd.DataFrame(man).to_csv(os.path.join(tmp, "manifest.csv"), index=False)
        mm = os.path.join(tmp, "mm.json")
        with open(mm, "w") as fh:
            json.dump({"CD4T": ["CD4"], "CD8T": ["CD8A", "CD8B"]}, fh)

        ris = {}
        for tr in ("exclude", "error", "correct"):
            args = argparse.Namespace(h5ad=h5, vdj_manifest=os.path.join(tmp, "manifest.csv"),
                                      marker_map=mm, reference_compartment="PBMC", tie_rule=tr)
            ris[tr] = esegui(args)

    e = ris["exclude"]
    verifica(e["n_cells_with_trb"] == 21, "cellule con TRB = 21 (barcode con suffisso abbinati)")
    verifica(e["n_reference_clones"] == 2, "cloni con riferimento = 2 (margine 0,20 esatto incluso)")
    ct = e["clone_tumor"]
    verifica((ct["n_units"], ct["n_ties"], ct["n_errors"], ct["n_patients"]) == (1, 1, 1, 1),
             "per clone, exclude: 1 unità, 1 pareggio, 1 errore, 1 paziente")
    verifica(_vicino(ct["rate"]["mean"], 1.0) and ct["rate"]["ci_low"] is None,
             "per clone, exclude: tasso 1, nessun intervallo")
    c = ris["error"]["clone_tumor"]
    verifica((c["n_units"], c["n_errors"], c["n_ties"]) == (2, 2, 1), "tie-rule error: 2 unità, 2 errori")
    c = ris["correct"]["clone_tumor"]
    verifica((c["n_units"], c["n_errors"]) == (2, 1), "tie-rule correct: 2 unità, 1 errore")
    ex = e["excess_PBMC_Tumor"]
    # X: pesato 11/54, uguale 1/6; Y: 1/3 e 1/3; Z: 0 e 0
    verifica(ex["n_pairs"] == 3 and ex["n_patients"] == 1, "eccesso: 3 coppie, 1 paziente")
    verifica(_vicino(ex["weighted"]["mean"], 29 / 162), "eccesso pesato medio = 29/162")
    verifica(_vicino(ex["equal"]["mean"], 1 / 6), "eccesso a pesi uguali medio = 1/6")
    verifica(ex["weighted"]["ci_low"] is None, "eccesso con 1 paziente: nessun intervallo")


# ---------------------------------------------------------------------------

def main(argv=None) -> int:
    p = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    p.add_argument("--selftest", action="store_true")
    p.add_argument("--h5ad")
    p.add_argument("--vdj-manifest")
    p.add_argument("--marker-map")
    p.add_argument("--reference-compartment", default="PBMC")
    p.add_argument("--tie-rule", choices=["exclude", "error", "correct"], default="exclude")
    p.add_argument("--out")
    args = p.parse_args(argv)
    if args.selftest:
        selftest()
        return 0
    for nome in ("h5ad", "vdj_manifest", "marker_map", "out"):
        if getattr(args, nome) is None:
            p.error(f"--{nome.replace('_', '-')} è obbligatorio")
    ris = esegui(args)
    with open(args.out, "w") as fh:
        json.dump(ris, fh, indent=2, ensure_ascii=False, default=_json_default)
    return 0


if __name__ == "__main__":
    sys.exit(main())
