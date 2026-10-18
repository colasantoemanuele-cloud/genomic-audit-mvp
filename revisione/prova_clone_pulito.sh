#!/bin/bash
# Prova di installazione su clone pulito del ramo revisione e suite completa con i dati reali.
# Uso: revisione/prova_clone_pulito.sh <cartella vuota di lavoro> <cartella dei log>
# I cloni stanno in una cartella senza pdac-ml accanto e con GSE278694_DIR non impostata:
# i test sui dati reali vengono saltati. La suite con i dati gira nella copia di lavoro.
set -u
REPO=$(cd "$(dirname "$0")/.." && pwd)
W=${1:?cartella di lavoro}; L=${2:?cartella dei log}
mkdir -p "$W" "$L"; L=$(cd "$L" && pwd)
unset GSE278694_DIR
BR=$(git -C "$REPO" rev-parse --abbrev-ref HEAD); SHA=$(git -C "$REPO" rev-parse --short HEAD)
log() { echo "[$(date '+%F %T')] $*"; }

prova_uv() {
  local d="$W/clone_uv"; rm -rf "$d"; git clone -q --branch "$BR" "$REPO" "$d"; cd "$d"
  { log "uv, commit $SHA"; uv --version; uv sync --extra dev 2>&1 | tail -3
    log "demo"; /usr/bin/time -f "demo: %e s, %M KB" uv run audit-sc demo --out results/demo_report.html > demo.txt 2>&1; tail -3 demo.txt
    log "suite completa senza i dati reali"
    uv run pytest -n auto -q -p no:cacheprovider -rsx 2>&1 | tail -15; } > "$L/clone_uv.txt" 2>&1
}
prova_pip() {
  local d="$W/clone_pip"; rm -rf "$d"; git clone -q --branch "$BR" "$REPO" "$d"; cd "$d"
  { log "pip con python3.13, commit $SHA"; python3.13 -m venv .venv && . .venv/bin/activate && python --version
    pip install -q -e ".[dev]" 2>&1 | tail -3
    log "demo"; audit-sc demo --out results/demo_report.html > demo.txt 2>&1; echo "demo exit $?"; tail -2 demo.txt
    log "pytest -m 'not slow'"; pytest -n auto -q -p no:cacheprovider -m "not slow" -rsx 2>&1 | tail -8; } > "$L/clone_pip.txt" 2>&1
}
prova_py312() {
  local d="$W/clone_py312"; rm -rf "$d"; git clone -q --branch "$BR" "$REPO" "$d"; cd "$d"
  # il venv di sistema di Python 3.12 non ha ensurepip: ambiente creato con uv sull'interprete 3.12
  { log "python3.12 (uv venv + uv pip), commit $SHA"; uv venv -q -p /usr/bin/python3.12 .venv && . .venv/bin/activate && python --version
    uv pip install -q -e ".[dev]" 2>&1 | tail -3
    log "pytest -m 'not slow'"; pytest -n auto -q -p no:cacheprovider -m "not slow" -rsx 2>&1 | tail -8; } > "$L/clone_py312.txt" 2>&1
}
prova_dati() {
  cd "$REPO"
  { log "copia di lavoro con i dati reali, commit $SHA"; uv run pytest -n auto -q -p no:cacheprovider -rsx 2>&1 | tail -15; } \
    > "$L/suite_con_dati.txt" 2>&1
}
prova_uv & prova_pip & prova_py312 & wait
prova_dati
log "fine" >> "$L/clone_uv.txt"
