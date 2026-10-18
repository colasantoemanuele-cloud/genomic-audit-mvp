#!/bin/bash
cd "$(dirname "$0")/../.."
start=$(date +%s)
uv run pytest -n auto -q -rfEsxX --durations=25 --junitxml=revisione/log/baseline_xdist.xml > revisione/log/baseline_xdist.txt 2>&1
echo "TEMPO_S $(( $(date +%s) - start ))" >> revisione/log/baseline_xdist.txt
start=$(date +%s)
uv run pytest -q -rfEsxX --durations=25 --junitxml=revisione/log/baseline_serial.xml > revisione/log/baseline_serial.txt 2>&1
echo "TEMPO_S $(( $(date +%s) - start ))" >> revisione/log/baseline_serial.txt
