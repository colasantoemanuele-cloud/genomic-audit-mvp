#!/bin/bash
# Riesegue le validazioni con il codice attuale. Uso: revisione/esegui_validazioni.sh [gse278694|esterna|tutto]
# Gli output pesanti (report HTML, copia .h5ad con i flag) vanno in $OUT, fuori dalla repo.
set -u
cd "$(dirname "$0")/.."
PART=${1:-tutto}
OUT=${OUT:-/tmp/genomic_audit_validazioni}
mkdir -p "$OUT"
PY="uv run python -W ignore"
PDAC=${GSE278694_DIR:-../pdac-ml}
H5AD=$PDAC/data/interim/sc_raw.h5ad
R=validation/results
log() { echo "[$(date +%H:%M:%S)] $*"; }

if [[ $PART == gse278694 || $PART == tutto ]]; then
  log "GSE278694: manifest VDJ"
  $PY -m validation.make_vdj_manifest --out $R/vdj_manifest_gse278694.csv
  echo '{"CD4T": ["CD4"], "CD8T": ["CD8A", "CD8B"]}' > "$OUT/markers.json"
  TCR="tcr --h5ad $H5AD --vdj-manifest $R/vdj_manifest_gse278694.csv --patient-col patients --compartment-col tissue
       --celltype-col all_celltype --barcode-col barcode --marker-map $OUT/markers.json --reference-compartment PBMC
       --compartment-order Tumor,Adjacent_normal"
  log "GSE278694: Modulo B, convenzioni della tesi"
  $PY cli.py $TCR --clone-error-labels NK --clone-marker-priority CD8T,CD4T --out "$OUT/tcr_tesi.html" \
      > $R/cli_tcr_real_conventions_tesi.txt 2>&1
  log "GSE278694: Modulo B, convenzioni di default"
  $PY cli.py $TCR --out "$OUT/tcr_default.html" > $R/cli_tcr_real_conventions_default.txt 2>&1
  log "GSE278694: Modulo B con export dei flag e frazione di CD8 sperimentale"
  /usr/bin/time -v $PY cli.py $TCR --export-flags --sperimentale-cd8 Tumor --out "$OUT/tcr_real.html" \
      > $R/cli_tcr_real.txt 2> "$OUT/tcr_real_time.txt"
  log "GSE278694: verifica di fedeltà del porting (scala delle definizioni)"
  $PY -m validation.gse278694_moduleB --vdj-manifest $R/vdj_manifest_gse278694.csv --out $R > "$OUT/moduleB.log" 2>&1
  log "GSE278694: quota delle etichette NK nella differenza fra convenzioni (B2)"
  $PY -m validation.gse278694_b2_nk --vdj-manifest $R/vdj_manifest_gse278694.csv --out $R > "$OUT/b2_nk.log" 2>&1
  log "GSE278694: scelte che dipendevano dall ordine delle righe (A8)"
  $PY -m validation.gse278694_pareggi --vdj-manifest $R/vdj_manifest_gse278694.csv --out $R > "$OUT/pareggi.log" 2>&1
  log "GSE278694: diagnostica della frazione di CD8"
  $PY -m validation.gse278694_cd8_diagnostic --vdj-manifest $R/vdj_manifest_gse278694.csv --out $R \
      > $R/pulizia4_cd8_diagnostic.txt 2>&1
  log "GSE278694: copertura dei flag e frazione di CD8 dai flag esportati"
  $PY -m validation.gse278694_cd8_flags --flags "$OUT/sc_raw_audit_flags.csv" --out $R > "$OUT/cd8_flags.log" 2>&1
  log "GSE278694: audit del disegno (metadati sc + sn, e solo sc da adata.obs)"
  $PY cli.py design --meta $R/gse278694_real_sheet.csv --patient-col patients --tissue-col tissue \
      --technical library=library --technical protocol=protocol \
      --compare tissue:Tumor:Adjacent_normal --compare protocol:scRNA:snRNA --out "$OUT/design_real.html" \
      > $R/cli_design_real.txt 2>&1
  $PY cli.py design --h5ad $H5AD --patient-col patients --tissue-col tissue --technical library=library \
      --compare tissue:Tumor:Adjacent_normal --compare tissue:Tumor:PBMC --out "$OUT/design_real_sc.html" \
      > $R/cli_design_real_sc.txt 2>&1
  sed -i "s#$OUT#<out>#g" $R/cli_*.txt
fi

if [[ $PART == esterna || $PART == tutto ]]; then
  E=validation_esterna/results
  log "Esterna: audit del disegno su tre dataset (D1-D4)"
  $PY -m validation_esterna.run_design > "$OUT/run_design.log" 2>&1
  log "Esterna: Modulo A su GSE125449 completo (base, permutazioni, formati)"
  $PY -m validation_esterna.run_moduleA base > $E/moduleA_base.txt 2>&1
  $PY -m validation_esterna.run_moduleA perm > $E/moduleA_perm.txt 2>&1
  $PY -m validation_esterna.run_moduleA formati > $E/moduleA_formati.txt 2>&1
  log "Esterna: A1 con la CLI su GSE125449 completo"
  $PY cli.py leakage --h5ad validation_esterna/data/GSE125449.h5ad --target-col Type --patient-col patient \
      --out "$OUT/leakage_GSE125449.html" > $E/A1_cli_stdout.txt 2> $E/A1_cli_stderr.txt
  log "Esterna: coorte ridotta (criteri originali e modifica M1)"
  for step in base perm formati; do
    $PY -m validation_esterna.run_moduleA $step data/demo/GSE125449_demo.h5ad _ridotta > $E/moduleA_${step}_ridotta.txt 2>&1
  done
  $PY cli.py leakage --h5ad data/demo/GSE125449_demo.h5ad --target-col Type --patient-col patient \
      --out "$OUT/leakage_ridotta.html" > $E/A1_ridotta_cli.txt 2> $E/A1_ridotta_cli.err; echo $? > $E/A1_ridotta_exit.txt
  log "Esterna: tempi del Modulo A (GSE125449 completo e 25.000 cellule sintetiche)"
  /usr/bin/time -v $PY -m validation_esterna.benchmark_moduleA reale > $E/benchmark_reale.log 2> $E/benchmark_reale.err
  /usr/bin/time -v $PY -m validation_esterna.benchmark_moduleA sintetico > $E/benchmark_sintetico.log 2> $E/benchmark_sintetico.err
  log "Esterna: verdetti con i criteri originali e modificati"
  $PY -m validation_esterna.criteri_originali > "$OUT/criteri.log" 2>&1
fi
log "fine"
