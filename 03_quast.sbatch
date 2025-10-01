#!/usr/bin/env bash

set -eo pipefail
module load miniconda 2>/dev/null || true
source "$HOME/miniconda3/etc/profile.d/conda.sh"
conda activate "$HOME/yeast_sv_project/env/ont_sv"

: "${SLURM_ARRAY_TASK_ID:?Run as job array, e.g. sbatch --array=1-9%4 scripts/03_quast.sbatch}"

SAMPLE="$(sed -n "${SLURM_ARRAY_TASK_ID}{s/\r$//;p;}" samples.txt)"

BASE="$HOME/yeast_sv_project"
ASM="$BASE/assemblies/${SAMPLE}/assembly.final.fasta"
REF="$BASE/ref/GCF_000146045.2_R64_genomic.fna"
GFF="$BASE/ref/genomic.gff"
OUT="$BASE/assemblies/${SAMPLE}/quast"

QUAST_BIN="$(command -v quast.py || command -v quast || true)"
: "${QUAST_BIN:?QUAST not found in PATH. Install quast in the env.}"

echo "[INFO] SAMPLE=${SAMPLE}"
echo "[INFO] Assembly: $ASM"
echo "[INFO] Ref: $REF"
echo "[INFO] GFF: $GFF"
echo "[INFO] QUAST: $QUAST_BIN"

[[ -s "$ASM" ]] || { echo "[ERROR] Final assembly not found: $ASM"; exit 2; }
[[ -s "$REF" ]] || { echo "[ERROR] Reference FASTA not found: $REF"; exit 2; }
[[ -s "$GFF" ]] || { echo "[WARN] GFF not found: $GFF (will run without genes)"; GFF=""; }

mkdir -p "$OUT"

if [[ -s "$OUT/report.tsv" || -s "$OUT/report.txt" || -s "$OUT/report.html" ]]; then
  echo "[SKIP] QUAST report already exists in: $OUT"
  exit 0
fi

QUAST_OPTS=(
  --threads "${SLURM_CPUS_PER_TASK}"
  -o "$OUT"
  -r "$REF"
  --est-ref-size 12000000
)
[[ -n "$GFF" ]] && QUAST_OPTS+=( -g "$GFF" )

echo "[STEP] Running QUAST..."
"$QUAST_BIN" "${QUAST_OPTS[@]}" "$ASM"

if [[ -s "$OUT/report.tsv" || -s "$OUT/report.txt" || -s "$OUT/report.html" ]]; then
  echo "[DONE] QUAST reports ready in: $OUT"
else
  echo "[FAIL] QUAST did not produce reports in: $OUT"
  exit 3
fi