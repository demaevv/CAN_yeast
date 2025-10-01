#!/usr/bin/env bash

set -eo pipefail
module load miniconda 2>/dev/null || true
source "$HOME/miniconda3/etc/profile.d/conda.sh"
conda activate "$HOME/yeast_sv_project/env/ont_sv"

: "${SLURM_ARRAY_TASK_ID:?Run as job array, e.g. sbatch --array=1-9%1 scripts/02_polish.sbatch}"
SAMPLE="$(sed -n "${SLURM_ARRAY_TASK_ID}{s/\r$//;p;}" samples.txt)"

BASE="$HOME/yeast_sv_project"
RAW="$BASE/data/raw/${SAMPLE}_nanopore.fastq.gz"
CLEAN="$BASE/data/clean/${SAMPLE}.fastq.gz"
READS="$RAW"; [[ -s "$CLEAN" ]] && READS="$CLEAN"

ASM_DIR="$BASE/assemblies/${SAMPLE}"
ASM0="$ASM_DIR/flye/assembly.fasta"          
ASM1="$ASM_DIR/assembly.racon1.fasta"
ASM2="$ASM_DIR/assembly.racon2.fasta"
MEDAKA_DIR="$ASM_DIR/medaka"
ASM_FINAL="$ASM_DIR/assembly.final.fasta"

PAF1="$ASM_DIR/r1.paf.gz"
PAF2="$ASM_DIR/r2.paf.gz"

THREADS="${SLURM_CPUS_PER_TASK}"
KEEP_PAF="${KEEP_PAF:-0}"
MEDAKA_MODEL="${MEDAKA_MODEL:-r941_min_sup_g507}"

echo "[INFO] SAMPLE=${SAMPLE}"
echo "[INFO] Using reads: $READS"
[[ -f "$ASM0" ]] || { echo "[ERROR] Flye assembly not found: $ASM0"; exit 2; }
[[ -f "$READS" ]] || { echo "[ERROR] Reads not found: $READS"; exit 2; }

if [[ -s "$ASM_FINAL" ]]; then
  echo "[SKIP] Final assembly exists: $ASM_FINAL"; exit 0
fi

echo "[INFO] Threads=$THREADS  Medaka model=$MEDAKA_MODEL"

# Racon — Round 1
echo "[STEP] minimap2 -> PAF (ASM0)"
minimap2 -x map-ont -t "$THREADS" "$ASM0" "$READS" | gzip > "$PAF1"
echo "[STEP] Racon round 1"
racon -t "$THREADS" "$READS" "$PAF1" "$ASM0" > "$ASM1"

# Racon — Round 2
echo "[STEP] minimap2 -> PAF (ASM1)"
minimap2 -x map-ont -t "$THREADS" "$ASM1" "$READS" | gzip > "$PAF2"
echo "[STEP] Racon round 2"
racon -t "$THREADS" "$READS" "$PAF2" "$ASM1" > "$ASM2"

# Medaka
if command -v medaka_consensus >/dev/null 2>&1; then
  echo "[STEP] Medaka consensus (model=$MEDAKA_MODEL)"
  mkdir -p "$MEDAKA_DIR"
  medaka_consensus -i "$READS" -d "$ASM2" -o "$MEDAKA_DIR" -t "$THREADS" -m "$MEDAKA_MODEL" \
    || echo "[WARN] Medaka failed; falling back to Racon2"
  if [[ -s "$MEDAKA_DIR/consensus.fasta" ]]; then
    ln -sf "$MEDAKA_DIR/consensus.fasta" "$ASM_FINAL"
    echo "[OK] Final assembly (medaka): $ASM_FINAL"
  else
    cp "$ASM2" "$ASM_FINAL"
    echo "[WARN] Using Racon2 as final: $ASM_FINAL"
  fi
else
  echo "[WARN] medaka_consensus not found; using Racon2 as final"
  cp "$ASM2" "$ASM_FINAL"
fi

if [[ "$KEEP_PAF" != "1" ]]; then
  rm -f "$PAF1" "$PAF2" 2>/dev/null || true
fi

[[ -s "$ASM_FINAL" ]] && echo "[DONE] Polishing completed: $ASM_FINAL" || { echo "[FAIL] No final assembly"; exit 3; }