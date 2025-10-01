#!/usr/bin/env bash

set -eo pipefail
module load miniconda 2>/dev/null || true
source "$HOME/miniconda3/etc/profile.d/conda.sh"
conda activate "$HOME/yeast_sv_project/env/ont_sv"

: "${SLURM_ARRAY_TASK_ID:?Run as job array, e.g. sbatch --array=1-9%1 scripts/01_assemble_flye.sbatch}"

SAMPLE="$(sed -n "${SLURM_ARRAY_TASK_ID}{s/\r$//;p;}" samples.txt)"

RAW="$HOME/yeast_sv_project/data/raw/${SAMPLE}_nanopore.fastq.gz"
CLEAN="$HOME/yeast_sv_project/data/clean/${SAMPLE}.fastq.gz"
READS="$RAW"
[[ -s "$CLEAN" ]] && READS="$CLEAN"

mkdir -p "$HOME/yeast_sv_project/logs" || true
OUTDIR="$HOME/yeast_sv_project/assemblies/${SAMPLE}"
FLYE_OUT="$OUTDIR/flye"
mkdir -p "$FLYE_OUT"

flye \
  --nano-raw "$READS" \
  --genome-size 12m \
  --threads "${SLURM_CPUS_PER_TASK}" \
  --out-dir "$FLYE_OUT"