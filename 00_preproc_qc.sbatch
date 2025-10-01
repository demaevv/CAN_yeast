#!/usr/bin/env bash

set -eo pipefail
module load miniconda 2>/dev/null || true
source "$HOME/miniconda3/etc/profile.d/conda.sh"
conda activate "$HOME/yeast_sv_project/env/ont_sv"

: "${SLURM_ARRAY_TASK_ID:?Run as job array, e.g. sbatch --array=1-9%2 scripts/00_preproc_qc.sbatch}"

SAMPLE="$(sed -n "${SLURM_ARRAY_TASK_ID}{s/\r$//;p;}" samples.txt)"

RAW="$HOME/yeast_sv_project/data/raw/${SAMPLE}_nanopore.fastq.gz"
CLEAN="$HOME/yeast_sv_project/data/clean/${SAMPLE}.fastq.gz"
READS="$RAW"
[[ -s "$CLEAN" ]] && READS="$CLEAN"

echo "[INFO] SLURM_ARRAY_TASK_ID=${SLURM_ARRAY_TASK_ID} SAMPLE=${SAMPLE}"
echo "[INFO] Using input: $READS"

if [[ ! -f "$READS" ]]; then
  echo "[ERROR] FASTQ not found: $READS"
  echo "[HINT] Available files in data/raw:"
  ls -1 "$HOME/yeast_sv_project/data/raw/" || true
  exit 1
fi

OUTDIR=~/yeast_sv_project/data/clean/${SAMPLE}
QC_RAW=~/yeast_sv_project/plots/nanoplot_raw/${SAMPLE}
QC_FILT=~/yeast_sv_project/plots/nanoplot_filt/${SAMPLE}

mkdir -p "$OUTDIR" "$QC_RAW" "$QC_FILT"

NanoPlot --fastq "$RAW" -t ${SLURM_CPUS_PER_TASK} -o "$QC_RAW" --raw

IN="$RAW"
TRIM="$OUTDIR/${SAMPLE}.trim.fastq.gz"
if command -v porechop_abi >/dev/null 2>&1; then
  porechop_abi -i "$IN" -o "$TRIM" -t ${SLURM_CPUS_PER_TASK}
  IN="$TRIM"
elif command -v porechop >/dev/null 2>&1; then
  porechop -i "$IN" -o "$TRIM" -t ${SLURM_CPUS_PER_TASK}
  IN="$TRIM"
fi

FILT="$OUTDIR/${SAMPLE}.filt.fastq.gz"
chopper -q 7 -l 1000 -i "$IN" | gzip > "$FILT"

# gunzip -c "$IN" | NanoFilt -q 7 --length 1000 | gzip > "$FILT"

NanoPlot --fastq "$FILT" -t ${SLURM_CPUS_PER_TASK} -o "$QC_FILT" --raw

ln -sf "$FILT" ~/yeast_sv_project/data/clean/${SAMPLE}.fastq.gz