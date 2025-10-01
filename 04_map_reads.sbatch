#!/usr/bin/env bash

set -eo pipefail
module load miniconda 2>/dev/null || true
source "$HOME/miniconda3/etc/profile.d/conda.sh"
conda activate "$HOME/yeast_sv_project/env/ont_sv"

: "${SLURM_ARRAY_TASK_ID:?Run as job array, e.g. --array=1-9%3}"

SAMPLE="$(sed -n "${SLURM_ARRAY_TASK_ID}{s/\r$//;p;}" samples.txt)"

BASE="$HOME/yeast_sv_project"
REF="$BASE/ref/GCF_000146045.2_R64_genomic.fna"
READS_RAW="$BASE/data/raw/${SAMPLE}_nanopore.fastq.gz"
READS_CLEAN="$BASE/data/clean/${SAMPLE}.fastq.gz"
READS="$READS_RAW"; [[ -s "$READS_CLEAN" ]] && READS="$READS_CLEAN"

OUTDIR="$BASE/align/${SAMPLE}"
mkdir -p "$OUTDIR"

BAM="$OUTDIR/${SAMPLE}.bam"
BAI="$BAM.bai"
PFX="$OUTDIR/${SAMPLE}"   

echo "[INFO] SAMPLE=$SAMPLE"
echo "[INFO] REF=$REF"
echo "[INFO] READS=$READS"
[[ -s "$REF" ]]   || { echo "[ERROR] Missing ref: $REF"; exit 2; }
[[ -s "$READS" ]] || { echo "[ERROR] Missing reads: $READS"; exit 2; }

if [[ -s "$BAM" && -s "$BAI" ]]; then
  echo "[SKIP] BAM already exists: $BAM"
else
  echo "[STEP] Mapping with minimap2 -> SAM | sort -> BAM"

  minimap2 -t "${SLURM_CPUS_PER_TASK}" -ax map-ont -Y --MD --secondary=no \
    -R "@RG\tID:${SAMPLE}\tSM:${SAMPLE}\tPL:ONT" \
    "$REF" "$READS" \
    | samtools sort -@ 8 -o "$BAM" -
  samtools index "$BAM"
fi

if [[ ! -s "$OUTDIR/flagstat.txt" ]]; then
  samtools flagstat   "$BAM" > "$OUTDIR/flagstat.txt"
fi
if [[ ! -s "$OUTDIR/stats.txt" ]]; then
  samtools stats      "$BAM" > "$OUTDIR/stats.txt"
fi

if [[ ! -s "${PFX}.regions.bed.gz" ]]; then
  mosdepth -t "${SLURM_CPUS_PER_TASK}" --by 1000 --mapq 20 --fast-mode "$PFX" "$BAM"
fi
if [[ ! -s "${PFX}.per-base.bed.gz" ]]; then
  mosdepth -t "${SLURM_CPUS_PER_TASK}" --mapq 20 --fast-mode --no-per-base 0 "$PFX" "$BAM" 2>/dev/null || true
fi
echo "[DONE] BAM+QC ready: $BAM"
