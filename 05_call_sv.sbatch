#!/usr/bin/env bash

set -eo pipefail
module load miniconda 2>/dev/null || true
source "$HOME/miniconda3/etc/profile.d/conda.sh"
conda activate "$HOME/yeast_sv_project/env/ont_sv"
export MPLBACKEND=Agg

: "${SLURM_ARRAY_TASK_ID:?Array required}"
SAMPLE="$(sed -n "${SLURM_ARRAY_TASK_ID}{s/\r$//;p;}" $HOME/yeast_sv_project/samples.txt)"

BASE="$HOME/yeast_sv_project"
REF="$BASE/ref/GCF_000146045.2_R64_genomic.fna"
BAM="$BASE/align/${SAMPLE}/${SAMPLE}.bam"
OUT="$BASE/sv/${SAMPLE}"
mkdir -p "$OUT" "$OUT/cutesv_work" "$OUT/svim"

[[ -s "$REF" ]] || { echo "[ERR] No ref: $REF"; exit 2; }
[[ -s "$BAM" ]] || { echo "[ERR] No BAM: $BAM"; exit 2; }

THREADS="${SLURM_CPUS_PER_TASK}"

# --- Sniffles2 ---
SNIFF_VCF="$OUT/${SAMPLE}.sniffles.vcf.gz"
if [[ -s "$SNIFF_VCF" ]]; then
  echo "[SKIP] Sniffles exists: $SNIFF_VCF"
else
  sniffles --input "$BAM" --vcf "$SNIFF_VCF" --reference "$REF" \
           --minsupport 10 --minsvlen 50 --mapq 20 --threads "$THREADS"
  tabix -f -p vcf "$SNIFF_VCF" || true
fi

# --- SVIM ---
SVIM_DIR="$OUT/svim"
SVIM_VCF="$SVIM_DIR/variants.vcf.gz"
if [[ -s "$SVIM_VCF" ]]; then
  echo "[SKIP] SVIM exists: $SVIM_VCF"
else
  set +e
  svim alignment "$SVIM_DIR" "$BAM" "$REF" --min_mapq 20 --min_sv_size 50 --max_sv_size 1000000
  rc=$?
  set -e
  if [[ -s "$SVIM_DIR/variants.vcf" ]]; then
    bcftools sort "$SVIM_DIR/variants.vcf" -Oz -o "$SVIM_VCF"
    tabix -f -p vcf "$SVIM_VCF"
  else
    echo "[WARN] SVIM produced no variants.vcf (exit=$rc)"
  fi
fi

# --- cuteSV ---
command -v cuteSV >/dev/null 2>&1 || { echo "[ERR] cuteSV not found in PATH"; exit 3; }
CUTESV_VCF="$OUT/${SAMPLE}.cutesv.vcf.gz"
if [[ -s "$CUTESV_VCF" ]]; then
  echo "[SKIP] cuteSV exists: $CUTESV_VCF"
else
  cuteSV "$BAM" "$REF" "$OUT/${SAMPLE}.cutesv.vcf" "$OUT/cutesv_work" \
        -t "$THREADS" -s 10 -l 50 -L 1000000 -q 20 --report_readid -S "$SAMPLE"
  bgzip -f "$OUT/${SAMPLE}.cutesv.vcf"
  tabix -f -p vcf "$CUTESV_VCF"
fi

echo "[DONE] SV calling for $SAMPLE"