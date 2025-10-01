#!/usr/bin/env bash

set -eo pipefail
module load miniconda 2>/dev/null || true
source "$HOME/miniconda3/etc/profile.d/conda.sh"
conda activate "$HOME/yeast_sv_project/env/ont_sv"
export MPLBACKEND=Agg

: "${SLURM_ARRAY_TASK_ID:?Run as job array (e.g. --array=1-9%3)}"
SAMPLE="$(sed -n "${SLURM_ARRAY_TASK_ID}{s/\r$//;p;}" $HOME/yeast_sv_project/samples.txt)"

BASE="$HOME/yeast_sv_project"
REF_FA="$BASE/ref/GCF_000146045.2_R64_genomic.fna"
REF_GFF="$BASE/ref/genomic.gff"
OUTDIR="$BASE/sv/${SAMPLE}"
MERGE_DIR="$OUTDIR/merge"
mkdir -p "$MERGE_DIR"

SNIFF="$OUTDIR/${SAMPLE}.sniffles.vcf.gz"
SVIM="$OUTDIR/svim/variants.vcf.gz"
CUTE="$OUTDIR/${SAMPLE}.cutesv.vcf.gz"

FOUND=()
for v in "$SNIFF" "$SVIM" "$CUTE"; do
  [[ -s "$v" ]] && FOUND+=("$v")
done
if [[ ${#FOUND[@]} -lt 1 ]]; then
  echo "[ERROR] No input VCFs found for $SAMPLE in $OUTDIR"
  exit 2
fi

if ! command -v SURVIVOR >/dev/null 2>&1; then
  echo "[ERROR] SURVIVOR not found. Install with: mamba install -y -c bioconda survivor"
  exit 3
fi

[[ -s "${REF_FA}.fai" ]] || samtools faidx "$REF_FA"

VCFLIST="$MERGE_DIR/${SAMPLE}.vcf.list"
printf "%s\n" "${FOUND[@]}" > "$VCFLIST"

fix_vcf() {
  local v="$1"
  if [[ ! -s "${v}.tbi" && ! -s "${v}.csi" ]]; then
    echo "[INFO] Sorting+indexing $v"
    local tmp="$v.tmp.gz"
    bcftools sort "$v" -Oz -o "$tmp"
    tabix -f -p vcf "$tmp"
    mv -f "$tmp" "$v"
  fi
}
for v in "${FOUND[@]}"; do fix_vcf "$v"; done

MAX_DIST="${MAX_DIST:-100}"     
MIN_SUPP="${MIN_SUPP:-2}"       
REQ_TYPE="${REQ_TYPE:-1}"       
REQ_STRAND="${REQ_STRAND:-1}"   
DIST_TYPE_REL="${DIST_TYPE_REL:-0}"  
MIN_SIZE="${MIN_SIZE:-50}"      

if [[ ${#FOUND[@]} -lt 2 ]]; then
  echo "[WARN] Only ${#FOUND[@]} VCF present; setting MIN_SUPP=1"
  MIN_SUPP=1
fi

MERGED_VCF="$MERGE_DIR/${SAMPLE}.merged.survivor.vcf"
MERGED_VCFGZ="${MERGED_VCF}.gz"

if [[ ! -s "$MERGED_VCFGZ" ]]; then
  echo "[STEP] SURVIVOR merge (${#FOUND[@]} files; MIN_SUPP=$MIN_SUPP; MAX_DIST=$MAX_DIST)"
  SURVIVOR merge "$VCFLIST" "$MAX_DIST" "$MIN_SUPP" "$REQ_TYPE" "$REQ_STRAND" "$DIST_TYPE_REL" "$MIN_SIZE" "$MERGED_VCF"
  bgzip -f "$MERGED_VCF"
  tabix -f -p vcf "$MERGED_VCFGZ"
else
  echo "[SKIP] Merged VCF exists: $MERGED_VCFGZ"
fi

read -r CAN1_CHR CAN1_START CAN1_END < <(
  awk -F'\t' '
    BEGIN{IGNORECASE=1}
    $0!~/^#/ && $3=="gene" && \
    ($9 ~ /(^|;)gene=CAN1(;|$)/ || $9 ~ /(^|;)(Name|gene_name)=CAN1(;|$)/ || $9 ~ /(^|;)locus_tag=YEL063C(;|$)/) {
      print $1, $4, $5; exit
    }' "$REF_GFF"
) || { echo "[ERROR] CAN1 gene not found in $REF_GFF"; exit 5; }

CHR_LEN=$(awk -v c="$CAN1_CHR" '$1==c{print $2}' "${REF_FA}.fai")
[[ -n "$CHR_LEN" ]] || { echo "[ERROR] Chromosome $CAN1_CHR not found in ${REF_FA}.fai"; exit 6; }

CAN1_FLANK="${CAN1_FLANK:-50000}"
REG_START=$(( CAN1_START - CAN1_FLANK ))
REG_END=$(( CAN1_END + CAN1_FLANK ))
(( REG_START < 1 )) && REG_START=1
(( REG_END > CHR_LEN )) && REG_END=$CHR_LEN
REGION="${CAN1_CHR}:${REG_START}-${REG_END}"
echo "[INFO] CAN1 window: $REGION  (flank=${CAN1_FLANK}bp)"

MERGED_CAN1_VCF="$MERGE_DIR/${SAMPLE}.merged.CAN1_${CAN1_FLANK}.vcf.gz"
if [[ ! -s "$MERGED_CAN1_VCF" ]]; then
  bcftools view -r "$REGION" -Oz -o "$MERGED_CAN1_VCF" "$MERGED_VCFGZ"
  tabix -f -p vcf "$MERGED_CAN1_VCF"
fi

CAN1_BED="$MERGE_DIR/CAN1_gene.bed"
awk -F'\t' '
  BEGIN{IGNORECASE=1}
  $0!~/^#/ && $3=="gene" && \
  ($9 ~ /(^|;)gene=CAN1(;|$)/ || $9 ~ /(^|;)(Name|gene_name)=CAN1(;|$)/ || $9 ~ /(^|;)locus_tag=YEL063C(;|$)/) {
    printf("%s\t%d\t%d\n", $1, $4-1, $5); exit
  }' "$REF_GFF" > "$CAN1_BED"

VCF_BED="$MERGE_DIR/${SAMPLE}.merged.CAN1_${CAN1_FLANK}.bed"
bcftools query -f '%CHROM\t%POS\t%INFO/END\t%ID\t%SVTYPE\t%SVLEN\n' "$MERGED_CAN1_VCF" \
  | awk 'BEGIN{OFS="\t"}{
      chr=$1; start=$2; end=$3; id=$4; t=$5; len=$6;
      if(end=="" || end=="." || end==0){ end=start+1 }   
      if(end<start){ tmp=start; start=end; end=tmp }
      print chr, start-1, end, id, t, len
    }' > "$VCF_BED"

OVR_IDS="$MERGE_DIR/${SAMPLE}.CAN1_overlap.ids.txt"
bedtools intersect -u -a "$VCF_BED" -b "$CAN1_BED" 2>/dev/null | cut -f4 > "$OVR_IDS" || true

TSV="$MERGE_DIR/${SAMPLE}.CAN1_${CAN1_FLANK}.summary.tsv"
{
  echo -e "SAMPLE\tCHROM\tPOS\tEND\tSVTYPE\tSVLEN\tQUAL\tSUPP\tSUPP_VEC\tOVERLAPS_CAN1"
  bcftools query -f "%CHROM\t%POS\t%INFO/END\t%SVTYPE\t%SVLEN\t%QUAL\t%INFO/SUPP\t%INFO/SUPP_VEC\t%ID\n" "$MERGED_CAN1_VCF" \
  | awk -v s="$SAMPLE" -v ids="$OVR_IDS" 'BEGIN{
      OFS="\t"; while((getline < ids)>0){hit[$1]=1} close(ids)
    }{
      chr=$1; pos=$2; end=$3; t=$4; len=$5; q=$6; supp=$7; vec=$8; id=$9;
      if(end=="" || end=="." || end==0){ end=pos }   # INS
      print s, chr, pos, end, t, len, q, supp, vec, (hit[id]?"YES":"NO")
    }'
} > "$TSV"

echo "[DONE] Merged genome-wide: $MERGED_VCFGZ"
echo "[DONE] Merged CAN1 region: $MERGED_CAN1_VCF"
echo "[DONE] CAN1 summary TSV:   $TSV"