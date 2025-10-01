#!/usr/bin/env bash

set -eo pipefail
module load miniconda 2>/dev/null || true
source "$HOME/miniconda3/etc/profile.d/conda.sh"
conda activate "$HOME/yeast_sv_project/env/ont_sv"
export MPLBACKEND=Agg

: "${SLURM_ARRAY_TASK_ID:?Run as job array, e.g. --array=1-9%4}"
SAMPLE="$(sed -n "${SLURM_ARRAY_TASK_ID}{s/\r$//;p;}" $HOME/yeast_sv_project/samples.txt)"

BASE="$HOME/yeast_sv_project"
REF_FA="$BASE/ref/GCF_000146045.2_R64_genomic.fna"
REF_GFF="$BASE/ref/genomic.gff"

CAN1_FLANK="${CAN1_FLANK:-50000}"   
YMAX="${YMAX:-auto}"                

ALIGN_DIR="$BASE/align/${SAMPLE}"
SV_DIR="$BASE/sv/${SAMPLE}"
MERGE_DIR="$SV_DIR/merge"
OUT_DIR="$BASE/reports/can1/${SAMPLE}"
mkdir -p "$OUT_DIR" "$MERGE_DIR"

[[ -s "${REF_FA}.fai" ]] || samtools faidx "$REF_FA"

BAM="$ALIGN_DIR/${SAMPLE}.bam"
[[ -s "$BAM" ]] || { echo "[ERROR] BAM not found: $BAM"; exit 2; }

pick_coverage() {
  local cov_file="" cov_mode=""
  local cands=(
    "$ALIGN_DIR/${SAMPLE}.covbase.per-base.bed.gz"
    "$ALIGN_DIR/${SAMPLE}.covbase.per-base.bed"
    "$ALIGN_DIR/${SAMPLE}.cov1k.regions.bed.gz"
    "$ALIGN_DIR/${SAMPLE}.per-base.bed.gz"
    "$ALIGN_DIR/${SAMPLE}.per-base.bed"
    "$ALIGN_DIR/${SAMPLE}.regions.bed.gz"
  )
  for f in "${cands[@]}"; do
    if [[ -s "$f" ]]; then
      cov_file="$f"
      if [[ "$f" == *"per-base"* ]]; then cov_mode="perbase"; else cov_mode="bins"; fi
      break
    fi
  done
  echo "$cov_file" "$cov_mode"
}

read -r COV_FILE COV_MODE < <(pick_coverage)

if [[ -z "${COV_FILE}" ]]; then
  echo "[INFO] Coverage files not found — computing with mosdepth..."
  mosdepth -t "${SLURM_CPUS_PER_TASK}" --mapq 20 "$ALIGN_DIR/${SAMPLE}.covbase" "$BAM"
  mosdepth -t "${SLURM_CPUS_PER_TASK}" --by 1000 --mapq 20 --fast-mode "$ALIGN_DIR/${SAMPLE}.cov1k" "$BAM"
  read -r COV_FILE COV_MODE < <(pick_coverage)
fi

[[ -n "$COV_FILE" ]] || { echo "[ERROR] Coverage not found and could not be generated."; exit 2; }
echo "[INFO] Coverage: $COV_FILE  (mode=$COV_MODE)"

MERGED_ALL="$MERGE_DIR/${SAMPLE}.merged.survivor.vcf.gz"
MERGED_CAN1="$MERGE_DIR/${SAMPLE}.merged.CAN1_${CAN1_FLANK}.vcf.gz"

read -r CAN1_CHR CAN1_START CAN1_END < <(
  awk -F'\t' '
    BEGIN{IGNORECASE=1}
    $0!~/^#/ && $3=="gene" && \
    ($9 ~ /(^|;)gene=CAN1(;|$)/ || $9 ~ /(^|;)(Name|gene_name)=CAN1(;|$)/ || $9 ~ /(^|;)locus_tag=YEL063C(;|$)/) {
      print $1, $4, $5; exit
    }' "$REF_GFF"
) || { echo "[ERROR] CAN1 gene not found in $REF_GFF"; exit 3; }

CHR_LEN=$(awk -v c="$CAN1_CHR" '$1==c{print $2}' "${REF_FA}.fai")
[[ -n "$CHR_LEN" ]] || { echo "[ERROR] Chromosome $CAN1_CHR not found in ${REF_FA}.fai"; exit 4; }

REG_START=$(( CAN1_START - CAN1_FLANK ))
REG_END=$(( CAN1_END + CAN1_FLANK ))
(( REG_START < 1 )) && REG_START=1
(( REG_END > CHR_LEN )) && REG_END=$CHR_LEN
REGION="${CAN1_CHR}:${REG_START}-${REG_END}"
echo "[INFO] REGION=${REGION}  (flank=${CAN1_FLANK})"

if [[ ! -s "$MERGED_CAN1" ]]; then
  if [[ ! -s "$MERGED_ALL" ]]; then
    echo "[ERROR] $MERGED_CAN1 not found and $MERGED_ALL not found. Run 06_merge_filter_CAN1.sbatch first."
    exit 5
  fi
  bcftools view -r "$REGION" -Oz -o "$MERGED_CAN1" "$MERGED_ALL"
  tabix -f -p vcf "$MERGED_CAN1"
fi

CAN1_BED="$OUT_DIR/CAN1_gene.bed"
awk -F'\t' '
  BEGIN{IGNORECASE=1}
  $0!~/^#/ && $3=="gene" && \
  ($9 ~ /(^|;)gene=CAN1(;|$)/ || $9 ~ /(^|;)(Name|gene_name)=CAN1(;|$)/ || $9 ~ /(^|;)locus_tag=YEL063C(;|$)/) {
    printf("%s\t%d\t%d\n", $1, $4-1, $5); exit
  }' "$REF_GFF" > "$CAN1_BED"

PYSCR="$OUT_DIR/plot_can1_${SAMPLE}.py"
cat > "$PYSCR" << 'PYCODE'
import argparse, gzip, re, statistics
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

def opener(path):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path, "rt")

def parse_region(r):
    m = re.match(r"^([^:]+):(\d+)-(\d+)$", r)
    if not m: raise ValueError("Bad region: "+r)
    return m.group(1), int(m.group(2)), int(m.group(3))

def load_gene(bed):
    with open(bed) as fh:
        for line in fh:
            if not line.strip(): continue
            chrom, s0, e = line.strip().split("\t")[:3]
            return chrom, int(s0), int(e)
    return None

def load_cov(path, mode, chrom, beg, end):
    xs, ys = [], []
    with opener(path) as fh:
        for line in fh:
            if not line.strip() or line[0]=='#': continue
            f = line.split()
            c = f[0]
            if c != chrom: continue
            s = int(f[1]); e = int(f[2])
            if e < beg or s > end: continue
            d = float(f[3])
            s2 = max(s, beg-1)
            e2 = min(e, end)
            xs.append(s2+1); ys.append(d)
    z = sorted(zip(xs,ys))
    return [p for p,_ in z], [v for _,v in z]

def load_vcf_events(vcf_path, chrom, beg, end):
    evts = []
    with opener(vcf_path) as fh:
        for line in fh:
            if line.startswith("#"): continue
            f = line.rstrip("\n").split("\t")
            c, pos, info, qual = f[0], int(f[1]), f[7], f[5]
            if c != chrom: continue
            d = {}
            for kv in info.split(";"):
                if "=" in kv:
                    k,v = kv.split("=",1); d[k]=v
                else:
                    d[kv]=True
            svt = d.get("SVTYPE","NA")
            endp = int(d.get("END", pos)) if d.get("END") else pos
            svlen = abs(int(d.get("SVLEN","0"))) if d.get("SVLEN") else 0
            if endp < beg or pos > end: 
                continue
            evts.append(dict(type=svt, pos=pos, end=endp, len=svlen, qual=qual, info=d))
    return evts

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--region", required=True)
    ap.add_argument("--cov", required=True)
    ap.add_argument("--cov_mode", choices=["perbase","bins"], required=True)
    ap.add_argument("--vcf", required=True)
    ap.add_argument("--gene_bed", required=True)
    ap.add_argument("--sample", required=True)
    ap.add_argument("--out_png", required=True)
    ap.add_argument("--out_pdf", required=True)
    ap.add_argument("--ymax", default="auto")
    a = ap.parse_args()

    chrom, beg, end = parse_region(a.region)
    gene = load_gene(a.gene_bed)
    xs, ys = load_cov(a.cov, a.cov_mode, chrom, beg, end)
    evs = load_vcf_events(a.vcf, chrom, beg, end)

    fig = plt.figure(figsize=(12,4))
    ax = plt.gca()

    if xs:
        ax.step(xs, ys, where='post', linewidth=1.0, label="Depth")
        med = statistics.median(ys) if ys else 0
        ax.axhline(med, linestyle="--", linewidth=0.8, alpha=0.6, label=f"median={med:.1f}")
        if a.ymax != "auto":
            ymax = float(a.ymax)
        else:
            p95 = sorted(ys)[int(0.95*len(ys))-1] if ys else 1.0
            ymax = max(5.0, p95*1.3)
        ax.set_ylim(0, ymax)
    else:
        ax.set_ylim(0, 1)
        ax.text(0.02, 0.9, "No coverage points in region", transform=ax.transAxes, color="red")

    if gene is not None and gene[0]==chrom:
        g_s = max(beg, gene[1]+1); g_e = min(end, gene[2])
        if g_s < g_e:
            ax.axvspan(g_s, g_e, color="grey", alpha=0.2, label="CAN1 gene")

    colors = {"DEL":"red","DUP":"blue","INV":"orange","INS":"green"}
    counts = {"DEL":0,"DUP":0,"INV":0,"INS":0}
    for ev in evs:
        t = ev["type"]
        if t not in colors: 
            continue
        counts[t]+=1
        if t=="INS":
            ax.axvline(ev["pos"], color=colors[t], alpha=0.7, linestyle="--")
        else:
            s = max(beg, ev["pos"]); e = min(end, ev["end"])
            if s<e:
                ax.axvspan(s, e, color=colors[t], alpha=0.25)

    ax.set_xlim(beg, end)
    ax.set_xlabel(f"{chrom} position (bp)")
    ax.set_ylabel("Depth (x)")
    ttl = f"{a.sample} — CAN1 locus ({chrom}:{beg}-{end})  |  SV: " + \
          ", ".join([f"{k}={v}" for k,v in counts.items()])
    ax.set_title(ttl, fontsize=10)
    ax.legend(loc="upper right", fontsize=8, frameon=False)
    plt.tight_layout()
    fig.savefig(a.out_png, dpi=200)
    fig.savefig(a.out_pdf)
    print(f"[OK] Saved: {a.out_png}\n[OK] Saved: {a.out_pdf}")

if __name__ == "__main__":
    main()
PYCODE

PNG="$OUT_DIR/${SAMPLE}.CAN1_${CAN1_FLANK}.png"
PDF="$OUT_DIR/${SAMPLE}.CAN1_${CAN1_FLANK}.pdf"

python "$PYSCR" \
  --region "$REGION" \
  --cov "$COV_FILE" \
  --cov_mode "$COV_MODE" \
  --vcf "$MERGED_CAN1" \
  --gene_bed "$CAN1_BED" \
  --sample "$SAMPLE" \
  --out_png "$PNG" \
  --out_pdf "$PDF" \
  --ymax "$YMAX"

echo "[DONE] Report: $PNG"