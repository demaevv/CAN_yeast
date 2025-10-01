#!/usr/bin/env bash

set -eo pipefail
module load miniconda 2>/dev/null || true
source "$HOME/miniconda3/etc/profile.d/conda.sh" 2>/dev/null || true

BASE="$HOME/yeast_sv_project"
OUTROOT="$BASE/reports/can1"

CAN1_FLANK="${CAN1_FLANK:-50000}"

SAMPLES=()
while IFS= read -r line; do
  [[ -n "$line" ]] || continue
  SAMPLES+=("$(echo "$line" | sed 's/\r$//')")
done < "$BASE/samples.txt"

mkdir -p "$OUTROOT"
HTML="$OUTROOT/CAN1_${CAN1_FLANK}_overview.html"
TS="$(date '+%Y-%m-%d %H:%M:%S %Z')"

cat > "$HTML" <<EOF
<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="utf-8">
<title>CAN1 locus overview (±${CAN1_FLANK} bp)</title>
<meta name="viewport" content="width=device-width, initial-scale=1">
<style>
  :root{--fg:#222;--muted:#666;--card:#fff;--bd:#ddd;--bd2:#eee;}
  body{font-family:system-ui,-apple-system,Segoe UI,Roboto,Arial,sans-serif;margin:16px;color:var(--fg);}
  h1{font-size:20px;margin:0 0 12px 0}
  .meta{color:var(--muted);font-size:12px;margin-bottom:16px}
  .grid{display:grid;grid-template-columns:repeat(auto-fill,minmax(360px,1fr));gap:16px}
  .card{border:1px solid var(--bd);border-radius:12px;padding:12px;box-shadow:0 1px 3px rgba(0,0,0,.04);background:var(--card)}
  .card h2{font-size:16px;margin:0 0 8px 0}
  .thumb{width:100%;height:auto;border-radius:8px;display:block;border:1px solid var(--bd2)}
  .links{margin-top:8px;font-size:13px}
  .links a{margin-right:10px;text-decoration:none;color:#0a58ca}
  .miss{color:#b00020;font-size:13px;margin:8px 0 0 0}
  .stats{margin-top:6px;font-size:12px;color:#333}
  .stats code{background:#f6f8fa;border:1px solid #e6e8eb;border-radius:6px;padding:2px 6px;margin-right:6px;display:inline-block}
  footer{margin-top:20px;color:var(--muted);font-size:12px}
</style>
</head>
<body>
<h1>CAN1 locus overview (±${CAN1_FLANK} bp) — ${#SAMPLES[@]} samples</h1>
<div class="meta">Generated: ${TS}. Paths are relative to <code>reports/can1/</code>.</div>
<div class="grid">
EOF

for s in "${SAMPLES[@]}"; do
  PNG="${s}/${s}.CAN1_${CAN1_FLANK}.png"
  PDF="${s}/${s}.CAN1_${CAN1_FLANK}.pdf"
  TSV_REL="../../sv/${s}/merge/${s}.CAN1_${CAN1_FLANK}.summary.tsv"
  VCF_REL="../../sv/${s}/merge/${s}.merged.CAN1_${CAN1_FLANK}.vcf.gz"

  echo '  <div class="card">' >> "$HTML"
  echo "    <h2>${s}</h2>" >> "$HTML"

  if [[ -s "$OUTROOT/$PNG" ]]; then
    echo "    <a href=\"$PNG\" target=\"_blank\"><img class=\"thumb\" src=\"$PNG\" alt=\"${s} CAN1\"></a>" >> "$HTML"
  else
    echo "    <div class=\"miss\">Image not found: <code>${PNG}</code>. Run 07_can_locus_report for this sample.</div>" >> "$HTML"
  fi

  STAT_HTML=""
  if [[ -s "$OUTROOT/$TSV_REL" ]]; then
    STAT_HTML=$(awk -F'\t' '
      BEGIN{c["DEL"]=c["INS"]=c["INV"]=c["DUP"]=0; g["DEL"]=g["INS"]=g["INV"]=g["DUP"]=0}
      NR>1{
        t=$5; if(t in c){ c[t]++ }
        if($10=="YES" && (t in g)){ g[t]++ }
      }
      END{
        for(t in c){ } # no-op for gawk ordering
        types[1]="DEL"; types[2]="INS"; types[3]="INV"; types[4]="DUP";
        for(i=1;i<=4;i++){
          t=types[i];
          ig=(g[t]>0)?sprintf(" (in gene=%d)", g[t]):"";
          printf("<code>%s=%d%s</code> ", t, c[t], ig);
        }
      }' "$OUTROOT/$TSV_REL")
  fi
  if [[ -n "$STAT_HTML" ]]; then
    echo "    <div class=\"stats\">$STAT_HTML</div>" >> "$HTML"
  fi

  echo '    <div class="links">' >> "$HTML"
  [[ -s "$OUTROOT/$PDF"     ]] && echo "      <a href=\"$PDF\" target=\"_blank\">PDF</a>" >> "$HTML"
  [[ -s "$OUTROOT/$PNG"     ]] && echo "      <a href=\"$PNG\" download>Download PNG</a>" >> "$HTML"
  [[ -s "$OUTROOT/$TSV_REL" ]] && echo "      <a href=\"$TSV_REL\" target=\"_blank\">Summary TSV</a>" >> "$HTML"
  [[ -s "$OUTROOT/$VCF_REL" ]] && echo "      <a href=\"$VCF_REL\" target=\"_blank\">VCF.gz</a>" >> "$HTML"
  echo '    </div>' >> "$HTML"

  if [[ ! -s "$OUTROOT/$TSV_REL" || ! -s "$OUTROOT/$VCF_REL" ]]; then
    echo "    <div class=\"miss\">" >> "$HTML"
    [[ ! -s "$OUTROOT/$TSV_REL" ]] && echo "Missing <code>Summary TSV</code> — run <code>06_merge_filter_CAN1.sbatch</code>." >> "$HTML"
    [[ ! -s "$OUTROOT/$VCF_REL" ]] && echo "Missing <code>VCF.gz</code> — run <code>06_merge_filter_CAN1.sbatch</code>." >> "$HTML"
    echo "    </div>" >> "$HTML"
  fi

  echo '  </div>' >> "$HTML"
done

cat >> "$HTML" <<'EOF'
</div>
<footer>
  Tip: click an image to open it full-size in a new tab. Use Summary TSV/VCF links for exact coordinates.
</footer>
</body>
</html>
EOF

echo "[DONE] HTML report: $HTML"