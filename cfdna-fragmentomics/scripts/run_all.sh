#!/usr/bin/env bash
# Run the whole pipeline for every sample in a sample sheet.
#
# Usage: scripts/run_all.sh samples.tsv REF.fa BLACKLIST.bed OUTDIR [THREADS] [BIN_SIZE]
#
# samples.tsv: tab-separated with header 'sample group r1 r2' (see samples.example.tsv).
# Steps already completed (output exists) are skipped, so the script can be re-run.
set -euo pipefail
if [[ $# -lt 4 ]]; then sed -n '2,8p' "$0"; exit 1; fi
SHEET=$1; REF=$2; BL=$3; OUT=$4; THREADS=${5:-8}; BIN=${6:-5000000}
HERE=$(cd "$(dirname "$0")/.." && pwd)
FRAG="python $HERE/fragmentomics/cfdna_frag.py"

tail -n +2 "$SHEET" | while IFS=$'\t' read -r S GROUP R1 R2 _; do
    [[ -z "$S" ]] && continue
    echo "== $S ($GROUP)"
    [[ -s "$OUT/bam/$S.bam" ]] || "$HERE/scripts/01_preprocess.sh" "$S" "$R1" "$R2" "$REF" "$OUT/bam" "$THREADS"
    [[ -s "$OUT/frags/$S.frags.tsv.gz" ]] || { mkdir -p "$OUT/frags"; \
        $FRAG extract -b "$OUT/bam/$S.bam" -r "$REF" --blacklist "$BL" -o "$OUT/frags/$S.frags.tsv.gz"; }
    $FRAG features -f "$OUT/frags/$S.frags.tsv.gz" -r "$REF" -o "$OUT/features/$S" -n "$S" --bin-size "$BIN"
    if [[ -n "${TSS_BED:-}" ]]; then
        $FRAG tss -f "$OUT/frags/$S.frags.tsv.gz" -t "$TSS_BED" -o "$OUT/features/$S"
    fi
done

$FRAG cohort -i "$OUT"/features/* -s "$SHEET" -o "$OUT/cohort"

# One-line QC per sample: mapping rate and duplication (BAL often has low human mapping).
{
  printf "sample\tmapped_pct\tdup_pct\n"
  for f in "$OUT"/bam/*.flagstat; do
      s=$(basename "$f" .flagstat)
      mapped=$(grep -m1 "mapped (" "$f" | sed -E 's/.*\(([0-9.]+)%.*/\1/')
      total=$(grep -m1 "in total" "$f" | cut -d' ' -f1)
      dups=$(grep -m1 "duplicates" "$f" | cut -d' ' -f1)
      printf "%s\t%s\t%.2f\n" "$s" "$mapped" "$(echo "100*$dups/$total" | bc -l)"
  done
} > "$OUT/cohort/alignment_qc.tsv"
echo "Done. See $OUT/cohort/"
