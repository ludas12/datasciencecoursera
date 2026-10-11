#!/usr/bin/env bash
# Step 11: BAL-versus-plasma checks for each subject/timepoint pair in
# downstream/samples.tsv: (1) imprinted promoters per gene with per-CpG coverage,
# CpG-resampling CIs and leave-one-CpG-out stability, sign test per timepoint;
# (2) BAL against plasma by methylation level (strata of reference methylation),
# testing whether BAL looks like uniform 5mC loss. Needs steps 04 and 05.
# Run with run_downstream.sh 11, or:  bash 11_compartment_check.sh
# shellcheck disable=SC1091
source "${PIPE_DIR:-$(dirname "$0")}/scripts/common.sh"
activate_env "$DOWNSTREAM_ENV"

C="$DOWN_DIR/compartment"
mkdir -p "$C"
IMPRINTED_BED="$OUTDIR/summary/imprinted_promoters.bed"
[[ -s "$IMPRINTED_BED" ]] || die "no $IMPRINTED_BED; run step 04 first"
[[ -d "$DOWN_DIR/regions/mapped" ]] || log "WARNING: no step 05 region counts; only the imprinted check will run"

# Per-CpG calls inside the imprinted promoters, one block per library
cpgs="$C/imprinted_cpgs.tsv"
: > "$cpgs.part"
while IFS=$'\t' read -r sample label _ _ include _; do
    [[ "$sample" == sample || "$include" != yes ]] && continue
    bg="$OUTDIR/samples/$sample/meth/${sample}_CpG.bedGraph.gz"
    [[ -s "$bg" ]] || { log "WARNING: no CpG calls for $sample; skipping"; continue; }
    zcat "$bg" | awk 'NR > 1' | bedtools intersect -wa -wb -a - -b "$IMPRINTED_BED" \
        | awk -v OFS='\t' -v l="$label" '{print l, $1, $2, $10, $5, $6}' | sort -u >> "$cpgs.part"
done < "$SAMPLE_INFO"
mv "$cpgs.part" "$cpgs"

python3 "$PIPE_DIR/downstream/compartment_check.py" --sample-info "$SAMPLE_INFO" \
    --imprinted-cpgs "$cpgs" --regions-dir "$DOWN_DIR/regions" \
    --min-calls "$MIN_REGION_CALLS" --max-cov-fold "$MAX_COV_FOLD" --out "$C" \
    | tee "$C/compartment_summary.txt"
log "results in $C"
