#!/usr/bin/env bash
# Step 11: BAL-versus-plasma checks for each subject/timepoint pair in
# downstream/samples.tsv: (1) imprinted promoters per gene with per-CpG coverage,
# CpG-resampling CIs and leave-one-CpG-out stability, sign test per timepoint;
# (2) BAL against plasma by methylation level (strata of reference methylation),
# testing whether BAL looks like uniform 5mC loss; (3) the same for short
# (< FRAG_SPLIT_BP) and long fragments separately, to see whether the BAL
# deficit is carried by its short fragments. Needs steps 04 and 05.
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

# Methylation by fragment length for plasma and BAL libraries, summed over the
# step 05 bins. Fragment length = |TLEN| + bases clipped before alignment, the
# same read-end exclusions as step 03 (MD_INCLUSION, CLIP_5P) are applied.
SZ="$C/size_split"
BINS="$DOWN_DIR/regions/beds/bins.bed"
if [[ -s "$BINS" ]]; then
    mkdir -p "$SZ/bins" "$SZ/tmp"
    while IFS=$'\t' read -r sample label type _ include _; do
        [[ "$sample" == sample || "$include" != yes ]] && continue
        [[ "${type,,}" == plasma || "${type,,}" == bal ]] || continue
        bam="$OUTDIR/samples/$sample/align/$sample.final.bam"
        [[ -s "$bam" ]] || { log "WARNING: no BAM for $sample; skipping size split"; continue; }
        [[ "$SZ/bins/$label.long.tsv" -nt "$bam" && "${FORCE:-0}" != 1 ]] && continue
        read -r c1 c2 < <(cat "$OUTDIR/samples/$sample/qc/$sample.clip5p.txt" 2>/dev/null || echo "$CLIP_5P_R1 $CLIP_5P_R2")
        ex1=$(( CLIP_5P_R1 - c1 )); ex2=$(( CLIP_5P_R2 - c2 ))
        (( ex1 < 0 )) && ex1=0; (( ex2 < 0 )) && ex2=0
        md_ends=(); (( ex1 > 0 || ex2 > 0 )) && md_ends=(--nOT "$ex1,0,$ex2,0" --nOB "$ex1,0,$ex2,0")
        t=$(( FRAG_SPLIT_BP - c1 - c2 ))      # TLEN threshold for the corrected length
        log "$label: methylation of fragments < and >= $FRAG_SPLIT_BP bp (|TLEN| threshold $t)"
        activate_env "$CONDA_ENV"
        for cls in short long; do
            if [[ $cls == short ]]; then expr="(tlen > 0 && tlen < $t) || (tlen < 0 && tlen > -$t)"
            else expr="tlen >= $t || tlen <= -$t"; fi
            sub="$SZ/tmp/$label.$cls.bam"
            samtools view -@ "$DOWNSTREAM_THREADS" -b -e "$expr" -o "$sub" "$bam"
            samtools index "$sub"
            # shellcheck disable=SC2086
            MethylDackel extract -@ "$DOWNSTREAM_THREADS" -q "$MIN_MAPQ" -p "$MIN_BASEQ" \
                ${md_ends[@]+"${md_ends[@]}"} $MD_INCLUSION --mergeContext \
                -o "$SZ/tmp/$label.$cls" "$REF_FA" "$sub"
        done
        activate_env "$DOWNSTREAM_ENV"
        for cls in short long; do
            awk -v OFS='\t' 'NR > 1 {print $1, $2, $3, $5, $6}' "$SZ/tmp/$label.${cls}_CpG.bedGraph" \
                | LC_ALL=C sort -k1,1 -k2,2n -S 2G -T "$SZ/tmp" \
                | bedtools map -a "$BINS" -b - -c 4,5 -o sum,sum > "$SZ/bins/$label.$cls.tsv.part"
            mv "$SZ/bins/$label.$cls.tsv.part" "$SZ/bins/$label.$cls.tsv"
        done
        rm -f "$SZ/tmp/$label".*
    done < "$SAMPLE_INFO"
    rmdir "$SZ/tmp" 2>/dev/null || true
else
    log "WARNING: no $BINS (step 05); size split skipped"
fi

python3 "$PIPE_DIR/downstream/compartment_check.py" --sample-info "$SAMPLE_INFO" \
    --imprinted-cpgs "$cpgs" --regions-dir "$DOWN_DIR/regions" --size-dir "$SZ" \
    --min-calls "$MIN_REGION_CALLS" --max-cov-fold "$MAX_COV_FOLD" --out "$C" \
    | tee "$C/compartment_summary.txt"
log "results in $C"
