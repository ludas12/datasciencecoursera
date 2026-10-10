#!/usr/bin/env bash
# Step 3: process ONE sample: trim -> align -> deduplicate -> QC -> M-bias ->
# methylation calls (genome-wide CpG + spike-in controls).
# Usage: bash 03_process_sample.sh <sample_index (1-based)>
#   (run_slurm.sh submits it as an array job using SLURM_ARRAY_TASK_ID).
# Finished sub-steps are skipped on re-run. To redo some, list them in REDO,
# e.g.  REDO=mbias,extract bash 03_process_sample.sh 1
# Under SLURM the script runs from a spool copy, so PIPE_DIR is passed in.
source "${PIPE_DIR:-$(dirname "$0")}/scripts/common.sh"
activate_env

idx="${1:-${SLURM_ARRAY_TASK_ID:-}}"
[[ -n "$idx" ]] || die "usage: $0 <sample_index>"
row="$(sample_row "$idx")"
[[ -n "$row" ]] || die "no sample number $idx in $SAMPLESHEET"
IFS=$'\t' read -r sample r1_list r2_list <<< "$row"
IFS=',' read -r -a r1_files <<< "$r1_list"
IFS=',' read -r -a r2_files <<< "${r2_list:-}"
paired=false; [[ ${#r2_files[@]} -gt 0 && -n "${r2_files[0]}" ]] && paired=true

S="$OUTDIR/samples/$sample"
mkdir -p "$S"/{trim,align,qc,mbias,meth}
[[ -n "${REDO:-}" ]] && for step in ${REDO//,/ }; do rm -f "$(done_marker "$S" "$step")"; done
[[ -s "$REF_FA.fai" ]] || die "reference not ready ($REF_FA.fai missing); run step 02 first"
log "sample $sample (#$idx), paired=$paired, aligner=$ALIGNER, threads=$THREADS"

FINAL_BAM="$S/align/$sample.final.bam"
T1="$S/trim/${sample}_R1.trimmed.fq.gz"
T2="$S/trim/${sample}_R2.trimmed.fq.gz"

# ---- 1. adapter/quality trimming ---------------------------------------------
if ! is_done "$S" trim && ! is_done "$S" align; then
    # Merge lanes if a sample has several files per read.
    in1="${r1_files[0]}"; in2="${r2_files[0]:-}"
    if [[ ${#r1_files[@]} -gt 1 ]]; then
        log "merging ${#r1_files[@]} lane files"
        in1="$S/trim/merged_R1.fq.gz"; cat "${r1_files[@]}" > "$in1"
        if $paired; then in2="$S/trim/merged_R2.fq.gz"; cat "${r2_files[@]}" > "$in2"; fi
    fi
    fastp_args=(-i "$in1" -o "$T1"
        --trim_front1 "$CLIP_5P_R1" --trim_tail1 "$CLIP_3P_R1"
        --length_required "$MIN_READ_LEN"
        --thread "$(( THREADS > 16 ? 16 : THREADS ))"
        --json "$S/qc/$sample.fastp.json" --html "$S/qc/$sample.fastp.html"
        --report_title "$sample")
    $paired && fastp_args+=(-I "$in2" -O "$T2" --detect_adapter_for_pe
        --trim_front2 "$CLIP_5P_R2" --trim_tail2 "$CLIP_3P_R2")
    if [[ "$TRIM_POLY_G" == true ]]; then fastp_args+=(--trim_poly_g)
    else fastp_args+=(--disable_trim_poly_g); fi
    log "fastp"
    fastp "${fastp_args[@]}" 2> "$S/qc/$sample.fastp.log"
    rm -f "$S/trim"/merged_R?.fq.gz
    mark_done "$S" trim
fi

# ---- 2. alignment + duplicate marking/removal --------------------------------
if ! is_done "$S" align; then
    reads=("$T1"); $paired && reads+=("$T2")
    tmp="$S/align/tmp"; mkdir -p "$tmp"
    mem_gb="${SAMPLE_MEM//[!0-9]/}"
    sort_mem="$(( mem_gb * 1024 / (2 * THREADS) ))M"   # half the job memory for sorting
    case "$ALIGNER" in
    bwameth)
        log "bwa-meth alignment"
        bwameth.py --reference "$REF_FA" -t "$THREADS" \
            --read-group "@RG\tID:$sample\tSM:$sample\tPL:ILLUMINA" "${reads[@]}" \
            2> "$S/align/$sample.bwameth.log" \
          | samtools fixmate -m -@ 2 - - \
          | samtools sort -@ "$THREADS" -m "$sort_mem" -T "$tmp/sort" -o "$tmp/$sample.sorted.bam" -
        log "marking duplicates"
        samtools markdup -@ "$THREADS" -T "$tmp/markdup" -f "$S/qc/$sample.markdup.txt" \
            "$tmp/$sample.sorted.bam" "$FINAL_BAM"
        ;;
    bismark)
        log "Bismark alignment"
        if $paired; then bm_in=(-1 "$T1" -2 "$T2"); else bm_in=("$T1"); fi
        bismark --genome "$REF_DIR" "${bm_in[@]}" --basename "$sample" \
            -p "$(( THREADS > 2 ? THREADS / 2 : 1 ))" --temp_dir "$tmp" -o "$tmp" \
            2> "$S/align/$sample.bismark.log"
        cp "$tmp/${sample}"_*E_report.txt "$S/qc/"
        log "deduplicate_bismark"
        bm_bam="$tmp/${sample}_pe.bam"; $paired || bm_bam="$tmp/$sample.bam"
        [[ -s "$bm_bam" ]] || die "Bismark output $bm_bam not found"
        deduplicate_bismark --bam --output_dir "$tmp" "$bm_bam" 2> "$S/align/$sample.dedup.log"
        cp "$tmp/"*deduplication_report.txt "$S/qc/"
        samtools sort -@ "$THREADS" -m "$sort_mem" -T "$tmp/sort" -o "$FINAL_BAM" "$tmp/"*.deduplicated.bam
        ;;
    *) die "unknown ALIGNER '$ALIGNER'" ;;
    esac
    samtools index -@ "$THREADS" "$FINAL_BAM"
    rm -rf "$tmp" "$S/trim"/*.fq.gz
    mark_done "$S" align
fi

# ---- 3. alignment QC -----------------------------------------------------------
if ! is_done "$S" qc; then
    log "alignment QC"
    samtools flagstat -@ "$THREADS" "$FINAL_BAM" > "$S/qc/$sample.flagstat.txt"
    samtools idxstats "$FINAL_BAM" > "$S/qc/$sample.idxstats.txt"
    samtools stats -@ "$THREADS" "$FINAL_BAM" > "$S/qc/$sample.samtools_stats.txt"
    mark_done "$S" qc
fi

# ---- 4. M-bias (position-dependent methylation along reads) ------------------
if ! is_done "$S" mbias; then
    log "MethylDackel mbias"
    MethylDackel mbias -@ "$THREADS" -q "$MIN_MAPQ" -p "$MIN_BASEQ" --txt \
        "$REF_FA" "$FINAL_BAM" "$S/mbias/$sample" \
        > "$S/mbias/$sample.mbias.txt" 2> "$S/mbias/$sample.mbias_suggestion.txt"
    mark_done "$S" mbias
fi

# ---- 5. methylation calls ------------------------------------------------------
if ! is_done "$S" extract; then
    # shellcheck disable=SC2206
    md=(-@ "$THREADS" -q "$MIN_MAPQ" -p "$MIN_BASEQ" $MD_INCLUSION)
    log "MethylDackel extract (CpG, strands merged)${MD_INCLUSION:+ with $MD_INCLUSION}"
    MethylDackel extract "${md[@]}" --mergeContext -o "$S/meth/$sample" "$REF_FA" "$FINAL_BAM"
    if [[ "$MD_METHYLKIT" == true ]]; then
        log "MethylDackel extract (methylKit format)"
        MethylDackel extract "${md[@]}" --methylKit -o "$S/meth/$sample" "$REF_FA" "$FINAL_BAM"
    fi
    # Spike-in controls, all contexts, for conversion efficiency.
    for ctrl in "$CONTROL_UNMETH" "$CONTROL_METH"; do
        [[ -n "$ctrl" ]] && cut -f1 "$REF_FA.fai" | grep -qx "$ctrl" || continue
        MethylDackel extract "${md[@]}" --CHG --CHH -r "$ctrl" \
            -o "$S/meth/${sample}.control_${ctrl}" "$REF_FA" "$FINAL_BAM"
    done
    gzip -f "$S/meth/$sample"_CpG*.bedGraph "$S/meth/$sample"*.methylKit 2>/dev/null || true
    mark_done "$S" extract
fi

log "sample $sample finished"
