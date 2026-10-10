#!/usr/bin/env bash
# Step 6: cell-type / tissue-of-origin deconvolution with UXM and the human
# methylation atlas of Loyfer et al. 2023 (Nature 613:355). UXM counts, at each
# atlas marker, the fraction of DNA fragments that are fully unmethylated and
# fits the cell-type mixture by non-negative least squares. It works on
# fragments, so it suits low-coverage plasma/BAL data.
# Needs downstream/setup_downstream.sh first. Submit with run_downstream.sh or:
#   bash 06_tissue_of_origin.sh
# shellcheck disable=SC1091
source "${PIPE_DIR:-$(dirname "$0")}/scripts/common.sh"
activate_env "$DOWNSTREAM_ENV"

WGBS="$TOOLS_DIR/wgbs_tools"
UXM="$TOOLS_DIR/UXM_deconv"
ATLAS="$UXM/supplemental/$UXM_ATLAS"
[[ -x "$WGBS/wgbstools" && -x "$UXM/uxm" ]] || die "wgbstools/UXM not installed; run downstream/setup_downstream.sh"
[[ -s "$ATLAS" ]] || die "atlas not found: $ATLAS"
export PATH="$WGBS:$UXM:$PATH"
T="$DOWN_DIR/tissue_of_origin"
mkdir -p "$T/pat" "$T/by_label" "$T/uxm_tmp"

# ---- 1. wgbstools genome (CpG index of hg38 chr1-22,X,Y,M from our FASTA) -----
if [[ ! -f "$TOOLS_DIR/.done.wgbstools_hg38" ]]; then
    log "wgbstools init_genome hg38 (once; ~30-60 min)"
    wgbstools init_genome hg38 --fasta_path "$REF_FA" -f
    touch "$TOOLS_DIR/.done.wgbstools_hg38"
fi

# ---- 2. BAM -> pat (only reads overlapping atlas markers, so this is quick) ---
tail -n +2 "$ATLAS" | cut -f1-3 | LC_ALL=C sort -k1,1 -k2,2n -u > "$T/atlas_markers.bed"
pats=()
while IFS=$'\t' read -r sample label _ _ include; do
    [[ "$sample" == sample || "$include" != yes ]] && continue
    bam="$OUTDIR/samples/$sample/align/$sample.final.bam"
    [[ -s "$bam" ]] || { log "WARNING: no BAM for $sample; skipping"; continue; }
    pat="$T/pat/$sample.final.pat.gz"
    if [[ ! -s "$pat" ]]; then
        log "$label: bam2pat"
        wgbstools bam2pat "$bam" --genome hg38 -L "$T/atlas_markers.bed" -q "$MIN_MAPQ" \
            -@ "$DOWNSTREAM_THREADS" -o "$T/pat" -f
    fi
    # UXM names samples after the pat file, so link it under the short label.
    ln -sf "../pat/$sample.final.pat.gz" "$T/by_label/$label.pat.gz"
    [[ -e "$pat.csi" ]] && ln -sf "../pat/$sample.final.pat.gz.csi" "$T/by_label/$label.pat.gz.csi"
    pats+=("$T/by_label/$label.pat.gz")
done < "$SAMPLE_INFO"
[[ ${#pats[@]} -gt 0 ]] || die "no samples to deconvolve"

# ---- 3. deconvolution + plot -----------------------------------------------------
log "UXM deconvolution of ${#pats[@]} samples with $(basename "$ATLAS")"
uxm deconv "${pats[@]}" --atlas "$ATLAS" -@ "$DOWNSTREAM_THREADS" -T "$T/uxm_tmp" \
    -o "$T/uxm_deconv.csv"
uxm plot "$T/uxm_deconv.csv" -o "$T/uxm_deconv.pdf" --min_rate 1 \
    || log "WARNING: uxm plot failed (the CSV is still valid)"
python3 "$PIPE_DIR/downstream/summarize_uxm.py" "$T/uxm_deconv.csv" "$T" | tee "$T/uxm_summary.txt"
log "tissue-of-origin results in $T"
