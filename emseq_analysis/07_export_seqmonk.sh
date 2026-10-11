#!/usr/bin/env bash
# Step 7: export per-CpG methylation as Bismark coverage files for SeqMonk.
# Writes OUTDIR/seqmonk/<label>.cov.gz (label from downstream/samples.tsv,
# otherwise the sample name), keeping chr1-22, X, Y and M only. In SeqMonk:
# File > Import Data > Bismark (cov). Run directly (a few minutes per sample),
# or it runs automatically at the end of step 04:  bash 07_export_seqmonk.sh
# shellcheck disable=SC1091
source "${PIPE_DIR:-$(dirname "$0")}/scripts/common.sh"
[[ -s "$SAMPLESHEET" ]] || die "no sample sheet at $SAMPLESHEET"
mkdir -p "$OUTDIR/seqmonk"

label_of() {
    local l=""
    [[ -s "$SAMPLE_INFO" ]] && l=$(awk -F'\t' -v s="$1" '$1 == s {print $2; exit}' "$SAMPLE_INFO")
    echo "${l:-$1}"
}

n=0
while IFS=$'\t' read -r sample _; do
    bg="$OUTDIR/samples/$sample/meth/${sample}_CpG.bedGraph.gz"
    [[ -s "$bg" ]] || { log "WARNING: no CpG calls for $sample; skipping"; continue; }
    out="$OUTDIR/seqmonk/$(label_of "$sample").cov.gz"
    [[ "$out" -nt "$bg" ]] && continue          # already exported and up to date
    log "$sample -> $(basename "$out")"
    # MethylDackel bedGraph (0-based start, track header) -> Bismark .cov
    # (1-based): chr, start, end, %meth, n_meth, n_unmeth
    zcat "$bg" | awk -v OFS='\t' 'NR > 1 && $1 ~ /^(chr)?([0-9]+|X|Y|M|MT)$/ {print $1, $2 + 1, $3, $4, $5, $6}' \
        | gzip > "$out.part"
    mv "$out.part" "$out"
    n=$((n + 1))
done < <(tail -n +2 "$SAMPLESHEET")
log "SeqMonk files in $OUTDIR/seqmonk ($n written this run)"
