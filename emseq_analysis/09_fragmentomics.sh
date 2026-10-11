#!/usr/bin/env bash
# Step 9: cfDNA fragmentomics from the aligned reads (no extra sequencing).
#   fragment sizes, genome-wide short/long ratio (5 Mb bins, GC-corrected),
#   nucleosome footprints at TSSs of marker gene sets, 4-mer end motifs.
# True fragment ends are restored for libraries whose reads were clipped before
# alignment (read from samples/<s>/qc/<s>.clip5p.txt, else CLIP_5P_R1/R2).
# Run with run_downstream.sh 09, or:  bash 09_fragmentomics.sh
# shellcheck disable=SC1091
source "${PIPE_DIR:-$(dirname "$0")}/scripts/common.sh"
activate_env "$DOWNSTREAM_ENV"

F="$DOWN_DIR/fragmentomics"
mkdir -p "$F/per_library"
REFSEQ="$ANNOT_DIR/ncbiRefSeqCurated.txt.gz"
[[ -s "$REFSEQ" ]] || log "WARNING: $REFSEQ missing; TSS footprints skipped (run downstream/setup_downstream.sh)"

jobs="$F/jobs.txt"; : > "$jobs"
while IFS=$'\t' read -r sample label _ _ include _; do
    [[ "$sample" == sample || "$include" != yes ]] && continue
    bam="$OUTDIR/samples/$sample/align/$sample.final.bam"
    [[ -s "$bam" ]] || { log "WARNING: no BAM for $sample; skipping"; continue; }
    out="$F/per_library/$label.frag.npz"
    [[ "$out" -nt "$bam" ]] && continue
    read -r c1 c2 < <(cat "$OUTDIR/samples/$sample/qc/$sample.clip5p.txt" 2>/dev/null || echo "$CLIP_5P_R1 $CLIP_5P_R2")
    printf '%s\t%s\t%s\t%s\t%s\n' "$bam" "$label" "$c1" "$c2" "$out" >> "$jobs"
done < "$SAMPLE_INFO"

n=$(wc -l < "$jobs")
if [[ $n -gt 0 ]]; then
    log "extracting fragments from $n libraries"
    export PIPE_DIR REF_FA REFSEQ
    # shellcheck disable=SC2016
    xargs -P "$DOWNSTREAM_THREADS" -L 1 bash -c 'python3 "$PIPE_DIR/downstream/fragmentomics.py" extract \
        --bam "$0" --label "$1" --clip "$2" "$3" --out "$4" \
        --fasta "$REF_FA" --fai "$REF_FA.fai" --refseq "$REFSEQ"' < "$jobs"
fi
python3 "$PIPE_DIR/downstream/fragmentomics.py" summarize --dir "$F/per_library" \
    --sample-info "$SAMPLE_INFO" --out "$F" | tee "$F/fragmentomics_summary.txt"
