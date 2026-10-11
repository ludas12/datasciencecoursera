#!/usr/bin/env bash
# Step 4: QC summary table across all samples + MultiQC report.
# Submitted by run_slurm.sh after all samples finish, or:  bash 04_summarize.sh
# Under SLURM the script runs from a spool copy, so PIPE_DIR is passed in.
source "${PIPE_DIR:-$(dirname "$0")}/scripts/common.sh"
activate_env
mkdir -p "$OUTDIR/summary"

# Imprinted-gene promoters (TSS +/- 1 kb, RefSeq) for the 5mC-protection check.
IMPRINTED_BED="$OUTDIR/summary/imprinted_promoters.bed"
REFSEQ="$OUTDIR/annotation/ncbiRefSeqCurated.txt.gz"
if [[ ! -s "$IMPRINTED_BED" ]]; then
    mkdir -p "$OUTDIR/annotation"
    [[ -s "$REFSEQ" ]] || curl -fsSL "$REFSEQ_URL" -o "$REFSEQ" 2>/dev/null || rm -f "$REFSEQ"
    if [[ -s "$REFSEQ" ]]; then
        zcat "$REFSEQ" | awk -v OFS='\t' -v genes="$IMPRINTED_GENES" '
            BEGIN { n = split(genes, g, " "); for (i = 1; i <= n; i++) want[g[i]] = 1 }
            ($13 in want) && $3 ~ /^chr[0-9XY]+$/ {
                tss = ($4 == "+") ? $5 : $6; s = tss - 1000; if (s < 0) s = 0
                print $3, s, tss + 1000, $13 }' | sort -k1,1 -k2,2n -u > "$IMPRINTED_BED"
    else
        log "WARNING: RefSeq annotation unavailable; imprinted-promoter check skipped"
    fi
fi

python3 "$PIPE_DIR/scripts/qc_summary.py" "$OUTDIR" \
    CONTROL_UNMETH="$CONTROL_UNMETH" CONTROL_METH="$CONTROL_METH" \
    MAX_LAMBDA_CPG_METH="$MAX_LAMBDA_CPG_METH" MIN_PUC19_CPG_METH="$MIN_PUC19_CPG_METH" \
    MIN_MAPPING_RATE="$MIN_MAPPING_RATE" MAX_DUP_RATE="$MAX_DUP_RATE" \
    MAX_CHH_METH="$MAX_CHH_METH" MIN_CONTROL_CALLS="$MIN_CONTROL_CALLS" \
    MIN_IMPRINTED_METH="$MIN_IMPRINTED_METH" MAX_IMPRINTED_METH="$MAX_IMPRINTED_METH" \
    MIN_GLOBAL_METH="$MIN_GLOBAL_METH" IMPRINTED_BED="$IMPRINTED_BED" \
    | tee "$OUTDIR/summary/qc_summary.txt"

if command -v multiqc >/dev/null 2>&1; then
    log "MultiQC"
    multiqc -f -q -o "$OUTDIR/summary/multiqc" -n emseq_multiqc_report "$OUTDIR/samples" \
        || log "WARNING: MultiQC failed (summary table above is unaffected)"
else
    log "multiqc not in env; skipping (pip/conda install multiqc to get the HTML report)"
fi
bash "$PIPE_DIR/07_export_seqmonk.sh" || log "WARNING: SeqMonk export failed (run 07_export_seqmonk.sh to retry)"
log "summary written to $OUTDIR/summary"
