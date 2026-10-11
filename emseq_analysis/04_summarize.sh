#!/usr/bin/env bash
# Step 4: QC summary table across all samples + MultiQC report.
# Submitted by run_slurm.sh after all samples finish, or:  bash 04_summarize.sh
# Under SLURM the script runs from a spool copy, so PIPE_DIR is passed in.
source "${PIPE_DIR:-$(dirname "$0")}/scripts/common.sh"
activate_env
mkdir -p "$OUTDIR/summary"

python3 "$PIPE_DIR/scripts/qc_summary.py" "$OUTDIR" \
    CONTROL_UNMETH="$CONTROL_UNMETH" CONTROL_METH="$CONTROL_METH" \
    MAX_LAMBDA_CPG_METH="$MAX_LAMBDA_CPG_METH" MIN_PUC19_CPG_METH="$MIN_PUC19_CPG_METH" \
    MIN_MAPPING_RATE="$MIN_MAPPING_RATE" MAX_DUP_RATE="$MAX_DUP_RATE" \
    MAX_CHH_METH="$MAX_CHH_METH" MIN_CONTROL_CALLS="$MIN_CONTROL_CALLS" \
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
