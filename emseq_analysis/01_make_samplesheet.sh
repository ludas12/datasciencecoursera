#!/usr/bin/env bash
# Step 1: create the fresh output directory and the sample sheet from the
# FASTQs of the new run. Run on the login node:  bash 01_make_samplesheet.sh
# Review OUTDIR/samplesheet.tsv afterwards (edit it to drop/rename samples).
# Under SLURM the script runs from a spool copy, so PIPE_DIR is passed in.
source "${PIPE_DIR:-$(dirname "$0")}/scripts/common.sh"

if [[ -e "$OUTDIR" && -n "$(ls -A "$OUTDIR" 2>/dev/null)" ]]; then
    if [[ ! -f "$OUTDIR/config.used.sh" ]]; then
        die "$OUTDIR already contains files not made by this pipeline. Set a new OUTDIR in config.sh for a fresh analysis."
    fi
    if [[ -f "$SAMPLESHEET" && "${FORCE:-0}" != 1 ]]; then
        die "$SAMPLESHEET already exists (keeping your edits). Re-run with FORCE=1 to regenerate it."
    fi
fi
mkdir -p "$OUTDIR" "$LOG_DIR"
cp "$EMSEQ_CONFIG" "$OUTDIR/config.used.sh"

mapfile -t fastqs < <(find -L "$FASTQ_DIR" -type f -name "$FASTQ_GLOB" | sort)
[[ ${#fastqs[@]} -gt 0 ]] || die "no files matching '$FASTQ_GLOB' in $FASTQ_DIR"
log "found ${#fastqs[@]} FASTQ files in $FASTQ_DIR"

status=0
python3 "$PIPE_DIR/scripts/make_samplesheet.py" "$SAMPLESHEET" "${fastqs[@]}" || status=$?
[[ $status -eq 1 ]] && exit 1
echo
echo "Review $SAMPLESHEET, then run: bash run_slurm.sh   (or bash run_local.sh)"
exit $status
