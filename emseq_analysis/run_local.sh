#!/usr/bin/env bash
# Run every step sequentially on the current machine (interactive node or
# workstation; no scheduler).  bash run_local.sh
# Under SLURM the script runs from a spool copy, so PIPE_DIR is passed in.
source "${PIPE_DIR:-$(dirname "$0")}/scripts/common.sh"
cd "$PIPE_DIR"
[[ -s "$SAMPLESHEET" ]] || die "no sample sheet; run 01_make_samplesheet.sh first"
mkdir -p "$LOG_DIR"

bash 02_prepare_reference.sh 2>&1 | tee "$LOG_DIR/02_reference.log"
failed=()
for i in $(seq 1 "$(n_samples)"); do
    bash 03_process_sample.sh "$i" > "$LOG_DIR/03_sample_$i.log" 2>&1 \
        && log "sample $i done" || { log "sample $i FAILED (see $LOG_DIR/03_sample_$i.log)"; failed+=("$i"); }
done
bash 04_summarize.sh 2>&1 | tee "$LOG_DIR/04_summary.log"
[[ ${#failed[@]} -eq 0 ]] || die "failed samples: ${failed[*]}"
