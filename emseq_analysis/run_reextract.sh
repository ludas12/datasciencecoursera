#!/usr/bin/env bash
# Re-call methylation for every sample without re-aligning (e.g. after setting
# MD_INCLUSION in config.sh from the M-bias suggestions), then rebuild what uses
# the calls: 04 summary (+ 07 SeqMonk export), 05 regions, 11 BAL-vs-plasma.
# Earlier 04/05/11 results are kept as *_prev_<date> for comparison.
# Run on the login node:  bash run_reextract.sh
source "${PIPE_DIR:-$(dirname "$0")}/scripts/common.sh"

[[ -n "$MD_INCLUSION" ]] || die "MD_INCLUSION is empty in config.sh; re-calling with unchanged settings changes nothing"
n=$(n_samples); [[ $n -gt 0 ]] || die "sample sheet is empty"
mkdir -p "$LOG_DIR"
stamp=$(date +%Y%m%d_%H%M)

# keep the previous results for a before/after comparison
[[ -s "$OUTDIR/summary/qc_summary.txt" ]] && cp "$OUTDIR/summary/qc_summary.txt" "$OUTDIR/summary/qc_summary_prev_$stamp.txt"
if [[ -d "$DOWN_DIR/regions/mapped" ]]; then
    mkdir -p "$DOWN_DIR/regions_prev_$stamp"
    mv "$DOWN_DIR/regions/mapped" "$DOWN_DIR/regions_prev_$stamp/"
    [[ -d "$DOWN_DIR/regions/results" ]] && mv "$DOWN_DIR/regions/results" "$DOWN_DIR/regions_prev_$stamp/"
fi
[[ -d "$DOWN_DIR/compartment" ]] && mv "$DOWN_DIR/compartment" "$DOWN_DIR/compartment_prev_$stamp"

export REDO="extract internal_conv"
common=(--parsable --export=ALL,EMSEQ_CONFIG="$EMSEQ_CONFIG",PIPE_DIR="$PIPE_DIR" --chdir="$PIPE_DIR")
[[ -n "$SLURM_ACCOUNT" ]]   && common+=(--account="$SLURM_ACCOUNT")
[[ -n "$SLURM_PARTITION" ]] && common+=(--partition="$SLURM_PARTITION")
# shellcheck disable=SC2206
[[ -n "$SLURM_EXTRA" ]]     && common+=($SLURM_EXTRA)
down=(--cpus-per-task="$DOWNSTREAM_THREADS" --mem="$DOWNSTREAM_MEM" --time="$DOWNSTREAM_TIME")

log "re-calling methylation with: $MD_INCLUSION"
j3=$(sbatch "${common[@]}" --job-name=emseq_reextract --array=1-"$n"%"$MAX_PARALLEL_SAMPLES" \
    --cpus-per-task="$THREADS" --mem="$SAMPLE_MEM" --time="$SAMPLE_TIME" \
    --output="$LOG_DIR/03_reextract_%a_%A.log" 03_process_sample.sh)
j4=$(sbatch "${common[@]}" --job-name=emseq_summary --dependency=afterok:"$j3" \
    --cpus-per-task=1 --mem="$SUMMARY_MEM" --time="$SUMMARY_TIME" \
    --output="$LOG_DIR/04_summary_%j.log" 04_summarize.sh)
j5=$(sbatch "${common[@]}" "${down[@]}" --job-name=emseq_regions --dependency=afterok:"$j4" \
    --output="$LOG_DIR/05_regions_%j.log" 05_regions.sh)
j11=$(sbatch "${common[@]}" "${down[@]}" --job-name=emseq_compartment --dependency=afterok:"$j5" \
    --output="$LOG_DIR/11_compartment_check_%j.log" 11_compartment_check.sh)
echo "03 re-extract: $j3 (array 1-$n)   04: $j4   05: $j5   11: $j11"
echo "Previous results: summary/qc_summary_prev_$stamp.txt, downstream/*_prev_$stamp"
echo "Monitor with: squeue -u \$USER     Logs: $LOG_DIR"
