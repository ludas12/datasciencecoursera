#!/usr/bin/env bash
# Submit the whole analysis to SLURM: reference (02) -> one array task per
# sample (03) -> summary (04), chained with job dependencies.
# Run on the login node after 00 and 01:  bash run_slurm.sh
# Under SLURM the script runs from a spool copy, so PIPE_DIR is passed in.
source "${PIPE_DIR:-$(dirname "$0")}/scripts/common.sh"

[[ -s "$SAMPLESHEET" ]] || die "no sample sheet; run 01_make_samplesheet.sh first"
n=$(n_samples); [[ $n -gt 0 ]] || die "sample sheet is empty"
mkdir -p "$LOG_DIR"

common=(--parsable --export=ALL,EMSEQ_CONFIG="$EMSEQ_CONFIG",PIPE_DIR="$PIPE_DIR" --chdir="$PIPE_DIR")
[[ -n "$SLURM_ACCOUNT" ]]   && common+=(--account="$SLURM_ACCOUNT")
[[ -n "$SLURM_PARTITION" ]] && common+=(--partition="$SLURM_PARTITION")
# shellcheck disable=SC2206
[[ -n "$SLURM_EXTRA" ]]     && common+=($SLURM_EXTRA)

ref_job=$(sbatch "${common[@]}" --job-name=emseq_ref \
    --cpus-per-task="$THREADS" --mem="$REF_MEM" --time="$REF_TIME" \
    --output="$LOG_DIR/02_reference_%j.log" 02_prepare_reference.sh)
echo "02 reference:  job $ref_job"

sample_job=$(sbatch "${common[@]}" --job-name=emseq_sample --dependency=afterok:"$ref_job" \
    --array=1-"$n"%"$MAX_PARALLEL_SAMPLES" \
    --cpus-per-task="$THREADS" --mem="$SAMPLE_MEM" --time="$SAMPLE_TIME" \
    --output="$LOG_DIR/03_sample_%a_%A.log" 03_process_sample.sh)
echo "03 samples:    job $sample_job (array 1-$n)"

# afterany: summarise even if some samples failed, so you can see which.
sum_job=$(sbatch "${common[@]}" --job-name=emseq_summary --dependency=afterany:"$sample_job" \
    --cpus-per-task=1 --mem="$SUMMARY_MEM" --time="$SUMMARY_TIME" \
    --output="$LOG_DIR/04_summary_%j.log" 04_summarize.sh)
echo "04 summary:    job $sum_job"
echo
echo "Monitor with: squeue -u \$USER     Logs: $LOG_DIR"
echo "Results: $OUTDIR/summary/qc_summary.txt and $OUTDIR/summary/multiqc/"
