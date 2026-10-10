#!/usr/bin/env bash
# Submit the downstream steps to SLURM (after the main pipeline and
# downstream/setup_downstream.sh):  bash run_downstream.sh [05] [06]
# With no arguments both steps are submitted; they run independently.
# shellcheck disable=SC1091
source "${PIPE_DIR:-$(dirname "$0")}/scripts/common.sh"
mkdir -p "$LOG_DIR"
steps=("$@"); [[ ${#steps[@]} -eq 0 ]] && steps=(05 06)

common=(--parsable --export=ALL,EMSEQ_CONFIG="$EMSEQ_CONFIG",PIPE_DIR="$PIPE_DIR" --chdir="$PIPE_DIR"
        --cpus-per-task="$DOWNSTREAM_THREADS" --mem="$DOWNSTREAM_MEM" --time="$DOWNSTREAM_TIME")
[[ -n "$SLURM_ACCOUNT" ]]   && common+=(--account="$SLURM_ACCOUNT")
[[ -n "$SLURM_PARTITION" ]] && common+=(--partition="$SLURM_PARTITION")
# shellcheck disable=SC2206
[[ -n "$SLURM_EXTRA" ]]     && common+=($SLURM_EXTRA)

for step in "${steps[@]}"; do
    case "$step" in
        05) script=05_regions.sh; name=emseq_regions ;;
        06) script=06_tissue_of_origin.sh; name=emseq_uxm ;;
        *) die "unknown step '$step' (use 05 and/or 06)" ;;
    esac
    job=$(sbatch "${common[@]}" --job-name="$name" --output="$LOG_DIR/${script%.sh}_%j.log" "$script")
    echo "$step $script: job $job"
done
echo "Logs: $LOG_DIR   Results: $DOWN_DIR"
