# Shared helpers, sourced by every step. Not meant to be run directly.
set -euo pipefail

PIPE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
export PIPE_DIR
EMSEQ_CONFIG="${EMSEQ_CONFIG:-$PIPE_DIR/config.sh}"
[[ -f "$EMSEQ_CONFIG" ]] || { echo "ERROR: config not found: $EMSEQ_CONFIG" >&2; exit 1; }
# shellcheck source=../config.sh
source "$EMSEQ_CONFIG"
export EMSEQ_CONFIG

REF_DIR="$OUTDIR/reference"
REF_FA="$REF_DIR/genome_with_controls.fa"
SAMPLESHEET="$OUTDIR/samplesheet.tsv"
LOG_DIR="$OUTDIR/logs"

log()  { echo "[$(date '+%F %T')] $*" >&2; }
die()  { log "ERROR: $*"; exit 1; }
done_marker() { echo "$1/.done.$2"; }
is_done()  { [[ -f "$(done_marker "$1" "$2")" ]]; }
mark_done() { touch "$(done_marker "$1" "$2")"; }

activate_env() {
    if [[ -n "${CONDA_SETUP:-}" ]]; then
        eval "$CONDA_SETUP"
    fi
    # conda activate scripts reference unset variables; relax -u while activating.
    set +u
    if [[ "${CONDA_DEFAULT_ENV:-}" == "$CONDA_ENV" ]]; then
        :
    elif command -v conda >/dev/null 2>&1; then
        source "$(conda info --base)/etc/profile.d/conda.sh"
        conda activate "$CONDA_ENV"
    elif command -v micromamba >/dev/null 2>&1; then
        eval "$(micromamba shell hook -s bash)"
        micromamba activate "$CONDA_ENV"
    else
        set -u
        die "neither conda nor micromamba found; set CONDA_SETUP in config.sh"
    fi
    set -u
}

# Number of samples in the sample sheet (excluding header).
n_samples() { tail -n +2 "$SAMPLESHEET" | grep -c . ; }

# Print the sample-sheet row (sample<TAB>r1<TAB>r2) for a 1-based index.
sample_row() { tail -n +2 "$SAMPLESHEET" | grep . | sed -n "${1}p"; }
