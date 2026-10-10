# Shared helpers, sourced by every step. Not meant to be run directly.
set -euo pipefail

PIPE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
export PIPE_DIR
EMSEQ_CONFIG="${EMSEQ_CONFIG:-$PIPE_DIR/config.sh}"
[[ -f "$EMSEQ_CONFIG" ]] || { echo "ERROR: config not found: $EMSEQ_CONFIG" >&2; exit 1; }
# shellcheck source=../config.sh
source "$EMSEQ_CONFIG"
export EMSEQ_CONFIG

# Settings added after the first release; override them in config.sh if needed.
: "${INTERNAL_CONV_REGION:=chr22}"   # region whose non-CpG (CHH) methylation estimates conversion
: "${MAX_CHH_METH:=1.0}"             # % CHH methylation above which conversion is flagged
: "${MIN_CONTROL_CALLS:=200}"        # spike-in CpG calls needed before judging lambda/pUC19
# Downstream analysis (steps 05/06)
: "${DOWNSTREAM_ENV:=emseq_downstream}"   # conda env made by downstream/setup_downstream.sh
: "${BIN_SIZE:=10000}"                    # genomic bin size for step 05
: "${MIN_REGION_CALLS:=10}"               # CpG calls a region needs in every sample (step 05)
: "${UXM_ATLAS:=Atlas.U250.l4.hg38.full.tsv}"  # Loyfer et al. 2023 atlas used by step 06
: "${DOWNSTREAM_MEM:=32G}"; : "${DOWNSTREAM_TIME:=12:00:00}"; : "${DOWNSTREAM_THREADS:=8}"

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
    local env="${1:-$CONDA_ENV}"
    if [[ -n "${CONDA_SETUP:-}" ]]; then
        eval "$CONDA_SETUP"
    fi
    # conda activate scripts reference unset variables; relax -u while activating.
    set +u
    if [[ "${CONDA_DEFAULT_ENV:-}" == "$env" ]]; then
        :
    elif command -v conda >/dev/null 2>&1; then
        source "$(conda info --base)/etc/profile.d/conda.sh"
        conda activate "$env"
    elif command -v micromamba >/dev/null 2>&1; then
        eval "$(micromamba shell hook -s bash)"
        micromamba activate "$env"
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

DOWN_DIR="$OUTDIR/downstream"
TOOLS_DIR="$OUTDIR/tools"
ANNOT_DIR="$OUTDIR/annotation"
: "${SAMPLE_INFO:=$PIPE_DIR/downstream/samples.tsv}"   # labels/types for steps 05-06
: "${COMPARISONS:=$PIPE_DIR/downstream/comparisons.tsv}"
