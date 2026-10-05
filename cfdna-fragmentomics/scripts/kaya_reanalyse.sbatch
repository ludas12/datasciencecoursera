#!/bin/bash
#SBATCH --job-name=cfdna-reanalyse
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=2:00:00
#SBATCH --output=logs/%x_%j.out
# Recompute per-sample features (and TSS profiles if TSS_BED is set) from existing fragment
# files, then refresh cohort/. Use after updating the code; no re-alignment.
# Submit from the pipeline directory:  sbatch scripts/kaya_reanalyse.sbatch
set -euo pipefail

WORK=/group/sbms001/hludewick/Fragmentomics
OUT=${OUT:-$WORK/results}
REF=${REF:-/group/sbms001/hludewick/G4_data/hg38/hg38.fa}   # its .fai defines the bins
TSS_BED=${TSS_BED:-$WORK/ref/gencode_v46_tss.bed}
BIN=${BIN:-5000000}

source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate "${CONDA_ENV:-cfdna-frag}"
FRAG="python fragmentomics/cfdna_frag.py"

tail -n +2 samples.tsv | while IFS=$'\t' read -r S _; do
    [[ -z "$S" ]] && continue
    echo "== $S"
    $FRAG features -f "$OUT/frags/$S.frags.tsv.gz" -r "$REF" -o "$OUT/features/$S" -n "$S" --bin-size "$BIN"
    if [[ -s "$TSS_BED" ]]; then
        $FRAG tss -f "$OUT/frags/$S.frags.tsv.gz" -t "$TSS_BED" -o "$OUT/features/$S" -r "$REF"
    fi
done
$FRAG cohort -i "$OUT"/features/* -s samples.tsv -o "$OUT/cohort"
echo "Done: $OUT/cohort/"
