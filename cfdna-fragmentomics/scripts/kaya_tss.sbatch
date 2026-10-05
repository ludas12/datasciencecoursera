#!/bin/bash
#SBATCH --job-name=cfdna-tss
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=2:00:00
#SBATCH --output=logs/%x_%j.out
# TSS coverage profiles + refreshed cohort plots, using fragment files from a finished run.
# Submit from the pipeline directory:  sbatch scripts/kaya_tss.sbatch
set -euo pipefail

WORK=/group/sbms001/hludewick/Fragmentomics
OUT=${OUT:-$WORK/results}
TSS_BED=${TSS_BED:-/group/sbms001/hludewick/G4_data/hg38/hg38_tss.bed}

source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate "${CONDA_ENV:-cfdna-frag}"
FRAG="python fragmentomics/cfdna_frag.py"

tail -n +2 samples.tsv | while IFS=$'\t' read -r S _; do
    [[ -z "$S" ]] && continue
    echo "== $S"
    $FRAG tss -f "$OUT/frags/$S.frags.tsv.gz" -t "$TSS_BED" -o "$OUT/features/$S"
done
$FRAG cohort -i "$OUT"/features/* -s samples.tsv -o "$OUT/cohort"
echo "Done: $OUT/cohort/tss_profiles.png, key_metrics.png, tss_central_coverage.tsv"
