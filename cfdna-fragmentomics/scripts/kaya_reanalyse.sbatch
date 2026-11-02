#!/bin/bash
#SBATCH --job-name=cfdna-reanalyse
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=12:00:00
#SBATCH --output=logs/%x_%j.out
# Recompute per-sample features (and TSS profiles if TSS_BED is set) from existing fragment
# files, then refresh cohort/. Use after updating the code; no re-alignment.
# Submit from the pipeline directory:  sbatch scripts/kaya_reanalyse.sbatch
# REEXTRACT=1 first rebuilds the fragment files from the BAMs (needed after the Nov 2026 change
# that stops requiring the proper-pair flag; ~10-20 min per sample). The old fragment files are
# moved to frags_proper_only/ so before and after can be compared.
set -euo pipefail

WORK=/group/sbms001/hludewick/Fragmentomics
OUT=${OUT:-$WORK/results}
REF=${REF:-/group/sbms001/hludewick/G4_data/hg38/hg38.fa}   # its .fai defines the bins
TSS_BED=${TSS_BED:-$WORK/ref/gencode_v46_tss.bed}
BIN=${BIN:-5000000}

source "$(conda info --base)/etc/profile.d/conda.sh"
DEFAULT_ENV=cfdna-frag
[[ -d /group/sbms001/hludewick/conda/cfdna-frag ]] && DEFAULT_ENV=/group/sbms001/hludewick/conda/cfdna-frag
conda activate "${CONDA_ENV:-$DEFAULT_ENV}"
FRAG="python fragmentomics/cfdna_frag.py"
ALIGN_REF=${ALIGN_REF:-$WORK/ref/hg38.fa}            # the FASTA the BAMs were aligned to
BLACKLIST=${BLACKLIST:-$WORK/ref/hg38-blacklist.v2.bed}
# GENESETS=path/to/genesets.tsv to use other gene sets (e.g. fragmentomics/genesets_hpa.tsv)

SHEET=${SHEET:-samples.tsv}
tail -n +2 "$SHEET" | while IFS=$'\t' read -r S _; do
    [[ -z "$S" ]] && continue
    echo "== $S"
    if [[ "${REEXTRACT:-0}" == 1 ]]; then
        mkdir -p "$OUT/frags_proper_only"
        F="$OUT/frags/$S.frags.tsv.gz"
        [[ -s "$F" && ! -e "$OUT/frags_proper_only/$S.frags.tsv.gz" ]] && mv "$F" "$OUT/frags_proper_only/"
        $FRAG extract -b "$OUT/bam/$S.bam" -r "$ALIGN_REF" --blacklist "$BLACKLIST" -o "$F"
    fi
    $FRAG features -f "$OUT/frags/$S.frags.tsv.gz" -r "$REF" -o "$OUT/features/$S" -n "$S" --bin-size "$BIN"
    if [[ -s "$TSS_BED" ]]; then
        $FRAG tss -f "$OUT/frags/$S.frags.tsv.gz" -t "$TSS_BED" -o "$OUT/features/$S" -r "$REF" ${GENESETS:+-g "$GENESETS"}
    fi
done
$FRAG cohort -i "$OUT"/features/* -s "$SHEET" -o "$OUT/cohort"
echo "Done: $OUT/cohort/"
