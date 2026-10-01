#!/usr/bin/env bash
# Optional: estimate tumour fraction / copy number with ichorCNA (designed for ~0.1x ULP-WGS).
# Useful alongside fragmentomics: tumour fraction explains much of the fragment-size shift.
#
# Usage: scripts/02_ichorcna.sh SAMPLE SAMPLE.bam OUTDIR
#
# Requires HMMcopy's readCounter and an ichorCNA install (R package + its extdata).
# The bundled wig files below are for hg38 at 1 Mb; edit ICHOR_EXT for your install.
set -euo pipefail
SAMPLE=$1; BAM=$2; OUT=$3
mkdir -p "$OUT"
ICHOR_EXT=${ICHOR_EXT:-$(Rscript -e 'cat(system.file("extdata", package="ichorCNA"))')}
RUN_ICHOR=${RUN_ICHOR:-runIchorCNA.R}   # scripts/runIchorCNA.R from the ichorCNA repo

readCounter --window 1000000 --quality 20 \
    --chromosome "chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22,chrX,chrY" \
    "$BAM" > "$OUT/$SAMPLE.wig"

Rscript "$RUN_ICHOR" --id "$SAMPLE" \
    --WIG "$OUT/$SAMPLE.wig" --ploidy "c(2,3)" \
    --normal "c(0.5,0.6,0.7,0.8,0.9,0.95,0.99)" --maxCN 5 \
    --gcWig "$ICHOR_EXT/gc_hg38_1000kb.wig" \
    --mapWig "$ICHOR_EXT/map_hg38_1000kb.wig" \
    --centromere "$ICHOR_EXT/GRCh38.GCA_000001405.2_centromere_acen.txt" \
    --normalPanel "$ICHOR_EXT/HD_ULP_PoN_hg38_1Mb_median_normAutosome_median.rds" \
    --includeHOMD False --chrs "c(1:22, 'X')" --chrTrain "c(1:22)" \
    --estimateNormal True --estimatePloidy True --estimateScPrevalence True \
    --scStates "c(1,3)" --txnE 0.9999 --txnStrength 10000 \
    --genomeBuild hg38 --genomeStyle UCSC --outDir "$OUT"
