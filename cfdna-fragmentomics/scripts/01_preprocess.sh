#!/usr/bin/env bash
# Preprocess paired-end cfDNA FASTQs into an analysis-ready BAM.
#
#   fastp (adapter/quality trim) -> bwa-mem2 -> sort -> markdup -> index
#
# Usage:
#   scripts/01_preprocess.sh SAMPLE R1.fastq.gz R2.fastq.gz REF.fa OUTDIR [THREADS]
#
# REF.fa must be indexed for bwa-mem2 (bwa-mem2 index REF.fa) and samtools (samtools faidx REF.fa).
# Use the SAME reference build (e.g. hg38 analysis set) for every sample and for the
# fragment extraction step.
set -euo pipefail

if [[ $# -lt 5 ]]; then
    sed -n '2,12p' "$0"; exit 1
fi

SAMPLE=$1; R1=$2; R2=$3; REF=$4; OUT=$5; THREADS=${6:-8}
mkdir -p "$OUT"

# 1. Trim adapters. Keep reads >= 30 bp: very short cfDNA fragments are real signal,
#    so don't use an aggressive length filter here.
fastp -i "$R1" -I "$R2" \
      -o "$OUT/$SAMPLE.trim.R1.fq.gz" -O "$OUT/$SAMPLE.trim.R2.fq.gz" \
      --detect_adapter_for_pe --trim_poly_g --length_required 30 --thread "$THREADS" \
      --json "$OUT/$SAMPLE.fastp.json" --html "$OUT/$SAMPLE.fastp.html"

# 2. Align, fixmate (adds MC/ms tags needed by markdup), sort.
bwa-mem2 mem -t "$THREADS" -R "@RG\tID:$SAMPLE\tSM:$SAMPLE\tPL:ILLUMINA" "$REF" \
         "$OUT/$SAMPLE.trim.R1.fq.gz" "$OUT/$SAMPLE.trim.R2.fq.gz" \
  | samtools fixmate -m -@ "$THREADS" - - \
  | samtools sort -@ "$THREADS" -o "$OUT/$SAMPLE.sorted.bam" -

# 3. Mark duplicates (do not remove; the extractor skips flagged duplicates).
samtools markdup -@ "$THREADS" -s -f "$OUT/$SAMPLE.markdup.stats" \
         "$OUT/$SAMPLE.sorted.bam" "$OUT/$SAMPLE.bam"
samtools index "$OUT/$SAMPLE.bam"
rm -f "$OUT/$SAMPLE.sorted.bam" "$OUT"/"$SAMPLE".trim.R?.fq.gz

# 4. QC.
samtools flagstat -@ "$THREADS" "$OUT/$SAMPLE.bam" > "$OUT/$SAMPLE.flagstat"
samtools stats -@ "$THREADS" "$OUT/$SAMPLE.bam" | grep '^SN' | cut -f2- > "$OUT/$SAMPLE.stats"
echo "Done: $OUT/$SAMPLE.bam"
