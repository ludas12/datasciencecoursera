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
#
# UMI_LEN (env, default 0): number of in-line bases at the start of BOTH reads that are not
# part of the cfDNA fragment (e.g. in-line UMIs). They are removed by fastp and kept in the read
# name. Read start = fragment end, so leaving them in shifts every fragment end and length.
# For the 2026-09-08 IDT xGen cfDNA & FFPE run, FastQC shows an 8 bp block: UMI_LEN=8.
set -euo pipefail

if [[ $# -lt 5 ]]; then
    sed -n '2,12p' "$0"; exit 1
fi

SAMPLE=$1; R1=$2; R2=$3; REF=$4; OUT=$5; THREADS=${6:-8}
UMI_LEN=${UMI_LEN:-0}
mkdir -p "$OUT"
UMI_OPTS=()
if (( UMI_LEN > 0 )); then
    UMI_OPTS=(--umi --umi_loc per_read --umi_len "$UMI_LEN")
fi

# 1. Trim adapters. Keep reads >= 30 bp: very short cfDNA fragments are real signal,
#    so don't use an aggressive length filter here.
fastp -i "$R1" -I "$R2" \
      -o "$OUT/$SAMPLE.trim.R1.fq.gz" -O "$OUT/$SAMPLE.trim.R2.fq.gz" \
      --detect_adapter_for_pe --trim_poly_g --length_required 30 --thread "$THREADS" \
      "${UMI_OPTS[@]}" \
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
# 5' soft-clipping check: if non-genomic bases remain at read starts, bwa soft-clips them and
# this rate is high. Expect a few % at most; >10% means UMI_LEN is probably wrong.
samtools view -F 0xF04 -q 30 "$OUT/$SAMPLE.bam" | head -n 2000000 | awk '
    { fwd = int($2 / 16) % 2 == 0; n++
      if ((fwd && $6 ~ /^[0-9]+S/) || (!fwd && $6 ~ /S$/)) clip++ }
    END { printf "reads\t%d\npct_5prime_softclipped\t%.2f\n", n, 100 * clip / n }' \
    > "$OUT/$SAMPLE.softclip.txt"
cat "$OUT/$SAMPLE.softclip.txt"
echo "Done: $OUT/$SAMPLE.bam"
