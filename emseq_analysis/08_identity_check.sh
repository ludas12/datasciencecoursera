#!/usr/bin/env bash
# Step 8: sample identity check -- catches swapped or mislabelled libraries.
#   * sex from chrX/chrY read depth (MAPQ >= 30) relative to autosomes
#   * genotype concordance between all pairs of libraries on IDENTITY_REGIONS
#     (default chr20-22), using only reference A/T positions: EM-seq conversion
#     changes C/G only, so any non-reference base seen there is a real allele
# Libraries with the same `subject` in downstream/samples.tsv should match each
# other and nobody else. Needs samtools + bcftools in the downstream env:
#   conda install -n emseq_downstream -c bioconda bcftools
# Run with run_downstream.sh 08, or:  bash 08_identity_check.sh
# shellcheck disable=SC1091
source "${PIPE_DIR:-$(dirname "$0")}/scripts/common.sh"
activate_env "$DOWNSTREAM_ENV"
bcftools --version >/dev/null 2>&1 || die "bcftools missing or broken (e.g. an old build needing libcrypto.so.1.0.0). Fix: conda install -n $DOWNSTREAM_ENV -c conda-forge -c bioconda 'bcftools>=1.17' 'htslib>=1.17'"

I="$OUTDIR/summary/identity"
mkdir -p "$I/vcf"
THREADS="$DOWNSTREAM_THREADS"

# ---- libraries to check ------------------------------------------------------------
: > "$I/bams.txt"
printf 'label\tchrX\tchrY\tautosomes\n' > "$I/sex_counts.tsv"
while IFS=$'\t' read -r sample label _ _ include _; do
    [[ "$sample" == sample || "$include" != yes ]] && continue
    bam="$OUTDIR/samples/$sample/align/$sample.final.bam"
    [[ -s "$bam" ]] || { log "WARNING: no BAM for $sample; skipping"; continue; }
    printf '%s\t%s\n' "$bam" "$label" >> "$I/bams.txt"
    # primary, mapped, non-duplicate reads with MAPQ >= 30
    cnt() { samtools view -c -q 30 -F 3844 -@ "$THREADS" "$bam" "$@"; }
    autos=$(cut -f1 "$REF_FA.fai" | grep -E '^chr([1-9]|1[0-9]|2[0-2])$' | tr '\n' ' ')
    # shellcheck disable=SC2086
    printf '%s\t%s\t%s\t%s\n' "$label" "$(cnt chrX)" "$(cnt chrY)" "$(cnt $autos)" >> "$I/sex_counts.tsv"
done < "$SAMPLE_INFO"
n=$(wc -l < "$I/bams.txt")
[[ $n -ge 2 ]] || die "need at least 2 libraries"
log "sex counts done for $n libraries"

# ---- joint variant calling on the identity regions (one job per region) -----------
cut -f1 "$I/bams.txt" > "$I/bam_list.txt"
cut -f2 "$I/bams.txt" > "$I/labels.txt"
call_region() {
    local r="$1" out="$I/vcf/$1.vcf.gz"
    [[ -s "$out" ]] && return 0
    # per-library allele depths at A/T reference positions where any library shows
    # a non-reference base
    bcftools mpileup -f "$REF_FA" -r "$r" -b "$I/bam_list.txt" -q 30 -Q 30 -B -I \
        --max-depth 2000 -a FORMAT/AD,FORMAT/DP -Ou 2>/dev/null \
      | bcftools view -i '(REF="A" || REF="T") && ALT!="<*>"' -Oz -o "$out.part" \
      && mv "$out.part" "$out"
}
export -f call_region; export I REF_FA
regions=$(tr ',' '\n' <<< "$IDENTITY_REGIONS" | while read -r r; do
    if cut -f1 "$REF_FA.fai" | grep -qx "${r%%:*}"; then echo "$r"; fi; done)
[[ -n "$regions" ]] || die "none of IDENTITY_REGIONS ($IDENTITY_REGIONS) are in the reference"
log "collecting allele counts on $(echo $regions)"
rm -f "$I"/vcf/*.part
xargs -P "$THREADS" -I{} bash -c 'call_region {}' <<< "$regions"
bcftools concat -Oz -o "$I/identity.vcf.gz" "$I"/vcf/*.vcf.gz
bcftools reheader -s "$I/labels.txt" "$I/identity.vcf.gz" -o "$I/identity.labelled.vcf.gz"

python3 "$PIPE_DIR/downstream/identity_check.py" \
    --vcf "$I/identity.labelled.vcf.gz" --sex "$I/sex_counts.tsv" --fai "$REF_FA.fai" \
    --sample-info "$SAMPLE_INFO" --out "$I" | tee "$I/identity_summary.txt"
