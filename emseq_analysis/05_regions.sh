#!/usr/bin/env bash
# Step 5: region-level methylation for all samples -- genomic bins, CpG islands
# and promoters -- plus correlation, PCA and the comparisons in
# downstream/comparisons.tsv. Low-coverage friendly: counts are summed per region.
# Needs downstream/setup_downstream.sh first. Submit with run_downstream.sh or:
#   bash 05_regions.sh
# shellcheck disable=SC1091
source "${PIPE_DIR:-$(dirname "$0")}/scripts/common.sh"
activate_env "$DOWNSTREAM_ENV"

R="$DOWN_DIR/regions"
B="$R/beds"
mkdir -p "$B" "$R/mapped"
THREADS="$DOWNSTREAM_THREADS"
sort_bed() { LC_ALL=C sort -k1,1 -k2,2n -S 2G --parallel="$THREADS" -T "$R" "$@"; }

# ---- region sets (autosomes only, so sex differences don't drive results) -----
grep -E '^chr([1-9]|1[0-9]|2[0-2])\s' "$REF_FA.fai" | cut -f1,2 > "$B/autosomes.genome"
if [[ ! -s "$B/bins.bed" ]]; then
    log "making ${BIN_SIZE} bp bins"
    bedtools makewindows -g "$B/autosomes.genome" -w "$BIN_SIZE" \
        | awk -v OFS='\t' '{print $1, $2, $3, $1":"$2"-"$3}' | sort_bed > "$B/bins.bed"
fi
if [[ ! -s "$B/cpg_islands.bed" && -s "$ANNOT_DIR/cpgIslandExt.txt.gz" ]]; then
    log "CpG islands"
    zcat "$ANNOT_DIR/cpgIslandExt.txt.gz" \
        | awk -v OFS='\t' 'NR==FNR {keep[$1]=1; next} ($2 in keep) {print $2, $3, $4, $2":"$3"-"$4}' \
            "$B/autosomes.genome" - | sort_bed -u > "$B/cpg_islands.bed"
fi
if [[ ! -s "$B/promoters.bed" && -s "$ANNOT_DIR/ncbiRefSeqCurated.txt.gz" ]]; then
    log "promoters (TSS +/- 1 kb, RefSeq curated)"
    zcat "$ANNOT_DIR/ncbiRefSeqCurated.txt.gz" \
        | awk -v OFS='\t' 'NR==FNR {keep[$1]=1; next} ($3 in keep) {
              tss = ($4 == "+") ? $5 : $6; s = tss - 1000; if (s < 0) s = 0
              print $3, s, tss + 1000, $13 }' "$B/autosomes.genome" - \
        | sort_bed | awk -v OFS='\t' '{k = $1":"$2"-"$3} !(k in seen) {seen[k]=1; print}' > "$B/promoters.bed"
fi
rtypes=()
for r in bins cpg_islands promoters; do [[ -s "$B/$r.bed" ]] && rtypes+=("$r"); done
log "region sets: ${rtypes[*]}"

# ---- per-sample counts per region ------------------------------------------------
while IFS=$'\t' read -r sample label _ _ include; do
    [[ "$sample" == sample || "$include" != yes ]] && continue
    bg="$OUTDIR/samples/$sample/meth/${sample}_CpG.bedGraph.gz"
    [[ -s "$bg" ]] || { log "WARNING: no CpG calls for $sample; skipping"; continue; }
    todo=()
    for r in "${rtypes[@]}"; do [[ -s "$R/mapped/$r/$label.tsv" ]] || todo+=("$r"); done
    [[ ${#todo[@]} -eq 0 ]] && continue
    log "$label: summing CpG calls over ${todo[*]}"
    tmp="$R/$label.cpg.tmp.bed"
    zcat "$bg" | awk -v OFS='\t' 'NR==FNR {keep[$1]=1; next} ($1 in keep) {print $1, $2, $3, $5, $6}' \
        "$B/autosomes.genome" - | sort_bed > "$tmp"
    for r in "${todo[@]}"; do
        mkdir -p "$R/mapped/$r"
        bedtools map -a "$B/$r.bed" -b "$tmp" -c 4,5 -o sum,sum > "$R/mapped/$r/$label.tsv.part"
        mv "$R/mapped/$r/$label.tsv.part" "$R/mapped/$r/$label.tsv"
    done
    rm -f "$tmp"
done < "$SAMPLE_INFO"

# ---- statistics and plots --------------------------------------------------------------
python3 "$PIPE_DIR/downstream/region_analysis.py" \
    --regions-dir "$R" --sample-info "$SAMPLE_INFO" \
    --comparisons "$COMPARISONS" \
    --min-calls "$MIN_REGION_CALLS" --out "$R/results"
log "region analysis written to $R/results"
