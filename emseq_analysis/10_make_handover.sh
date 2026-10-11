#!/usr/bin/env bash
# Step 10: assemble a handover folder for an independent QC review (no re-alignment).
# Collects pipeline code + versions, reference contents, per-library QC tables,
# spike-in and endogenous conversion metrics, M-bias files, identity and
# fragmentomics outputs, a 1 Mb methylation matrix, chr20 per-CpG files and all
# downstream results/figures into OUTDIR/handover_<date>/ (+ a .tar.gz).
# Run on the login node:  bash 10_make_handover.sh
# shellcheck disable=SC1091
source "${PIPE_DIR:-$(dirname "$0")}/scripts/common.sh"
activate_env
H="$OUTDIR/handover_$(date +%Y%m%d)"
mkdir -p "$H"/{A_provenance,C_alignment,D_conversion/mbias,E_identity,F_compact_data,G_downstream,work}

# ---- A. provenance --------------------------------------------------------------------
A="$H/A_provenance"
if git -C "$PIPE_DIR" rev-parse HEAD >/dev/null 2>&1; then
    git -C "$PIPE_DIR" rev-parse HEAD > "$A/pipeline_git_commit.txt"
    git -C "$PIPE_DIR" archive --format=tar.gz -o "$A/pipeline_scripts.tar.gz" HEAD -- .
else
    tar czf "$A/pipeline_scripts.tar.gz" -C "$PIPE_DIR" --exclude=.git .
fi
cp "$OUTDIR/config.used.sh" "$A/config.used.sh" 2>/dev/null || true
cp "$SAMPLESHEET" "$A/samplesheet.tsv"; cp "$SAMPLE_INFO" "$A/samples.tsv" 2>/dev/null || true
{
    for t in fastp bwameth.py bwa samtools MethylDackel multiqc python3; do
        printf '%-14s ' "$t"
        if [[ $t == bwa ]]; then { bwa 2>&1 || true; } | grep -m1 Version || echo NA
        else { "$t" --version 2>&1 || true; } | grep -m1 . || echo NA; fi
    done
} > "$A/tool_versions.txt"
if command -v conda >/dev/null; then
    conda env export -n "$CONDA_ENV" > "$A/env_${CONDA_ENV##*/}.yml" 2>/dev/null || true
    conda env export -n "$DOWNSTREAM_ENV" > "$A/env_${DOWNSTREAM_ENV##*/}.yml" 2>/dev/null || true
fi
{
    echo "Reference FASTA used for alignment: $REF_FA (built from GENOME_FA=$GENOME_FA)"
    echo "Contigs: $(wc -l < "$REF_FA.fai")"
    echo; echo "Non-standard contigs (not chr1-22/X/Y/M):"
    grep -vE '^chr([0-9]+|X|Y|M)\s' "$REF_FA.fai" | cut -f1,2 || true
    echo; for c in phage_lambda plasmid_puc19c phiX174 chrEBV chrM; do
        cut -f1 "$REF_FA.fai" | grep -qix "$c" && echo "$c: present" || echo "$c: absent"; done
} > "$A/reference_contents.txt"
grep -h "bwameth.py\|\[M::process\]" "$OUTDIR"/samples/*/align/*.bwameth.log 2>/dev/null | grep -m6 "CMD\|bwameth" > "$A/aligner_command_lines.txt" || true

# ---- C/D. per-library QC ------------------------------------------------------------------
for f in "$OUTDIR"/summary/qc_summary.{tsv,txt} "$OUTDIR"/summary/imprinted_promoters.tsv; do
    [[ -e "$f" ]] && cp "$f" "$H/C_alignment/"; done
cp -r "$OUTDIR/summary/multiqc" "$H/C_alignment/" 2>/dev/null || true
cp "$OUTDIR"/preseq*.txt "$OUTDIR"/dup_*.txt "$H/C_alignment/" 2>/dev/null || true
while IFS=$'\t' read -r sample _; do
    [[ "$sample" == sample ]] && continue
    q="$OUTDIR/samples/$sample/qc"; m="$OUTDIR/samples/$sample/mbias"
    mkdir -p "$H/C_alignment/per_library/$sample"
    cp "$q"/*.{flagstat,idxstats,samtools_stats,markdup,fastp}.* "$H/C_alignment/per_library/$sample/" 2>/dev/null || true
    cp "$m"/* "$H/D_conversion/mbias/" 2>/dev/null || true
    # chrM CpG methylation (mitochondrial DNA is ~unmethylated)
    bam="$OUTDIR/samples/$sample/align/$sample.final.bam"
    if [[ -s "$bam" ]] && cut -f1 "$REF_FA.fai" | grep -qx chrM; then
        MethylDackel extract -q "$MIN_MAPQ" -p "$MIN_BASEQ" -r chrM -o "$H/work/$sample.chrM" "$REF_FA" "$bam" 2>/dev/null \
          && awk 'NR>1{m+=$5;u+=$6} END{if(m+u) printf "%.3f\n", 100*m/(m+u)}' "$H/work/$sample.chrM_CpG.bedGraph" \
             > "$H/work/$sample.chrM.txt"
    fi
done < "$SAMPLESHEET"
python3 "$PIPE_DIR/scripts/handover_tables.py" "$OUTDIR" "$H/work" "$SAMPLE_INFO" "$CLIP_5P_R1" "$CLIP_5P_R2"
mv "$H/work/alignment_qc.csv" "$H/C_alignment/"
mv "$H/work/spikein_controls.csv" "$H/work/conversion_proxies.csv" "$H/D_conversion/"
mkdir -p "$H/F_compact_data/fragment_lengths"
mv "$H/work"/fragment_lengths_*.tsv "$H/F_compact_data/fragment_lengths/"
mv "$H/work/methylation_1Mb_bins.csv" "$H/work/chr20_cpg" "$H/F_compact_data/"

# ---- E/F/G. identity, fragmentomics, downstream ------------------------------------------
cp -r "$OUTDIR/summary/identity/"*.{tsv,png,txt} "$H/E_identity/" 2>/dev/null || true
FR="$DOWN_DIR/fragmentomics"
[[ -d "$FR" ]] && { mkdir -p "$H/F_compact_data/fragmentomics"; cp "$FR"/*.{tsv,png,txt} "$H/F_compact_data/fragmentomics/" 2>/dev/null || true; }
TO="$DOWN_DIR/tissue_of_origin"
[[ -d "$TO" ]] && { mkdir -p "$H/G_downstream/tissue_of_origin"; cp "$TO"/*.{csv,tsv,txt,png,pdf,bed} "$H/G_downstream/tissue_of_origin/" 2>/dev/null || true; }
RR="$DOWN_DIR/regions/results"
[[ -d "$RR" ]] && { mkdir -p "$H/G_downstream/regions"; cp "$RR"/*.{tsv,png} "$H/G_downstream/regions/" 2>/dev/null || true; }
rm -rf "$H/work"

# ---- manifest -------------------------------------------------------------------------------
( cd "$H" && find . -type f | sort | while read -r f; do printf '%8s  %s\n' "$(du -h "$f" | cut -f1)" "$f"; done ) > "$H/MANIFEST.txt"
tar czf "$H.tar.gz" -C "$(dirname "$H")" "$(basename "$H")"
log "handover folder: $H  ($(du -sh "$H.tar.gz" | cut -f1) as $H.tar.gz)"
