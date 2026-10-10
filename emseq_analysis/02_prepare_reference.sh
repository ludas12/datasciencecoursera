#!/usr/bin/env bash
# Step 2: build genome + spike-in control reference and the aligner index.
# Submitted by run_slurm.sh, or run directly:  bash 02_prepare_reference.sh
# Under SLURM the script runs from a spool copy, so PIPE_DIR is passed in.
source "${PIPE_DIR:-$(dirname "$0")}/scripts/common.sh"
activate_env

CONTROLS_URL="https://raw.githubusercontent.com/nebiolabs/EM-seq/master/assets/methylation_controls.fa"
mkdir -p "$REF_DIR/src"

if ! is_done "$REF_DIR" fasta; then
    log "building $REF_FA"
    zcat -f "$GENOME_FA" > "$REF_FA.tmp"
    if [[ -n "$CONTROLS_FA" ]]; then
        ctrl="$CONTROLS_FA"
        if [[ "$CONTROLS_FA" == "download" ]]; then
            ctrl="$REF_DIR/src/methylation_controls.fa"
            [[ -s "$ctrl" ]] || curl -fsSL "$CONTROLS_URL" -o "$ctrl" \
                || die "could not download controls; download $CONTROLS_URL yourself and set CONTROLS_FA to it"
        fi
        # Append each control contig unless the genome already contains it.
        grep '^>' "$REF_FA.tmp" | sed 's/^>//; s/[[:space:]].*//' > "$REF_DIR/src/genome_contigs.txt"
        zcat -f "$ctrl" > "$REF_DIR/src/controls.tmp.fa"
        awk 'FNR == NR { seen[$1] = 1; next }
             /^>/ { name = substr($1, 2); keep = !(name in seen)
                    if (!keep) print "  control " name " already in genome, not appended" > "/dev/stderr" }
             keep' "$REF_DIR/src/genome_contigs.txt" "$REF_DIR/src/controls.tmp.fa" >> "$REF_FA.tmp"
        rm -f "$REF_DIR/src/controls.tmp.fa"
    fi
    mv "$REF_FA.tmp" "$REF_FA"
    samtools faidx "$REF_FA"
    mark_done "$REF_DIR" fasta
fi

for ctrl_name in "$CONTROL_UNMETH" "$CONTROL_METH"; do
    [[ -z "$ctrl_name" ]] && continue
    cut -f1 "$REF_FA.fai" | grep -qx "$ctrl_name" \
        || log "WARNING: control contig '$ctrl_name' not in reference; conversion QC for it will be empty"
done

if ! is_done "$REF_DIR" "index_$ALIGNER"; then
    log "indexing reference for $ALIGNER (this can take a few hours for a mammalian genome)"
    case "$ALIGNER" in
        bwameth) bwameth.py index "$REF_FA" ;;
        bismark) bismark_genome_preparation --parallel "$(( THREADS > 1 ? THREADS / 2 : 1 ))" "$REF_DIR" ;;
        *) die "unknown ALIGNER '$ALIGNER'" ;;
    esac
    mark_done "$REF_DIR" "index_$ALIGNER"
fi
log "reference ready: $REF_FA"
