#!/usr/bin/env bash
# Step 0: check the conda env, inputs and config before anything is submitted.
# Run on the login node:  bash 00_check_env.sh
# Under SLURM the script runs from a spool copy, so PIPE_DIR is passed in.
source "${PIPE_DIR:-$(dirname "$0")}/scripts/common.sh"
activate_env

fail=0
need() {
    if command -v "$1" >/dev/null 2>&1; then
        local v
        if [[ "$1" == bwa ]]; then v=$(bwa 2>&1 | grep -m1 Version || true)
        else v=$("$1" --version 2>&1 | grep -m1 . || true); fi
        printf "  %-16s OK   %s\n" "$1" "$v"
    else
        printf "  %-16s MISSING\n" "$1"; fail=1
    fi
}

echo "Conda env: ${CONDA_DEFAULT_ENV:-?}  (python: $(command -v python3 || echo none))"
echo "Required tools:"
for t in fastp samtools MethylDackel python3; do need "$t"; done
case "$ALIGNER" in
    bwameth) need bwameth.py; need bwa ;;
    bismark) need bismark; need bowtie2; need deduplicate_bismark; need bismark_genome_preparation ;;
    *) echo "  ALIGNER must be 'bwameth' or 'bismark' (got '$ALIGNER')"; fail=1 ;;
esac
echo "Optional tools:"
for t in multiqc fastqc; do
    command -v "$t" >/dev/null 2>&1 && printf "  %-16s OK\n" "$t" || printf "  %-16s not found (step 04 skips it)\n" "$t"
done

echo "Inputs:"
if [[ -d "$FASTQ_DIR" ]]; then
    n=$(find -L "$FASTQ_DIR" -type f -name "$FASTQ_GLOB" | wc -l)
    echo "  FASTQ_DIR        $FASTQ_DIR ($n files matching '$FASTQ_GLOB')"
    [[ $n -gt 0 ]] || fail=1
else
    echo "  FASTQ_DIR        NOT FOUND: $FASTQ_DIR"; fail=1
fi
if [[ -f "$GENOME_FA" ]]; then echo "  GENOME_FA        $GENOME_FA"
else echo "  GENOME_FA        NOT FOUND: $GENOME_FA"; fail=1; fi
if [[ "$CONTROLS_FA" == "download" ]]; then echo "  CONTROLS_FA      will be downloaded from NEB GitHub"
elif [[ -z "$CONTROLS_FA" ]]; then
    in_genome=""
    if [[ -n "$CONTROL_UNMETH" && -f "$GENOME_FA.fai" ]] && cut -f1 "$GENOME_FA.fai" | grep -qx "$CONTROL_UNMETH"; then
        in_genome=yes
    fi
    if [[ -n "$in_genome" ]]; then echo "  CONTROLS_FA      none needed: $CONTROL_UNMETH/$CONTROL_METH already in GENOME_FA"
    else echo "  CONTROLS_FA      none (conversion QC disabled unless the controls are in GENOME_FA)"; fi
elif [[ -f "$CONTROLS_FA" ]]; then echo "  CONTROLS_FA      $CONTROLS_FA"
else echo "  CONTROLS_FA      NOT FOUND: $CONTROLS_FA"; fail=1; fi

echo "Output:"
if [[ -e "$OUTDIR" && -n "$(ls -A "$OUTDIR" 2>/dev/null)" ]]; then
    if [[ -f "$OUTDIR/config.used.sh" ]]; then
        echo "  OUTDIR           $OUTDIR (existing run of this pipeline: finished steps are skipped)"
    else
        echo "  OUTDIR           $OUTDIR is NOT empty and was not created by this pipeline."
        echo "                   Choose a new OUTDIR for a fresh analysis."; fail=1
    fi
else
    echo "  OUTDIR           $OUTDIR (new)"
fi

if [[ $fail -ne 0 ]]; then echo; echo "Problems found -- fix them before running."; exit 1; fi
echo; echo "All checks passed. Next: bash 01_make_samplesheet.sh"
