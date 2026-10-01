#!/bin/bash
#SBATCH --job-name=cfdna-frag
#SBATCH --cpus-per-task=16
#SBATCH --mem=48G
#SBATCH --time=12:00:00
#SBATCH --output=logs/%x_%j.out
# Submit from the pipeline directory:  cd $WORK/cfdna-fragmentomics && sbatch scripts/kaya_run.sbatch
set -euo pipefail

WORK=/group/sbms001/hludewick/Fragmentomics
REF=$WORK/ref/hg38.fa                      # needs bwa-mem2 index + .fai
BLACKLIST=$WORK/ref/hg38-blacklist.v2.bed
OUT=$WORK/results
export UMI_LEN=8      # 8 bp in-line block at the start of R1 and R2 (see FastQC per-base content)

source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate "${CONDA_ENV:-cfdna-frag}"   # env name or path; override: sbatch --export=ALL,CONDA_ENV=fragtools ...
scripts/run_all.sh samples.tsv "$REF" "$BLACKLIST" "$OUT" "${SLURM_CPUS_PER_TASK:-16}"
