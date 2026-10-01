#!/bin/bash
#SBATCH --job-name=bwa-index
#SBATCH --cpus-per-task=4
#SBATCH --mem=80G
#SBATCH --time=6:00:00
#SBATCH --output=logs/%x_%j.out
# One-off: download hg38 analysis set + ENCODE blacklist and build indexes.
set -euo pipefail
REFDIR=/group/sbms001/hludewick/Fragmentomics/ref
mkdir -p "$REFDIR" && cd "$REFDIR"
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate cfdna-frag

[[ -s hg38.fa ]] || { wget -q https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/analysisSet/hg38.analysisSet.fa.gz
                      gunzip -c hg38.analysisSet.fa.gz > hg38.fa && rm hg38.analysisSet.fa.gz; }
[[ -s hg38-blacklist.v2.bed ]] || { wget -q https://github.com/Boyle-Lab/Blacklist/raw/master/lists/hg38-blacklist.v2.bed.gz
                                    gunzip hg38-blacklist.v2.bed.gz; }
samtools faidx hg38.fa
bwa-mem2 index hg38.fa
