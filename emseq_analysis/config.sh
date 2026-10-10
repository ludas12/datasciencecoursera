# =============================================================================
# EM-seq analysis configuration -- edit this file, then run run_slurm.sh
# (or run_local.sh). Everything the pipeline writes goes under OUTDIR, so a
# fresh OUTDIR is a fresh analysis: nothing from earlier attempts is read.
# =============================================================================

# ---- Input ------------------------------------------------------------------
# Directory holding the NEW sequencing run's FASTQ files (searched recursively).
FASTQ_DIR="/path/to/new_run/fastq"
# Glob used to find FASTQs inside FASTQ_DIR.
FASTQ_GLOB="*.f*q.gz"

# ---- Output -----------------------------------------------------------------
# Must be a new/empty directory (the pipeline refuses to reuse a directory that
# holds results from a different config). Do NOT point this at old attempts.
OUTDIR="/path/to/emseq_fresh_run2"

# ---- Software ---------------------------------------------------------------
CONDA_ENV="emseq"
# Command(s) needed before `conda` is available on your HPC, e.g.
#   CONDA_SETUP="module load Anaconda3/2024.02"
# Leave empty if conda/mamba/micromamba is already on PATH in batch jobs.
CONDA_SETUP=""

# ---- Reference --------------------------------------------------------------
# Genome FASTA (plain or .gz). It is copied into OUTDIR/reference and indexed.
GENOME_FA="/path/to/genome.fa"
# EM-seq spike-in controls (unmethylated lambda, CpG-methylated pUC19) used to
# measure conversion efficiency. "download" fetches NEB's control FASTA from
# GitHub (needs internet on the node running 02); or give a local FASTA path;
# or "" to skip (no conversion QC -- not recommended).
CONTROLS_FA="download"
CONTROL_UNMETH="phage_lambda"     # contig name of the unmethylated control
CONTROL_METH="plasmid_puc19c"     # contig name of the CpG-methylated control

# ---- Aligner ----------------------------------------------------------------
# "bwameth" (fast, NEB's choice) or "bismark" (Bowtie2-based).
ALIGNER="bwameth"

# ---- Trimming (fastp) -------------------------------------------------------
# Adapters are auto-detected. EM-seq libraries have end-repair bias at read
# ends; NEB/Bismark recommend clipping ~10 bp from the 5' end of both reads.
CLIP_5P_R1=10
CLIP_5P_R2=10
CLIP_3P_R1=0
CLIP_3P_R2=0
MIN_READ_LEN=30
TRIM_POLY_G=true                  # true for NovaSeq/NextSeq 2-colour chemistry

# ---- Methylation calling (MethylDackel) -------------------------------------
MIN_MAPQ=10
MIN_BASEQ=5
# Extra per-strand inclusion bounds, e.g. "--nOT 0,0,0,5 --nOB 0,0,5,0".
# Leave empty first; check the M-bias plots / suggestions in qc/mbias and
# re-run step 03 with SKIP_TO=extract if needed.
MD_INCLUSION=""
# Also write per-strand methylKit-format files (an extra extraction pass;
# the default merged CpG bedGraph already works with methylKit/DSS in R).
MD_METHYLKIT=false

# ---- QC thresholds (used to flag samples in the summary table) --------------
MAX_LAMBDA_CPG_METH=1.0     # % methylated CpG on lambda (should be ~<0.5%)
MIN_PUC19_CPG_METH=90.0     # % methylated CpG on pUC19 (should be ~>95%)
MIN_MAPPING_RATE=70.0       # % reads mapped
MAX_DUP_RATE=40.0           # % duplicates

# ---- Compute ----------------------------------------------------------------
THREADS=16
# SLURM options (ignored by run_local.sh). Leave ACCOUNT/PARTITION empty if
# your cluster does not need them.
SLURM_ACCOUNT=""
SLURM_PARTITION=""
SLURM_EXTRA=""                    # e.g. "--qos=normal"
REF_MEM="64G";    REF_TIME="12:00:00"
SAMPLE_MEM="64G"; SAMPLE_TIME="48:00:00"
SUMMARY_MEM="8G"; SUMMARY_TIME="2:00:00"
MAX_PARALLEL_SAMPLES=8            # array throttle
