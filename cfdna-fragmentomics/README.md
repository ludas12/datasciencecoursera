# cfDNA fragmentomics for low-coverage WGS

A pipeline for paired-end, shallow (~0.1–5x) cfDNA whole-genome sequencing. It takes FASTQs
and produces per-sample fragmentomic features plus cohort-level matrices ready for modelling.

```
FASTQ ─► 01_preprocess.sh ─► BAM ─► cfdna_frag.py extract ─► fragments.tsv.gz
                              │                                   │
                              └► 02_ichorcna.sh (tumour fraction) └► features ─► cohort
```

## Setup
```bash
conda env create -f environment.yml && conda activate cfdna-frag
# Reference: use one build everywhere, e.g. the hg38 analysis set (no alt contigs)
bwa-mem2 index hg38.fa && samtools faidx hg38.fa
# Blacklist: https://github.com/Boyle-Lab/Blacklist (hg38-blacklist.v2.bed.gz, gunzip it)
```

## Run
```bash
# 1. FASTQ -> duplicate-marked BAM
scripts/01_preprocess.sh S1 S1_R1.fq.gz S1_R2.fq.gz hg38.fa bam/ 16

# 2. BAM -> one row per fragment (MAPQ>=30, proper pairs, non-duplicate, autosomes, blacklist removed)
python fragmentomics/cfdna_frag.py extract -b bam/S1.bam -r hg38.fa \
    --blacklist hg38-blacklist.v2.bed -o frags/S1.frags.tsv.gz

# 3. Per-sample features
python fragmentomics/cfdna_frag.py features -f frags/S1.frags.tsv.gz -r hg38.fa -o features/S1

# 4. Merge all samples
python fragmentomics/cfdna_frag.py cohort -i features/* -o cohort/

# Optional: tumour fraction / copy number
scripts/02_ichorcna.sh S1 bam/S1.bam ichor/
```

## Features
| Feature | Output | Notes |
|---|---|---|
| Fragment size distribution | `size_hist.tsv`, `size_distribution.png` | Modal ~167 bp in healthy plasma; tumour DNA is shifted shorter |
| Size summary | `summary.json` | fraction <150 bp, 90–150 bp, short (100–150)/long (151–220) ratio, di-nucleosomal fraction, ~10 bp periodicity power |
| 5′ end motifs | `end_motifs.tsv` | 256 4-mer frequencies over both fragment ends + motif diversity score (MDS; Jiang et al. 2020) |
| DELFI-style profile | `bins.tsv`, `delfi_short_long_ratio.tsv` | short/long counts per 5 Mb bin, lowess GC-corrected (Cristiano et al. 2019) |
| Normalised coverage | `coverage_norm.tsv` | GC-corrected bin coverage (crude CNA signal; use ichorCNA for calls) |

`cohort/` contains sample × feature matrices (`size_motif_summary.tsv`, `delfi_short_long_ratio.tsv`,
`end_motif_freq.tsv`, `coverage_norm.tsv`) — only bins passing QC in every sample are kept.

## Low-coverage notes
- **Depth needed:** size distribution and end motifs are stable from ~1M fragments. DELFI ratios at
  5 Mb bins work from ~0.1x (~1–2M fragments → a few thousand per bin). Use `--bin-size 10000000` below that.
- **Not feasible at low depth per locus:** windowed protection score / nucleosome positioning at
  individual genes needs >30x. At low depth you can only do *aggregate* profiles (e.g. coverage
  summed over thousands of TSSs or TF binding sites — Griffin, Ulz et al.).
- **Batch effects dominate:** fragment sizes are sensitive to extraction kit, library prep
  (single- vs double-stranded), size selection and pre-analytics (tube type, time to spin).
  Process cases and controls together, and include healthy controls from the *same* workflow.
- **Tumour fraction:** most size/ratio differences scale with tumour fraction; report ichorCNA TF
  alongside fragment features, and be cautious with samples below ~3% TF.
- End motifs are read from the reference at the fragment ends; soft-clipped ends are not trimmed back.

## Modelling suggestion
With small cohorts, keep it simple: penalised logistic regression (or gradient boosting) on
PCA of `delfi_short_long_ratio.tsv` + `size_motif_summary.tsv` + `end_motif_freq.tsv`, evaluated with
repeated stratified cross-validation; standardise features inside the CV folds.

## Test
`tests/run_smoke_test.sh` simulates a small genome and three samples (one with shortened fragments)
and runs `extract → features → cohort` end to end.
