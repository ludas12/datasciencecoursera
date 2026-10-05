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
python fragmentomics/cfdna_frag.py cohort -i features/* -s samples.tsv -o cohort/

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

## This study: lung transplant plasma + BAL (6 samples, IDT xGen cfDNA & FFPE)

`samples.tsv` is filled in for the current run (HC7/HC8 healthy plasma; CF0020A and CF0007C
transplant plasma + BAL from the same patients). Replace `FASTQ_DIR` with the real path:

```bash
sed -i "s|FASTQ_DIR|/path/to/fastq|g" samples.tsv
```

**1. Check the FastQC reports first**
```bash
python scripts/00_fastqc_summary.py /path/to/fastqc > fastqc_summary.tsv
```
Look at `max_cycle_bias_1to10` and `first10_*`. A strong composition bias in the first few cycles
(especially on R2) means the library prep added non-templated bases at read starts. Read starts
are the fragment ends, so those bases would distort fragment lengths and end motifs.
Check the IDT xGen cfDNA & FFPE analysis guide for any recommended trimming, and trim with fastp
`--trim_front1/--trim_front2` if it applies. Also check the guide for where the UMIs are: if they
are in the index read, normal duplicate marking is fine at this depth.

**2. Run everything** (on Kaya, inside an `sbatch` job or `salloc` session):
```bash
#!/bin/bash
#SBATCH --job-name=cfdna-frag --cpus-per-task=16 --mem=48G --time=12:00:00
conda activate cfdna-frag
scripts/run_all.sh samples.tsv /ref/hg38.fa /ref/hg38-blacklist.v2.bed results 16
```
Outputs in `results/cohort/`: `size_distribution.png`, `pca.png`, `delfi_ratio_heatmap.png`,
`size_motif_summary.tsv` (one row per sample), `group_summary.tsv`, `alignment_qc.tsv`.

**What to expect / how to interpret**
- **n = 2 per group:** treat this as descriptive / pilot data. Show distributions and per-sample
  values; don't run classifiers or p-values across groups.
- **BAL vs plasma is a different matrix, not just a different group.** BAL cfDNA comes from local
  lung cells, immune cells and microbes. It often has more long fragments and a lower human
  mapping rate (check `alignment_qc.tsv`). Compare BAL to plasma within a patient
  (CF0020A_P vs CF0020A_BAL), not as two independent groups.
- **Donor-derived DNA:** if a donor and recipient are sex-mismatched (e.g. male donor, female
  recipient), `frac_chrY` in `size_motif_summary.tsv` estimates donor fraction:
  `donor_fraction ≈ frac_chrY(sample) / frac_chrY(male control)`. Calibrate with a known male
  sample processed the same way. HC7/HC8 are useful here if one of them is male.
- **ichorCNA** estimates *tumour* fraction; it is not meaningful for transplant samples.

## TSS coverage profiles (tissue-of-origin signal)
Active promoters are nucleosome-depleted, so cfDNA from cells expressing a gene shows a coverage
dip at its transcription start site. `cfdna_frag.py tss` aggregates fragment coverage (100–220 bp
fragments) over the TSSs of each gene set in `fragmentomics/genesets.tsv`: lung epithelium,
neutrophil, lymphocyte, megakaryocyte/erythroid, liver, housekeeping (positive control) and
inactive genes (negative control). Profiles are strand-oriented and normalised to the
±2–3 kb flanks; `central_coverage` is the mean over ±150 bp (lower = more open/active).

```bash
python fragmentomics/cfdna_frag.py tss -f frags/S1.frags.tsv.gz -t hg38_tss.bed -o features/S1
```
The TSS BED needs gene symbols in one column (detected automatically) and ideally strand.
On Kaya, `sbatch scripts/kaya_tss.sbatch` runs it for all samples of a finished run and refreshes
`cohort/` (`tss_profiles.png`, `tss_central_coverage.tsv`, `key_metrics.png`).

Robustness: TSS windows with mean coverage >5× or <0.1× the median window are dropped
(repeats, CNVs, unmappable), and each TSS's per-position coverage is capped at the pooled 99.9th
percentile so a single artefact locus cannot dominate a set (`n_tss_excluded` in
`tss_summary.tsv`). With `-r REF` each fragment is weighted by expected/observed density of its GC
content (expected from 500k random genomic fragments with the sample's length distribution);
weights are saved to `tss_gc_weights.tsv`. The Kaya scripts pass `-r` by default.

Windows: `ndr_coverage` (−200 to −20 bp, the nucleosome-depleted region) is the primary measure
and drives `key_metrics.png`; `central_coverage` (±150 bp) and `plus1_coverage` (+60 to +200 bp,
first downstream nucleosome) are also reported. Core-trimmed fragments (BAL) mark the +1
nucleosome sharply, which can cancel the dip in the central window. TSSs of the same gene within
500 bp are merged so genes with many transcripts are counted once per promoter.

Interpretation: compare each set with `inactive_control` (TSSs are GC-rich, so some dip can be
technical) and with `housekeeping`. At low coverage only set-level differences between samples are
meaningful; a lung-epithelium dip is expected in BAL but is usually too small to see in plasma.

## BAL-oriented size metrics
BAL cfDNA lacks the 167 bp peak and mixes nuclease-trimmed nucleosomal DNA with randomly cut DNA.
Three metrics in `summary.json` / `size_motif_summary.tsv` capture this:

| Metric | What it measures | Plasma | BAL |
|---|---|---|---|
| `peak167_prominence` | mean count at 163–171 bp / shoulders (140–149, 185–194 bp) | high (>2) | ~1 or below |
| `ladder_amp_80_150` | amplitude of a 10–11 bp sinusoid fitted to the detrended 80–150 bp histogram (0.10 = peaks 10% above trend) | low | high |
| `frac_400_600` | fraction of fragments 400–600 bp (structureless long tail) | low | high |

The ladder amplitude has a noise floor that falls with fragment count (~0.1–0.2 at 60k
fragments, roughly 10× lower at several million), so compare samples of similar depth.
Extra columns in `samples.tsv` (e.g. `patient`, `diagnosis`) are carried into the cohort tables.
To recompute features for a finished run after a code update: `sbatch scripts/kaya_reanalyse.sbatch`.

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
