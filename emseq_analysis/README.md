# EM-seq analysis: fresh run

A small, self-contained pipeline that runs a fresh EM-seq analysis of a new sequencing run on an HPC.
It uses your existing conda env (`emseq`) and needs no Nextflow.
All results go into one new output directory. Nothing from earlier attempts is read or reused, and the pipeline refuses to write into a directory that it did not create.

```
FASTQ ─► fastp (adapter + 5′ end-repair-bias trimming)
      ─► bwa-meth  (or Bismark) ─► duplicate marking/removal
      ─► samtools QC ─► MethylDackel M-bias ─► MethylDackel CpG calls
      ─► spike-in conversion QC (lambda / pUC19) ─► summary table + MultiQC
```

## Quick start (on the HPC)

```bash
git clone -b claude/nice-fermi-cp8li3 https://github.com/ludas12/datasciencecoursera.git
cd datasciencecoursera/emseq_analysis

nano config.sh                 # set FASTQ_DIR, OUTDIR (new!), GENOME_FA, SLURM_ACCOUNT ...
bash 00_check_env.sh           # checks tools in the emseq env, inputs, output dir
bash 01_make_samplesheet.sh    # finds/pairs FASTQs, writes OUTDIR/samplesheet.tsv -> review it
bash run_slurm.sh              # submits 02 -> 03 (array, one task per sample) -> 04
```

There's no SLURM on an interactive node or workstation, so use `bash run_local.sh` there instead of `run_slurm.sh`.

The genome is indexed once inside `OUTDIR/reference`. For a mammalian genome this takes a few hours, and the sample jobs wait for it automatically.

### Your conda env must contain
`fastp samtools methyldackel python>=3.8` plus the aligner's tools:
- bwa-meth: `bwameth bwa`
- Bismark: `bismark bowtie2`

`multiqc` is optional. `00_check_env.sh` tells you what is missing. To add a tool:
`conda install -n emseq -c conda-forge -c bioconda <tool>`

## Files

| File | What it does |
|---|---|
| `config.sh` | The only file you edit: paths, aligner, trimming, QC thresholds, SLURM resources |
| `00_check_env.sh` | Pre-flight check (run on the login node) |
| `01_make_samplesheet.sh` | Creates `OUTDIR` and `samplesheet.tsv`. It merges lanes (`_L001`, `_L002`, …), recognises the `_R1_001`, `_R1`, `_1` and `.1` naming styles, and skips I1/I2 index-read files |
| `02_prepare_reference.sh` | Genome + NEB spike-in controls (lambda, pUC19, T4, Xp12) → `samtools faidx` → aligner index |
| `03_process_sample.sh` | One sample: trim → align → dedup → QC → M-bias → methylation calls |
| `04_summarize.sh` | `summary/qc_summary.tsv` with PASS/WARN flags + MultiQC report |
| `run_slurm.sh` / `run_local.sh` | Run everything |

Each step writes `.done.*` markers, so re-running skips finished work. This lets you just resubmit after a timeout or failure.

## Outputs (`OUTDIR/`)

```
samplesheet.tsv  config.used.sh  logs/
reference/genome_with_controls.fa(+ index)
samples/<sample>/
    align/<sample>.final.bam(.bai)       deduplicated/marked, sorted, indexed
    qc/      fastp, flagstat, idxstats, samtools stats, dup reports
    mbias/   <sample>_OT.svg/_OB.svg, mbias.txt, mbias_suggestion.txt
    meth/    <sample>_CpG.bedGraph.gz    per-CpG (both strands merged):
                                         chr start end %meth n_meth n_unmeth
             <sample>.control_*_{CpG,CHG,CHH}.bedGraph   spike-in calls
summary/qc_summary.tsv, qc_summary.txt, multiqc/emseq_multiqc_report.html
```

## Judging whether the new data are good

Check `summary/qc_summary.txt` first:

| Metric | Good EM-seq | Problem it reveals |
|---|---|---|
| `lambda_pct_meth_CpG` / `conversion_efficiency_pct` | ≲ 0.5 % / ≥ 99.5 % | incomplete enzymatic conversion (APOBEC step); judged only with ≥ `MIN_CONTROL_CALLS` (200) lambda CpG calls |
| `genome_pct_meth_CHH` / `internal_conversion_pct` | ≲ 1 % / ≥ 99 % | same, measured on the sample's own DNA (CHH on `INTERNAL_CONV_REGION`, default chr22); works when spike-in reads are scarce |
| `puc19_pct_meth_CpG` | ≳ 95 % | 5mC not protected (TET2/oxidation step failed) |
| `pct_mapped` | ≳ 70–80 % (bwa-meth) | adapter dimers, contamination, wrong genome |
| `pct_duplicates` | low (depends on input / depth) | low library complexity / too little input DNA |
| `genome_pct_meth_CpG` | ~70–80 % (mammalian somatic tissue) | global conversion problems |
| `imprinted_pct_meth` | ~50 % (30–70 % accepted) | **5mC protection failure** (TET2 step) when there is no pUC19 spike-in. This is the median methylation of promoters with an imprinting DMR (SNRPN, KCNQ1OT1, PEG3, MEST, PEG10, NNAT, …; set in `IMPRINTED_GENES`), which are ~50 % methylated in every normal tissue. Per-gene values are in `summary/imprinted_promoters.tsv`. Genes need ≥ 20 calls, so at ~1× depth only a few genes contribute |
| `CpGs_cov10` | as many as possible | sequencing depth |
| M-bias plots | flat lines along the read | end-repair bias → trim more |

If the M-bias plots (`samples/*/mbias/*.svg`) show methylation rising or falling at the read ends, put the bounds into `MD_INCLUSION` in `config.sh`. Use the bounds MethylDackel suggests (`mbias_suggestion` column) or set your own, e.g. `--nOT 0,0,0,5 --nOB 0,0,5,0`. Then re-call methylation without re-aligning:

```bash
for i in $(seq 1 $(($(wc -l < OUTDIR/samplesheet.tsv) - 1))); do
    REDO=extract bash 03_process_sample.sh $i
done
bash 04_summarize.sh
```

(or submit the same via `sbatch --array`).

## Notes
- **Aligner:** `ALIGNER="bwameth"` (default) is fast and is what NEB uses for EM-seq. `ALIGNER="bismark"` uses Bismark/Bowtie2 alignment and deduplication. With either aligner, methylation is called with MethylDackel, so the outputs are identical in format.
- **Spike-in controls:** with `CONTROLS_FA="download"`, NEB's control FASTA is downloaded from GitHub, which needs internet on the node that runs step 02. If compute nodes are offline, run `bash 02_prepare_reference.sh` on a login/data-mover node, or download the file and set `CONTROLS_FA` to its path.
- **Library type:** EM-seq libraries are directional, which both aligners assume by default.
- **Downstream analysis:** `meth/*_CpG.bedGraph.gz` files load directly into R for differential methylation (methylKit, DSS, dmrseq) once the QC looks good.

## Downstream analysis (steps 05–06)

These steps run after the main pipeline. Sample labels, types and timepoints are set in `downstream/samples.tsv`, and the pairs to compare in `downstream/comparisons.tsv`. Edit both for a new project.

```bash
bash downstream/setup_downstream.sh   # once, on the LOGIN node (internet):
                                      #  conda env emseq_downstream, wgbstools + UXM, UCSC annotation
bash run_downstream.sh                # submits 05 and 06 (or: bash run_downstream.sh 05)
```

### 05: region-level methylation (`OUTDIR/downstream/regions/results/`)
At ~1–2× CpG depth, single CpGs are too noisy, so CpG calls are summed over larger regions: 10 kb bins, CpG islands (UCSC) and promoters (RefSeq TSS ± 1 kb). Only autosomes are used, and a region must have at least `MIN_REGION_CALLS` (10) calls in every sample. Outputs for each region set:
- methylation matrix
- per-sample summary
- correlation heatmap
- PCA
- distributions
- one table and scatter plot per comparison, with the difference (`group_b − group_a`), a two-proportion z-test and BH FDR

Comparisons only use regions that pass three filters:
- **Coverage:** at least `MIN_COMPARE_CALLS` (50) calls in both samples.
- **Blacklist:** no overlap with the ENCODE hg38 blacklist (problem and repeat regions).
- **No extreme coverage:** fewer than `MAX_COV_FOLD` (10) × the median calls, since higher counts typically mean collapsed repeats.

CpG calls are not independent: one read covers several CpGs, especially in CpG islands. The z-statistics are therefore corrected for over-dispersion. The factor φ is estimated from each comparison and reported in `*.comparisons_summary.tsv`, together with how many regions each filter excluded.

A genome-wide shift between two samples also inflates φ. Read that kind of shift from `mean_delta` and the distributions; the test then picks out regions that differ beyond the typical difference.

With one sample per group, the p-values still only measure technical noise, not biological variation. Use them to rank regions, not as formal statistical evidence.

### 06: tissue / cell-of-origin (`OUTDIR/downstream/tissue_of_origin/`)
This step uses [UXM](https://github.com/nloyfer/UXM_deconv) with the human methylation atlas of [Loyfer et al. 2023, *Nature*](https://www.nature.com/articles/s41586-022-05580-6), which has 39 cell types and the `U250` marker set by default. It estimates what fraction of the DNA comes from each cell type, for example lung alveolar and bronchial epithelium, neutrophils, monocytes/macrophages, lymphocytes, endothelium and liver.

Outputs:
- `uxm_deconv.csv`: all cell types
- `uxm_grouped.tsv/.png`: grouped by compartment
- `uxm_deconv.pdf`

In a lung transplant recipient, plasma cfDNA from lung epithelium is expected to come largely from the donor lung. This method cannot tell donor DNA from recipient DNA; that needs genotype-based donor-derived cfDNA assays.

UXM and wgbstools are under the authors' academic/research licences. They are downloaded by the setup script, not redistributed here. If you use them, cite Loyfer et al. 2023.

## SeqMonk export (step 07)

`07_export_seqmonk.sh` writes `OUTDIR/seqmonk/<label>.cov.gz`, one per sample. These are Bismark coverage files (1-based: chr, start, end, % methylated, methylated count, unmethylated count) for chromosomes 1–22, X, Y and M. Labels come from `downstream/samples.tsv`; any sample not listed there keeps its sample name.

The export runs automatically at the end of step 04. To run it on its own, use `bash 07_export_seqmonk.sh` or `bash run_downstream.sh 07`. Files that are already up to date are skipped.

In SeqMonk:
1. File → New Project → Homo sapiens GRCh38.
2. File → Import Data → **Bismark (cov)**.
3. Make probes: running windows, or a fixed number of CpGs per probe, suit low coverage.
4. Data → Quantitation → **Bisulphite methylation over features**, with a minimum call count per probe.

Give SeqMonk ≥ 8 GB of memory (Edit → Preferences).

## Capture panel draft (`panel/`)

`panel/build_capture_panel.py` drafts a targeted methylation capture panel for lung-derived DNA from the Loyfer et al. 2023 atlas: about 2,200 markers in about 0.9 Mb. See `panel/README.md`.

## Sample identity check (step 08)

`bash run_downstream.sh 08` (or `bash 08_identity_check.sh`) checks that every library belongs to the person it is labelled as. It catches swapped or mislabelled libraries, for example from mirrored i7/i5 indexes entered the wrong way round. Libraries from the same person must have the same `subject` in `downstream/samples.tsv`.

- **Sex:** chrX and chrY read density relative to autosomes (MAPQ ≥ 30).
  - Female: `rel_X` ≈ 1 and `rel_Y` ≈ 0. Male: `rel_X` ≈ 0.5 and `rel_Y` clearly > 0.
  - A small chrY signal in a female recipient can be DNA from a male donor lung.
- **Genotype relatedness between every pair of libraries,** on `IDENTITY_REGIONS` (default chr20–22).
  - Only reference A/T positions are used. EM-seq conversion changes C/G only, so any non-reference base seen there is a real allele.
  - Sites must be common variants among the other libraries.
  - Scale: about 0 for unrelated people (up to ~0.6 when there are few libraries), about 0.5 for first-degree relatives, and about 1 for the same person. Pairs above 0.8 count as the same person.
- **Flags:**
  - A library that doesn't match its own subject.
  - A library that matches a different subject.
  - Sex that differs within a subject.

Outputs are in `OUTDIR/summary/identity/`: `identity_check.tsv`, `identity_matrix.tsv` and `.png`.

This needs `bcftools` in the downstream env. `downstream/setup_downstream.sh` installs it, or install it with `conda install -n emseq_downstream -c bioconda bcftools`. I tested the step on simulated EM-seq libraries at ~1× and ~7× depth, including a deliberately swapped label.

## Fragmentomics (step 09)

`bash run_downstream.sh 09` (or `bash 09_fragmentomics.sh`) analyses cfDNA fragmentation from the same BAMs, with no extra sequencing. Results go to `OUTDIR/downstream/fragmentomics/`:

| Output | Content |
|---|---|
| `fragment_summary.tsv`, `fragment_sizes.png/.tsv` | Size distribution per library: median, mode, % < 150 bp, short (100–150) / long (151–220) ratio, % di-nucleosomal (250–400 bp), end-motif diversity, CCCA % |
| `short_long_bins.tsv`, `short_long_profile.png` | Genome-wide short/long ratio in 5 Mb bins (DELFI-style), GC-corrected, shown relative to healthy controls (`type` = healthy) or to all libraries |
| `tss_footprints.tsv/.png` | Fragment coverage ±2 kb around the TSSs of marker gene sets (lung epithelium, neutrophil, megakaryocyte, liver, housekeeping). A dip at the TSS means the promoter was open in the source cells. These sets have 9–24 genes, so at ~1× coverage the curves are noisy; they become informative with depth or capture |
| `end_motifs.tsv` | 4-mer frequencies at fragment ends, read from the reference, so EM-seq conversion doesn't affect them |

**Fragment ends.** Fragmentomics needs the true fragment ends. Since this change, `KEEP_FRAGMENT_ENDS=true` (the default) stops fastp from clipping the 5′ read ends. Instead, MethylDackel ignores those `CLIP_5P_R1/R2` bases at methylation calling (`--nOT/--nOB`), so methylation results are equivalent. Each library records what was clipped before alignment in `qc/<sample>.clip5p.txt`.

Libraries aligned earlier, such as the pilot, were clipped. Step 09 restores their ends from the `CLIP_5P` settings, so their sizes are correct as well.
