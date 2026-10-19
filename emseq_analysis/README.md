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

If the M-bias plots (`samples/*/mbias/*.svg`) show methylation rising or falling at the read ends, put inclusion bounds into `MD_INCLUSION` in `config.sh`. Use the bounds MethylDackel suggests (`mbias_suggestion` column), ideally one common setting for all libraries so they stay comparable, e.g. `--OT 3,0,12,40 --OB 3,0,12,40` (read 1 start, read 1 end, read 2 start, read 2 end; 0 = no bound). Then re-call methylation without re-aligning:

```bash
bash run_reextract.sh
```

This submits the re-extraction for every sample, then steps 04 (with the SeqMonk export), 05 and 11. The previous summary and step 05/11 results are kept as `*_prev_<date>` for a before/after comparison.

Bounds are positions in the aligned reads. Libraries aligned with `KEEP_FRAGMENT_ENDS=true` keep 10 more bases at the 5′ end than the pilot, so check their own M-bias before reusing pilot bounds.

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
| `tss_summary.tsv`, `tss_ndr_relative.tsv`, `tss_footprints.png`, `tss_profiles.tsv.gz` | Promoter nucleosome footprints, using the same method as the cfdna-fragmentomics pipeline. Fragments are 100–220 bp. Each fragment is GC-weighted (expected/observed fragment-GC density). There is one TSS per promoter (same-gene starts within 500 bp are merged). Outlier windows are dropped and spikes capped. Coverage is normalised to the 2–3 kb flanks. The primary measure is NDR coverage (−200 to −20 bp) relative to the library's own `inactive_control` set (< 1 means more open than background). Gene sets are in `downstream/genesets.tsv` (curated). For the larger Human Protein Atlas sets, run `python3 downstream/build_genesets_hpa.py` once with internet access; step 09 then uses `downstream/genesets_hpa.tsv` automatically, or whatever `FRAG_GENESETS` points to |
| `end_motifs.tsv` | 4-mer frequencies at fragment ends, read from the reference, so EM-seq conversion doesn't affect them |

**Fragment ends.** Fragmentomics needs the true fragment ends. Since this change, `KEEP_FRAGMENT_ENDS=true` (the default) stops fastp from clipping the 5′ read ends. Instead, MethylDackel ignores those `CLIP_5P_R1/R2` bases at methylation calling (`--nOT/--nOB`), so methylation results are equivalent. Each library records what was clipped before alignment in `qc/<sample>.clip5p.txt`.

Libraries aligned earlier, such as the pilot, were clipped. Step 09 restores their ends from the `CLIP_5P` settings, so their sizes are correct as well.

## BAL vs plasma checks (step 11)

`bash run_downstream.sh 11` (or `bash 11_compartment_check.sh`) needs steps 04 and 05. Pairs are one plasma and one BAL library with the same `subject` and `timepoint` in `downstream/samples.tsv`. Results go to `OUTDIR/downstream/compartment/`; `compartment_summary.txt` has the headline numbers.

| Output | Content |
|---|---|
| `imprinted_genes.tsv` | Per imprinted promoter and library: methylation, calls, CpGs covered, share of calls on the top 3 CpGs |
| `imprinted_pairs.tsv`, `imprinted_pairs.png` | Per pair and gene: BAL − plasma difference, 95 % CI from resampling CpGs, and whether leaving out any single CpG flips the sign. Shows whether a difference rests on a few positions |
| `imprinted_pairs_summary.tsv` | Per pair (each timepoint separately): genes lower/higher in BAL, two-sided sign test, genes whose CI excludes 0 |
| `bal_vs_plasma_fit.tsv`, `bal_vs_plasma_strata.tsv`, `bal_vs_plasma.png` | BAL against plasma in 10 kb bins, CpG islands and promoters. Regions are grouped by the reference libraries' methylation (independent of both samples, so the grouping doesn't bias the comparison), and calls are pooled per group. Uniform 5mC loss in BAL predicts intercept ≈ 0 and the same BAL/plasma ratio in highly (reference ≥ 0.7) and moderately (0.2–0.7) methylated regions. A cell-composition difference gives a non-zero intercept or ratios that differ |
| `size_split_fit.tsv`, `size_split_strata.tsv` | The same comparison by methylation level for fragments shorter than `FRAG_SPLIT_BP` (default 150 bp) and longer, in 10 kb bins: plasma short / plasma long, BAL long / plasma long, BAL short / plasma long, BAL short / BAL long. BAL long ≈ plasma long means the BAL deficit is carried by its short fragments; a deficit in long fragments too points to another source |
| `size_classes.tsv` | Finer size classes (edges in `FRAG_SIZE_EDGES`, default 100 150 220 400 bp): each class of plasma and BAL against plasma long fragments, with its share of CpG calls. Mono-nucleosomal fragments are 150–220 bp, di-nucleosomal 220–400 bp; DNA above 400 bp usually comes from lysed cells rather than apoptosis. Patient BAMs are split by fragment length (absolute TLEN plus any bases clipped before alignment) and methylation is re-called per class with the step 03 settings, so this takes about as long as one methylation call per library; it is skipped when up to date (`FORCE=1` recomputes) |
| `size_profile.tsv` | Within every library, each size class against the library's own 150–220 bp class. The reference DNAs are included: they are sheared genomic DNA with spike-in-verified conversion and protection, so a flat profile there means fragment-length effects in plasma or BAL are not technical. A reference is stratified by the other reference only, so its own reads don't define the strata |

By default MethylDackel counts only read pairs the aligner flagged as properly paired. bwa-meth sets that flag from each library's own insert-size distribution, so in plasma (tightly centred on ~167 bp) most fragments longer than ~300 bp lose it and the long size classes are a selected subset. `FRAG_ALL_PAIRS=true bash run_downstream.sh 11` repeats the size analysis counting all pairs up to `FRAG_MAX_BP` (default 1000 bp; MethylDackel `--keepDiscordant`) and writes to `compartment/all_pairs/`, leaving the default results in place.

The imprinted regions are RefSeq promoters (TSS ± 1 kb) of the genes in `IMPRINTED_GENES`, not curated imprinting-control coordinates. With one library per compartment, all of this describes one patient's samples, not BAL versus plasma in general.

## Comparison with the xGen fragmentomics pipeline

`downstream/compare_xgen.py` compares step 09 with the xGen `cfdna-fragmentomics` pipeline for the same cfDNA samples. Size and end-motif metrics are recomputed from both pipelines' raw tables with the same code and size range (50–600 bp), so differences reflect the library methods. Pairs are listed in `downstream/xgen_pairs.tsv`; xGen samples not in a pair (healthy plasma) are reported as context.

```bash
conda activate emseq_downstream
python3 downstream/compare_xgen.py \
    --emseq OUTDIR/downstream/fragmentomics \
    --xgen /group/sbms001/hludewick/Fragmentomics/cfdna-fragmentomics/results \
    --pairs downstream/xgen_pairs.tsv --out OUTDIR/downstream/xgen_comparison
```

Outputs: `metrics.tsv` (sizes, end motifs per library), `agreement.tsv` (per pair: correlation of size distributions, of the 256 end-motif frequencies and of 5 Mb short/long ratios), `tss_ndr.tsv` (promoter NDR per gene set, both methods) and `size_overlay.png`.

## Handover for independent QC review (step 10)

Run `bash 10_make_handover.sh` on the login node. It needs no re-alignment and builds `OUTDIR/handover_<date>/` plus a `.tar.gz`:

| Folder | Contents |
|---|---|
| `A_provenance/` | Pipeline code at its git commit, config used, conda env exports, tool versions, reference contents (lambda, pUC19, phiX, chrEBV, chrM present?), sample sheets |
| `C_alignment/` | `alignment_qc.csv` (one row per library: read pairs, mapping, duplicate/optical/PCR %, unique pairs, approximate coverage, median insert raw and clip-corrected, % chrM, % unmapped), MultiQC, per-library samtools/fastp reports, Preseq and optical-duplicate outputs |
| `D_conversion/` | `spikein_controls.csv`, `conversion_proxies.csv` (global mCpG, CHH on chr22, chrM mCpG, imprinted promoters, CpG coverage counts), M-bias files |
| `E_identity/` | Step 08 outputs |
| `F_compact_data/` | 1 Mb binned methylation matrix, chr20 per-CpG bedGraphs, fragment-length tables, step 09 outputs |
| `G_downstream/` | Tissue-of-origin and region results with figures |

`MANIFEST.txt` lists every file.
