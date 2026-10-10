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

With one sample per group, the p-values only measure read-sampling noise, not biological variation. Use them to rank regions, not as formal statistical evidence.

### 06: tissue / cell-of-origin (`OUTDIR/downstream/tissue_of_origin/`)
This step uses [UXM](https://github.com/nloyfer/UXM_deconv) with the human methylation atlas of [Loyfer et al. 2023, *Nature*](https://www.nature.com/articles/s41586-022-05580-6), which has 39 cell types and the `U250` marker set by default. It estimates what fraction of the DNA comes from each cell type, for example lung alveolar and bronchial epithelium, neutrophils, monocytes/macrophages, lymphocytes, endothelium and liver.

Outputs:
- `uxm_deconv.csv`: all cell types
- `uxm_grouped.tsv/.png`: grouped by compartment
- `uxm_deconv.pdf`

In a lung transplant recipient, plasma cfDNA from lung epithelium is expected to come largely from the donor lung. This method cannot tell donor DNA from recipient DNA; that needs genotype-based donor-derived cfDNA assays.

UXM and wgbstools are under the authors' academic/research licences. They are downloaded by the setup script, not redistributed here. If you use them, cite Loyfer et al. 2023.
