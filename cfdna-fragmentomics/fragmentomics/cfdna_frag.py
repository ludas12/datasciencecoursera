#!/usr/bin/env python3
"""cfDNA fragmentomics for low-coverage paired-end WGS.

Subcommands
-----------
extract   BAM -> per-fragment table (coords, length, GC, 5' end motifs of both ends)
features  fragment table -> size distribution, size summary, end-motif frequencies + MDS,
          DELFI-style GC-corrected short/long ratios in large genomic bins
cohort    merge the per-sample feature outputs into model-ready matrices and QC plots

Typical run
-----------
  python cfdna_frag.py extract  -b S1.bam -r hg38.fa -o S1.frags.tsv.gz --blacklist hg38-blacklist.v2.bed
  python cfdna_frag.py features -f S1.frags.tsv.gz -r hg38.fa -o features/S1
  python cfdna_frag.py cohort   -i features/* -o cohort/
"""
import argparse
import gzip
import itertools
import json
import os
import re
import sys

import numpy as np
import pandas as pd

AUTOSOMES = r"^(chr)?([1-9]|1[0-9]|2[0-2])$"
AUTOSOMES_XY = r"^(chr)?([1-9]|1[0-9]|2[0-2]|X|Y)$"
MOTIFS = ["".join(p) for p in itertools.product("ACGT", repeat=4)]
COMP = str.maketrans("ACGTN", "TGCAN")

# Size windows (bp). Short/long follow Cristiano et al. 2019 (DELFI).
SHORT = (100, 150)
LONG = (151, 220)
MAX_LEN = 1000


def revcomp(s):
    return s.translate(COMP)[::-1]


def log(msg):
    print(msg, file=sys.stderr, flush=True)


# --------------------------------------------------------------------------- blacklist
def load_blacklist(path):
    """BED -> {chrom: (starts, ends)} of merged, sorted intervals."""
    if not path:
        return {}
    bl = pd.read_csv(path, sep="\t", header=None, usecols=[0, 1, 2],
                     names=["chrom", "start", "end"], comment="#")
    out = {}
    for chrom, d in bl.sort_values(["chrom", "start"]).groupby("chrom"):
        starts, ends = [], []
        for s, e in zip(d.start, d.end):
            if starts and s <= ends[-1]:
                ends[-1] = max(ends[-1], e)
            else:
                starts.append(s)
                ends.append(e)
        out[chrom] = (np.array(starts), np.array(ends))
    return out


def overlaps_blacklist(chrom, start, end, bl):
    """Vectorised: boolean mask of fragments [start, end) overlapping any blacklist interval."""
    if chrom not in bl:
        return np.zeros(len(start), dtype=bool)
    bs, be = bl[chrom]
    # last blacklist interval starting before the fragment end
    idx = np.searchsorted(bs, end, side="left") - 1
    valid = idx >= 0
    hit = np.zeros(len(start), dtype=bool)
    hit[valid] = be[idx[valid]] > start[valid]
    return hit


# --------------------------------------------------------------------------- extract
def extract(args):
    import pysam

    bam = pysam.AlignmentFile(args.bam)
    fasta = pysam.FastaFile(args.reference)
    chrom_re = re.compile(args.chroms)
    bl = load_blacklist(args.blacklist)
    chroms = [c for c in bam.references if chrom_re.match(c)]
    if not chroms:
        sys.exit(f"No BAM contigs match --chroms {args.chroms!r}")

    stats = dict(pairs_seen=0, kept=0, blacklisted=0)
    first = True
    with gzip.open(args.out, "wt") as fh:
        for chrom in chroms:
            seq = fasta.fetch(chrom).upper()
            rows = []
            for r in bam.fetch(chrom):
                if (r.is_unmapped or r.mate_is_unmapped or not r.is_paired or not r.is_proper_pair
                        or r.is_secondary or r.is_supplementary or r.is_duplicate or r.is_qcfail):
                    continue
                tlen = r.template_length
                # one record per fragment: the leftmost read (tie -> read1)
                if tlen <= 0 or (r.reference_start == r.next_reference_start and not r.is_read1):
                    continue
                stats["pairs_seen"] += 1
                mate_mq = r.get_tag("MQ") if r.has_tag("MQ") else r.mapping_quality
                if min(r.mapping_quality, mate_mq) < args.min_mapq:
                    continue
                if not (args.min_len <= tlen <= MAX_LEN):
                    continue
                start = r.reference_start
                end = start + tlen
                if end > len(seq):
                    continue
                frag = seq[start:end]
                gc = (frag.count("G") + frag.count("C")) / max(1, tlen - frag.count("N"))
                rows.append((start, end, tlen, min(r.mapping_quality, mate_mq), round(gc, 4),
                             frag[:4], revcomp(frag[-4:])))
            if not rows:
                continue
            df = pd.DataFrame(rows, columns=["start", "end", "length", "mapq", "gc",
                                             "motif_5p", "motif_3p"])
            hit = overlaps_blacklist(chrom, df.start.values, df.end.values, bl)
            stats["blacklisted"] += int(hit.sum())
            df = df[~hit]
            df.insert(0, "chrom", chrom)
            stats["kept"] += len(df)
            df.to_csv(fh, sep="\t", index=False, header=first)
            first = False
            log(f"{chrom}: {len(df):,} fragments")
    log(json.dumps(stats))


# --------------------------------------------------------------------------- features
def size_features(lengths):
    hist = np.bincount(lengths, minlength=MAX_LEN + 1)[: MAX_LEN + 1]
    n = hist.sum()
    in_rng = lambda lo, hi: hist[lo:hi + 1].sum()  # noqa: E731
    short, long_ = in_rng(*SHORT), in_rng(*LONG)

    # ~10 bp periodicity (nucleosome/DNA-helix footprint) in the sub-nucleosomal range:
    # fraction of detrended spectral power at period 10-11 bp.
    seg = hist[60:151].astype(float) / max(n, 1)
    trend = np.convolve(seg, np.ones(11) / 11, mode="same")
    resid = (seg - trend)[5:-5]
    power = np.abs(np.fft.rfft(resid - resid.mean())) ** 2
    freqs = np.fft.rfftfreq(len(resid))
    band = (freqs >= 1 / 11) & (freqs <= 1 / 9.5)
    period10 = float(power[band].sum() / power[1:].sum()) if power[1:].sum() > 0 else np.nan

    summary = {
        "n_fragments": int(n),
        "median_len": float(np.median(lengths)) if n else np.nan,
        "mode_len": int(np.argmax(hist[50:400]) + 50),
        "frac_lt150": in_rng(0, 149) / n,
        "frac_90_150": in_rng(90, 150) / n,
        "frac_short_100_150": short / n,
        "frac_long_151_220": long_ / n,
        "short_long_ratio": short / long_ if long_ else np.nan,
        "frac_dinucleosome_250_450": in_rng(250, 450) / n,
        "period10_power": period10,
    }
    return hist, summary


def motif_features(frags):
    counts = pd.concat([frags.motif_5p, frags.motif_3p]).value_counts()
    counts = counts.reindex(MOTIFS, fill_value=0)  # drops motifs containing N
    freq = counts / counts.sum()
    p = freq[freq > 0].values
    mds = float(-(p * np.log(p)).sum() / np.log(len(MOTIFS)))  # Jiang et al. 2020
    return freq, mds


def lowess_correct(counts, gc, frac=0.75):
    """DELFI-style GC correction: residual of a lowess fit of counts ~ GC, plus the median."""
    from statsmodels.nonparametric.smoothers_lowess import lowess

    ok = np.isfinite(counts) & np.isfinite(gc)
    out = np.full(len(counts), np.nan)
    if ok.sum() < 10:
        return out
    fit = lowess(counts[ok], gc[ok], frac=frac, return_sorted=False)
    out[ok] = counts[ok] - fit + np.median(counts[ok])
    return out


def bin_features(frags, bin_size, fai, chrom_re):
    sizes = pd.read_csv(fai, sep="\t", header=None, usecols=[0, 1], names=["chrom", "len"])
    sizes = sizes[sizes.chrom.str.match(chrom_re)]
    bins = pd.concat([pd.DataFrame({"chrom": c, "start": np.arange(0, L, bin_size)})
                      for c, L in zip(sizes.chrom, sizes.len)], ignore_index=True)
    bins["end"] = bins.start + bin_size

    f = frags.assign(start=(frags.start // bin_size) * bin_size)
    f["is_short"] = f.length.between(*SHORT)
    f["is_long"] = f.length.between(*LONG)
    agg = f.groupby(["chrom", "start"]).agg(
        n=("length", "size"), short=("is_short", "sum"), long=("is_long", "sum"),
        gc=("gc", "mean")).reset_index()
    b = bins.merge(agg, on=["chrom", "start"], how="left")
    b[["n", "short", "long"]] = b[["n", "short", "long"]].fillna(0)

    # drop bins with poor coverage (gaps, centromeres, heavily blacklisted)
    med = b.n[b.n > 0].median()
    b["pass"] = b.n >= 0.25 * med
    p = b["pass"]
    for col in ("short", "long", "n"):
        b[f"{col}_corr"] = np.nan
        b.loc[p, f"{col}_corr"] = lowess_correct(b.loc[p, col].values.astype(float),
                                                b.loc[p, "gc"].values)
    b["ratio"] = b.short_corr / b.long_corr
    b["ratio_centered"] = b.ratio - b.ratio[p].mean()
    b["coverage_norm"] = b.n_corr / b.n_corr[p].median()
    return b


def features(args):
    os.makedirs(args.out, exist_ok=True)
    name = args.name or os.path.basename(args.out.rstrip("/"))
    frags = pd.read_csv(args.fragments, sep="\t",
                        dtype={"chrom": str, "motif_5p": str, "motif_3p": str})
    frags = frags[frags.length.between(args.min_len, args.max_len)]

    # Sex-chromosome fractions. In a sex-mismatched transplant (male donor -> female recipient)
    # the chrY fraction relative to male controls approximates the donor-derived fraction.
    chrom = frags.chrom.str.replace("^chr", "", regex=True)
    sex = {"frac_chrX": float((chrom == "X").mean()), "frac_chrY": float((chrom == "Y").mean())}
    frags = frags[frags.chrom.str.match(args.chroms)]
    log(f"{name}: {len(frags):,} fragments on autosomes")

    hist, summary = size_features(frags.length.values)
    motif_freq, mds = motif_features(frags)
    summary["motif_diversity_score"] = mds
    summary.update(sex)
    for m in ("CCCA", "CCAG", "CCTG", "TGGG"):  # motifs often reported shifted in cancer
        summary[f"motif_{m}"] = float(motif_freq[m])
    bins = bin_features(frags, args.bin_size, args.reference + ".fai", args.chroms)
    summary["n_bins_pass"] = int(bins["pass"].sum())
    summary["sample"] = name

    pd.DataFrame({"length": np.arange(len(hist)), "count": hist}).to_csv(
        f"{args.out}/size_hist.tsv", sep="\t", index=False)
    motif_freq.rename("freq").rename_axis("motif").to_csv(f"{args.out}/end_motifs.tsv", sep="\t")
    bins.to_csv(f"{args.out}/bins.tsv", sep="\t", index=False, float_format="%.5g")
    with open(f"{args.out}/summary.json", "w") as fh:
        json.dump(summary, fh, indent=2, default=float)
    log(json.dumps(summary, default=float))


# --------------------------------------------------------------------------- cohort
def cohort(args):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    os.makedirs(args.out, exist_ok=True)
    summaries, ratios, covs, motifs, hists = [], {}, {}, {}, {}
    for d in args.inputs:
        with open(f"{d}/summary.json") as fh:
            s = json.load(fh)
        name = s["sample"]
        summaries.append(s)
        b = pd.read_csv(f"{d}/bins.tsv", sep="\t", dtype={"chrom": str})
        key = b.chrom + ":" + b.start.astype(str)
        ratios[name] = pd.Series(b.ratio_centered.values, index=key)
        covs[name] = pd.Series(b.coverage_norm.values, index=key)
        motifs[name] = pd.read_csv(f"{d}/end_motifs.tsv", sep="\t", index_col=0)["freq"]
        h = pd.read_csv(f"{d}/size_hist.tsv", sep="\t")["count"]
        hists[name] = h / h.sum()

    summary = pd.DataFrame(summaries).set_index("sample")
    groups = pd.Series("all", index=summary.index, name="group")
    if args.samples:
        sheet = pd.read_csv(args.samples, sep="\t", dtype=str).set_index("sample")
        missing = set(summary.index) - set(sheet.index)
        if missing:
            sys.exit(f"Samples missing from {args.samples}: {sorted(missing)}")
        groups = sheet.loc[summary.index, "group"]
    order = groups.sort_values(kind="stable").index
    summary = summary.loc[order]
    summary.insert(0, "group", groups.loc[order])
    # keep bins that pass in every sample so the matrix is complete for modelling
    ratio_m = pd.DataFrame(ratios).dropna().T.loc[order]
    cov_m = pd.DataFrame(covs).loc[ratio_m.columns].T.loc[order]
    motif_m = pd.DataFrame(motifs).T.loc[order]

    summary.to_csv(f"{args.out}/size_motif_summary.tsv", sep="\t")
    ratio_m.to_csv(f"{args.out}/delfi_short_long_ratio.tsv", sep="\t", float_format="%.5g")
    cov_m.to_csv(f"{args.out}/coverage_norm.tsv", sep="\t", float_format="%.5g")
    motif_m.to_csv(f"{args.out}/end_motif_freq.tsv", sep="\t", float_format="%.6g")
    num = summary.drop(columns="group").select_dtypes("number")
    num.groupby(summary.group).agg(["mean", "min", "max"]).T.to_csv(
        f"{args.out}/group_summary.tsv", sep="\t", float_format="%.5g")
    log(f"{len(summary)} samples, {ratio_m.shape[1]} bins shared")

    # -- plots
    palette = ["#1b6ca8", "#d1495b", "#2e933c", "#edae49", "#7b4b94", "#5c5c5c"]
    gcol = {g: palette[i % len(palette)] for i, g in enumerate(groups.loc[order].unique())}
    styles = ["-", "--", ":", "-."]

    fig, axes = plt.subplots(1, 2, figsize=(12, 4))
    rep = groups.loc[order].groupby(groups.loc[order]).cumcount()  # replicate index within group
    for name in order:
        g = groups[name]
        ls = styles[rep[name] % len(styles)]
        h = hists[name]
        for ax in axes:
            ax.plot(h.index, h.values, lw=1.1, ls=ls, color=gcol[g], label=f"{name} ({g})")
    axes[0].set(xlim=(50, 450), xlabel="Fragment length (bp)", ylabel="Proportion",
                title="Fragment size distribution")
    axes[1].set(xlim=(50, 600), yscale="log", xlabel="Fragment length (bp)",
                title="Log scale (long / multi-nucleosomal fragments)")
    for ax in axes:
        ax.axvline(167, color="grey", ls=":", lw=0.8)
    if len(hists) <= 12:
        axes[0].legend(fontsize=7, frameon=False)
    fig.tight_layout()
    fig.savefig(f"{args.out}/size_distribution.png", dpi=150)

    # PCA of the genome-wide ratio profile and of end-motif frequencies
    fig, axes = plt.subplots(1, 2, figsize=(10, 4.2))
    for ax, (title, m) in zip(axes, [("Short/long ratio profile", ratio_m),
                                     ("End-motif frequencies", motif_m)]):
        x = m.values - m.values.mean(axis=0)
        x = x / np.where(x.std(axis=0) > 0, x.std(axis=0), 1)
        if len(m) < 3:
            ax.set_axis_off()
            continue
        u, sv, _ = np.linalg.svd(x, full_matrices=False)
        pcs, var = u * sv, sv ** 2 / (sv ** 2).sum()
        for i, name in enumerate(m.index):
            ax.scatter(pcs[i, 0], pcs[i, 1], color=gcol[groups[name]], s=40)
            ax.annotate(name, (pcs[i, 0], pcs[i, 1]), fontsize=7, xytext=(3, 3),
                        textcoords="offset points")
        ax.set(title=title, xlabel=f"PC1 ({var[0]:.0%})", ylabel=f"PC2 ({var[1]:.0%})")
    handles = [plt.Line2D([], [], marker="o", ls="", color=c, label=g) for g, c in gcol.items()]
    fig.legend(handles=handles, loc="lower center", ncol=len(gcol), frameon=False, fontsize=8)
    fig.tight_layout(rect=(0, 0.06, 1, 1))
    fig.savefig(f"{args.out}/pca.png", dpi=150)

    chroms = [c.split(":")[0] for c in ratio_m.columns]
    fig, ax = plt.subplots(figsize=(12, max(2, 0.3 * len(ratio_m) + 1)))
    lim = np.nanpercentile(np.abs(ratio_m.values), 99) if ratio_m.size else 1
    ax.imshow(ratio_m.values, aspect="auto", cmap="RdBu_r", vmin=-lim, vmax=lim,
              interpolation="nearest")
    ax.set_yticks(range(len(ratio_m)), ratio_m.index, fontsize=7)
    breaks = [i for i in range(1, len(chroms)) if chroms[i] != chroms[i - 1]]
    for x in breaks:
        ax.axvline(x - 0.5, color="k", lw=0.4)
    ticks = [0] + breaks
    ax.set_xticks([(a + b) / 2 for a, b in zip(ticks, ticks[1:] + [len(chroms)])],
                  [chroms[t].replace("chr", "") for t in ticks], fontsize=6)
    ax.set_title("Centered short/long fragment ratio (GC-corrected)")
    fig.tight_layout()
    fig.savefig(f"{args.out}/delfi_ratio_heatmap.png", dpi=150)


# --------------------------------------------------------------------------- cli
def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawTextHelpFormatter)
    sub = p.add_subparsers(dest="cmd", required=True)

    e = sub.add_parser("extract", help="BAM -> fragment table")
    e.add_argument("-b", "--bam", required=True)
    e.add_argument("-r", "--reference", required=True, help="FASTA (faidx-indexed), same build as BAM")
    e.add_argument("-o", "--out", required=True, help="output .tsv.gz")
    e.add_argument("--blacklist", help="BED of regions to exclude (e.g. ENCODE hg38 blacklist v2)")
    e.add_argument("--min-mapq", type=int, default=30)
    e.add_argument("--min-len", type=int, default=20)
    e.add_argument("--chroms", default=AUTOSOMES_XY,
                   help="regex of contigs to keep (default autosomes + X + Y)")
    e.set_defaults(func=extract)

    f = sub.add_parser("features", help="fragment table -> per-sample features")
    f.add_argument("-f", "--fragments", required=True)
    f.add_argument("-r", "--reference", required=True, help="FASTA; its .fai defines the bins")
    f.add_argument("-o", "--out", required=True, help="output directory")
    f.add_argument("-n", "--name", help="sample name (default: output dir name)")
    f.add_argument("--bin-size", type=int, default=5_000_000)
    f.add_argument("--min-len", type=int, default=50)
    f.add_argument("--max-len", type=int, default=600)
    f.add_argument("--chroms", default=AUTOSOMES,
                   help="regex of contigs used for size/motif/bin features (default autosomes)")
    f.set_defaults(func=features)

    c = sub.add_parser("cohort", help="merge per-sample features")
    c.add_argument("-i", "--inputs", nargs="+", required=True, help="feature directories")
    c.add_argument("-o", "--out", required=True)
    c.add_argument("-s", "--samples", help="TSV with columns 'sample' and 'group' (colours plots)")
    c.set_defaults(func=cohort)

    args = p.parse_args(argv)
    args.func(args)


if __name__ == "__main__":
    main()
