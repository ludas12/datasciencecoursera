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
def ladder_amplitude(hist, lo=80, hi=150):
    """Relative amplitude of the ~10 bp ladder in [lo, hi] bp.

    The histogram is divided by its 11 bp moving average, and a sinusoid with period 10-11 bp
    is fitted by least squares; the best amplitude is returned (0.10 = peaks ~10% above trend).
    Fitting one frequency averages out counting noise, unlike a peak-to-trough measure."""
    seg = hist[lo - 5:hi + 6].astype(float)
    if seg.sum() < 1000:
        return np.nan
    trend = np.convolve(seg, np.ones(11) / 11, mode="same")[5:-5]
    rel = seg[5:-5] / np.where(trend > 0, trend, np.nan) - 1
    x = np.arange(lo, hi + 1)
    ok = np.isfinite(rel)
    best = 0.0
    for period in np.arange(10.0, 11.01, 0.05):
        w = 2 * np.pi * x[ok] / period
        design = np.column_stack([np.sin(w), np.cos(w), np.ones(ok.sum())])
        coef, *_ = np.linalg.lstsq(design, rel[ok], rcond=None)
        best = max(best, float(np.hypot(coef[0], coef[1])))
    return best


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
        # BAL-oriented shape metrics
        # mono-nucleosome peak: mean count at 163-171 bp over the shoulders either side
        "peak167_prominence": float(hist[163:172].mean() / np.concatenate(
            [hist[140:150], hist[185:195]]).mean()) if n else np.nan,
        # 10 bp ladder strength in the sub-nucleosomal range (nuclease trimming of nucleosomes)
        "ladder_amp_80_150": ladder_amplitude(hist),
        # long, structureless tail: randomly cut DNA (no di/tri-nucleosome peaks)
        "frac_400_600": in_rng(400, 600) / n,
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


# --------------------------------------------------------------------------- tss
GENESETS = os.path.join(os.path.dirname(os.path.abspath(__file__)), "genesets.tsv")
TSS_FLANK = 3000          # profile spans TSS +/- this many bp
TSS_NORM = (2000, 3000)   # |distance| range used as the local baseline
TSS_CENTRAL = (-150, 150) # nucleosome-depleted region; mean coverage here / baseline
# finer windows: the nucleosome-depleted region just upstream of the TSS, and the strongly
# positioned +1 nucleosome downstream. Core-trimmed fragments (e.g. BAL) mark the +1 nucleosome
# sharply, which can cancel the NDR dip inside the +/-150 bp central window.
TSS_NDR = (-200, -20)
TSS_PLUS1 = (60, 200)


def load_tss(bed, genes):
    """TSS BED -> DataFrame(chrom, pos, strand, gene) for genes in `genes`.

    The gene-name column is detected as the column with the most matches to `genes`;
    the strand column as one holding only +/-/. values. A 1-bp interval is taken as the TSS;
    otherwise start (+ strand) or end-1 (- strand)."""
    df = pd.read_csv(bed, sep="\t", header=None, comment="#", dtype=str)
    df = df[~df[0].str.startswith(("track", "browser"))]
    gene_set = set(genes)
    hits = {c: df[c].isin(gene_set).sum() for c in df.columns[3:]}
    if not hits or max(hits.values()) == 0:
        sys.exit(f"No gene-set symbols found in {bed}; first line: {df.iloc[0].tolist()}")
    name_col = max(hits, key=hits.get)
    strand_col = next((c for c in df.columns[3:]
                       if df[c].isin(["+", "-", "."]).all()), None)
    out = pd.DataFrame({"chrom": df[0], "start": df[1].astype(int), "end": df[2].astype(int),
                        "gene": df[name_col],
                        "strand": df[strand_col] if strand_col is not None else "+"})
    out = out[out.gene.isin(gene_set)]
    one_bp = (out.end - out.start) == 1
    out["pos"] = np.where(one_bp | (out.strand != "-"), out.start, out.end - 1)
    return collapse_tss(out[["chrom", "pos", "strand", "gene"]])


TSS_MERGE = 500  # TSSs of the same gene within this distance are one promoter


def collapse_tss(sites, merge=TSS_MERGE):
    """One TSS per promoter: start sites of the same gene (and strand) within `merge` bp are
    collapsed to the position used by most transcripts (ties: nearest the cluster median).

    Without this, a gene with many near-identical transcript starts gets many overlapping
    windows, so a single high-coverage locus is counted many times at almost the same offset."""
    counts = sites.groupby(["chrom", "strand", "gene", "pos"]).size().rename("n").reset_index()
    keep = []
    for _, g in counts.groupby(["chrom", "strand", "gene"], sort=False):
        g = g.sort_values("pos")
        pos, n = g.pos.values, g.n.values
        start = 0
        for i in range(1, len(pos) + 1):
            if i == len(pos) or pos[i] - pos[start] > merge:
                cp, cn = pos[start:i], n[start:i]
                best = np.flatnonzero(cn == cn.max())
                med = np.median(cp)
                keep.append(g.iloc[start + best[np.argmin(np.abs(cp[best] - med))]])
                start = i
    return pd.DataFrame(keep)[["chrom", "pos", "strand", "gene"]].reset_index(drop=True)


TSS_OUTLIER_HIGH = 5.0    # drop TSS windows with mean coverage > this x the median window
TSS_OUTLIER_LOW = 0.1     # ... or < this x the median (unmappable / deleted)
TSS_CAP_QUANTILE = 0.999  # per-position coverage of one TSS is capped at this pooled quantile
GC_BINS = np.linspace(0, 1, 51)


def tss_coverage(frags, tss, flank=TSS_FLANK, weights=None):
    """Strand-oriented (optionally weighted) fragment coverage around each TSS
    -> array (n_tss, 2*flank+1)."""
    width = 2 * flank + 1
    out = np.zeros((len(tss), width))
    max_len = int(frags.length.max()) if len(frags) else 0
    w_all = np.ones(len(frags)) if weights is None else np.asarray(weights, float)
    chrom_vals = frags.chrom.values
    for chrom, t in tss.groupby("chrom"):
        sel = chrom_vals == chrom
        if not sel.any():
            continue
        st0, en0, w0_ = frags.start.values[sel], frags.end.values[sel], w_all[sel]
        order = np.argsort(st0, kind="stable")
        st, en, wt = st0[order], en0[order], w0_[order]
        for i, (pos, strand) in zip(t.index, zip(t.pos.values, t.strand.values)):
            w0, w1 = pos - flank, pos + flank + 1
            lo = np.searchsorted(st, w0 - max_len, side="left")
            hi = np.searchsorted(st, w1, side="left")
            s_, e_, ww = st[lo:hi], en[lo:hi], wt[lo:hi]
            keep = e_ > w0
            s_, e_, ww = np.clip(s_[keep], w0, w1) - w0, np.clip(e_[keep], w0, w1) - w0, ww[keep]
            diff = np.zeros(width + 1)
            np.add.at(diff, s_, ww)
            np.add.at(diff, e_, -ww)
            cov = np.cumsum(diff)[:width]
            out[tss.index.get_loc(i)] = cov[::-1] if strand == "-" else cov
    return out


def gc_weights(frags, reference, n_sample=500_000, seed=0):
    """Per-fragment GC-bias weights: expected / observed density of fragment GC.

    Expected = GC of random genomic fragments (same chromosomes and length distribution as the
    sample, N-free). Weights are clipped to [0.25, 4], sparse GC bins take the nearest
    well-sampled bin's weight, and weights are scaled
    to mean 1 over the sample's fragments. Returns (weights, table for QC)."""
    import pysam

    rng = np.random.default_rng(seed)
    fa = pysam.FastaFile(reference)
    chroms = [c for c in frags.chrom.unique() if c in fa.references]
    lens = np.array([fa.get_reference_length(c) for c in chroms], float)
    pick = rng.choice(len(chroms), n_sample, p=lens / lens.sum())
    flen = rng.choice(frags.length.values, n_sample)
    exp_gc = []
    for ci, c in enumerate(chroms):
        m = pick == ci
        L = int(lens[ci])
        starts = rng.integers(0, L - 1000, m.sum())
        for st, ln in zip(starts, flen[m]):
            seq = fa.fetch(c, int(st), int(st + ln)).upper()
            if "N" in seq or not seq:
                continue
            exp_gc.append((seq.count("G") + seq.count("C")) / len(seq))
    exp_n, _ = np.histogram(exp_gc, GC_BINS)
    obs_n, _ = np.histogram(frags.gc.values, GC_BINS)
    exp_p, obs_p = exp_n / max(exp_n.sum(), 1), obs_n / max(obs_n.sum(), 1)
    ok = (exp_n >= 50) & (obs_n >= 50)
    centers = (GC_BINS[:-1] + GC_BINS[1:]) / 2
    with np.errstate(divide="ignore", invalid="ignore"):
        ratio = exp_p / obs_p
    # sparse bins take the weight of the nearest well-sampled bin (edges held constant)
    w = np.interp(centers, centers[ok], ratio[ok]) if ok.sum() >= 2 else np.ones(len(centers))
    w = np.clip(w, 0.25, 4.0)
    b = np.clip(np.digitize(frags.gc.values, GC_BINS) - 1, 0, len(w) - 1)
    weights = w[b]
    weights /= weights.mean()
    table = pd.DataFrame({"gc_lo": GC_BINS[:-1], "gc_hi": GC_BINS[1:], "observed": obs_n,
                          "expected": exp_n, "weight": w / np.mean(w[b])})
    return weights, table


def tss(args):
    os.makedirs(args.out, exist_ok=True)
    sets = pd.read_csv(args.genesets, sep="\t", comment="#")
    sites = load_tss(args.tss_bed, sets.gene).reset_index(drop=True)
    frags = pd.read_csv(args.fragments, sep="\t", usecols=["chrom", "start", "end", "length", "gc"],
                        dtype={"chrom": str})
    frags = frags[frags.length.between(args.min_len, args.max_len)].reset_index(drop=True)
    sites = sites[sites.chrom.isin(set(frags.chrom))].reset_index(drop=True)

    weights = None
    if args.reference:
        weights, gct = gc_weights(frags, args.reference)
        gct.to_csv(f"{args.out}/tss_gc_weights.tsv", sep="\t", index=False, float_format="%.4g")
        log(f"GC weights: range {weights.min():.2f}-{weights.max():.2f}")
    cov = tss_coverage(frags, sites, weights=weights)

    # outlier windows: abnormally high (repeats, CNVs, artefacts) or near-empty (unmappable)
    wmean = cov.mean(axis=1)
    med = np.median(wmean[wmean > 0]) if (wmean > 0).any() else 0
    good = (wmean <= TSS_OUTLIER_HIGH * med) & (wmean >= TSS_OUTLIER_LOW * med)
    # single-position spikes inside otherwise normal windows
    cap = max(3.0, float(np.quantile(cov[good], TSS_CAP_QUANTILE))) if good.any() else np.inf
    cov = np.minimum(cov, cap)

    x = np.arange(-TSS_FLANK, TSS_FLANK + 1)
    base = (np.abs(x) >= TSS_NORM[0]) & (np.abs(x) <= TSS_NORM[1])
    central = (x >= TSS_CENTRAL[0]) & (x <= TSS_CENTRAL[1])
    ndr = (x >= TSS_NDR[0]) & (x <= TSS_NDR[1])
    plus1 = (x >= TSS_PLUS1[0]) & (x <= TSS_PLUS1[1])

    profiles, rows = {}, []
    for name, genes in sets.groupby("set").gene:
        in_set = sites.gene.isin(set(genes)).values
        idx = np.flatnonzero(in_set & good)
        if len(idx) == 0:
            log(f"{name}: no usable TSSs")
            continue
        agg = cov[idx].sum(axis=0)
        baseline = agg[base].mean()
        prof = agg / baseline if baseline > 0 else np.full_like(agg, np.nan)
        profiles[name] = prof
        rows.append({"set": name, "n_genes_found": sites.loc[in_set, "gene"].nunique(),
                     "n_genes_in_set": genes.nunique(), "n_tss": len(idx),
                     "n_tss_excluded": int((in_set & ~good).sum()),
                     "baseline_coverage": baseline / len(idx),
                     "central_coverage": float(prof[central].mean()),
                     "ndr_coverage": float(prof[ndr].mean()),
                     "plus1_coverage": float(prof[plus1].mean())})
    summary = pd.DataFrame(rows)
    pd.DataFrame(profiles, index=pd.Index(x, name="position")).to_csv(
        f"{args.out}/tss_profiles.tsv", sep="\t", float_format="%.5g")
    summary.to_csv(f"{args.out}/tss_summary.tsv", sep="\t", index=False, float_format="%.5g")
    log(f"coverage cap per TSS position: {cap:.2f}; GC correction: {'on' if weights is not None else 'off'}")
    log(summary.to_string(index=False))


# --------------------------------------------------------------------------- cohort
def cohort(args):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    os.makedirs(args.out, exist_ok=True)
    summaries, ratios, covs, motifs, hists = [], {}, {}, {}, {}
    tss_sum, tss_prof = {}, {}
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
        if os.path.exists(f"{d}/tss_summary.tsv"):
            tss_sum[name] = pd.read_csv(f"{d}/tss_summary.tsv", sep="\t").set_index("set")
            tss_prof[name] = pd.read_csv(f"{d}/tss_profiles.tsv", sep="\t", index_col=0)

    summary = pd.DataFrame(summaries).set_index("sample")
    groups = pd.Series("all", index=summary.index, name="group")
    patients = pd.Series(summary.index, index=summary.index, name="patient")
    if args.samples:
        sheet = pd.read_csv(args.samples, sep="\t", dtype=str).set_index("sample")
        missing = set(summary.index) - set(sheet.index)
        if missing:
            sys.exit(f"Samples missing from {args.samples}: {sorted(missing)}")
        groups = sheet.loc[summary.index, "group"]
        if "patient" in sheet.columns:
            patients = sheet.loc[summary.index, "patient"].fillna(pd.Series(summary.index, index=summary.index))
    order = groups.sort_values(kind="stable").index
    summary = summary.loc[order]
    summary.insert(0, "group", groups.loc[order])
    if args.samples:
        extra = [c for c in sheet.columns if c not in ("group", "r1", "r2")]
        for i, c in enumerate(extra, start=1):
            summary.insert(i, c, sheet.loc[order, c])
    # keep bins that pass in every sample so the matrix is complete for modelling
    ratio_m = pd.DataFrame(ratios).dropna().T.loc[order]
    cov_m = pd.DataFrame(covs).loc[ratio_m.columns].T.loc[order]
    motif_m = pd.DataFrame(motifs).T.loc[order]

    summary.to_csv(f"{args.out}/size_motif_summary.tsv", sep="\t")
    ratio_m.to_csv(f"{args.out}/delfi_short_long_ratio.tsv", sep="\t", float_format="%.5g")
    cov_m.to_csv(f"{args.out}/coverage_norm.tsv", sep="\t", float_format="%.5g")
    motif_m.to_csv(f"{args.out}/end_motif_freq.tsv", sep="\t", float_format="%.6g")
    num = summary.select_dtypes("number")
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

    # TSS coverage profiles (if `tss` was run)
    if tss_sum:
        central = pd.DataFrame({n: t["central_coverage"] for n, t in tss_sum.items()}).T.loc[
            [n for n in order if n in tss_sum]]
        central.insert(0, "group", groups.loc[central.index])
        central.to_csv(f"{args.out}/tss_central_coverage.tsv", sep="\t", float_format="%.4g")
        # primary TSS measure = the upstream nucleosome-depleted window when available
        primary, primary_label = central, "Central TSS dip"
        for col in ("ndr_coverage", "plus1_coverage"):
            if all(col in t.columns for t in tss_sum.values()):
                m = pd.DataFrame({n: t[col] for n, t in tss_sum.items()}).T.loc[central.index]
                m.insert(0, "group", groups.loc[m.index])
                m.to_csv(f"{args.out}/tss_{col}.tsv", sep="\t", float_format="%.4g")
                if col == "ndr_coverage":
                    primary = m
                    primary_label = "Upstream TSS dip"
        sets = list(next(iter(tss_prof.values())).columns)
        ncol = min(4, len(sets))
        nrow = int(np.ceil(len(sets) / ncol))
        fig, axes = plt.subplots(nrow, ncol, figsize=(3.6 * ncol, 3 * nrow), squeeze=False,
                                 sharex=True)
        for ax, st in zip(axes.flat, sets):
            for name in central.index:
                prof = tss_prof[name][st].rolling(101, center=True, min_periods=1).mean()
                ax.plot(prof.index, prof.values, lw=1, color=gcol[groups[name]],
                        ls=styles[rep[name] % len(styles)], label=name)
            n_tss = int(tss_sum[central.index[0]].loc[st, "n_tss"])
            ax.set_title(f"{st} ({n_tss} TSSs)", fontsize=9)
            ax.axvspan(*TSS_NDR, color="grey", alpha=0.15, lw=0)   # primary (upstream) window
            ax.axvline(0, color="grey", lw=0.5, ls=":")
            ax.axhline(1, color="grey", lw=0.5)
        for ax in axes.flat[len(sets):]:
            ax.set_axis_off()
        for ax in axes[-1]:
            ax.set_xlabel("Distance to TSS (bp)")
        for ax in axes[:, 0]:
            ax.set_ylabel("Relative coverage")
        axes.flat[0].legend(fontsize=6, frameon=False)
        fig.suptitle(f"Shaded: upstream window used for the TSS metric ({TSS_NDR[0]} to {TSS_NDR[1]} bp)",
                     fontsize=8)
        fig.tight_layout()
        fig.savefig(f"{args.out}/tss_profiles.png", dpi=150)

    # Key metrics per sample; samples from the same patient joined by a line
    metrics = {"median_len": "Median length (bp)", "frac_lt150": "Fraction < 150 bp",
               "short_long_ratio": "Short/long ratio",
               "frac_dinucleosome_250_450": "Di-nucleosomal fraction",
               "motif_CCCA": "CCCA end-motif freq.", "motif_diversity_score": "Motif diversity",
               "peak167_prominence": "167 bp peak prominence",
               "ladder_amp_80_150": "10 bp ladder amplitude (80-150 bp)",
               "frac_400_600": "Fraction 400-600 bp (random cuts)"}
    panels = [(summary[k], v) for k, v in metrics.items() if k in summary]
    if tss_sum:
        for st in ("housekeeping", "inactive_control", "lung_epithelium", "neutrophil"):
            if st in primary.columns:
                panels.append((primary[st], f"{primary_label}: {st}"))
    glist = list(gcol)
    ncol = 4
    nrow = int(np.ceil(len(panels) / ncol))
    fig, axes = plt.subplots(nrow, ncol, figsize=(3.4 * ncol, 3 * nrow), squeeze=False)
    for ax, (vals, title) in zip(axes.flat, panels):
        xpos = {n: glist.index(groups[n]) for n in vals.index}
        for pt, members in vals.groupby(patients.loc[vals.index]).groups.items():
            if len(members) > 1:
                ax.plot([xpos[m] for m in members], vals.loc[members], color="grey", lw=0.8, zorder=1)
        for n, v in vals.items():
            ax.scatter(xpos[n], v, color=gcol[groups[n]], s=36, zorder=2)
        ax.set_xticks(range(len(glist)), glist, fontsize=7, rotation=20)
        ax.set_title(title, fontsize=9)
        ax.set_xlim(-0.5, len(glist) - 0.5)
        ax.ticklabel_format(axis="y", useOffset=False)
    for ax in axes.flat[len(panels):]:
        ax.set_axis_off()
    fig.suptitle("Per-sample values (lines join samples from the same patient)", fontsize=9)
    fig.tight_layout()
    fig.savefig(f"{args.out}/key_metrics.png", dpi=150)


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

    t = sub.add_parser("tss", help="fragment table -> aggregate TSS coverage per gene set")
    t.add_argument("-f", "--fragments", required=True)
    t.add_argument("-t", "--tss-bed", required=True,
                   help="BED of TSSs with gene symbols (e.g. hg38_tss.bed / RefSeq TSS BED)")
    t.add_argument("-o", "--out", required=True, help="sample feature directory")
    t.add_argument("-g", "--genesets", default=GENESETS, help="TSV: set<TAB>gene")
    t.add_argument("-r", "--reference", help="FASTA (faidx-indexed); enables GC-bias correction")
    t.add_argument("--min-len", type=int, default=100)
    t.add_argument("--max-len", type=int, default=220)
    t.set_defaults(func=tss)

    c = sub.add_parser("cohort", help="merge per-sample features")
    c.add_argument("-i", "--inputs", nargs="+", required=True, help="feature directories")
    c.add_argument("-o", "--out", required=True)
    c.add_argument("-s", "--samples", help="TSV with columns 'sample' and 'group' (colours plots)")
    c.set_defaults(func=cohort)

    args = p.parse_args(argv)
    args.func(args)


if __name__ == "__main__":
    main()
