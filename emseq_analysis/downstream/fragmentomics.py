#!/usr/bin/env python3
"""cfDNA fragmentomics from paired-end EM-seq alignments.

extract   one library -> OUT/<label>.frag.npz
          fragments = primary, non-duplicate pairs, both mates MAPQ >= 30, same chromosome,
          FR orientation, <= 1,000 bp (the aligner's proper-pair flag is not required).
          True fragment ends are restored when reads were clipped before alignment
          (--clip R1 R2) and from soft clips (own CIGAR and the mate's MC tag).
          Collects: size histogram; short (100-150 bp) / long (151-220 bp) counts per
          5 Mb bin with bin GC; fragments (100-220 bp) overlapping +/-3 kb of the TSSs of
          the gene sets in --genesets, with their GC; observed and expected fragment-GC
          histograms; 4-mer end motifs read from the reference (EM-seq conversion does not
          affect them).
          TSS method as in the cfdna-fragmentomics pipeline (cfdna_frag.py tss): one TSS per
          promoter (same-gene starts within 500 bp merged), GC-bias weights per fragment
          (expected/observed fragment-GC density), outlier windows dropped (>5x or <0.1x the
          median window), per-position cap at the pooled 99.9th percentile, coverage
          normalised to the 2-3 kb flanks; NDR (-200..-20 bp) reported relative to each
          library's inactive_control set.
summarize all libraries -> tables and plots in OUT/
"""
import argparse
import bisect
import gzip
import json
import subprocess
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np

MAX_LEN = 1000
BIN = 5_000_000
SHORT, LONG = (100, 150), (151, 220)
TSS_FLANK = 3000          # profile spans TSS +/- this many bp
TSS_NORM = (2000, 3000)   # |distance| range used as the local baseline
TSS_CENTRAL = (-150, 150)
TSS_NDR = (-200, -20)     # nucleosome-depleted region (primary measure)
TSS_PLUS1 = (60, 200)     # +1 nucleosome
TSS_MERGE = 500           # TSSs of the same gene within this distance are one promoter
TSS_OUTLIER_HIGH = 5.0    # drop TSS windows with mean coverage > this x the median window
TSS_OUTLIER_LOW = 0.1     # ... or < this x the median (unmappable / deleted)
TSS_CAP_QUANTILE = 0.999  # per-position coverage is capped at this pooled quantile
TSS_LEN = (100, 220)      # fragment lengths used for TSS coverage
GC_BINS = np.linspace(0, 1, 51)
N_GC_SAMPLE = 500_000     # random genomic fragments for the expected GC distribution
COMP = str.maketrans("ACGTN", "TGCAN")
BASES = "ACGT"
KMERS = [a + b + c + d for a in BASES for b in BASES for c in BASES for d in BASES]
KIDX = {k: i for i, k in enumerate(KMERS)}
DEFAULT_GENESETS = Path(__file__).resolve().parent / "genesets.tsv"


def read_genesets(path):
    """set <tab> gene file (comments allowed) -> {set: set(genes)}."""
    sets = defaultdict(set)
    for line in open(path):
        if line.startswith("#") or not line.strip():
            continue
        f = line.rstrip("\n").split("\t")
        if len(f) >= 2 and f[0] != "set":
            sets[f[0]].add(f[1])
    return dict(sets)


def load_fai(fai):
    out = {}
    for line in open(fai):
        f = line.split("\t")
        out[f[0]] = (int(f[1]), int(f[2]), int(f[3]), int(f[4]))
    return out


def chrom_seq(fasta, fai, chrom):
    length, offset, lb, lbytes = fai[chrom]
    nlines = length // lb + 1
    with open(fasta, "rb") as fh:
        fh.seek(offset)
        raw = fh.read(nlines * lbytes)
    return raw.replace(b"\n", b"").replace(b"\r", b"")[:length].decode().upper()


def load_tss(refseq, genes):
    """RefSeq table (UCSC ncbiRefSeqCurated) -> one TSS per promoter for `genes`.

    Start sites of the same gene and strand within TSS_MERGE bp are collapsed to the position
    used by most transcripts (ties: nearest the cluster median). Returns a list of
    (chrom, pos, strand, gene) sorted by chrom and pos."""
    counts = defaultdict(lambda: defaultdict(int))
    if refseq and Path(refseq).exists():
        with gzip.open(refseq, "rt") as fh:
            for line in fh:
                f = line.rstrip("\n").split("\t")
                if len(f) > 12 and f[12] in genes:
                    pos = int(f[4]) if f[3] == "+" else int(f[5]) - 1
                    counts[(f[2], f[3], f[12])][pos] += 1
    out = []
    for (chrom, strand, gene), pc in counts.items():
        pos = sorted(pc)
        start = 0
        for i in range(1, len(pos) + 1):
            if i == len(pos) or pos[i] - pos[start] > TSS_MERGE:
                cl = pos[start:i]
                top = max(pc[p] for p in cl)
                best = [p for p in cl if pc[p] == top]
                med = float(np.median(cl))
                out.append((chrom, min(best, key=lambda p: abs(p - med)), strand, gene))
                start = i
    return sorted(out)


def frag_gc(sq):
    acgt = len(sq) - sq.count("N")
    return (sq.count("G") + sq.count("C")) / acgt if acgt else -1.0


def softclip(cigar, end):
    """Length of the soft clip at the start ('L') or end ('R') of a CIGAR string."""
    if cigar in ("*", ""):
        return 0
    num = ""
    ops = []
    for ch in cigar:
        if ch.isdigit():
            num += ch
        else:
            ops.append((int(num), ch))
            num = ""
    if not ops:
        return 0
    n, op = ops[0] if end == "L" else ops[-1]
    return n if op == "S" else 0


def extract(a):
    fai = load_fai(a.fai)
    gsets = read_genesets(a.genesets)
    all_genes = set().union(*gsets.values()) if gsets else set()
    tss_all = load_tss(a.refseq, all_genes)
    tss = defaultdict(list)                     # chrom -> [(pos, strand, global index)]
    for k, (c, pos, strand, _) in enumerate(tss_all):
        tss[c].append((pos, strand, k))
    rec_tss, rec_lo, rec_hi, rec_gc = [], [], [], []
    obs_gc = np.zeros(len(GC_BINS) - 1, dtype=np.int64)
    sizes = np.zeros(MAX_LEN + 1, dtype=np.int64)
    autos = [c for c in fai if c.replace("chr", "").isdigit() and 1 <= int(c.replace("chr", "")) <= 22]
    nbins = {c: fai[c][0] // BIN + 1 for c in autos}
    bins_short = {c: np.zeros(nbins[c], dtype=np.int64) for c in autos}
    bins_long = {c: np.zeros(nbins[c], dtype=np.int64) for c in autos}
    bins_gc, bins_n = {}, {}

    def bin_gc(chrom, sq):
        gc, nn = np.zeros(nbins[chrom]), np.zeros(nbins[chrom])
        for b in range(nbins[chrom]):
            w = sq[b * BIN:(b + 1) * BIN]
            acgt = len(w) - w.count("N")
            gc[b] = (w.count("G") + w.count("C")) / acgt if acgt else np.nan
            nn[b] = w.count("N") / max(len(w), 1)
        bins_gc[chrom], bins_n[chrom] = gc, nn
    motifs = np.zeros(len(KMERS), dtype=np.int64)
    clip = {64: a.clip[0], 128: a.clip[1]}            # read1 / read2 flag bits

    # Not "-f 0x2": bwa-meth sets the proper-pair flag from each library's own insert-size
    # distribution, which can drop di-nucleosomal fragments in sharply peaked plasma
    # libraries. Instead: both mates mapped to the same chromosome in FR orientation,
    # fragment <= MAX_LEN, MAPQ >= mapq for both mates (MQ tag from samtools fixmate).
    cmd = ["samtools", "view", "-f", "0x1", "-F", "0xF0C", "-q", str(a.mapq), a.bam]
    proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, text=True, bufsize=1 << 20)
    cur, seq, tss_pos, tss_list = None, "", [], []
    n_frag = 0
    for line in proc.stdout:
        f = line.split("\t", 12)
        tlen = int(f[8])
        if tlen <= 0 or f[6] != "=":        # one record per pair (the leftmost mate), same chrom
            continue
        flag0 = int(f[1])
        if flag0 & 0x10 or not flag0 & 0x20:  # leftmost mate forward, its mate reverse (FR)
            continue
        chrom = f[2]
        if chrom != cur:
            cur = chrom
            seq = chrom_seq(a.fasta, fai, chrom) if chrom in fai else ""
            tss_list = tss.get(chrom, [])
            tss_pos = [t[0] for t in tss_list]
            if chrom in nbins and chrom not in bins_gc:
                bin_gc(chrom, seq)
        flag = int(f[1])
        this_clip = clip[flag & 64] if flag & 64 else clip[128]
        mate_clip = clip[128] if flag & 64 else clip[64]
        tags = "\t".join(f[11:])
        mc = next((t[5:] for t in tags.split("\t") if t.startswith("MC:Z:")), "")
        mq = next((t[5:] for t in tags.split("\t") if t.startswith("MQ:i:")), "")
        if mq and int(mq) < a.mapq:
            continue
        start = int(f[3]) - 1 - softclip(f[5], "L") - this_clip          # 0-based
        end = int(f[3]) - 1 + tlen + softclip(mc, "R") + mate_clip      # exclusive
        length = end - start
        if length < 30 or length > MAX_LEN or start < 0:
            continue
        n_frag += 1
        sizes[length] += 1
        if seq:
            b = (start + length // 2) // BIN
            if chrom not in bins_short:
                b = None
            if b is not None and SHORT[0] <= length <= SHORT[1]:
                bins_short[chrom][b] += 1
            elif b is not None and LONG[0] <= length <= LONG[1]:
                bins_long[chrom][b] += 1
            m5 = seq[start:start + 4]
            m3 = seq[end - 4:end].translate(COMP)[::-1]
            for m in (m5, m3):
                k = KIDX.get(m)
                if k is not None:
                    motifs[k] += 1
        if seq and TSS_LEN[0] <= length <= TSS_LEN[1]:
            g = frag_gc(seq[start:end])
            gbin = min(int(g * (len(GC_BINS) - 1)), len(GC_BINS) - 2) if g >= 0 else -1
            if gbin >= 0 and chrom in nbins:
                obs_gc[gbin] += 1
            if tss_pos and gbin >= 0:
                i = bisect.bisect_left(tss_pos, start - TSS_FLANK)
                while i < len(tss_pos) and tss_pos[i] <= end + TSS_FLANK:
                    t, strand, k = tss_list[i]
                    lo, hi = max(start - t, -TSS_FLANK), min(end - 1 - t, TSS_FLANK)
                    if lo <= hi:
                        if strand == "-":
                            lo, hi = -hi, -lo
                        rec_tss.append(k); rec_lo.append(lo); rec_hi.append(hi); rec_gc.append(gbin)
                    i += 1
    proc.wait()
    if proc.returncode:
        sys.exit(f"samtools view failed on {a.bam}")
    # Expected fragment GC: random N-free genomic fragments on the autosomes, lengths drawn
    # from this library's 100-220 bp size distribution. Also fills bin GC for unvisited chroms.
    rng = np.random.default_rng(0)
    lens = np.array([fai[c][0] for c in autos], float)
    sz = sizes[TSS_LEN[0]:TSS_LEN[1] + 1].astype(float)
    p_len = sz / sz.sum() if sz.sum() else np.full(len(sz), 1 / len(sz))
    exp_gc = np.zeros(len(GC_BINS) - 1, dtype=np.int64)
    n_per = rng.multinomial(N_GC_SAMPLE, lens / lens.sum()) if len(autos) else []
    for c, n_c in zip(autos, n_per):
        sq = chrom_seq(a.fasta, fai, c)
        if c not in bins_gc:
            bin_gc(c, sq)
        if n_c == 0 or len(sq) < 2 * TSS_LEN[1]:
            continue
        flen = rng.choice(np.arange(TSS_LEN[0], TSS_LEN[1] + 1), n_c, p=p_len)
        st = rng.integers(0, len(sq) - TSS_LEN[1] - 1, n_c)
        for s0, ln in zip(st, flen):
            w = sq[s0:s0 + ln]
            if "N" in w:
                continue
            exp_gc[min(int(frag_gc(w) * (len(GC_BINS) - 1)), len(GC_BINS) - 2)] += 1
    chroms = autos
    np.savez_compressed(
        a.out, label=a.label, n_fragments=n_frag, sizes=sizes, motifs=motifs,
        bin_chrom=np.array([c for c in chroms for _ in range(len(bins_short[c]))]),
        bin_index=np.concatenate([np.arange(len(bins_short[c])) for c in chroms]) if chroms else np.array([]),
        bin_short=np.concatenate([bins_short[c] for c in chroms]) if chroms else np.array([]),
        bin_long=np.concatenate([bins_long[c] for c in chroms]) if chroms else np.array([]),
        bin_gc=np.concatenate([bins_gc[c] for c in chroms]) if chroms else np.array([]),
        bin_n=np.concatenate([bins_n[c] for c in chroms]) if chroms else np.array([]),
        obs_gc=obs_gc, exp_gc=exp_gc,
        tss_chrom=np.array([t[0] for t in tss_all], dtype=str),
        tss_pos=np.array([t[1] for t in tss_all], dtype=np.int64),
        tss_strand=np.array([t[2] for t in tss_all], dtype=str),
        tss_gene=np.array([t[3] for t in tss_all], dtype=str),
        rec_tss=np.array(rec_tss, dtype=np.int32), rec_lo=np.array(rec_lo, dtype=np.int16),
        rec_hi=np.array(rec_hi, dtype=np.int16), rec_gc=np.array(rec_gc, dtype=np.int8))
    print(f"{a.label}: {n_frag:,} fragments", file=sys.stderr)


def summarize(a):
    import pandas as pd
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    out = Path(a.out)
    info = pd.read_csv(a.sample_info, sep="\t", dtype=str)
    typ = dict(zip(info["label"], info["type"]))
    libs = []
    for p in sorted(Path(a.dir).glob("*.frag.npz")):
        d = np.load(p, allow_pickle=False)
        libs.append({k: d[k] for k in d.files})
    if not libs:
        sys.exit("no *.frag.npz files")
    labels = [str(d["label"]) for d in libs]

    # ---- sizes and end motifs ------------------------------------------------------
    L = np.arange(MAX_LEN + 1)
    rows = []
    for d, lab in zip(libs, labels):
        s = d["sizes"].astype(float)
        tot = s.sum()
        cdf = np.cumsum(s) / tot
        mot = d["motifs"].astype(float)
        fr = mot / mot.sum() if mot.sum() else mot
        mds = float(-(fr[fr > 0] * np.log(fr[fr > 0])).sum() / np.log(len(KMERS))) if mot.sum() else np.nan
        top = [KMERS[i] for i in np.argsort(-fr)[:3]]
        rows.append({
            "label": lab, "type": typ.get(lab, "NA"), "fragments": int(tot),
            "median_bp": int(np.searchsorted(cdf, 0.5)), "mode_bp": int(np.argmax(s[50:]) + 50),
            "pct_lt150": 100 * s[:150].sum() / tot,
            "short_long_ratio": s[SHORT[0]:SHORT[1] + 1].sum() / max(s[LONG[0]:LONG[1] + 1].sum(), 1),
            "pct_dinucleosome_250_400": 100 * s[250:401].sum() / tot,
            "motif_diversity": mds, "CCCA_pct": 100 * fr[KIDX["CCCA"]] if mot.sum() else np.nan,
            "top_motifs": ",".join(top)})
    summ = pd.DataFrame(rows)
    summ.to_csv(out / "fragment_summary.tsv", sep="\t", index=False, float_format="%.4f")

    pd.DataFrame({lab: d["sizes"] for d, lab in zip(libs, labels)}, index=L).to_csv(
        out / "fragment_sizes.tsv", sep="\t", index_label="length_bp")
    pd.DataFrame({lab: d["motifs"] / max(d["motifs"].sum(), 1) for d, lab in zip(libs, labels)},
                 index=KMERS).to_csv(out / "end_motifs.tsv", sep="\t", index_label="motif", float_format="%.6f")

    colors = {t: c for t, c in zip(sorted(set(typ.get(l, "NA") for l in labels)),
                                   plt.rcParams["axes.prop_cycle"].by_key()["color"])}
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.5))
    for d, lab in zip(libs, labels):
        s = d["sizes"].astype(float)
        for ax, rng in zip(axes, ((50, 500), (50, 200))):
            ax.plot(L[rng[0]:rng[1]], s[rng[0]:rng[1]] / s.sum(), lw=1.2,
                    color=colors[typ.get(lab, "NA")], label=lab)
    axes[0].set_title("Fragment size distribution")
    axes[1].set_title("Short fragments (10 bp periodicity)")
    for ax in axes:
        ax.set_xlabel("fragment length (bp)")
        ax.set_ylabel("fraction of fragments")
        ax.axvline(167, color="grey", lw=0.6, ls="--")
    axes[0].legend(fontsize=7)
    fig.tight_layout()
    fig.savefig(out / "fragment_sizes.png", dpi=150)
    plt.close(fig)

    # ---- genome-wide short/long ratio (DELFI-style), GC-corrected --------------------
    d0 = libs[0]
    if len(d0["bin_short"]):
        keep = (d0["bin_n"] < 0.1) & ~np.isnan(d0["bin_gc"])
        ratios = {}
        for d, lab in zip(libs, labels):
            sh, lo = d["bin_short"].astype(float), d["bin_long"].astype(float)
            k = keep & (sh + lo > 0)
            r = np.full(len(sh), np.nan)
            r[k] = sh[k] / np.maximum(lo[k], 1)
            if k.sum() > 10:                       # remove GC trend (quadratic fit)
                co = np.polyfit(d["bin_gc"][k], r[k], 2)
                r[k] = r[k] - np.polyval(co, d["bin_gc"][k]) + np.nanmedian(r[k])
            ratios[lab] = r
        bins = pd.DataFrame({"chrom": d0["bin_chrom"], "bin_start": d0["bin_index"] * BIN})
        rat = pd.concat([bins, pd.DataFrame(ratios)], axis=1)
        rat.to_csv(out / "short_long_bins.tsv", sep="\t", index=False, float_format="%.4f")
        m = rat[labels].dropna()
        if len(m) > 5:
            # profile relative to healthy controls if present, else all libraries
            ref = [l for l in labels if str(typ.get(l, "")).lower() in ("healthy", "healthy_control", "control")]
            base = m[ref].median(axis=1) if ref else m.median(axis=1)
            dev = m.sub(base, axis=0)
            fig, ax = plt.subplots(figsize=(13, 0.6 * len(labels) + 2))
            for i, lab in enumerate(labels):
                ax.plot(np.arange(len(dev)), dev[lab].to_numpy() / max(dev.abs().max().max(), 1e-9) * 0.45 + i,
                        lw=0.8, color=colors[typ.get(lab, "NA")])
            ax.set_yticks(range(len(labels)), labels, fontsize=8)
            starts = rat.loc[m.index].reset_index(drop=True).groupby("chrom").head(1)
            ax.set_xticks(starts.index, [c.replace("chr", "") for c in starts["chrom"]], fontsize=6)
            ax.set_xlabel("chromosome (5 Mb bins)")
            ax.set_title("Short/long fragment ratio relative to " + ("healthy controls" if ref else "all libraries"))
            fig.tight_layout()
            fig.savefig(out / "short_long_profile.png", dpi=150)
            plt.close(fig)
            m.corr().to_csv(out / "short_long_correlation.tsv", sep="\t", float_format="%.3f")

    # ---- TSS nucleosome footprints (ported from cfdna_frag.py tss) ---------------------
    gsets = read_genesets(a.genesets)
    x = np.arange(-TSS_FLANK, TSS_FLANK + 1)
    W = len(x)
    base = (np.abs(x) >= TSS_NORM[0]) & (np.abs(x) <= TSS_NORM[1])
    central = (x >= TSS_CENTRAL[0]) & (x <= TSS_CENTRAL[1])
    ndr = (x >= TSS_NDR[0]) & (x <= TSS_NDR[1])
    plus1 = (x >= TSS_PLUS1[0]) & (x <= TSS_PLUS1[1])
    centers = (GC_BINS[:-1] + GC_BINS[1:]) / 2
    tss_rows, profiles = [], {}
    for d, lab in zip(libs, labels):
        if "rec_tss" not in d or len(d["tss_gene"]) == 0:
            continue
        # GC-bias weight per fragment: expected / observed fragment-GC density
        obs, exp = d["obs_gc"].astype(float), d["exp_gc"].astype(float)
        ok = (obs >= 50) & (exp >= 50)
        with np.errstate(divide="ignore", invalid="ignore"):
            ratio = (exp / exp.sum()) / (obs / obs.sum())
        wbin = np.interp(centers, centers[ok], ratio[ok]) if ok.sum() >= 2 else np.ones(len(centers))
        wbin = np.clip(wbin, 0.25, 4.0)
        wbin /= (wbin * obs).sum() / max(obs.sum(), 1)        # mean weight 1 over fragments
        rt, lo, hi = d["rec_tss"].astype(np.int64), d["rec_lo"].astype(np.int64), d["rec_hi"].astype(np.int64)
        w = wbin[d["rec_gc"].astype(np.int64)]
        n_tss = len(d["tss_gene"])
        diff = np.zeros((n_tss, W + 1), dtype=np.float64)
        np.add.at(diff, (rt, lo + TSS_FLANK), w)
        np.add.at(diff, (rt, hi + TSS_FLANK + 1), -w)
        cov = np.cumsum(diff, axis=1)[:, :W]
        del diff
        wmean = cov.mean(axis=1)
        med = np.median(wmean[wmean > 0]) if (wmean > 0).any() else 0
        good = (wmean <= TSS_OUTLIER_HIGH * med) & (wmean >= TSS_OUTLIER_LOW * med)
        cap = max(3.0, float(np.quantile(cov[good], TSS_CAP_QUANTILE))) if good.any() else np.inf
        cov = np.minimum(cov, cap)
        genes = d["tss_gene"]
        for sname, sgenes in gsets.items():
            in_set = np.isin(genes, list(sgenes))
            idx = np.flatnonzero(in_set & good)
            if len(idx) == 0:
                continue
            agg = cov[idx].sum(axis=0)
            bl = agg[base].mean()
            prof = agg / bl if bl > 0 else np.full_like(agg, np.nan)
            profiles[(lab, sname)] = prof
            tss_rows.append({"label": lab, "type": typ.get(lab, "NA"), "set": sname,
                             "n_genes_in_set": len(sgenes), "n_genes_found": len(set(genes[in_set])),
                             "n_tss": len(idx), "n_tss_excluded": int((in_set & ~good).sum()),
                             "baseline_coverage_per_tss": bl / len(idx),
                             "central_coverage": float(np.nanmean(prof[central])),
                             "ndr_coverage": float(np.nanmean(prof[ndr])),
                             "plus1_coverage": float(np.nanmean(prof[plus1]))})
        del cov
    if tss_rows:
        tdf = pd.DataFrame(tss_rows)
        inact = tdf[tdf["set"] == "inactive_control"].set_index("label")["ndr_coverage"]
        tdf["ndr_relative_to_inactive"] = [r.ndr_coverage / inact[r.label] if r.label in inact.index
                                           and inact[r.label] > 0 else np.nan for r in tdf.itertuples()]
        tdf.to_csv(out / "tss_summary.tsv", sep="\t", index=False, float_format="%.4f")
        wide = tdf.pivot(index="label", columns="set", values="ndr_relative_to_inactive")
        if "inactive_control" in wide.columns:
            wide["inactive_control_raw"] = tdf[tdf["set"] == "inactive_control"].set_index("label")["ndr_coverage"]
            wide = wide.drop(columns="inactive_control")
        wide.to_csv(out / "tss_ndr_relative.tsv", sep="\t", float_format="%.4f")
        pd.DataFrame({f"{l}|{s_}": v for (l, s_), v in profiles.items()},
                     index=pd.Index(x, name="position")).to_csv(out / "tss_profiles.tsv.gz", sep="\t",
                                                                float_format="%.4g")
        set_names = list(gsets)
        fig, axes = plt.subplots(1, len(set_names), figsize=(3.2 * len(set_names), 3.4), sharey=True)
        axes = np.atleast_1d(axes)
        for ax, sname in zip(axes, set_names):
            for lab in labels:
                if (lab, sname) in profiles:
                    ax.plot(x, profiles[(lab, sname)], lw=0.9, color=colors[typ.get(lab, "NA")], label=lab)
            n = tdf[tdf["set"] == sname]["n_tss"].max() if (tdf["set"] == sname).any() else 0
            ax.set_title(f"{sname.replace('_', ' ')} ({n} TSS)", fontsize=9)
            ax.axvline(0, color="grey", lw=0.5)
            ax.set_xlim(-1500, 1500)
            ax.set_xlabel("distance to TSS (bp)")
        axes[0].set_ylabel("coverage / flank (GC-corrected)")
        axes[0].legend(fontsize=6)
        fig.tight_layout()
        fig.savefig(out / "tss_footprints.png", dpi=150)
        plt.close(fig)

    with pd.option_context("display.width", 200, "display.max_columns", 20):
        print(summ.drop(columns=["top_motifs"]).round(3).to_string(index=False))
        if tss_rows:
            print("\nPromoter NDR coverage (-200..-20 bp) relative to each library's inactive-control set")
            print("(<1 = more open than background; last column = raw inactive-control NDR coverage):")
            print(wide.round(2).to_string())
            print("\nTSSs used per set:", dict(tdf.groupby("set")["n_tss"].max()))
    print(f"\nResults in {out}")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    e = sub.add_parser("extract")
    e.add_argument("--bam", required=True)
    e.add_argument("--label", required=True)
    e.add_argument("--fasta", required=True)
    e.add_argument("--fai", required=True)
    e.add_argument("--refseq")
    e.add_argument("--genesets", default=str(DEFAULT_GENESETS))
    e.add_argument("--clip", type=int, nargs=2, default=[0, 0], metavar=("R1", "R2"))
    e.add_argument("--mapq", type=int, default=30)
    e.add_argument("--out", required=True)
    s = sub.add_parser("summarize")
    s.add_argument("--dir", required=True)
    s.add_argument("--sample-info", required=True)
    s.add_argument("--genesets", default=str(DEFAULT_GENESETS))
    s.add_argument("--out", required=True)
    a = ap.parse_args()
    extract(a) if a.cmd == "extract" else summarize(a)


if __name__ == "__main__":
    main()
