#!/usr/bin/env python3
"""cfDNA fragmentomics from paired-end EM-seq alignments.

extract   one library -> OUT/<label>.frag.npz
          fragments = primary, non-duplicate pairs, both mates MAPQ >= 30, same chromosome,
          FR orientation, <= 1,000 bp (the aligner's proper-pair flag is not required).
          True fragment ends are restored when reads were clipped before alignment
          (--clip R1 R2) and from soft clips (own CIGAR and the mate's MC tag).
          Collects: size histogram; short (100-150 bp) / long (151-220 bp) counts per
          5 Mb bin with bin GC; fragment coverage around TSSs of marker gene sets;
          4-mer end motifs read from the reference (EM-seq conversion does not affect
          them).
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
TSS_FLANK = 2000
COMP = str.maketrans("ACGTN", "TGCAN")
BASES = "ACGT"
KMERS = [a + b + c + d for a in BASES for b in BASES for c in BASES for d in BASES]
KIDX = {k: i for i, k in enumerate(KMERS)}

# Marker gene sets for promoter nucleosome footprints. A deep, uniform dip at the TSS
# means the promoter was open (gene active) in the cells the cfDNA came from.
GENE_SETS = {
    "lung_epithelium": ["SFTPC", "SFTPB", "SFTPA1", "SFTPA2", "SFTPD", "NAPSA", "NKX2-1", "LAMP3",
                        "ABCA3", "SLC34A2", "AGER", "HOPX", "EMP2", "SCGB1A1", "SCGB3A2", "SCGB3A1",
                        "FOXJ1", "KRT5", "TP63", "MUC5B", "BPIFA1", "CAPS", "RSPH1", "CLDN18"],
    "neutrophil": ["ELANE", "MPO", "CSF3R", "FCGR3B", "CXCR2", "CXCR1", "S100A8", "S100A9",
                   "S100A12", "FPR1", "MME", "CEACAM8", "LTF", "CAMP", "MMP9", "AQP9", "PADI4"],
    "megakaryocyte": ["PF4", "PPBP", "ITGA2B", "GP9", "GP1BA", "GP6", "TUBB1", "VWF", "SELP"],
    "liver": ["ALB", "APOA1", "APOA2", "APOB", "APOC3", "FGA", "FGB", "FGG", "TF", "SERPINA1",
              "TTR", "HP", "AHSG", "CYP3A4", "CYP2E1"],
    "housekeeping": ["ACTB", "GAPDH", "RPLP0", "PGK1", "TBP", "HPRT1", "B2M", "PPIA", "YWHAZ",
                     "UBC", "EEF1A1", "TUBB", "RPL13A", "SDHA", "HMBS", "GUSB"],
}


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


def load_tss(refseq):
    """{chrom: sorted [(tss, strand, set)]} for genes in GENE_SETS, one TSS per transcript."""
    want = {g: s for s, genes in GENE_SETS.items() for g in genes}
    tss = defaultdict(set)
    if refseq and Path(refseq).exists():
        with gzip.open(refseq, "rt") as fh:
            for line in fh:
                f = line.rstrip("\n").split("\t")
                if len(f) > 12 and f[12] in want:
                    pos = int(f[4]) if f[3] == "+" else int(f[5])
                    tss[f[2]].add((pos, f[3], want[f[12]]))
    return {c: sorted(v) for c, v in tss.items()}


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
    tss = load_tss(a.refseq)
    sets = list(GENE_SETS)
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
    tss_cov = {s: np.zeros(2 * TSS_FLANK + 1, dtype=np.float64) for s in sets}
    tss_n = {s: 0 for s in sets}
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
        if tss_pos and 100 <= length <= 220:
            i = bisect.bisect_left(tss_pos, start - TSS_FLANK)
            while i < len(tss_pos) and tss_pos[i] <= end + TSS_FLANK:
                t, strand, gset = tss_list[i]
                lo, hi = max(start - t, -TSS_FLANK), min(end - 1 - t, TSS_FLANK)
                if lo <= hi:
                    if strand == "-":
                        lo, hi = -hi, -lo
                    tss_cov[gset][lo + TSS_FLANK:hi + TSS_FLANK + 1] += 1
                i += 1
    proc.wait()
    if proc.returncode:
        sys.exit(f"samtools view failed on {a.bam}")
    for c in autos:                         # GC for chromosomes without fragments
        if c not in bins_gc:
            bin_gc(c, chrom_seq(a.fasta, fai, c))
    for gset in sets:
        tss_n[gset] = sum(1 for v in tss.values() for t in v if t[2] == gset)
    chroms = autos
    np.savez_compressed(
        a.out, label=a.label, n_fragments=n_frag, sizes=sizes, motifs=motifs,
        bin_chrom=np.array([c for c in chroms for _ in range(len(bins_short[c]))]),
        bin_index=np.concatenate([np.arange(len(bins_short[c])) for c in chroms]) if chroms else np.array([]),
        bin_short=np.concatenate([bins_short[c] for c in chroms]) if chroms else np.array([]),
        bin_long=np.concatenate([bins_long[c] for c in chroms]) if chroms else np.array([]),
        bin_gc=np.concatenate([bins_gc[c] for c in chroms]) if chroms else np.array([]),
        bin_n=np.concatenate([bins_n[c] for c in chroms]) if chroms else np.array([]),
        tss_sets=np.array(sets), tss_cov=np.stack([tss_cov[s] for s in sets]),
        tss_n=np.array([tss_n[s] for s in sets]))
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

    # ---- TSS nucleosome footprints --------------------------------------------------
    sets = [str(s) for s in d0["tss_sets"]]
    x = np.arange(-TSS_FLANK, TSS_FLANK + 1)
    flank = (np.abs(x) >= 1500)
    center = (np.abs(x) <= 150)
    dip_rows = []
    fig, axes = plt.subplots(1, len(sets), figsize=(3.2 * len(sets), 3.4), sharey=True)
    for si, sname in enumerate(sets):
        ax = axes[si] if len(sets) > 1 else axes
        for d, lab in zip(libs, labels):
            cov = d["tss_cov"][si]
            fl = cov[flank].mean()
            norm = cov / fl if fl > 0 else cov * np.nan
            ax.plot(x, norm, lw=1, color=colors[typ.get(lab, "NA")], label=lab)
            dip_rows.append({"label": lab, "gene_set": sname, "genes_tss": int(d["tss_n"][si]),
                             "flank_coverage": fl,
                             "center_over_flank": cov[center].mean() / fl if fl > 0 else np.nan})
        ax.set_title(sname.replace("_", " "), fontsize=9)
        ax.axvline(0, color="grey", lw=0.5)
        ax.set_xlabel("distance to TSS (bp)")
    (axes[0] if len(sets) > 1 else axes).set_ylabel("coverage / flank")
    (axes[0] if len(sets) > 1 else axes).legend(fontsize=6)
    fig.tight_layout()
    fig.savefig(out / "tss_footprints.png", dpi=150)
    plt.close(fig)
    pd.DataFrame(dip_rows).to_csv(out / "tss_footprints.tsv", sep="\t", index=False, float_format="%.4f")

    with pd.option_context("display.width", 200, "display.max_columns", 20):
        print(summ.drop(columns=["top_motifs"]).round(3).to_string(index=False))
        dip = pd.DataFrame(dip_rows).pivot(index="label", columns="gene_set", values="center_over_flank")
        print("\nTSS coverage dip (centre/flank; lower = more open promoters in the source cells):")
        print(dip.round(2).to_string())
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
    e.add_argument("--clip", type=int, nargs=2, default=[0, 0], metavar=("R1", "R2"))
    e.add_argument("--mapq", type=int, default=30)
    e.add_argument("--out", required=True)
    s = sub.add_parser("summarize")
    s.add_argument("--dir", required=True)
    s.add_argument("--sample-info", required=True)
    s.add_argument("--out", required=True)
    a = ap.parse_args()
    extract(a) if a.cmd == "extract" else summarize(a)


if __name__ == "__main__":
    main()
