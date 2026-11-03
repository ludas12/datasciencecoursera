#!/usr/bin/env python3
"""Rough chromosome-arm copy-number screen from the cohort's binned coverage (no ichorCNA needed).

Usage:  python scripts/arm_cna.py RESULTS/cohort [--plot CF0009A_BAL CF0009A_P] [-o RESULTS/cohort/arm_cna]

Reads cohort/coverage_norm.tsv (GC-corrected coverage per bin, each sample scaled to its median)
and the group column of cohort/size_motif_summary.tsv. For each autosome arm and sample:
  arm log2 = mean log2(coverage) over the arm's bins, minus the median of the comparison panel;
  z = that difference / the panel's spread for that arm.
The panel is every sample of the same compartment (plasma or BAL) from OTHER recipients, so a
sample is never compared with itself or with its own other visits. Healthy plasma is in the
plasma panel. A tumour contributing a fraction f of the DNA shifts a one-copy gain by about
log2(1 + f/2) and a one-copy loss by about log2(1 - f/2): at f = 10%, +0.07 and -0.07.

This is a screen at 5 Mb resolution, not a tumour-fraction estimate. Arms with |z| >= 3 and a shift of at least --min-shift in a
sample of interest, with few or none in the other samples, are the signal; every sample is
scored so the false-positive rate is visible. With only 6 other BAL samples (5 from one
recipient), the BAL panel's spread is rough.
Outputs: <out>_z.tsv (z per arm and sample), <out>_log2.tsv (arm log2 differences),
<out>_summary.tsv (arms called per sample) and <out>.png when --plot is given.
"""
import argparse

import numpy as np
import pandas as pd

# Approximate hg38 centromere midpoints (Mb); acrocentric arms (13, 14, 15, 21, 22) keep q only.
CENTROMERE_MB = {"chr1": 123.4, "chr2": 93.9, "chr3": 90.9, "chr4": 50.0, "chr5": 48.8, "chr6": 59.8,
                 "chr7": 60.1, "chr8": 45.2, "chr9": 43.0, "chr10": 39.8, "chr11": 53.4, "chr12": 35.5,
                 "chr13": 17.7, "chr14": 17.2, "chr15": 19.0, "chr16": 36.8, "chr17": 25.1, "chr18": 18.5,
                 "chr19": 26.2, "chr20": 28.1, "chr21": 12.0, "chr22": 15.0}
ACROCENTRIC = {"chr13", "chr14", "chr15", "chr21", "chr22"}
MIN_BINS = 3
SD_FLOOR = 0.01  # log2; small panels can understate the spread


def arm_of(key):
    chrom, start = key.rsplit(":", 1)
    if chrom not in CENTROMERE_MB:
        return None
    arm = "p" if int(start) / 1e6 < CENTROMERE_MB[chrom] else "q"
    if arm == "p" and chrom in ACROCENTRIC:
        return None
    return f"{chrom[3:]}{arm}"


def recipient(sample):
    return sample.split("_")[0][:6]


def compartment(sample, group):
    return "bal" if group.get(sample, "").endswith("bal") or sample.endswith("_BAL") else "plasma"


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawTextHelpFormatter)
    ap.add_argument("cohort")
    ap.add_argument("--plot", nargs="*", default=[], help="samples to draw as genome plots")
    ap.add_argument("-o", "--out", help="output prefix (default: <cohort>/arm_cna)")
    ap.add_argument("--min-shift", type=float, default=0.03,
                    help="smallest arm log2 shift reported (0.03 ~ a one-copy change at 4%% tumour DNA)")
    a = ap.parse_args()
    out = a.out or f"{a.cohort.rstrip('/')}/arm_cna"
    cov = pd.read_csv(f"{a.cohort}/coverage_norm.tsv", sep="\t", index_col=0)
    grp = pd.read_csv(f"{a.cohort}/size_motif_summary.tsv", sep="\t", index_col=0)["group"].to_dict()
    arms = pd.Series({k: arm_of(k) for k in cov.columns}).dropna()
    lg = np.log2(cov[arms.index].where(cov[arms.index] > 0))
    arm_vals = lg.T.groupby(arms).agg(lambda x: x.mean() if x.notna().sum() >= MIN_BINS else np.nan).T
    order = sorted(arm_vals.columns, key=lambda s: (int(s[:-1]), s[-1]))
    arm_vals = arm_vals[order]

    z = pd.DataFrame(index=arm_vals.index, columns=order, dtype=float)
    d = pd.DataFrame(index=arm_vals.index, columns=order, dtype=float)
    panel_size = {}
    for s in arm_vals.index:
        comp = compartment(s, grp)
        panel = [t for t in arm_vals.index if t != s and recipient(t) != recipient(s) and compartment(t, grp) == comp]
        panel_size[s] = len(panel)
        if len(panel) < 3:
            continue
        ref = arm_vals.loc[panel]
        med = ref.median()
        sd = (ref.quantile(0.84) - ref.quantile(0.16)) / 2  # robust spread
        sd = sd.where(sd > 0, ref.std()).clip(lower=SD_FLOOR)
        d.loc[s] = arm_vals.loc[s] - med
        z.loc[s] = (arm_vals.loc[s] - med) / sd

    rows = []
    for s in z.index:
        if panel_size[s] < 3:
            rows.append({"sample": s, "compartment": compartment(s, grp), "panel_n": panel_size[s],
                         "arms_called": np.nan, "arms": "not scored: fewer than 3 panel samples"})
            continue
        keep = (z.loc[s].abs() >= 3) & (d.loc[s].abs() >= a.min_shift)
        hits = z.loc[s][keep].dropna()
        rows.append({"sample": s, "compartment": compartment(s, grp), "panel_n": panel_size[s],
                     "arms_called": len(hits),
                     "arms": ", ".join(f"{k} {'+' if d.loc[s, k] > 0 else '-'}{abs(d.loc[s, k]):.3f} (z {z.loc[s, k]:.1f})"
                                       for k in hits.index)})
    summ = pd.DataFrame(rows).sort_values(["compartment", "arms_called"], ascending=[True, False])
    z.to_csv(f"{out}_z.tsv", sep="\t", float_format="%.2f")
    d.to_csv(f"{out}_log2.tsv", sep="\t", float_format="%.4f")
    summ.to_csv(f"{out}_summary.tsv", sep="\t", index=False)
    pd.set_option("display.width", 200, "display.max_colwidth", 120)
    print(summ.to_string(index=False))
    print(f"\nWrote {out}_z.tsv, {out}_log2.tsv, {out}_summary.tsv")

    if a.plot:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        keys = list(arms.index)
        chroms = [k.rsplit(":", 1)[0] for k in keys]
        x = np.arange(len(keys))
        fig, axes = plt.subplots(len(a.plot), 1, figsize=(10, 2.2 * len(a.plot)), sharex=True, squeeze=False)
        for ax, s in zip(axes[:, 0], a.plot):
            if s not in lg.index:
                ax.set_title(f"{s}: not found")
                continue
            comp = compartment(s, grp)
            panel = [t for t in lg.index if t != s and recipient(t) != recipient(s) and compartment(t, grp) == comp]
            dev = lg.loc[s, keys] - lg.loc[panel, keys].median()
            ax.scatter(x, dev, s=3, color="#2f6db3" if comp == "bal" else "#555555")
            for arm in order:
                m = (arms == arm).values
                ax.hlines(d.loc[s, arm], x[m].min(), x[m].max(), color="#d62728", lw=1.5)
            for i in range(1, len(keys)):
                if chroms[i] != chroms[i - 1]:
                    ax.axvline(i - 0.5, color="#dddddd", lw=0.5)
            ax.axhline(0, color="#999999", lw=0.5)
            ax.set_ylim(-0.4, 0.4)
            ax.set_ylabel("log2 vs panel", fontsize=8)
            ax.set_title(f"{s} (panel: {len(panel)} {comp} samples from other recipients; red = arm mean)", fontsize=8)
        starts = [i for i in range(len(keys)) if i == 0 or chroms[i] != chroms[i - 1]]
        axes[-1, 0].set_xticks(starts, [chroms[i][3:] for i in starts], fontsize=7)
        fig.tight_layout()
        fig.savefig(f"{out}.png", dpi=200)
        print(f"Wrote {out}.png")


if __name__ == "__main__":
    main()
