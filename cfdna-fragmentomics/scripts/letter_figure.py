#!/usr/bin/env python3
"""Draw the research letter's figure (panels A-D) and write the numbers behind it.

Usage:
    python scripts/letter_figure.py RESULTS/cohort [--emseq EMSEQ/downstream/compartment/all_pairs/size_profile.tsv]
                                    [-o letter_figure]

Inputs (all written by the pipeline):
  RESULTS/cohort/size_motif_summary.tsv   per-sample features and group
  RESULTS/cohort/motif_scale_pairs.tsv    whole-sample end-motif divergences (scripts/motif_scale.py)
  RESULTS/features/<sample>/size_hist.tsv fragment-length histogram per sample
  --emseq size_profile.tsv                EM-seq methylation by fragment-length class (all pairs);
                                          panel D is left out if not given
Outputs: <out>.pdf, <out>.png (300 dpi) and <out>_stats.tsv (the values and P values shown).

Panels:
  A  fragment-length density, 50-600 bp: each BAL sample over its same-visit plasma
  B  paired features (167 bp peak, nucleosome ladder, CCCA ends, motif diversity), lines join
     each visit's BAL and plasma; exact two-sided Wilcoxon signed-rank P
  C  end-motif divergence: BAL vs same-visit plasma, plasma vs plasma, BAL vs BAL
  D  EM-seq methylation of each length class relative to 150-220 bp, in regions highly
     methylated in reference DNA (ratio_high), by sample type
Recipients are labelled R1, R2, R3 in the order of RECIPIENT_LABELS, as in the letter.
"""
import argparse
import os

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from scipy.stats import wilcoxon  # noqa: E402

RECIPIENT_LABELS = {"CF0020": "R1", "CF0007": "R2", "CF0009": "R3"}
BAL_C, PL_C, REF_C, INK = "#2f6db3", "#8a8a8a", "#c4c4c4", "#222222"
MARKERS = {"R1": "o", "R2": "s", "R3": "^"}
FEATURES = [("peak167_prominence", "167 bp peak\nprominence"),
            ("nucleosome_ladder_score", "Nucleosome ladder\n(di-, tri-)"),
            ("motif_CCCA", "CCCA end\nfrequency"),
            ("motif_diversity_score", "End-motif\ndiversity")]
CLASS_ORDER = ["0-100", "100-150", "150-220", "220-400", "400-inf"]
CLASS_LABEL = ["<100", "100-150", "150-220", "220-400", ">400"]


def pairs_from(summary):
    bal = [s for s in summary.index if s.endswith("_BAL") and s[:-4] + "_P" in summary.index]
    rows = []
    for b in bal:
        visit = b[:-4]
        r = RECIPIENT_LABELS.get(visit[:6], visit[:6])
        rows.append({"visit": visit, "recipient": r, "label": f"{r}-{visit[6:]}", "bal": b, "plasma": visit + "_P"})
    return pd.DataFrame(rows).sort_values(["recipient", "visit"]).reset_index(drop=True)


def density(path, lo=50, hi=600):
    h = pd.read_csv(path, sep="\t").set_index("length")["count"].reindex(range(lo, hi + 1), fill_value=0)
    return h / h.sum()


def style(ax):
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    ax.tick_params(labelsize=7, length=2.5)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawTextHelpFormatter)
    ap.add_argument("cohort")
    ap.add_argument("--features", help="per-sample feature folder (default: ../features next to cohort)")
    ap.add_argument("--emseq", help="EM-seq size_profile.tsv (all pairs) for panel D")
    ap.add_argument("-o", "--out", default="letter_figure")
    a = ap.parse_args()
    feat_dir = a.features or os.path.join(os.path.dirname(os.path.abspath(a.cohort.rstrip("/"))), "features")
    s = pd.read_csv(f"{a.cohort}/size_motif_summary.tsv", sep="\t", index_col=0)
    pr = pairs_from(s)
    if pr.empty:
        raise SystemExit("no BAL/plasma pairs (<visit>_BAL with <visit>_P) in size_motif_summary.tsv")
    stats = []

    plt.rcParams.update({"font.size": 7, "font.family": "sans-serif", "axes.linewidth": 0.6})
    n = len(pr)
    fig = plt.figure(figsize=(7.2, 7.6))
    gs = fig.add_gridspec(3, 1, height_ratios=[1.0, 1.15, 1.25], hspace=0.55)

    # A: size distributions, one small plot per pair
    ga = gs[0].subgridspec(1, n, wspace=0.12)
    for i, p in pr.iterrows():
        ax = fig.add_subplot(ga[0, i])
        for smp, c, lw in ((p.plasma, PL_C, 0.9), (p.bal, BAL_C, 1.1)):
            f = f"{feat_dir}/{smp}/size_hist.tsv"
            if os.path.exists(f):
                d = density(f)
                ax.plot(d.index, d.values * 100, color=c, lw=lw)
        ax.set_xlim(50, 600)
        ax.set_xticks([167, 334, 500])
        ax.set_title(p.label, fontsize=7, pad=2)
        style(ax)
        if i == 0:
            ax.set_ylabel("Fragments (% per bp)", fontsize=7)
            ax.text(-0.55, 1.18, "A", transform=ax.transAxes, fontsize=10, fontweight="bold")
        else:
            ax.set_yticklabels([])
        if i == n // 2:
            ax.set_xlabel("Fragment length (bp)", fontsize=7)
    fig.legend(handles=[plt.Line2D([], [], color=BAL_C, lw=1.2, label="BAL"),
                        plt.Line2D([], [], color=PL_C, lw=1.0, label="Plasma, same visit")],
               loc="upper right", bbox_to_anchor=(0.98, 0.995), ncol=2, frameon=False, fontsize=7)

    # B: paired features
    gb = gs[1].subgridspec(1, len(FEATURES), wspace=0.6)
    for j, (col, label) in enumerate(FEATURES):
        ax = fig.add_subplot(gb[0, j])
        b, pl = s.loc[pr.bal, col].values, s.loc[pr.plasma, col].values
        for k, p in pr.iterrows():
            ax.plot([0, 1], [pl[k], b[k]], color="#b0b0b0", lw=0.6, zorder=1)
            ax.scatter([0], [pl[k]], color=PL_C, marker=MARKERS.get(p.recipient, "o"), s=14, zorder=2)
            ax.scatter([1], [b[k]], color=BAL_C, marker=MARKERS.get(p.recipient, "o"), s=14, zorder=2)
        pval = wilcoxon(b, pl, alternative="two-sided", method="exact").pvalue if n >= 2 else np.nan
        stats.append({"panel": "B", "measure": col, "BAL_min": b.min(), "BAL_max": b.max(),
                      "plasma_min": pl.min(), "plasma_max": pl.max(), "n_pairs": n,
                      "BAL_gt_plasma": int((b > pl).sum()), "P_wilcoxon": pval})
        ax.set_xticks([0, 1], ["Plasma", "BAL"])
        ax.set_xlim(-0.4, 1.4)
        ax.set_title(f"{label}\nP = {pval:.3f}", fontsize=7, pad=3)
        style(ax)
        if j == 0:
            ax.text(-0.75, 1.25, "B", transform=ax.transAxes, fontsize=10, fontweight="bold")
    ax.legend(handles=[plt.Line2D([], [], color="#555555", marker=m, ls="", ms=4, label=r)
                       for r, m in MARKERS.items() if r in set(pr.recipient)],
              loc="center left", bbox_to_anchor=(1.02, 0.5), frameon=False, fontsize=7, title="Recipient",
              title_fontsize=7)

    # C and D side by side
    gcd = gs[2].subgridspec(1, 2, wspace=0.35, width_ratios=[1, 1.25])
    ax = fig.add_subplot(gcd[0, 0])
    mp = f"{a.cohort}/motif_scale_pairs.tsv"
    if os.path.exists(mp):
        m = pd.read_csv(mp, sep="\t")
        order = [("BAL vs plasma, same visit", "BAL vs plasma\nsame visit", BAL_C),
                 ("plasma vs plasma", "Plasma vs\nplasma", PL_C), ("BAL vs BAL", "BAL vs\nBAL", "#7fa6d6")]
        rng = np.random.default_rng(0)
        for x, (key, lab, c) in enumerate(order):
            v = m.loc[m.comparison == key, "js"].values
            if not len(v):
                continue
            ax.scatter(x + rng.uniform(-0.18, 0.18, len(v)), v, s=6, color=c, alpha=0.75, lw=0)
            ax.plot([x - 0.28, x + 0.28], [np.median(v)] * 2, color=INK, lw=1.0)
            stats.append({"panel": "C", "measure": key, "n_pairs": len(v), "min": v.min(),
                          "median": np.median(v), "max": v.max()})
        ax.set_yscale("log")
        ax.set_xticks(range(len(order)), [o[1] for o in order])
        ax.set_ylabel("End-motif divergence (JSD, log scale)", fontsize=7)
    else:
        ax.text(0.5, 0.5, "run scripts/motif_scale.py first", ha="center", transform=ax.transAxes)
    style(ax)
    ax.text(-0.3, 1.06, "C", transform=ax.transAxes, fontsize=10, fontweight="bold")

    ax = fig.add_subplot(gcd[0, 1])
    if a.emseq and os.path.exists(a.emseq):
        e = pd.read_csv(a.emseq, sep="\t")
        e["ratio_high"] = np.where(e.size_class == "150-220", 1.0, e.ratio_high)
        e = e[e.size_class.isin(CLASS_ORDER) & e.ratio_high.notna()]
        x = {c: i for i, c in enumerate(CLASS_ORDER)}
        for typ, lab, c, ls in (("reference", "Reference DNA", REF_C, "--"), ("plasma", "Plasma", PL_C, "-"),
                                ("bal", "BAL", BAL_C, "-")):
            t = e[e.type.str.lower() == typ].groupby("size_class").ratio_high.agg(["mean", "min", "max"])
            t = t.reindex([k for k in CLASS_ORDER if k in t.index])
            xs = [x[k] for k in t.index]
            ax.plot(xs, t["mean"], color=c, ls=ls, lw=1.4, marker="o", ms=3.5, label=lab)
            ax.vlines(xs, t["min"], t["max"], color=c, lw=1.2)
            for k, r in t.iterrows():
                stats.append({"panel": "D", "measure": f"{typ} {k}", "mean": r["mean"], "min": r["min"],
                              "max": r["max"]})
        ax.axhline(1, color="#999999", lw=0.6, ls=":")
        ax.set_xticks(range(len(CLASS_ORDER)), CLASS_LABEL)
        ax.set_xlabel("Fragment length (bp)", fontsize=7)
        ax.set_ylabel("Methylation vs 150-220 bp\n(highly methylated regions)", fontsize=7)
        ax.legend(frameon=False, fontsize=7, loc="lower left")
    else:
        ax.text(0.5, 0.5, "EM-seq size_profile.tsv not given (--emseq)", ha="center", transform=ax.transAxes,
                fontsize=7)
        ax.set_xticks([])
        ax.set_yticks([])
    style(ax)
    ax.text(-0.28, 1.06, "D", transform=ax.transAxes, fontsize=10, fontweight="bold")

    fig.savefig(f"{a.out}.pdf", bbox_inches="tight")
    fig.savefig(f"{a.out}.png", dpi=300, bbox_inches="tight")
    pd.DataFrame(stats).to_csv(f"{a.out}_stats.tsv", sep="\t", index=False, float_format="%.5g")
    print(f"Wrote {a.out}.pdf, {a.out}.png and {a.out}_stats.tsv ({n} pairs)")


if __name__ == "__main__":
    main()
