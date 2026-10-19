#!/usr/bin/env python3
"""Compare EM-seq fragmentomics (step 09) with the xGen cfdna-fragmentomics pipeline.

Pairs are the same cfDNA sample sequenced with both library methods. Every size and
end-motif metric is recomputed here from each pipeline's raw tables with the same code
and the same size range (default 50-600 bp, the xGen pipeline's), so differences come
from the libraries, not from how each pipeline summarises them.

Inputs
  --emseq   step 09 output: fragment_sizes.tsv, end_motifs.tsv, short_long_bins.tsv,
            tss_ndr_relative.tsv
  --xgen    xGen results dir: features/<sample>/{size_hist.tsv,end_motifs.tsv,bins.tsv},
            cohort/tss_ndr_relative.tsv
  --pairs   TSV: emseq_label, xgen_sample (header). xGen samples not in a pair (e.g.
            healthy plasma) are reported as context.

Outputs (--out)
  metrics.tsv          size and end-motif metrics per library, both methods
  agreement.tsv        per pair: correlation of size distributions (50-600 bp), of the
                       256 end-motif frequencies, and of GC-corrected short/long ratios in
                       5 Mb bins
  tss_ndr.tsv          promoter NDR (relative to inactive control) per gene set, both methods
  size_overlay.png     size distributions per pair
"""
import argparse
import itertools
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

MOTIFS = ["".join(p) for p in itertools.product("ACGT", repeat=4)]
SHORT, LONG = (100, 150), (151, 220)


def size_metrics(h, lo, hi):
    """h: counts indexed by length (bp)."""
    s = h.reindex(range(lo, hi + 1), fill_value=0).astype(float)
    tot = s.sum()
    cdf = s.cumsum() / tot
    return {"fragments": int(tot),
            "median_bp": int(cdf.index[np.searchsorted(cdf.values, 0.5)]),
            "mode_bp": int(s.idxmax()),
            "pct_lt150": 100 * s.loc[:149].sum() / tot,
            "short_long_ratio": s.loc[SHORT[0]:SHORT[1]].sum() / max(s.loc[LONG[0]:LONG[1]].sum(), 1),
            "pct_250_400": 100 * s.loc[250:400].sum() / tot}


def motif_metrics(f):
    f = f.reindex(MOTIFS, fill_value=0).astype(float)
    f = f / f.sum()
    p = f[f > 0].values
    return {"motif_diversity": float(-(p * np.log(p)).sum() / np.log(len(MOTIFS))),
            "CCCA_pct": 100 * f["CCCA"], "CCAG_pct": 100 * f["CCAG"], "TGGG_pct": 100 * f["TGGG"]}, f


def load_xgen(xdir, name):
    d = Path(xdir) / "features" / name
    h = pd.read_csv(d / "size_hist.tsv", sep="\t").set_index("length")["count"]
    m = pd.read_csv(d / "end_motifs.tsv", sep="\t", index_col=0)["freq"]
    b = pd.read_csv(d / "bins.tsv", sep="\t", dtype={"chrom": str})
    b = b[b["pass"].astype(str).str.lower().isin(["true", "1"])]
    r = pd.Series(b["ratio"].values, index=b["chrom"] + ":" + b["start"].astype(str))
    return h, m, r


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--emseq", required=True)
    ap.add_argument("--xgen", required=True)
    ap.add_argument("--pairs", required=True)
    ap.add_argument("--min-len", type=int, default=50)
    ap.add_argument("--max-len", type=int, default=600)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    em, xg = Path(a.emseq), Path(a.xgen)
    pairs = pd.read_csv(a.pairs, sep="\t", dtype=str, comment="#")

    em_sizes = pd.read_csv(em / "fragment_sizes.tsv", sep="\t", index_col=0)
    em_motifs = pd.read_csv(em / "end_motifs.tsv", sep="\t", index_col=0)
    em_bins = None
    if (em / "short_long_bins.tsv").exists():
        b = pd.read_csv(em / "short_long_bins.tsv", sep="\t", dtype={"chrom": str})
        b.index = b["chrom"] + ":" + b["bin_start"].astype(str)
        em_bins = b
    xg_samples = sorted(p.name for p in (xg / "features").iterdir()
                        if (p / "size_hist.tsv").exists())

    rows, agree, motif_freqs, hists = [], [], {}, {}
    paired_x = set(pairs["xgen_sample"])
    for lab in em_sizes.columns:
        mm, f = motif_metrics(em_motifs[lab])
        rows.append({"library": lab, "method": "EM-seq", **size_metrics(em_sizes[lab], a.min_len, a.max_len), **mm})
        motif_freqs[("EM-seq", lab)], hists[("EM-seq", lab)] = f, em_sizes[lab]
    xg_bins = {}
    for s in xg_samples:
        h, m, r = load_xgen(xg, s)
        mm, f = motif_metrics(m)
        rows.append({"library": s, "method": "xGen" + ("" if s in paired_x else " (context)"),
                     **size_metrics(h, a.min_len, a.max_len), **mm})
        motif_freqs[("xGen", s)], hists[("xGen", s)], xg_bins[s] = f, h, r
    met = pd.DataFrame(rows)
    met.to_csv(out / "metrics.tsv", sep="\t", index=False, float_format="%.4f")

    rng = range(a.min_len, a.max_len + 1)
    for p in pairs.itertuples():
        if p.emseq_label not in em_sizes.columns or p.xgen_sample not in xg_bins:
            print(f"skipping {p.emseq_label} / {p.xgen_sample}: missing on one side")
            continue
        he = em_sizes[p.emseq_label].reindex(rng, fill_value=0).astype(float)
        hx = hists[("xGen", p.xgen_sample)].reindex(rng, fill_value=0).astype(float)
        row = {"emseq_label": p.emseq_label, "xgen_sample": p.xgen_sample,
               "size_dist_r": np.corrcoef(he / he.sum(), hx / hx.sum())[0, 1],
               "end_motif_r": np.corrcoef(motif_freqs[("EM-seq", p.emseq_label)],
                                          motif_freqs[("xGen", p.xgen_sample)])[0, 1]}
        if em_bins is not None and p.emseq_label in em_bins.columns:
            j = pd.concat([em_bins[p.emseq_label], xg_bins[p.xgen_sample]], axis=1, join="inner").dropna()
            row["short_long_bins_r"] = j.corr().iloc[0, 1] if len(j) > 10 else np.nan
            row["bins_compared"] = len(j)
        agree.append(row)
    ag = pd.DataFrame(agree)
    ag.to_csv(out / "agreement.tsv", sep="\t", index=False, float_format="%.3f")

    tss = []
    for path, method in ((em / "tss_ndr_relative.tsv", "EM-seq"), (xg / "cohort" / "tss_ndr_relative.tsv", "xGen")):
        if path.exists():
            t = pd.read_csv(path, sep="\t", index_col=0).drop(columns=["group", "inactive_control_raw"], errors="ignore")
            t = t.stack().rename("ndr_rel").reset_index()
            t.columns = ["library", "gene_set", "ndr_rel"]
            t["method"] = method
            tss.append(t)
    if tss:
        tss = pd.concat(tss)
        lab2pair = {**{r.emseq_label: r.emseq_label for r in pairs.itertuples()},
                    **{r.xgen_sample: r.emseq_label for r in pairs.itertuples()}}
        tss["pair"] = tss["library"].map(lab2pair)
        wide = (tss.dropna(subset=["pair"]).pivot_table(index=["pair", "gene_set"], columns="method",
                                                        values="ndr_rel").reset_index())
        wide.to_csv(out / "tss_ndr.tsv", sep="\t", index=False, float_format="%.3f")

    n = len(ag)
    if n:
        fig, axes = plt.subplots(1, n, figsize=(5 * n, 3.8), squeeze=False)
        for ax, p in zip(axes[0], ag.itertuples()):
            for (method, lab), col in ((("EM-seq", p.emseq_label), "#1f77b4"), (("xGen", p.xgen_sample), "#d62728")):
                h = hists[(method, lab)].reindex(rng, fill_value=0).astype(float)
                ax.plot(list(rng), h / h.sum(), color=col, lw=1.2, label=f"{method} {lab}")
            ax.axvline(167, color="grey", lw=0.6, ls="--")
            ax.set_xlim(a.min_len, 400)
            ax.set_xlabel("fragment length (bp)"); ax.set_ylabel("fraction")
            ax.set_title(f"{p.emseq_label}  (r = {p.size_dist_r:.3f})", fontsize=9)
            ax.legend(fontsize=7)
        fig.tight_layout()
        fig.savefig(out / "size_overlay.png", dpi=150)
        plt.close(fig)

    pd.set_option("display.width", 200)
    cols = ["library", "method", "fragments", "median_bp", "mode_bp", "pct_lt150", "short_long_ratio",
            "pct_250_400", "CCCA_pct", "motif_diversity"]
    print(f"Size and end-motif metrics ({a.min_len}-{a.max_len} bp, same code for both methods)")
    print(met[cols].to_string(index=False, float_format=lambda v: f"{v:.3f}"))
    if n:
        print("\nAgreement per pair (Pearson r)")
        print(ag.to_string(index=False, float_format=lambda v: f"{v:.3f}"))
    if isinstance(tss, pd.DataFrame) and len(tss):
        print("\nPromoter NDR relative to inactive control (< 1 = more open)")
        print(wide.to_string(index=False, float_format=lambda v: f"{v:.3f}"))


if __name__ == "__main__":
    main()
