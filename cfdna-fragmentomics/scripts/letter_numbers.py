#!/usr/bin/env python3
"""Print every number the research letter quotes, from a cohort output folder.

Usage:  python scripts/letter_numbers.py RESULTS/cohort

Reads size_motif_summary.tsv (and size_class_motifs.tsv, motif_scale_*.tsv when present).
Pairs are BAL and plasma from the same visit: <visit>_BAL with <visit>_P. The recipient is the
visit ID without its final letter (CF0020A -> CF0020).

Sections:
 1. Paired BAL vs plasma: ranges, exact Wilcoxon signed-rank P, and whether every pair and
    every recipient goes the same way.
 2. Healthy vs transplant plasma: ranges only (3 controls; descriptive).
 3. chrY fraction per sample, and a donor-DNA estimate for CF0012 (male donor, female recipient).
 4. Long size-class end-motif excess divergence (needs size_class_motifs.tsv).
"""
import os
import sys

import pandas as pd
from scipy.stats import wilcoxon

FEATURES = [  # column, label, decimals
    ("peak167_prominence", "167 bp peak prominence", 2),
    ("nucleosome_ladder_score", "Nucleosome ladder (di, tri)", 2),
    ("motif_CCCA", "CCCA end frequency", 4),
    ("motif_diversity_score", "Motif diversity", 3),
    ("median_len", "Median length (bp)", 0),
    ("frac_lt150", "Fraction < 150 bp", 3),
    ("frac_400_600", "Fraction 400-600 bp", 3),
]
MALE_REFERENCE = ["CF0003"]  # male recipient, male donor: chrY level of an all-male sample
FEMALE_SAME_SEX = ["CF0007", "CF0009", "CF0020"]  # female recipient, female donor: background
SEX_MISMATCHED = {"CF0012": "male donor, female recipient"}


def rng(s, d):
    return f"{s.min():.{d}f}-{s.max():.{d}f}"


def main():
    d = sys.argv[1] if len(sys.argv) > 1 else "cohort"
    s = pd.read_csv(f"{d}/size_motif_summary.tsv", sep="\t", index_col=0)
    bal = [x for x in s.index if x.endswith("_BAL") and x[:-4] + "_P" in s.index]
    pla = [x[:-4] + "_P" for x in bal]
    recip = pd.Series([x[:6] for x in bal], index=bal)
    print(f"== 1. Paired BAL vs plasma: {len(bal)} pairs, {recip.nunique()} recipients "
          f"({', '.join(sorted(recip.unique()))})\n")
    rows = []
    for col, label, dp in FEATURES:
        if col not in s.columns:
            continue
        b, p = s.loc[bal, col].values, s.loc[pla, col].values
        diff = pd.Series(b - p, index=bal)
        stat = wilcoxon(b, p, alternative="two-sided", method="exact") if len(bal) >= 2 else None
        up = (diff > 0).sum()
        per_recip = diff.groupby(recip).apply(lambda x: "+" if (x > 0).all() else ("-" if (x < 0).all() else "mixed"))
        rows.append({"feature": label, "BAL": rng(s.loc[bal, col], dp), "plasma": rng(s.loc[pla, col], dp),
                     "BAL>plasma": f"{up}/{len(bal)}", "P (Wilcoxon)": f"{stat.pvalue:.3f}" if stat else "",
                     "per recipient": " ".join(f"{k}{v}" for k, v in per_recip.items())})
    print(pd.DataFrame(rows).to_string(index=False))

    print("\n== 2. Plasma: healthy vs transplant (descriptive)\n")
    grp = s["group"] if "group" in s.columns else pd.Series("", index=s.index)
    hp, tp = s.index[grp == "healthy_plasma"], s.index[grp == "ltx_plasma"]
    rows = [{"feature": label, f"healthy (n={len(hp)})": rng(s.loc[hp, col], dp),
             f"transplant (n={len(tp)})": rng(s.loc[tp, col], dp)}
            for col, label, dp in FEATURES if col in s.columns and len(hp) and len(tp)]
    print(pd.DataFrame(rows).to_string(index=False))

    if "frac_chrY" in s.columns:
        print("\n== 3. chrY fraction (fragments 50-600 bp)\n")
        y = s["frac_chrY"]
        print(y.to_string(float_format="%.6f"))
        male = y[[i for i in y.index if i[:6] in MALE_REFERENCE and i.endswith("_P")]]
        bg = y[[i for i in y.index if i[:6] in FEMALE_SAME_SEX]]
        if len(male) and len(bg):
            print(f"\nmale reference (CF0003 plasma) median {male.median():.6f}; "
                  f"female background (same-sex pairs) median {bg.median():.6f}, max {bg.max():.6f}")
            for r, note in SEX_MISMATCHED.items():
                for i in [i for i in y.index if i.startswith(r)]:
                    est = (y[i] - bg.median()) / (male.median() - bg.median())
                    flag = "above background max" if y[i] > bg.max() else "WITHIN background range"
                    print(f"{i} ({note}): chrY {y[i]:.6f} -> donor DNA about {100 * est:.1f}% ({flag})")

    sc_path = f"{d}/size_class_motifs.tsv"
    if os.path.exists(sc_path):
        print("\n== 4. Long-fragment end-motif excess divergence (vs 150-220 bp)\n")
        sc = pd.read_csv(sc_path, sep="\t", index_col=0)
        sc = sc[sc.size_class.isin(["401_600", "601_1000"]) & sc.js_excess.notna()]
        print(sc.groupby(["group", "size_class"]).js_excess.agg(["count", "min", "median", "max"])
              .to_string(float_format="%.4f"))
        if "js_p_vs_null" in sc.columns:
            print(f"\nlargest js_p_vs_null: {sc.js_p_vs_null.max():.4f} (floor = 1/(n_boot+1))")
    for f in ("motif_scale_pairs.tsv", "motif_scale_long.tsv"):
        if not os.path.exists(f"{d}/{f}"):
            print(f"\n(no {f}: run  python scripts/motif_scale.py {d}  for the divergence scale)")


if __name__ == "__main__":
    main()
