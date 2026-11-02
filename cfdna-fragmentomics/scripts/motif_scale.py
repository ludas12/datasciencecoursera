#!/usr/bin/env python3
"""Put the long-fragment end-motif divergences on a biological scale.

js_p_vs_null bottoms out at 1/(n_boot+1) for every sample once classes hold tens of thousands of
fragments, so "different from noise" says nothing about how different. This compares each
sample's js_excess for its long classes against end-motif divergences between whole samples:

  BAL vs plasma, same visit       (two nuclease environments: the scale of a real process change)
  plasma vs plasma, between people (biological spread within one compartment)
  BAL vs BAL, between samples

Usage (from the cohort output folder):
    python scripts/motif_scale.py cohort/            # reads end_motif_freq.tsv, size_class_motifs.tsv
"""
import itertools
import os
import sys

import pandas as pd

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
from fragmentomics.cfdna_frag import js_divergence  # noqa: E402


def main():
    d = sys.argv[1] if len(sys.argv) > 1 else "cohort"
    motifs = pd.read_csv(f"{d}/end_motif_freq.tsv", sep="\t", index_col=0).fillna(0)
    sc = pd.read_csv(f"{d}/size_class_motifs.tsv", sep="\t", index_col=0)
    group = sc.groupby(level=0)["group"].first()
    motifs = motifs.loc[[s for s in motifs.index if s in group.index]]

    def js(a, b):
        return js_divergence(motifs.loc[a].values, motifs.loc[b].values)

    rows = []
    bal = [s for s in motifs.index if group[s] == "ltx_bal"]
    plasma = [s for s in motifs.index if group[s].endswith("plasma")]
    for b in bal:
        p = b.replace("_BAL", "_P")
        if p in motifs.index:
            rows.append(("BAL vs plasma, same visit", f"{b} / {p}", js(b, p)))
    for a, b in itertools.combinations(plasma, 2):
        rows.append(("plasma vs plasma", f"{a} / {b}", js(a, b)))
    for a, b in itertools.combinations(bal, 2):
        rows.append(("BAL vs BAL", f"{a} / {b}", js(a, b)))
    ref = pd.DataFrame(rows, columns=["comparison", "pair", "js"])
    print("Whole-sample end-motif divergences (all fragments):")
    print(ref.groupby("comparison")["js"].describe()[["count", "min", "50%", "max"]]
          .to_string(float_format="%.4f"))

    scale = ref.loc[ref.comparison == "BAL vs plasma, same visit", "js"].median()
    long_ = sc[sc.size_class.isin(["401_600", "601_1000"])][["group", "size_class", "n", "js_excess"]]
    long_ = long_.assign(frac_of_BAL_vs_plasma=long_.js_excess / scale)
    print(f"\nLong-class js_excess as a fraction of the median BAL-vs-plasma divergence ({scale:.4f}):")
    print(long_.to_string(float_format="%.4f"))
    ref.to_csv(f"{d}/motif_scale_pairs.tsv", sep="\t", index=False, float_format="%.5g")
    long_.to_csv(f"{d}/motif_scale_long.tsv", sep="\t", float_format="%.5g")


if __name__ == "__main__":
    main()
