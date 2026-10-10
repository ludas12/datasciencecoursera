#!/usr/bin/env python3
"""Summarise a UXM deconvolution CSV (rows = cell types, columns = samples).

usage: summarize_uxm.py uxm_deconv.csv OUTDIR
Writes OUTDIR/uxm_grouped.tsv and OUTDIR/uxm_grouped.png and prints
percentages per sample: grouped compartments plus the top cell types.
"""
import sys
from pathlib import Path

import pandas as pd
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

GROUPS = {
    "Lung epithelium": ["Lung-Ep-Alveo", "Lung-Ep-Bron"],
    "Neutrophils/granulocytes": ["Blood-Granul"],
    "Monocytes/macrophages": ["Blood-Mono+Macro"],
    "Lymphocytes (T/B/NK)": ["Blood-T", "Blood-B", "Blood-NK"],
    "Megakaryocytes/erythroid": ["Megakaryocytes", "Eryth-prog"],
    "Endothelium": ["Endothel"],
    "Liver": ["Liver-Hep"],
    "Stromal (fibro/muscle)": ["Colon-Fibro", "Dermal-Fibro", "Heart-Fibro", "Smooth-Musc",
                               "Skeletal-Musc", "Heart-Cardio"],
}


def main(csv, outdir):
    df = pd.read_csv(csv).set_index("CellType") * 100
    grouped = {g: df.loc[[c for c in cells if c in df.index]].sum() for g, cells in GROUPS.items()}
    grouped = pd.DataFrame(grouped).T
    grouped.loc["Other"] = df.sum() - grouped.sum()
    grouped.round(2).to_csv(Path(outdir) / "uxm_grouped.tsv", sep="\t")

    print("Estimated contribution (%) by compartment:")
    print(grouped.round(1).to_string())
    print("\nTop cell types per sample (%):")
    for s in df.columns:
        top = df[s].sort_values(ascending=False).head(6)
        print(f"  {s}: " + ", ".join(f"{c} {v:.1f}" for c, v in top.items()))

    fig, ax = plt.subplots(figsize=(1.3 * len(df.columns) + 3, 5))
    bottom = pd.Series(0.0, index=grouped.columns)
    for g in grouped.index:
        ax.bar(grouped.columns, grouped.loc[g], bottom=bottom, label=g)
        bottom += grouped.loc[g]
    ax.set_ylabel("% of DNA")
    ax.set_title("UXM tissue/cell-of-origin (Loyfer et al. 2023 atlas)")
    ax.tick_params(axis="x", rotation=45)
    ax.legend(fontsize=8, bbox_to_anchor=(1.01, 1), loc="upper left")
    fig.tight_layout()
    fig.savefig(Path(outdir) / "uxm_grouped.png", dpi=150)


if __name__ == "__main__":
    if len(sys.argv) != 3:
        sys.exit(__doc__)
    main(sys.argv[1], sys.argv[2])
