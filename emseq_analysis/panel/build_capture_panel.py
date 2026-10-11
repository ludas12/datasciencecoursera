#!/usr/bin/env python3
"""Draft a targeted methylation capture panel from the Loyfer et al. 2023 atlas.

Selects cell-type-specific unmethylated markers (UXM "U250" set), weighted
towards the question in lung transplant recipients: lung epithelium
(alveolar, bronchial) against the cell types that make up plasma and BAL
cfDNA (blood/immune cells, endothelium, liver), plus a smaller set for every
other atlas cell type so deconvolution does not force unrelated signal into
the main compartments.

Lung markers are ranked by their contrast against the *highest* unmethylated
fraction among plasma/BAL background cell types: a ~1 % lung signal is only
detectable at markers that are almost never unmethylated in blood. Erythroid
progenitors are globally hypomethylated and would dominate that maximum, so
they are left out of the ranking and reported separately (eryth_U) for review.

usage: build_capture_panel.py MARKERS.U250.hg38.tsv ATLAS.U250.l4.hg38.full.tsv OUTDIR
writes OUTDIR/panel_markers.tsv, panel_targets.bed, panel_summary.tsv
"""
import sys
from pathlib import Path

import pandas as pd

PAD = 60                         # bp added each side so overlapping fragments are captured
LUNG = ["Lung-Ep-Alveo", "Lung-Ep-Bron"]
BACKGROUND = ["Blood-Granul", "Blood-Mono+Macro", "Blood-T", "Blood-B", "Blood-NK",
              "Megakaryocytes", "Eryth-prog", "Endothel", "Liver-Hep"]
STROMAL = ["Smooth-Musc", "Dermal-Fibro", "Heart-Fibro", "Colon-Fibro", "Skeletal-Musc",
           "Heart-Cardio", "Adipocytes"]
N_PER_TYPE = {"lung": 250, "background": 100, "stromal": 40, "other": 25}
MIN_CPGS = 4                     # UXM counts fragments with >= 4 CpGs in the marker


def tier(cell):
    if cell in LUNG:
        return "lung"
    if cell in BACKGROUND:
        return "background"
    if cell in STROMAL:
        return "stromal"
    return "other"


def main(markers_path, atlas_path, outdir):
    out = Path(outdir)
    out.mkdir(parents=True, exist_ok=True)
    m = pd.read_csv(markers_path, sep="\t")
    m = m.rename(columns={"#chr": "chr"})
    m["n_cpgs"] = m["lenCpG"].str.replace("CpGs", "").astype(int)
    m["length_bp"] = m["bp"].str.replace("bp", "").astype(int)
    m = m[(m["direction"] == "U") & (m["n_cpgs"] >= MIN_CPGS)]

    # Per-cell-type unmethylated-fragment fractions at each marker.
    atlas = pd.read_csv(atlas_path, sep="\t")
    bg_cols = [c for c in BACKGROUND if c in atlas.columns and c != "Eryth-prog"]
    # highest unmethylated fraction among background cell types other than the target
    atlas["max_bg_U"] = [max(row[c] for c in bg_cols if c != row["target"]) for _, row in atlas.iterrows()]
    atlas["eryth_U"] = atlas["Eryth-prog"]
    atlas["target_U"] = [row[row["target"]] if row["target"] in atlas.columns else float("nan")
                         for _, row in atlas.iterrows()]
    m = m.merge(atlas[["name", "target_U", "max_bg_U", "eryth_U"]], left_on="region", right_on="name", how="left")
    m["contrast_vs_bg"] = m["target_U"] - m["max_bg_U"]
    m["tier"] = m["target"].map(tier)

    picks = []
    for cell, g in m.groupby("target"):
        t = tier(cell)
        if t == "lung":
            # low blood/endothelium/liver background first, then effect size
            g = g.sort_values(["contrast_vs_bg", "delta_means"], ascending=False)
        else:
            g = g.sort_values(["delta_means", "n_cpgs"], ascending=False)
        picks.append(g.head(N_PER_TYPE[t]))
    panel = pd.concat(picks)
    order = {"lung": 0, "background": 1, "stromal": 2, "other": 3}
    panel = panel.sort_values(["tier", "target", "chr", "start"], key=lambda s: s.map(order) if s.name == "tier" else s)

    cols = ["chr", "start", "end", "region", "target", "tier", "n_cpgs", "length_bp",
            "tg_mean", "bg_mean", "delta_means", "target_U", "max_bg_U", "eryth_U", "contrast_vs_bg", "ttest"]
    panel[cols].to_csv(out / "panel_markers.tsv", sep="\t", index=False, float_format="%.4g")

    # Padded, merged capture targets.
    bed = panel[["chr", "start", "end", "target"]].copy()
    bed["start"] = (bed["start"] - PAD).clip(lower=0)
    bed["end"] = bed["end"] + PAD
    bed = bed.sort_values(["chr", "start"])
    merged = []
    for chrom, g in bed.groupby("chr", sort=True):
        cur = None
        for s, e, t in g[["start", "end", "target"]].itertuples(index=False):
            if cur and s <= cur[2]:
                cur[2] = max(cur[2], e)
                cur[3].add(t)
            else:
                if cur:
                    merged.append(cur)
                cur = [chrom, s, e, {t}]
        if cur:
            merged.append(cur)
    with open(out / "panel_targets.bed", "w") as fh:
        for c, s, e, t in merged:
            fh.write(f"{c}\t{s}\t{e}\t{','.join(sorted(t))}\n")

    summ = (panel.groupby(["tier", "target"])
            .agg(markers=("region", "size"), total_bp=("length_bp", "sum"),
                 median_cpgs=("n_cpgs", "median"), median_delta=("delta_means", "median"),
                 median_max_bg_U=("max_bg_U", "median"),
                 markers_bg_le_2pct=("max_bg_U", lambda x: int((x <= 0.02).sum())),
                 markers_bg_le_5pct=("max_bg_U", lambda x: int((x <= 0.05).sum())))
            .reset_index())
    summ["tier"] = pd.Categorical(summ["tier"], ["lung", "background", "stromal", "other"])
    summ = summ.sort_values(["tier", "target"])
    summ.to_csv(out / "panel_summary.tsv", sep="\t", index=False, float_format="%.3f")

    target_bp = sum(e - s for _, s, e, _ in merged)
    print(f"{len(panel)} markers from {panel['target'].nunique()} cell types; "
          f"{len(merged)} capture targets, {target_bp / 1e6:.2f} Mb with ±{PAD} bp padding")
    print(summ.to_string(index=False))


if __name__ == "__main__":
    if len(sys.argv) != 4:
        sys.exit(__doc__)
    main(*sys.argv[1:])
