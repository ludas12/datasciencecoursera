#!/usr/bin/env python3
"""Region-level methylation analysis for low-coverage EM-seq.

For every region set (bins, CpG islands, promoters) found in
REGIONS_DIR/mapped/<set>/<label>.tsv (bedtools map output: chr start end name
n_meth n_unmeth) this writes to OUT/:

  <set>.methylation.tsv.gz   per-region methylation fraction + calls per sample
  <set>.sample_summary.tsv   regions used, mean/median methylation per sample
  <set>.correlation.tsv/.png Pearson correlation between samples
  <set>.pca.png              PCA of the most variable regions
  <set>.distribution.png     methylation distribution per sample
  <set>.<comparison>.tsv.gz  per-region difference (group_b - group_a) with a
                             two-proportion z-test on read counts + BH FDR
  <set>.<comparison>.png     scatter of group_a vs group_b

With one sample per group the p-values only reflect read-sampling noise, not
biological variability: use them to rank regions, not as formal evidence.
"""
import argparse
import gzip
import math
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

N_VARIABLE = 5000          # regions used for PCA
MIN_DELTA = 0.2            # |difference| reported as "large"
FDR = 0.05


def load_set(set_dir, labels):
    meth, calls = {}, {}
    keys = None
    for lab in labels:
        p = set_dir / f"{lab}.tsv"
        if not p.exists():
            continue
        df = pd.read_csv(p, sep="\t", header=None, na_values=".",
                         names=["chr", "start", "end", "name", "m", "u"])
        if keys is None:
            keys = df[["chr", "start", "end", "name"]]
        m, u = df["m"].fillna(0).to_numpy(), df["u"].fillna(0).to_numpy()
        n = m + u
        with np.errstate(invalid="ignore", divide="ignore"):
            meth[lab] = np.where(n > 0, m / n, np.nan)
        calls[lab] = n
    if keys is None:
        return None, None, None
    return keys, pd.DataFrame(meth), pd.DataFrame(calls)


def bh(p):
    p = np.asarray(p, dtype=float)
    order = np.argsort(p)
    ranked = p[order] * len(p) / np.arange(1, len(p) + 1)
    q = np.minimum.accumulate(ranked[::-1])[::-1]
    out = np.empty_like(q)
    out[order] = np.minimum(q, 1)
    return out


def two_prop_z(m1, n1, m2, n2):
    p1, p2 = m1 / n1, m2 / n2
    p = (m1 + m2) / (n1 + n2)
    se = np.sqrt(p * (1 - p) * (1 / n1 + 1 / n2))
    with np.errstate(invalid="ignore", divide="ignore"):
        z = np.where(se > 0, (p2 - p1) / se, 0.0)
    pval = np.array([math.erfc(abs(x) / math.sqrt(2)) for x in z])
    return z, pval


def heatmap(mat, labels, title, path):
    fig, ax = plt.subplots(figsize=(1.1 * len(labels) + 2, 1.1 * len(labels) + 1))
    im = ax.imshow(mat, cmap="viridis", vmin=np.nanmin(mat), vmax=1)
    ax.set_xticks(range(len(labels)), labels, rotation=45, ha="right")
    ax.set_yticks(range(len(labels)), labels)
    for i in range(len(labels)):
        for j in range(len(labels)):
            ax.text(j, i, f"{mat[i, j]:.2f}", ha="center", va="center",
                    color="white" if mat[i, j] < (np.nanmin(mat) + 1) / 2 else "black", fontsize=8)
    fig.colorbar(im, ax=ax, label="Pearson r")
    ax.set_title(title)
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)


def analyse_set(name, keys, meth, calls, info, comparisons, min_calls, out):
    labels = list(meth.columns)
    ok = (calls >= min_calls).all(axis=1).to_numpy()
    print(f"[{name}] {len(keys):,} regions, {ok.sum():,} with >= {min_calls} calls in all "
          f"{len(labels)} samples")

    table = keys.copy()
    for lab in labels:
        table[f"{lab}_meth"] = meth[lab].round(4)
        table[f"{lab}_calls"] = calls[lab].astype(int)
    with gzip.open(out / f"{name}.methylation.tsv.gz", "wt") as fh:
        table.to_csv(fh, sep="\t", index=False)

    m_ok = meth[ok]
    summary = pd.DataFrame({
        "sample": labels,
        "regions_used": [int(ok.sum())] * len(labels),
        "mean_meth": [m_ok[l].mean() for l in labels],
        "median_meth": [m_ok[l].median() for l in labels],
        "median_calls_per_region": [calls[l][ok].median() for l in labels],
    })
    summary.to_csv(out / f"{name}.sample_summary.tsv", sep="\t", index=False, float_format="%.4f")
    if ok.sum() < 10:
        print(f"[{name}] too few regions pass the coverage filter; skipping plots")
        return

    corr = m_ok.corr()
    corr.to_csv(out / f"{name}.correlation.tsv", sep="\t", float_format="%.4f")
    heatmap(corr.to_numpy(), labels, f"{name}: methylation correlation", out / f"{name}.correlation.png")

    # PCA on the most variable regions
    var = m_ok.var(axis=1).sort_values(ascending=False)
    x = m_ok.loc[var.index[:N_VARIABLE]].to_numpy().T
    x = x - x.mean(axis=0)
    u, s, _ = np.linalg.svd(x, full_matrices=False)
    pcs = u * s
    expl = s ** 2 / np.sum(s ** 2)
    fig, ax = plt.subplots(figsize=(6, 5))
    types = info.set_index("label")["type"].reindex(labels).fillna("NA")
    for t in sorted(types.unique()):
        idx = [i for i, l in enumerate(labels) if types[l] == t]
        ax.scatter(pcs[idx, 0], pcs[idx, 1], s=60, label=t)
    for i, l in enumerate(labels):
        ax.annotate(l, (pcs[i, 0], pcs[i, 1]), fontsize=8, xytext=(4, 4), textcoords="offset points")
    ax.set_xlabel(f"PC1 ({100 * expl[0]:.1f}%)")
    ax.set_ylabel(f"PC2 ({100 * expl[1]:.1f}%)" if len(expl) > 1 else "PC2")
    ax.set_title(f"{name}: PCA of {x.shape[1]:,} most variable regions")
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(out / f"{name}.pca.png", dpi=150)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(1.2 * len(labels) + 2, 4))
    ax.violinplot([m_ok[l].dropna() for l in labels], showmedians=True)
    ax.set_xticks(range(1, len(labels) + 1), labels, rotation=45, ha="right")
    ax.set_ylabel("methylation fraction")
    ax.set_title(f"{name}: per-region methylation")
    fig.tight_layout()
    fig.savefig(out / f"{name}.distribution.png", dpi=150)
    plt.close(fig)

    rows = []
    for _, c in comparisons.iterrows():
        a, b = c["group_a"], c["group_b"]
        if a not in labels or b not in labels:
            print(f"[{name}] comparison {c['name']}: {a} or {b} missing; skipped")
            continue
        keep = (calls[a] >= min_calls) & (calls[b] >= min_calls)
        na, nb = calls[a][keep].to_numpy(), calls[b][keep].to_numpy()
        ma, mb = meth[a][keep].to_numpy() * na, meth[b][keep].to_numpy() * nb
        z, p = two_prop_z(ma, na, mb, nb)
        q = bh(p)
        res = keys[keep.to_numpy()].copy()
        res[f"{a}_meth"], res[f"{b}_meth"] = (ma / na).round(4), (mb / nb).round(4)
        res[f"{a}_calls"], res[f"{b}_calls"] = na.astype(int), nb.astype(int)
        res["delta"] = (mb / nb - ma / na).round(4)
        res["z"], res["p"], res["fdr"] = z.round(3), p, q
        res = res.sort_values("p")
        with gzip.open(out / f"{name}.{c['name']}.tsv.gz", "wt") as fh:
            res.to_csv(fh, sep="\t", index=False)
        big = (res["fdr"] < FDR) & (res["delta"].abs() >= MIN_DELTA)
        rows.append({"comparison": c["name"], "a": a, "b": b, "regions_tested": len(res),
                     "pearson_r": np.corrcoef(ma / na, mb / nb)[0, 1],
                     "mean_delta": res["delta"].mean(),
                     f"hyper_in_b(fdr<{FDR},delta>={MIN_DELTA})": int((big & (res["delta"] > 0)).sum()),
                     f"hypo_in_b(fdr<{FDR},delta<=-{MIN_DELTA})": int((big & (res["delta"] < 0)).sum())})

        fig, ax = plt.subplots(figsize=(5, 5))
        ax.scatter(ma / na, mb / nb, s=2, alpha=0.2, color="grey", rasterized=True)
        sig = res[big]
        ax.scatter(sig[f"{a}_meth"], sig[f"{b}_meth"], s=4, color="crimson", label=f"FDR<{FDR}, |Δ|≥{MIN_DELTA}")
        ax.plot([0, 1], [0, 1], "k--", lw=0.8)
        ax.set_xlabel(f"{a} methylation")
        ax.set_ylabel(f"{b} methylation")
        ax.set_title(f"{name}: {c['name']} ({len(res):,} regions)")
        ax.legend(fontsize=8, loc="upper left")
        fig.tight_layout()
        fig.savefig(out / f"{name}.{c['name']}.png", dpi=150)
        plt.close(fig)
    if rows:
        pd.DataFrame(rows).to_csv(out / f"{name}.comparisons_summary.tsv", sep="\t",
                                  index=False, float_format="%.4f")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--regions-dir", required=True, type=Path)
    ap.add_argument("--sample-info", required=True, type=Path)
    ap.add_argument("--comparisons", type=Path)
    ap.add_argument("--min-calls", type=int, default=10)
    ap.add_argument("--out", required=True, type=Path)
    args = ap.parse_args()

    info = pd.read_csv(args.sample_info, sep="\t", dtype=str)
    info = info[info["include"].str.lower() == "yes"]
    labels = list(info["label"])
    comparisons = (pd.read_csv(args.comparisons, sep="\t", dtype=str)
                   if args.comparisons and args.comparisons.exists()
                   else pd.DataFrame(columns=["name", "group_a", "group_b"]))
    args.out.mkdir(parents=True, exist_ok=True)

    for set_dir in sorted((args.regions_dir / "mapped").iterdir()):
        keys, meth, calls = load_set(set_dir, labels)
        if keys is None:
            continue
        analyse_set(set_dir.name, keys, meth, calls, info, comparisons, args.min_calls, args.out)


if __name__ == "__main__":
    main()
