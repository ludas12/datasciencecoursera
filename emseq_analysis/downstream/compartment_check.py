#!/usr/bin/env python3
"""BAL-versus-plasma checks of the 5mC-protection question (step 11).

Pairs are one plasma and one BAL library from the same subject and timepoint
(downstream/samples.tsv, include = yes). Writes to OUT/:

Imprinted promoters (from IMPRINTED_CPGS: label chr pos gene n_meth n_unmeth)
  imprinted_genes.tsv          per gene and library: methylation, calls, CpGs
                               covered, share of calls on the top 3 CpGs
  imprinted_pairs.tsv          per pair and gene: BAL - plasma difference with a
                               95 % CI from resampling CpGs, and whether leaving
                               out any single CpG flips its sign
  imprinted_pairs_summary.tsv  per pair: genes lower/higher in BAL, two-sided
                               sign test, genes whose CI excludes 0
  imprinted_pairs.png

BAL against plasma by methylation level (from step 05 region counts)
  Regions are split into strata by the methylation of the reference libraries
  (independent of both samples, so stratifying does not bias the comparison).
  Within each stratum the BAL and plasma calls are pooled. If BAL lost 5mC
  uniformly, BAL = k x plasma in every stratum: intercept ~0 and the same
  BAL/plasma ratio in highly and moderately methylated regions. A difference
  in cell composition gives a non-zero intercept or a ratio that changes
  with methylation level.
  bal_vs_plasma_strata.tsv     per pair, region set and stratum
  bal_vs_plasma_fit.tsv        intercept, slope, ratios, bootstrap 95 % CIs
  bal_vs_plasma.png

Fragment size (from --size-dir: bins/<label>.<class>.tsv, classes from --size-classes)
  size_split_fit.tsv, size_split_strata.tsv: the same comparison for plasma
  short / plasma long, BAL long / plasma long, BAL short / plasma long and
  BAL short / BAL long. BAL long ~ plasma long means the BAL deficit is
  carried by its short fragments; a deficit in long fragments too points to
  cell composition.
  size_classes.tsv: each size class of plasma and BAL against plasma long
  fragments, with its share of CpG calls.

With one library per compartment these describe this patient's samples; they
are not tests of BAL versus plasma in general.
"""
import argparse
import math
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

MIN_GENE_CALLS = 20        # calls a gene needs in both libraries (as in step 04)
N_BOOT = 1000              # CpG resamples per gene
N_BOOT_BINS = 200          # region resamples per fit
N_STRATA = 20              # equal-width strata of reference methylation
MIN_STRATUM_REGIONS = 50
HIGH = (0.7, 1.0)          # reference methylation of "high" regions
MID = (0.2, 0.7)
SEED = 1


def read_info(path):
    info = pd.read_csv(path, sep="\t", dtype=str).fillna("NA")
    info["type"] = info["type"].str.lower()
    return info[info["include"].str.lower() == "yes"]


def find_pairs(info):
    pairs = []
    pat = info[info["type"].isin(["plasma", "bal"])]
    for (subj, tp), g in pat.groupby(["subject", "timepoint"], sort=True):
        pl = g.loc[g["type"] == "plasma", "label"].tolist()
        bal = g.loc[g["type"] == "bal", "label"].tolist()
        if len(pl) == 1 and len(bal) == 1:
            pairs.append((f"{subj}_{tp}", pl[0], bal[0]))
    return pairs


def sign_test(k, n):
    """Two-sided exact sign test: k of n differences in one direction."""
    if n == 0:
        return float("nan")
    lo = min(k, n - k)
    tail = sum(math.comb(n, i) for i in range(lo + 1)) / 2 ** n
    return min(1.0, 2 * tail)


def fmt(v, d=3):
    return "NA" if v is None or (isinstance(v, float) and math.isnan(v)) else f"{v:.{d}f}"


# ---------------------------------------------------------------- imprinted --
def gene_table(cpg):
    rows = []
    for (lab, gene), g in cpg.groupby(["label", "gene"], sort=True):
        n = (g["m"] + g["u"]).to_numpy()
        tot = n.sum()
        top3 = np.sort(n)[::-1][:3].sum() / tot if tot else float("nan")
        rows.append({"label": lab, "gene": gene, "pct_meth": 100 * g["m"].sum() / tot if tot else float("nan"),
                     "calls": int(tot), "cpgs_covered": int((n > 0).sum()), "top3_cpg_share": top3})
    return pd.DataFrame(rows)


def pair_gene(gp, gb, rng):
    """BAL - plasma methylation for one gene, resampling and leaving out CpGs."""
    j = pd.merge(gp, gb, on=["chr", "pos"], how="outer", suffixes=("_p", "_b")).fillna(0)
    mp, np_ = j["m_p"].to_numpy(), (j["m_p"] + j["u_p"]).to_numpy()
    mb, nb = j["m_b"].to_numpy(), (j["m_b"] + j["u_b"]).to_numpy()
    if np_.sum() < MIN_GENE_CALLS or nb.sum() < MIN_GENE_CALLS:
        return None
    d = mb.sum() / nb.sum() - mp.sum() / np_.sum()
    k = len(j)
    idx = rng.integers(0, k, size=(N_BOOT, k))
    sp, sb = np_[idx].sum(1), nb[idx].sum(1)
    ok = (sp > 0) & (sb > 0)
    boot = mb[idx].sum(1)[ok] / sb[ok] - mp[idx].sum(1)[ok] / sp[ok]
    lo, hi = np.percentile(boot, [2.5, 97.5]) if ok.sum() > 10 else (float("nan"), float("nan"))
    # leave one CpG out
    lp, lb = np_.sum() - np_, nb.sum() - nb
    okl = (lp > 0) & (lb > 0)
    loo = (mb.sum() - mb)[okl] / lb[okl] - (mp.sum() - mp)[okl] / lp[okl]
    stable = bool(d != 0 and np.all(np.sign(loo) == np.sign(d)))
    return {"pl_pct": 100 * mp.sum() / np_.sum(), "bal_pct": 100 * mb.sum() / nb.sum(),
            "diff_pct": 100 * d, "ci_lo": 100 * lo, "ci_hi": 100 * hi,
            "pl_calls": int(np_.sum()), "bal_calls": int(nb.sum()),
            "cpgs_union": k, "max_loo_shift_pct": 100 * float(np.max(np.abs(loo - d))) if okl.any() else float("nan"),
            "sign_stable_loo": stable}


def imprinted(cpg_path, pairs, info, out):
    cpg = pd.read_csv(cpg_path, sep="\t", header=None, names=["label", "chr", "pos", "gene", "m", "u"])
    cpg = cpg[cpg["label"].isin(info["label"])]
    genes = gene_table(cpg)
    genes.to_csv(out / "imprinted_genes.tsv", sep="\t", index=False, float_format="%.3f")
    rng = np.random.default_rng(SEED)
    rows, summ = [], []
    for name, pl, bal in pairs:
        cp, cb = cpg[cpg["label"] == pl], cpg[cpg["label"] == bal]
        for gene in sorted(set(cp["gene"]) | set(cb["gene"])):
            r = pair_gene(cp[cp["gene"] == gene][["chr", "pos", "m", "u"]],
                          cb[cb["gene"] == gene][["chr", "pos", "m", "u"]], rng)
            if r:
                rows.append({"pair": name, "plasma": pl, "bal": bal, "gene": gene, **r})
    pg = pd.DataFrame(rows)
    pg.to_csv(out / "imprinted_pairs.tsv", sep="\t", index=False, float_format="%.3f")
    if pg.empty:
        print("imprinted promoters: no gene had enough calls in both libraries of any pair")
        return
    print("Imprinted promoters, BAL - plasma (each timepoint separately; genes with "
          f">= {MIN_GENE_CALLS} calls in both)")
    print("pair\tgenes\tlower_in_BAL\thigher_in_BAL\tsign_p_2sided\tmedian_diff_pct\t"
          "CI_below_0\tCI_above_0\tlower_and_LOO_stable")
    for name, g in pg.groupby("pair", sort=True):
        nz = g[g["diff_pct"] != 0]
        lower = int((nz["diff_pct"] < 0).sum())
        s = {"pair": name, "genes": len(g), "lower_in_BAL": lower, "higher_in_BAL": len(nz) - lower,
             "sign_p_2sided": sign_test(lower, len(nz)), "median_diff_pct": float(g["diff_pct"].median()),
             "ci_below_0": int((g["ci_hi"] < 0).sum()), "ci_above_0": int((g["ci_lo"] > 0).sum()),
             "lower_and_loo_stable": int(((g["diff_pct"] < 0) & g["sign_stable_loo"]).sum())}
        summ.append(s)
        print("\t".join([name] + [f"{v:.2g}" if k == "sign_p_2sided" else fmt(v, 1) if isinstance(v, float)
                                   else str(v) for k, v in s.items() if k != "pair"]))
    pd.DataFrame(summ).to_csv(out / "imprinted_pairs_summary.tsv", sep="\t", index=False, float_format="%.4f")

    names = sorted(pg["pair"].unique())
    fig, axes = plt.subplots(1, len(names), figsize=(4.2 * len(names), 0.32 * pg["gene"].nunique() + 1.5),
                             sharey=True, squeeze=False)
    order = sorted(pg["gene"].unique(), reverse=True)
    for ax, name in zip(axes[0], names):
        g = pg[pg["pair"] == name].set_index("gene")
        for i, gene in enumerate(order):
            if gene not in g.index:
                continue
            r = g.loc[gene]
            col = "#c0392b" if r["ci_hi"] < 0 else "#2471a3" if r["ci_lo"] > 0 else "#7f8c8d"
            ax.plot([r["ci_lo"], r["ci_hi"]], [i, i], color=col, lw=1.5)
            ax.plot(r["diff_pct"], i, "o", color=col, mfc=col if r["sign_stable_loo"] else "white")
        ax.axvline(0, color="black", lw=0.8)
        ax.set_yticks(range(len(order)), order)
        ax.set_xlabel("BAL - plasma methylation (%)")
        ax.set_title(name)
    fig.suptitle("Imprinted promoters: 95 % CI from resampling CpGs; open = sign flips when one CpG is left out",
                 fontsize=9)
    fig.tight_layout()
    fig.savefig(out / "imprinted_pairs.png", dpi=150)
    plt.close(fig)


# -------------------------------------------------------- stratified scatter --
def load_counts(path):
    df = pd.read_csv(path, sep="\t", header=None, na_values=".",
                     names=["chr", "start", "end", "name", "m", "u"])
    m, u = df["m"].fillna(0).to_numpy(), df["u"].fillna(0).to_numpy()
    return df["name"].to_numpy(), m, m + u


def fit_strata(s, mp, np_, mb, nb, keep_strata):
    """Pooled plasma and BAL methylation per stratum, line fit and ratios."""
    xs = np.bincount(s, mp, N_STRATA) / np.maximum(np.bincount(s, np_, N_STRATA), 1)
    ys = np.bincount(s, mb, N_STRATA) / np.maximum(np.bincount(s, nb, N_STRATA), 1)
    if keep_strata.sum() < 3:
        return xs, ys, float("nan"), float("nan")
    b, a = np.polyfit(xs[keep_strata], ys[keep_strata], 1)
    return xs, ys, a, b


def ratio(sel, mp, np_, mb, nb):
    if not sel.any() or np_[sel].sum() == 0 or nb[sel].sum() == 0 or mp[sel].sum() == 0:
        return float("nan")
    return (mb[sel].sum() / nb[sel].sum()) / (mp[sel].sum() / np_[sel].sum())


def region_flags(d, flag_dir, rset, refs, min_calls):
    """Reference counts and usable-region mask (not blacklisted, enough reference calls)."""
    ref = [load_counts(d / f"{r}.tsv") for r in refs if (d / f"{r}.tsv").exists()]
    if not ref:
        return None
    names = ref[0][0]
    rm, rn = sum(r[1] for r in ref), sum(r[2] for r in ref)
    black = set()
    bl = flag_dir / f"{rset}.blacklisted.bed"
    if bl.exists() and bl.stat().st_size:
        black = set(pd.read_csv(bl, sep="\t", header=None, usecols=[3])[3])
    return rm, rn, ~np.isin(names, list(black)) & (rn >= min_calls)


def compare(rm, rn, base_ok, x, y, min_calls, max_fold, rng):
    """y against x pooled within strata of reference methylation, with bootstrap CIs."""
    (_, mp, np_), (_, mb, nb) = x, y
    ok = base_ok & (np_ >= min_calls) & (nb >= min_calls)
    for n in (np_, nb):
        med = np.median(n[n > 0]) if (n > 0).any() else 0
        ok &= n <= max_fold * med
    if ok.sum() < N_STRATA * MIN_STRATUM_REGIONS:
        return None
    edges = np.linspace(0, 1, N_STRATA + 1)
    refm = rm[ok] / rn[ok]
    s = np.clip(np.digitize(refm, edges) - 1, 0, N_STRATA - 1)
    a_mp, a_np, a_mb, a_nb = mp[ok], np_[ok], mb[ok], nb[ok]
    nreg = np.bincount(s, minlength=N_STRATA)
    keep = nreg >= MIN_STRATUM_REGIONS
    xs, ys, a, b = fit_strata(s, a_mp, a_np, a_mb, a_nb, keep)
    hi = (refm >= HIGH[0]) & (refm <= HIGH[1])
    mid = (refm >= MID[0]) & (refm < MID[1])
    boot = []
    k = len(s)
    for _ in range(N_BOOT_BINS):
        i = rng.integers(0, k, k)
        kp = keep & (np.bincount(s[i], minlength=N_STRATA) >= MIN_STRATUM_REGIONS)
        _, _, ba, bb = fit_strata(s[i], a_mp[i], a_np[i], a_mb[i], a_nb[i], kp)
        boot.append((ba, bb, ratio(hi[i], a_mp[i], a_np[i], a_mb[i], a_nb[i]),
                     ratio(mid[i], a_mp[i], a_np[i], a_mb[i], a_nb[i])))
    ci = np.nanpercentile(np.array(boot), [2.5, 97.5], axis=0)
    mids = (edges[:-1] + edges[1:]) / 2
    strata = [{"ref_meth_mid": mids[j], "regions": int(nreg[j]), "x_meth": xs[j], "y_meth": ys[j],
               "y_over_x": ys[j] / xs[j] if xs[j] > 0 else float("nan")} for j in np.flatnonzero(keep)]
    fit = {"regions": int(ok.sum()), "x_meth": a_mp.sum() / a_np.sum(), "y_meth": a_mb.sum() / a_nb.sum(),
           "intercept": a, "intercept_lo": ci[0, 0], "intercept_hi": ci[1, 0],
           "slope": b, "slope_lo": ci[0, 1], "slope_hi": ci[1, 1],
           "ratio_high": ratio(hi, a_mp, a_np, a_mb, a_nb), "ratio_high_lo": ci[0, 2], "ratio_high_hi": ci[1, 2],
           "ratio_mid": ratio(mid, a_mp, a_np, a_mb, a_nb), "ratio_mid_lo": ci[0, 3], "ratio_mid_hi": ci[1, 3]}
    plot = (a_mp / a_np, a_mb / a_nb, xs[keep], ys[keep], a, b)
    return fit, strata, plot


def ci_txt(r, key):
    return f"{fmt(r[key])} [{fmt(r[key + '_lo'])}, {fmt(r[key + '_hi'])}]"


def stratified(regions_dir, flag_dir, pairs, info, min_calls, max_fold, out):
    mapped = regions_dir / "mapped"
    refs = info.loc[info["type"] == "reference", "label"].tolist()
    if not refs:
        print("\nBAL vs plasma by methylation level: no reference libraries in samples.tsv; skipped")
        return
    rng = np.random.default_rng(SEED)
    srows, frows, plots = [], [], []
    for rset in ("bins", "cpg_islands", "promoters"):
        d = mapped / rset
        flags = region_flags(d, flag_dir, rset, refs, min_calls) if d.is_dir() else None
        if flags is None:
            continue
        rm, rn, base_ok = flags
        for name, pl, bal in pairs:
            if not ((d / f"{pl}.tsv").exists() and (d / f"{bal}.tsv").exists()):
                continue
            res = compare(rm, rn, base_ok, load_counts(d / f"{pl}.tsv"), load_counts(d / f"{bal}.tsv"),
                          min_calls, max_fold, rng)
            if res is None:
                print(f"\n{rset} {name}: too few regions with >= {min_calls} calls; skipped")
                continue
            fit, strata, plot = res
            for st in strata:
                srows.append({"pair": name, "set": rset, "ref_meth_mid": st["ref_meth_mid"],
                              "regions": st["regions"], "plasma_meth": st["x_meth"], "bal_meth": st["y_meth"],
                              "bal_over_plasma": st["y_over_x"]})
            fit = {"pair": name, "set": rset, **fit}
            fit["plasma_meth"], fit["bal_meth"] = fit.pop("x_meth"), fit.pop("y_meth")
            frows.append(fit)
            if rset == "bins":
                plots.append((name, pl, bal, *plot))
    if not frows:
        print("\nBAL vs plasma by methylation level: no step 05 region counts (run step 05 first)")
        return
    pd.DataFrame(srows).to_csv(out / "bal_vs_plasma_strata.tsv", sep="\t", index=False, float_format="%.4f")
    fit = pd.DataFrame(frows)
    fit.to_csv(out / "bal_vs_plasma_fit.tsv", sep="\t", index=False, float_format="%.4f")
    print(f"\nBAL vs plasma by methylation level (strata of reference methylation: {', '.join(refs)})")
    print("Uniform 5mC loss predicts intercept ~0 and ratio_high ~ ratio_mid (both < 1).")
    print("pair\tset\tregions\tplasma\tBAL\tintercept [95% CI]\tslope [95% CI]\t"
          "BAL/PL high-meth [CI]\tBAL/PL mid-meth [CI]")
    for _, r in fit.iterrows():
        print(f"{r['pair']}\t{r['set']}\t{r['regions']}\t{fmt(r['plasma_meth'])}\t{fmt(r['bal_meth'])}\t"
              f"{ci_txt(r, 'intercept')}\t{ci_txt(r, 'slope')}\t{ci_txt(r, 'ratio_high')}\t{ci_txt(r, 'ratio_mid')}")

    fig, axes = plt.subplots(1, len(plots), figsize=(4.6 * len(plots), 4.4), squeeze=False)
    for ax, (name, pl, bal, x, y, sx, sy, a, b) in zip(axes[0], plots):
        ax.hexbin(x, y, gridsize=60, bins="log", cmap="Greys", mincnt=1, extent=(0, 1, 0, 1))
        ax.plot([0, 1], [0, 1], color="black", lw=0.8, ls=":")
        ax.plot(sx, sy, "o", color="#c0392b", ms=4, label="strata (pooled)")
        ax.plot([0, 1], [a, a + b], color="#c0392b", lw=1, label=f"fit: {a:.3f} + {b:.2f}x")
        ax.set_xlim(0, 1); ax.set_ylim(0, 1)
        ax.set_xlabel(f"{pl} methylation"); ax.set_ylabel(f"{bal} methylation")
        ax.set_title(f"{name}: genomic bins")
        ax.legend(fontsize=7, loc="upper left")
    fig.tight_layout()
    fig.savefig(out / "bal_vs_plasma.png", dpi=150)
    plt.close(fig)


def sum_counts(parts):
    names = parts[0][0]
    return names, sum(p[1] for p in parts), sum(p[2] for p in parts)


def size_split(size_dir, flag_dir, regions_dir, pairs, info, classes, split_bp, min_calls, max_fold, out):
    """Methylation by fragment size: which sizes carry the BAL deficit?"""
    refs = info.loc[info["type"] == "reference", "label"].tolist()
    rd, sd = regions_dir / "mapped" / "bins", size_dir / "bins"
    flags = region_flags(rd, flag_dir, "bins", refs, min_calls) if refs and rd.is_dir() else None
    if flags is None or not sd.is_dir() or not classes:
        return
    rm, rn, base_ok = flags
    lo = {c: float(c.split("-")[0]) for c in classes}
    short_cls = [c for c in classes if lo[c] < split_bp]
    long_cls = [c for c in classes if lo[c] >= split_bp]
    rng = np.random.default_rng(SEED)
    rows, strows, crows = [], [], []
    for name, pl, bal in pairs:
        files = {(lab, c): sd / f"{lab}.{c}.tsv" for lab in (pl, bal) for c in classes}
        if not all(p.exists() for p in files.values()):
            continue
        cnt = {k: load_counts(p) for k, p in files.items()}
        c = {"PL_short": sum_counts([cnt[(pl, k)] for k in short_cls]),
             "PL_long": sum_counts([cnt[(pl, k)] for k in long_cls]),
             "BAL_short": sum_counts([cnt[(bal, k)] for k in short_cls]),
             "BAL_long": sum_counts([cnt[(bal, k)] for k in long_cls])}
        calls = {k: v[2].sum() for k, v in c.items()}
        for y, x in (("PL_short", "PL_long"), ("BAL_long", "PL_long"), ("BAL_short", "PL_long"),
                     ("BAL_short", "BAL_long")):
            res = compare(rm, rn, base_ok, c[x], c[y], min_calls, max_fold, rng)
            if res is None:
                continue
            fit, strata, _ = res
            comp = y.split("_")[0]
            rows.append({"pair": name, "comparison": f"{y} / {x}",
                         "share_short_calls": calls[comp + "_short"] / (calls[comp + "_short"] + calls[comp + "_long"]),
                         **fit})
            for st in strata:
                strows.append({"pair": name, "comparison": f"{y} / {x}", **st})
        # each size class against plasma long fragments; classes are sparse, so any
        # region with calls in both is used (selection on coverage, not methylation)
        tot = {lab: sum(cnt[(lab, k)][2].sum() for k in classes) for lab in (pl, bal)}
        for comp, lab in (("PL", pl), ("BAL", bal)):
            for k in classes:
                res = compare(rm, rn, base_ok, c["PL_long"], cnt[(lab, k)], 1, max_fold, rng)
                row = {"pair": name, "library": comp, "size_class": k,
                       "share_calls": cnt[(lab, k)][2].sum() / tot[lab] if tot[lab] else float("nan")}
                if res is not None:
                    fit, _, _ = res
                    row.update({"regions": fit["regions"], "meth": fit["y_meth"],
                                **{f: fit[f] for f in ("ratio_high", "ratio_high_lo", "ratio_high_hi",
                                                       "ratio_mid", "ratio_mid_lo", "ratio_mid_hi")}})
                crows.append(row)
    if not rows:
        print("\nShort vs long fragments: too few regions with enough calls per size class; skipped")
        return
    df = pd.DataFrame(rows)
    df.to_csv(out / "size_split_fit.tsv", sep="\t", index=False, float_format="%.4f")
    pd.DataFrame(strows).to_csv(out / "size_split_strata.tsv", sep="\t", index=False, float_format="%.4f")
    print(f"\nShort (< {split_bp:g} bp) vs long fragments (genomic bins, ratios of pooled methylation; "
          "strata of reference methylation)")
    print("If BAL_long / PL_long is ~1 in both columns, the BAL deficit is carried by short fragments;")
    print("if it keeps the high > mid pattern, BAL long fragments differ too (cell composition).")
    print("pair\tcomparison\tregions\tx_meth\ty_meth\tratio high-meth [CI]\tratio mid-meth [CI]")
    for _, r in df.iterrows():
        print(f"{r['pair']}\t{r['comparison']}\t{r['regions']}\t{fmt(r['x_meth'])}\t{fmt(r['y_meth'])}\t"
              f"{ci_txt(r, 'ratio_high')}\t{ci_txt(r, 'ratio_mid')}")
    for name in df["pair"].unique():
        sub = df[df["pair"] == name].set_index("comparison")["share_short_calls"]
        pl_share = sub.get("PL_short / PL_long", float("nan"))
        bal_share = sub.get("BAL_short / BAL_long", float("nan"))
        print(f"{name}: share of CpG calls from short fragments: plasma {fmt(pl_share, 2)}, BAL {fmt(bal_share, 2)}")
    if crows:
        cdf = pd.DataFrame(crows)
        cdf.to_csv(out / "size_classes.tsv", sep="\t", index=False, float_format="%.4f")
        print("\nBy size class (bp), each against the same timepoint's plasma long fragments")
        print("pair\tlibrary\tsize_class\tshare_calls\tregions\tmeth\tratio high-meth [CI]\tratio mid-meth [CI]")
        for _, r in cdf.iterrows():
            if "ratio_high" in r and pd.notna(r.get("ratio_high")):
                print(f"{r['pair']}\t{r['library']}\t{r['size_class']}\t{fmt(r['share_calls'], 3)}\t"
                      f"{int(r['regions'])}\t{fmt(r['meth'])}\t{ci_txt(r, 'ratio_high')}\t{ci_txt(r, 'ratio_mid')}")
            else:
                print(f"{r['pair']}\t{r['library']}\t{r['size_class']}\t{fmt(r['share_calls'], 3)}\t"
                      "too few regions")



def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--sample-info", required=True)
    ap.add_argument("--imprinted-cpgs", help="label chr pos gene n_meth n_unmeth (no header)")
    ap.add_argument("--regions-dir", help="step 05 output (OUTDIR/downstream/regions)")
    ap.add_argument("--size-dir", help="region counts by fragment size: <set>/<label>.<short|long>.tsv")
    ap.add_argument("--size-classes", default="", help="comma-separated size classes, e.g. 0-100,100-150,...,400-inf")
    ap.add_argument("--split-bp", type=float, default=150, help="short/long boundary (a class edge)")
    ap.add_argument("--min-calls", type=int, default=10, help="calls per region in every library")
    ap.add_argument("--max-cov-fold", type=float, default=10)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    info = read_info(a.sample_info)
    pairs = find_pairs(info)
    if not pairs:
        raise SystemExit("no subject/timepoint with exactly one plasma and one BAL library in samples.tsv")
    print("pairs: " + ", ".join(f"{n} ({p} / {b})" for n, p, b in pairs))
    if a.imprinted_cpgs and Path(a.imprinted_cpgs).exists():
        imprinted(a.imprinted_cpgs, pairs, info, out)
    if a.regions_dir:
        rd = Path(a.regions_dir)
        stratified(rd, rd / "beds", pairs, info, a.min_calls, a.max_cov_fold, out)
        if a.size_dir:
            classes = [c for c in a.size_classes.split(",") if c]
            size_split(Path(a.size_dir), rd / "beds", rd, pairs, info, classes, a.split_bp,
                       a.min_calls, a.max_cov_fold, out)


if __name__ == "__main__":
    main()
