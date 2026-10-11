#!/usr/bin/env python3
"""Sample identity check for EM-seq libraries: sex and genotype concordance.

Genotype score for a pair (i, j): over A/T-reference sites where library i shows
the site's main non-reference allele, and the *other* libraries show it in >= 2
reads at >= 5 % (a real, common variant rather than a sequencing error), count
how often j shows that allele too, relative to what a heterozygous carrier would
give at j's depth, and scale it between what an unrelated person (Hardy-Weinberg,
allele frequency from the other libraries) and the same person would give:
0 ~ unrelated, ~0.5 ~ first-degree relative, 1 ~ same person.

Sex: chrX and chrY read density relative to autosomes (MAPQ >= 30).
Female: X ~1, Y ~0. Male: X ~0.5, Y clearly > 0. A small chrY signal in a
female recipient's plasma or BAL can be donor-derived DNA from a male donor lung.
"""
import argparse
import gzip
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

BASES = set("ACGT")


def read_vcf(path):
    """Per-library depth and read count of the dominant non-reference allele per site."""
    samples, alt_rows, dep_rows = None, [], []
    with gzip.open(path, "rt") as fh:
        for line in fh:
            if line.startswith("##"):
                continue
            f = line.rstrip("\n").split("\t")
            if line.startswith("#"):
                samples = f[9:]
                continue
            alts = f[4].split(",")
            real = [k + 1 for k, a in enumerate(alts) if a in BASES]
            if not real:
                continue
            ad_i = f[8].split(":").index("AD")
            ads = []
            for s in f[9:]:
                ad = [int(x) if x != "." else 0 for x in s.split(":")[ad_i].split(",")]
                ads.append(ad + [0] * (len(alts) + 1 - len(ad)))
            best = max(real, key=lambda k: sum(ad[k] for ad in ads))
            alt_rows.append([ad[best] for ad in ads])
            dep_rows.append([sum(ad) for ad in ads])
    return samples, np.array(alt_rows, dtype=np.int32), np.array(dep_rows, dtype=np.int32)


MIN_ALT_OTHERS = 2       # allele reads required in libraries other than the pair compared
MIN_AF_OTHERS = 0.05     # ... at this allele fraction, so the site is a real, common variant


def concordance(alt, dep):
    n = alt.shape[1]
    has = alt > 0
    alt_tot, dep_tot = alt.sum(axis=1), dep.sum(axis=1)
    score = np.full((n, n), np.nan)
    sites = np.zeros((n, n), dtype=int)
    for i in range(n):
        for j in range(i + 1, n):
            ao = alt_tot - alt[:, i] - alt[:, j]
            do = dep_tot - dep[:, i] - dep[:, j]
            with np.errstate(invalid="ignore", divide="ignore"):
                common = (ao >= MIN_ALT_OTHERS) & (ao >= MIN_AF_OTHERS * do)
            # allele frequency from the other libraries, leaving out the reads that were
            # needed to select the site (otherwise it is biased upwards at low depth)
            p = np.clip((ao - MIN_ALT_OTHERS + 0.5) / np.maximum(do - MIN_ALT_OTHERS + 1, 1), 0.01, 0.99)
            obs = exp_unrel = exp_same = 0.0
            nsite = 0
            for a, b in ((i, j), (j, i)):
                m = common & has[:, a] & (dep[:, b] > 0)
                if not m.any():
                    continue
                pm, see_het = p[m], 1 - 0.5 ** dep[m, b]
                # chance b shows the allele: unrelated person (Hardy-Weinberg) vs the same
                # person as a, who is known to carry it (het or hom-alt)
                exp_unrel += (2 * pm * (1 - pm) * see_het + pm ** 2).sum()
                hom_given_carrier = pm / (2 - pm)
                exp_same += ((1 - hom_given_carrier) * see_het + hom_given_carrier).sum()
                obs += has[m, b].sum()
                nsite += int(m.sum())
            if nsite and exp_same > exp_unrel:
                score[i, j] = score[j, i] = (obs - exp_unrel) / (exp_same - exp_unrel)
                sites[i, j] = sites[j, i] = nsite
    return score, sites


def sex_calls(sex_path, fai_path):
    sx = pd.read_csv(sex_path, sep="\t")
    fai = pd.read_csv(fai_path, sep="\t", header=None, usecols=[0, 1], names=["chr", "len"])
    L = dict(zip(fai["chr"], fai["len"]))
    auto_len = sum(v for k, v in L.items() if k.replace("chr", "").isdigit() and 1 <= int(k.replace("chr", "")) <= 22)
    dens = sx["autosomes"] / auto_len
    sx["rel_X"] = sx["chrX"] / L.get("chrX", np.nan) / dens
    sx["rel_Y"] = sx["chrY"] / L.get("chrY", np.nan) / dens

    def call(r):
        if r["rel_Y"] > 0.05 and r["rel_X"] < 0.75:
            return "male"
        if r["rel_Y"] < 0.02 and r["rel_X"] > 0.8:
            return "female"
        return "unclear"
    sx["sex"] = sx.apply(call, axis=1)
    return sx[["label", "rel_X", "rel_Y", "sex"]]


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--vcf", required=True)
    ap.add_argument("--sex", required=True)
    ap.add_argument("--fai", required=True)
    ap.add_argument("--sample-info", required=True)
    ap.add_argument("--out", required=True, type=Path)
    a = ap.parse_args()

    info = pd.read_csv(a.sample_info, sep="\t", dtype=str)
    if "subject" not in info.columns:
        info["subject"] = info["label"]
    subj = dict(zip(info["label"], info["subject"].fillna(info["label"])))

    labels, alt, dep = read_vcf(a.vcf)
    score, sites = concordance(alt, dep)
    n = len(labels)
    mat = pd.DataFrame(score, index=labels, columns=labels)
    mat.round(3).to_csv(a.out / "identity_matrix.tsv", sep="\t")
    pd.DataFrame(sites, index=labels, columns=labels).to_csv(a.out / "identity_sites.tsv", sep="\t")

    same = [score[i, j] for i in range(n) for j in range(i + 1, n)
            if subj.get(labels[i]) == subj.get(labels[j]) and not np.isnan(score[i, j])]
    diff = [score[i, j] for i in range(n) for j in range(i + 1, n)
            if subj.get(labels[i]) != subj.get(labels[j]) and not np.isnan(score[i, j])]
    cut = 0.8    # above: same person. Unrelated people score ~0-0.6 (higher with few libraries)
    sx = sex_calls(a.sex, a.fai).set_index("label")
    rows = []
    for i, lab in enumerate(labels):
        s = subj.get(lab, lab)
        mates = [j for j in range(n) if j != i and subj.get(labels[j]) == s]
        others = [j for j in range(n) if j != i and subj.get(labels[j]) != s]
        best = int(np.nanargmax(np.where(np.isnan(score[i]), -1, score[i])))
        flags = []
        if mates:
            worst_mate = np.nanmin(score[i, mates])
            if worst_mate < cut:
                flags.append(f"does not match own subject (min score {worst_mate:.2f})")
        if others:
            j = others[int(np.nanargmax(np.where(np.isnan(score[i, others]), -1, score[i, others])))]
            if score[i, j] >= cut:
                flags.append(f"matches {labels[j]} ({subj.get(labels[j])}) score {score[i, j]:.2f}")
        sexes = {sx.loc[labels[j], "sex"] for j in mates + [i] if labels[j] in sx.index} - {"unclear"}
        if len(sexes) > 1:
            flags.append("sex differs within subject")
        rows.append({"label": lab, "subject": s,
                     "sex": sx.loc[lab, "sex"] if lab in sx.index else "NA",
                     "rel_X": sx.loc[lab, "rel_X"] if lab in sx.index else np.nan,
                     "rel_Y": sx.loc[lab, "rel_Y"] if lab in sx.index else np.nan,
                     "best_match": labels[best], "best_score": score[i, best],
                     "mean_same_subject": np.nanmean(score[i, mates]) if mates else np.nan,
                     "max_other_subject": np.nanmax(score[i, others]) if others else np.nan,
                     "status": "CHECK" if flags else "OK", "flags": "; ".join(flags)})
    res = pd.DataFrame(rows)
    res.to_csv(a.out / "identity_check.tsv", sep="\t", index=False, float_format="%.3f")

    fig, ax = plt.subplots(figsize=(0.45 * n + 3, 0.45 * n + 2))
    im = ax.imshow(score, cmap="magma", vmin=0, vmax=1.1)
    ax.set_xticks(range(n), labels, rotation=90, fontsize=7)
    ax.set_yticks(range(n), [f"{l} [{subj.get(l)}]" for l in labels], fontsize=7)
    fig.colorbar(im, ax=ax, label="relatedness (0 unrelated, 1 same person)")
    ax.set_title("Sample identity")
    fig.tight_layout()
    fig.savefig(a.out / "identity_matrix.png", dpi=150)

    print(f"{alt.shape[0]:,} informative A/T sites; same-subject median score "
          f"{np.median(same) if same else float('nan'):.2f}, different-subject median "
          f"{np.median(diff) if diff else float('nan'):.2f} (same person if > {cut})")
    with pd.option_context("display.width", 200, "display.max_columns", 20):
        print(res[["label", "subject", "sex", "rel_X", "rel_Y", "best_match", "best_score",
                   "max_other_subject", "status"]].round(3).to_string(index=False))
    for r in rows:
        if r["flags"]:
            print(f"  ! {r['label']}: {r['flags']}")
    print(f"\nFull results: {a.out}/identity_check.tsv, identity_matrix.tsv/.png")


if __name__ == "__main__":
    main()
