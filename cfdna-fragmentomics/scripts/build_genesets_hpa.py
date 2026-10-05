#!/usr/bin/env python3
"""Build large tissue / cell-type gene sets for TSS profiles from the Human Protein Atlas.

Usage (needs internet once; e.g. the Kaya login node):
    python scripts/build_genesets_hpa.py -o fragmentomics/genesets_hpa.tsv
    python scripts/build_genesets_hpa.py --hpa proteinatlas.tsv.zip -o ...   # use a local copy

Sets are hundreds of genes each instead of the 15-50 curated ones, which reduces the noise of
aggregate TSS profiles at low coverage roughly with the square root of the set size:

  lung_epithelium     single-cell type enriched in alveolar type 1/2, club, ciliated, basal or
                      other respiratory epithelial cells, and not blood-cell specific
  lung_tissue         tissue enriched in lung (bulk), not blood-cell specific
  neutrophil          blood-cell enriched in neutrophils
  monocyte_macrophage blood-cell enriched in monocytes, or single-cell enriched in macrophages
  lymphocyte          blood-cell enriched in T, B or NK cells
  liver               tissue enriched in liver
  housekeeping        low tissue specificity and detected in all tissues (random 2,000)
  inactive_control    tissue enriched (single tissue) in testis, brain, retina, skeletal muscle
                      or pancreas, and not blood-cell specific (random 1,000)
  megakaryocyte_erythroid   copied from the curated fragmentomics/genesets.tsv (no HPA equivalent
                      among mature blood cells)

Column names and category labels differ slightly between HPA releases; the script finds them by
pattern, prints what it used and the size of every set, and stops if a set is too small.
"""
import argparse
import io
import os
import re
import sys
import urllib.request
import zipfile

import numpy as np
import pandas as pd

HPA_URLS = [  # the download path has changed between HPA releases; tried in order
    "https://www.proteinatlas.org/download/proteinatlas.tsv.zip",
    "https://www.proteinatlas.org/download/tsv/proteinatlas.tsv.zip",
    "https://v24.proteinatlas.org/download/proteinatlas.tsv.zip",
    "https://v23.proteinatlas.org/download/proteinatlas.tsv.zip",
]
HERE = os.path.dirname(os.path.abspath(__file__))
CURATED = os.path.join(HERE, "..", "fragmentomics", "genesets.tsv")
ENRICHED = re.compile(r"enriched", re.I)          # "Tissue enriched", "Group enriched", ...
ENR_OR_ENH = re.compile(r"enriched|enhanced", re.I)
MIN_SET = 20


def load_hpa(path):
    if path is None:
        raw, errors = None, []
        for url in HPA_URLS:
            print(f"Downloading {url}", file=sys.stderr)
            try:
                req = urllib.request.Request(url, headers={"User-Agent": "Mozilla/5.0"})
                with urllib.request.urlopen(req, timeout=300) as r:
                    raw = r.read()
                with open("proteinatlas.tsv.zip", "wb") as fh:   # reuse later with --hpa
                    fh.write(raw)
                print("  saved a copy as proteinatlas.tsv.zip", file=sys.stderr)
                break
            except Exception as e:  # noqa: BLE001 - report and try the next mirror
                errors.append(f"{url}: {e}")
        if raw is None:
            sys.exit("Could not download the HPA table:\n  " + "\n  ".join(errors) +
                     "\nDownload proteinatlas.tsv.zip from https://www.proteinatlas.org/about/download"
                     " and rerun with --hpa proteinatlas.tsv.zip")
    else:
        with open(path, "rb") as fh:
            raw = fh.read()
    if raw[:2] == b"PK":
        with zipfile.ZipFile(io.BytesIO(raw)) as z:
            name = next(n for n in z.namelist() if n.endswith(".tsv"))
            raw = z.read(name)
    return pd.read_csv(io.BytesIO(raw), sep="\t", dtype=str, low_memory=False)


def find_col(df, pattern, what):
    hits = [c for c in df.columns if re.fullmatch(pattern, c.strip(), re.I)]
    if not hits:
        rna = [c for c in df.columns if c.lower().startswith("rna")]
        sys.exit(f"Could not find the {what} column (pattern {pattern!r}).\n"
                 f"RNA columns in this file: {rna}")
    print(f"  {what:32s} <- {hits[0]!r}", file=sys.stderr)
    return hits[0]


def parse_specific(cell):
    """'lung: 123.4;liver: 5' -> {'lung': 123.4, 'liver': 5.0}"""
    out = {}
    if isinstance(cell, str):
        for part in cell.split(";"):
            if ":" in part:
                k, v = part.rsplit(":", 1)
                try:
                    out[k.strip().lower()] = float(v)
                except ValueError:
                    pass
    return out


def select(df, cat_col, spec_col, key_pattern, category=ENRICHED, single_only=False):
    keys = re.compile(key_pattern, re.I)
    cat = df[cat_col].fillna("")
    ok = cat.str.contains(category)
    if single_only:   # exclude "Group enriched": one tissue only
        ok &= ~cat.str.contains("group", case=False)
    spec = df[spec_col].map(parse_specific)
    hit = spec.map(lambda d: any(keys.search(k) for k in d))
    return set(df.loc[ok & hit, "Gene"])


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawTextHelpFormatter)
    ap.add_argument("--hpa", help="local proteinatlas.tsv(.zip); default: download")
    ap.add_argument("-o", "--out", default=os.path.join(HERE, "..", "fragmentomics", "genesets_hpa.tsv"))
    ap.add_argument("--n-housekeeping", type=int, default=2000)
    ap.add_argument("--n-inactive", type=int, default=1000)
    ap.add_argument("--seed", type=int, default=1)
    args = ap.parse_args()

    df = load_hpa(args.hpa)
    print(f"HPA table: {len(df):,} genes, {df.shape[1]} columns. Columns used:", file=sys.stderr)
    gene = find_col(df, r"gene", "gene symbol")
    df = df.rename(columns={gene: "Gene"})
    t_cat = find_col(df, r"RNA tissue specificity", "tissue specificity")
    t_dist = find_col(df, r"RNA tissue distribution", "tissue distribution")
    t_spec = find_col(df, r"RNA tissue specific [np]?[CT]PM", "tissue specific nTPM")
    c_cat = find_col(df, r"RNA single cell type specificity", "single-cell type specificity")
    c_spec = find_col(df, r"RNA single cell type specific [np]?[CT]PM", "single-cell type specific nTPM")
    b_cat = find_col(df, r"RNA blood cell specificity", "blood cell specificity")
    b_spec = find_col(df, r"RNA blood cell specific [np]?[CT]PM", "blood cell specific nTPM")

    names = sorted({k for d in df[c_spec].map(parse_specific) for k in d})
    lungish = [n for n in names if re.search(r"alveol|club|cilia|respir|basal|ionocyte|secretory|macrophage", n, re.I)]
    print(f"\nSingle-cell types in this release: {len(names)}; lung/macrophage-related: {lungish}",
          file=sys.stderr)
    blood = sorted({k for d in df[b_spec].map(parse_specific) for k in d})
    print(f"Blood cell types: {blood}", file=sys.stderr)

    blood_specific = set(df.loc[df[b_cat].fillna("").str.contains(ENR_OR_ENH), "Gene"])
    rng = np.random.default_rng(args.seed)

    def sample(genes, n):
        genes = sorted(genes)
        return set(rng.choice(genes, n, replace=False)) if len(genes) > n else set(genes)

    sets = {
        "lung_epithelium": select(df, c_cat, c_spec,
                                  r"alveolar cells type|club|ciliated|basal respiratory|"
                                  r"respiratory|ionocyte|secretory cells") - blood_specific,
        "lung_tissue": select(df, t_cat, t_spec, r"^lung$") - blood_specific,
        "neutrophil": select(df, b_cat, b_spec, r"neutrophil"),
        "monocyte_macrophage": select(df, b_cat, b_spec, r"monocyte")
                               | select(df, c_cat, c_spec, r"macrophage"),
        "lymphocyte": select(df, b_cat, b_spec, r"t-cell|b-cell|nk-cell|\bt cell|\bb cell|nk cell"),
        "liver": select(df, t_cat, t_spec, r"^liver$"),
        "housekeeping": sample(set(df.loc[df[t_cat].fillna("").str.contains("low tissue specificity", case=False)
                                          & df[t_dist].fillna("").str.contains("detected in all", case=False),
                                          "Gene"]), args.n_housekeeping),
        "inactive_control": sample(select(df, t_cat, t_spec,
                                          r"^(testis|brain|cerebral cortex|retina|skeletal muscle|pancreas)$",
                                          single_only=True) - blood_specific, args.n_inactive),
    }
    curated = pd.read_csv(CURATED, sep="\t", comment="#")
    sets["megakaryocyte_erythroid"] = set(curated.gene[curated.set == "megakaryocyte_erythroid"])

    print("\nSet sizes:", file=sys.stderr)
    too_small = []
    for name, genes in sets.items():
        print(f"  {name:24s} {len(genes):6,d}", file=sys.stderr)
        if len(genes) < MIN_SET and name != "megakaryocyte_erythroid":
            too_small.append(name)
    if too_small:
        sys.exit(f"Sets with < {MIN_SET} genes: {too_small}. The category labels or cell-type names "
                 f"may differ in this HPA release; check the columns printed above.")

    rows = [(s, g) for s, genes in sets.items() for g in sorted(genes)]
    with open(args.out, "w") as fh:
        fh.write("# Gene sets from the Human Protein Atlas (proteinatlas.tsv), built by\n"
                 "# scripts/build_genesets_hpa.py. megakaryocyte_erythroid is the curated set.\n")
        pd.DataFrame(rows, columns=["set", "gene"]).to_csv(fh, sep="\t", index=False)
    print(f"\nWrote {len(rows):,} rows to {args.out}", file=sys.stderr)


if __name__ == "__main__":
    main()
