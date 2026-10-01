#!/usr/bin/env python3
"""Summarise a directory of FastQC *_fastqc.zip files into one table.

Usage: python scripts/00_fastqc_summary.py FASTQC_DIR [> fastqc_summary.tsv]

Besides the usual totals it reports base composition over the first 10 cycles, which shows
whether the library prep leaves non-templated / biased bases at read starts. That matters
for fragmentomics: read start = fragment end, so those bases would shift fragment
lengths and end motifs.
"""
import glob
import os
import re
import sys
import zipfile


def modules(text):
    out, name, rows = {}, None, []
    for line in text.splitlines():
        if line.startswith(">>END_MODULE"):
            out[name] = rows
            name, rows = None, []
        elif line.startswith(">>"):
            name, rows = line[2:].split("\t")[0], []
        elif name and line:
            rows.append(line.split("\t"))
    return out


def main(d):
    cols = ["file", "total_reads", "read_length", "pct_gc", "pct_dup_remaining",
            "first10_A", "first10_C", "first10_G", "first10_T", "max_cycle_bias_1to10",
            "adapter_max_pct", "top_overrepresented"]
    if not os.path.isdir(d):
        sys.exit(f"Not a directory: {d}")
    zips = sorted(glob.glob(os.path.join(d, "**", "*_fastqc.zip"), recursive=True))
    if not zips:
        sys.exit(f"No *_fastqc.zip files found in {os.path.abspath(d)} (searched subfolders too)")
    print(f"Found {len(zips)} FastQC zips", file=sys.stderr)
    print("\t".join(cols))
    for z in zips:
        with zipfile.ZipFile(z) as zf:
            data = zf.read(next(n for n in zf.namelist() if n.endswith("fastqc_data.txt"))).decode()
        m = modules(data)
        basic = {r[0]: r[1] for r in m["Basic Statistics"] if not r[0].startswith("#")}
        dup = next((r[1] for r in m.get("Sequence Duplication Levels", [])
                    if r[0].startswith("#Total Deduplicated")), "NA")
        # per-base content rows: Base, G, A, T, C  (Base can be "1" or "10-14")
        early = [r for r in m["Per base sequence content"][1:] if int(re.split("-", r[0])[0]) <= 10]
        g, a, t, c = (sum(float(r[i]) for r in early) / len(early) for i in (1, 2, 3, 4))
        bias = max(max(float(x) for x in r[1:5]) - min(float(x) for x in r[1:5]) for r in early)
        ad = m.get("Adapter Content", [])
        ad_max = max((float(x) for r in ad[1:] for x in r[1:]), default=0.0)
        over = m.get("Overrepresented sequences", [])
        top = f"{over[1][0][:30]}({float(over[1][2]):.2f}%)" if len(over) > 1 else "none"
        print("\t".join([os.path.basename(z).replace("_fastqc.zip", ""),
                         basic["Total Sequences"], basic["Sequence length"], basic["%GC"],
                         dup, f"{a:.1f}", f"{c:.1f}", f"{g:.1f}", f"{t:.1f}", f"{bias:.1f}",
                         f"{ad_max:.1f}", top]))


if __name__ == "__main__":
    if len(sys.argv) != 2:
        sys.exit(__doc__)
    main(sys.argv[1])
