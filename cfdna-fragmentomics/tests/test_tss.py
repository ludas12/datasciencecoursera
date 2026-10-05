"""TSS profile test with a known answer: fragments are depleted around 'housekeeping' TSSs
(both strands) and untouched around 'inactive_control' TSSs."""
import os
import subprocess
import sys
import tempfile

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
TOOL = os.path.join(HERE, "..", "fragmentomics", "cfdna_frag.py")


def main():
    rng = np.random.default_rng(7)
    d = tempfile.mkdtemp()
    L, n = 20_000_000, 2_000_000
    sets = pd.read_csv(os.path.join(HERE, "..", "fragmentomics", "genesets.tsv"), sep="\t", comment="#")
    active = sets.gene[sets.set == "housekeeping"].tolist()
    inactive = sets.gene[sets.set == "inactive_control"].tolist()
    genes = active + inactive
    pos = np.sort(rng.choice(np.arange(10_000, L - 10_000, 12_000), len(genes), replace=False))
    strand = rng.choice(["+", "-"], len(genes))
    # 6-column BED, 1-bp TSS intervals, name in col 4, strand in col 6
    pd.DataFrame({"c": "chr1", "s": pos, "e": pos + 1, "g": genes, "sc": 0, "st": strand}).to_csv(
        f"{d}/tss.bed", sep="\t", header=False, index=False)

    start = rng.integers(0, L - 400, n)
    length = rng.normal(167, 15, n).astype(int).clip(100, 220)
    end = start + length
    keep = np.ones(n, bool)
    # deplete only DOWNSTREAM of each active TSS (in transcript orientation), so the test also
    # checks that minus-strand profiles are flipped
    for p, st in zip(pos[: len(active)], strand[: len(active)]):
        lo, hi = (p, p + 300) if st == "+" else (p - 300, p)
        hit = (start < hi) & (end > lo)
        keep[hit & (rng.random(n) < 0.7)] = False
    f = pd.DataFrame({"chrom": "chr1", "start": start[keep], "end": end[keep], "length": length[keep],
                      "mapq": 60, "gc": 0.4, "motif_5p": "ACGT", "motif_3p": "ACGT"})
    f.sort_values("start").to_csv(f"{d}/frags.tsv.gz", sep="\t", index=False)

    subprocess.run([sys.executable, TOOL, "tss", "-f", f"{d}/frags.tsv.gz", "-t", f"{d}/tss.bed",
                    "-o", f"{d}/out"], check=True)
    s = pd.read_csv(f"{d}/out/tss_summary.tsv", sep="\t").set_index("set")
    print(s)
    hk, ctrl = s.loc["housekeeping", "central_coverage"], s.loc["inactive_control", "central_coverage"]
    assert s.loc["housekeeping", "n_genes_found"] == len(active)
    assert hk < 0.85, hk
    assert 0.9 < ctrl < 1.1, ctrl
    prof = pd.read_csv(f"{d}/out/tss_profiles.tsv", sep="\t", index_col=0)["housekeeping"]
    down, up = prof.loc[100:250].mean(), prof.loc[-250:-100].mean()
    print(f"downstream {down:.2f}  upstream {up:.2f}")
    assert down < 0.6 and up > 0.8, (down, up)   # dip is downstream for both strands
    print("OK")


if __name__ == "__main__":
    main()
