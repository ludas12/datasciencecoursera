"""GC-correction test with a known answer.

A synthetic genome has GC-rich (CpG-island-like) promoters around every TSS. Fragments are
sampled uniformly, then kept with a probability that rises with fragment GC (library GC bias),
and active (housekeeping) TSSs are depleted. Without correction the GC-rich inactive TSSs show a
false coverage peak; with correction (-r REF) they return to ~1 and the housekeeping dip remains."""
import os
import subprocess
import sys
import tempfile

import numpy as np
import pandas as pd
import pysam

HERE = os.path.dirname(os.path.abspath(__file__))
TOOL = os.path.join(HERE, "..", "fragmentomics", "cfdna_frag.py")


def main():
    rng = np.random.default_rng(3)
    d = tempfile.mkdtemp()
    L = 20_000_000
    sets = pd.read_csv(os.path.join(HERE, "..", "fragmentomics", "genesets.tsv"), sep="\t", comment="#")
    active = sets.gene[sets.set == "housekeeping"].tolist()
    inactive = sets.gene[sets.set == "inactive_control"].tolist()
    genes = active + inactive
    pos = np.sort(rng.choice(np.arange(20_000, L - 20_000, 15_000), len(genes), replace=False))
    strand = rng.choice(["+", "-"], len(genes))

    gc_p = np.full(L, 0.40)
    # GC-rich islands across the genome (most GC-rich DNA is not at the tested promoters)
    for c in rng.integers(10_000, L - 10_000, 3000):
        gc_p[c - 600:c + 600] = 0.68
    for p in pos:
        gc_p[p - 600:p + 600] = 0.68                      # GC-rich promoter
    is_gc = rng.random(L) < gc_p
    seq = np.where(is_gc, rng.choice(np.array(list("GC")), L), rng.choice(np.array(list("AT")), L))
    with open(f"{d}/ref.fa", "w") as fh:
        fh.write(">chr1\n")
        s = "".join(seq)
        fh.write("\n".join(s[i:i + 80] for i in range(0, L, 80)) + "\n")
    pysam.faidx(f"{d}/ref.fa")
    pd.DataFrame({"c": "chr1", "s": pos, "e": pos + 1, "g": genes, "sc": 0, "st": strand}).to_csv(
        f"{d}/tss.bed", sep="\t", header=False, index=False)

    n = 4_000_000
    start = rng.integers(0, L - 400, n)
    length = rng.normal(167, 15, n).astype(int).clip(100, 220)
    end = start + length
    cum = np.concatenate([[0], np.cumsum(is_gc)])
    gc = (cum[end] - cum[start]) / length
    keep = rng.random(n) < np.clip(0.15 + 1.6 * (gc - 0.3), 0.02, 1)   # GC-biased library
    for p in pos[: len(active)]:
        hit = (start < p + 150) & (end > p - 150)
        keep[hit & (rng.random(n) < 0.6)] = False
    pd.DataFrame({"chrom": "chr1", "start": start[keep], "end": end[keep], "length": length[keep],
                  "mapq": 60, "gc": gc[keep].round(4), "motif_5p": "ACGT", "motif_3p": "ACGT"}
                 ).sort_values("start").to_csv(f"{d}/frags.tsv.gz", sep="\t", index=False)

    res = {}
    for label, extra in [("raw", []), ("gc", ["-r", f"{d}/ref.fa"])]:
        subprocess.run([sys.executable, TOOL, "tss", "-f", f"{d}/frags.tsv.gz", "-t", f"{d}/tss.bed",
                        "-o", f"{d}/{label}"] + extra, check=True, capture_output=True)
        res[label] = pd.read_csv(f"{d}/{label}/tss_summary.tsv", sep="\t").set_index("set")["central_coverage"]
    print(pd.DataFrame(res))
    assert res["raw"]["inactive_control"] > 1.3, "test setup: GC bias should create a false peak"
    assert 0.85 < res["gc"]["inactive_control"] < 1.15, res["gc"]["inactive_control"]
    assert res["gc"]["housekeeping"] < 0.85 * res["gc"]["inactive_control"]
    print("OK")


if __name__ == "__main__":
    main()
