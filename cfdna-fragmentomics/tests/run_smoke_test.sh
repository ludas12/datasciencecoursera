#!/usr/bin/env bash
# End-to-end smoke test on simulated data: two "healthy" and one "shortened" sample.
set -euo pipefail
cd "$(dirname "$0")/.."
T=$(mktemp -d)
python tests/make_test_data.py "$T" 0  healthy1 1
python tests/make_test_data.py "$T" 0  healthy2 2   # same reference (seeded from the outdir's ref)
python tests/make_test_data.py "$T" 15 cancer1  3
# TSS BED for every gene-set gene at random positions on the simulated genome
python - "$T" <<'PY'
import sys, numpy as np, pandas as pd
g = pd.read_csv("fragmentomics/genesets.tsv", sep="\t", comment="#").gene.unique()
rng = np.random.default_rng(0)
pd.DataFrame({"c": rng.choice(["chr1", "chr2"], len(g)), "s": rng.integers(10000, 1900000, len(g))}
            ).assign(e=lambda d: d.s + 1, g=g, sc=0, st=rng.choice(["+", "-"], len(g))
            ).to_csv(f"{sys.argv[1]}/tss.bed", sep="\t", header=False, index=False)
PY
for s in healthy1 healthy2 cancer1; do
  python fragmentomics/cfdna_frag.py extract -b "$T/$s.bam" -r "$T/ref.fa" -o "$T/$s.frags.tsv.gz" --blacklist "$T/blacklist.bed"
  python fragmentomics/cfdna_frag.py features -f "$T/$s.frags.tsv.gz" -r "$T/ref.fa" -o "$T/features/$s" --bin-size 500000
  python fragmentomics/cfdna_frag.py tss -f "$T/$s.frags.tsv.gz" -t "$T/tss.bed" -o "$T/features/$s"
done
printf "sample\tgroup\tpatient\nhealthy1\thealthy\tP1\nhealthy2\thealthy\tP2\ncancer1\tcase\tP1\n" > "$T/samples.tsv"
python fragmentomics/cfdna_frag.py cohort -i "$T"/features/* -s "$T/samples.tsv" -o "$T/cohort"
python - "$T" <<'PY'
import sys, pandas as pd
s = pd.read_csv(f"{sys.argv[1]}/cohort/size_motif_summary.tsv", sep="\t", index_col=0)
print(s[["n_fragments", "mode_len", "short_long_ratio", "motif_diversity_score", "period10_power"]])
assert s.loc["cancer1", "short_long_ratio"] > s.loc["healthy1", "short_long_ratio"]
print("OK")
PY
echo "outputs in $T"
