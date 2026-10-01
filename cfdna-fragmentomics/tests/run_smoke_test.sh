#!/usr/bin/env bash
# End-to-end smoke test on simulated data: two "healthy" and one "shortened" sample.
set -euo pipefail
cd "$(dirname "$0")/.."
T=$(mktemp -d)
python tests/make_test_data.py "$T" 0  healthy1 1
python tests/make_test_data.py "$T" 0  healthy2 2   # same reference (seeded from the outdir's ref)
python tests/make_test_data.py "$T" 15 cancer1  3
for s in healthy1 healthy2 cancer1; do
  python fragmentomics/cfdna_frag.py extract -b "$T/$s.bam" -r "$T/ref.fa" -o "$T/$s.frags.tsv.gz" --blacklist "$T/blacklist.bed"
  python fragmentomics/cfdna_frag.py features -f "$T/$s.frags.tsv.gz" -r "$T/ref.fa" -o "$T/features/$s" --bin-size 500000
done
python fragmentomics/cfdna_frag.py cohort -i "$T"/features/* -o "$T/cohort"
python - "$T" <<'PY'
import sys, pandas as pd
s = pd.read_csv(f"{sys.argv[1]}/cohort/size_motif_summary.tsv", sep="\t", index_col=0)
print(s[["n_fragments", "mode_len", "short_long_ratio", "motif_diversity_score", "period10_power"]])
assert s.loc["cancer1", "short_long_ratio"] > s.loc["healthy1", "short_long_ratio"]
print("OK")
PY
echo "outputs in $T"
