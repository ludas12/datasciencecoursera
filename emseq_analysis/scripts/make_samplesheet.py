#!/usr/bin/env python3
"""Group FASTQ files into samples and read pairs.

Understands the common naming schemes, e.g.
  SAMPLE_S1_L001_R1_001.fastq.gz   (Illumina bcl2fastq/BCLConvert; lanes merged)
  SAMPLE_R1.fq.gz  SAMPLE_1.fastq.gz  SAMPLE.1.fastq.gz  SAMPLE.R1.fastq.gz
Files with no read number are treated as single-end.
Writes a TSV: sample, r1 (comma-separated), r2 (comma-separated or empty).
"""
import re
import sys
from collections import defaultdict
from pathlib import Path

PATTERN = re.compile(
    r"^(?P<sample>.+?)(?P<sid>_S\d+)?(?P<lane>_L\d{3})?"
    r"(?P<sep>[._])(?P<rtag>R?)(?P<read>[12])(?P<suffix>_\d{3})?$"
)
EXT = re.compile(r"\.(fastq|fq)(\.gz)?$")
INDEX_READ = re.compile(r"_I[12](_\d{3})?$")   # Illumina index reads: not used


def main(out_path, files):
    pairs = defaultdict(dict)            # (sample, pairkey) -> {read: path}
    single = defaultdict(list)           # sample -> [path]
    for f in sorted(files):
        stem = EXT.sub("", Path(f).name)
        if INDEX_READ.search(stem):
            print(f"  skipping index-read file {f}")
            continue
        m = PATTERN.match(stem)
        if m:
            key = stem[: m.start("sep")] + (m.group("suffix") or "")
            pairs[(m.group("sample"), key)][m.group("read")] = f
        else:
            single[stem].append(f)

    samples = defaultdict(lambda: ([], []))
    problems = []
    for (sample, key), reads in sorted(pairs.items()):
        if "1" in reads and "2" in reads:
            samples[sample][0].append(reads["1"])
            samples[sample][1].append(reads["2"])
        elif "1" in reads:
            samples[sample][0].append(reads["1"])
        else:
            problems.append(f"R2 without R1: {reads['2']}")
    for sample, fs in single.items():
        if sample in samples:
            problems.append(f"{sample}: has both paired and unpaired files: {fs}")
        samples[sample][0].extend(fs)

    for sample, (r1, r2) in samples.items():
        if r2 and len(r1) != len(r2):
            problems.append(f"{sample}: {len(r1)} R1 files but {len(r2)} R2 files")
    layouts = {bool(r2) for r1, r2 in samples.values()}
    if len(layouts) > 1:
        problems.append("mix of paired-end and single-end samples")

    with open(out_path, "w") as out:
        out.write("sample\tr1\tr2\n")
        for sample in sorted(samples):
            r1, r2 = samples[sample]
            out.write(f"{sample}\t{','.join(r1)}\t{','.join(r2)}\n")

    for sample in sorted(samples):
        r1, r2 = samples[sample]
        kind = "paired-end" if r2 else "single-end"
        print(f"  {sample:30s} {kind:10s} {len(r1)} file(s) per read")
    print(f"{len(samples)} samples written to {out_path}")
    if problems:
        print("\nWARNING -- please check the sample sheet by hand:")
        for p in problems:
            print("  - " + p)
        return 2
    return 0


if __name__ == "__main__":
    if len(sys.argv) < 3:
        sys.exit("usage: make_samplesheet.py OUT.tsv FASTQ [FASTQ ...]")
    sys.exit(main(sys.argv[1], sys.argv[2:]))
