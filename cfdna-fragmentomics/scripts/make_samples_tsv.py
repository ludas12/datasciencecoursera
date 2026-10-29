#!/usr/bin/env python3
"""Build a sample sheet (samples.tsv) from a folder of paired FASTQs.

Usage:
    python scripts/make_samples_tsv.py FASTQ_DIR [-o samples_batch2.tsv] [--batch batch2]

Expects facility names like
    RL8957_2026_09_08_HL_LC017_HC7_Plasma_S15_R1_001.fastq.gz
    <run>_HL_LC001_LR_CF0020A_Plasma_S1_L001_R1_001.fastq.gz
and derives, per library:
    sample     CF0020A_P / CF0020A_BAL / HC7_P   (visit ID + sample type)
    group      healthy_plasma | ltx_plasma | ltx_bal
    patient    visit ID (CF0020A): joins a visit's plasma and BAL in the plots
    subject    recipient (CF0020), timepoint (A), library (LC001), batch
    diagnosis  left empty: fill it in from the clinical records before running cohort
Files split over several lanes (_L001_, _L002_ ...) are listed in a merge script, because
run_all.sh takes one R1 and one R2 per sample. Names it cannot parse are reported, not guessed.
"""
import argparse
import os
import re
import sys
from collections import defaultdict

FASTQ = re.compile(r"^(?P<stem>.+?)(?:_L(?P<lane>\d{3}))?_R(?P<read>[12])_\d{3}\.f(?:ast)?q\.gz$")
LIBRARY = re.compile(r"LC_?0*(\d+)", re.I)
RECIPIENT = re.compile(r"CF_?(\d+)([A-Z])(?=_|$)", re.I)
HEALTHY = re.compile(r"(?<![A-Z])(HC_?\d+)(?=_|$)", re.I)


def parse(stem):
    """Return the sheet fields for one library, or None if the name is not recognised."""
    kind = "BAL" if re.search(r"(?:^|_)BAL(?:_|$)", stem, re.I) else (
        "P" if re.search(r"(?:^|_)Plasma(?:_|$)", stem, re.I) else None)
    lib = LIBRARY.search(stem)
    rec, hc = RECIPIENT.search(stem), HEALTHY.search(stem)
    if kind is None or not (rec or hc):
        return None
    if rec:
        subject = "CF" + rec.group(1).zfill(4)
        timepoint = rec.group(2).upper()
        visit = subject + timepoint
        group = "ltx_bal" if kind == "BAL" else "ltx_plasma"
    else:
        visit = subject = hc.group(1).upper().replace("_", "")
        timepoint = ""
        group = "healthy_bal" if kind == "BAL" else "healthy_plasma"
    return {"sample": f"{visit}_{kind}", "group": group, "patient": visit, "subject": subject,
            "timepoint": timepoint, "library": f"LC{int(lib.group(1)):03d}" if lib else ""}


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawTextHelpFormatter)
    ap.add_argument("fastq_dir")
    ap.add_argument("-o", "--out", default="samples_batch2.tsv")
    ap.add_argument("--batch", default="batch2")
    ap.add_argument("--merged-dir", default=None, help="where merged multi-lane FASTQs go "
                    "(default: <fastq_dir>/merged)")
    args = ap.parse_args()

    fq_dir = os.path.abspath(args.fastq_dir)
    if not os.path.isdir(fq_dir):
        sys.exit(f"{fq_dir} is not a folder: give the folder that holds the FASTQ files")
    files = defaultdict(lambda: {"1": [], "2": []})
    unparsed = []
    for name in sorted(os.listdir(fq_dir)):
        m = FASTQ.match(name)
        if not m:
            continue
        files[m.group("stem")][m.group("read")].append(os.path.join(fq_dir, name))
    if not files:
        sys.exit(f"No *_R1_001.fastq.gz / *_R2_001.fastq.gz files in {fq_dir}")

    merged_dir = args.merged_dir or os.path.join(fq_dir, "merged")
    rows, merges, seen = [], [], {}
    for stem, reads in sorted(files.items()):
        info = parse(stem)
        if info is None:
            unparsed.append(stem)
            continue
        if len(reads["1"]) != len(reads["2"]) or not reads["1"]:
            sys.exit(f"{stem}: {len(reads['1'])} R1 vs {len(reads['2'])} R2 files")
        if info["sample"] in seen:
            sys.exit(f"Two libraries map to sample {info['sample']}: {seen[info['sample']]} and {stem}. "
                     "Rename one in the sheet by hand.")
        seen[info["sample"]] = stem
        if len(reads["1"]) == 1:
            r1, r2 = reads["1"][0], reads["2"][0]
        else:  # several lanes: concatenate (gzip streams concatenate safely)
            r1 = os.path.join(merged_dir, f"{info['sample']}_R1.fastq.gz")
            r2 = os.path.join(merged_dir, f"{info['sample']}_R2.fastq.gz")
            merges.append(f"cat {' '.join(sorted(reads['1']))} > {r1}")
            merges.append(f"cat {' '.join(sorted(reads['2']))} > {r2}")
        rows.append([info["sample"], info["group"], r1, r2, info["patient"], "",
                     info["subject"], info["timepoint"], info["library"], args.batch])

    with open(args.out, "w") as fh:
        fh.write("sample\tgroup\tr1\tr2\tpatient\tdiagnosis\tsubject\ttimepoint\tlibrary\tbatch\n")
        for r in rows:
            fh.write("\t".join(r) + "\n")
    print(f"Wrote {len(rows)} samples to {args.out}. Fill in the diagnosis column before running.",
          file=sys.stderr)
    if merges:
        script = os.path.splitext(args.out)[0] + "_merge_lanes.sh"
        with open(script, "w") as fh:
            fh.write(f"#!/bin/bash\nset -euo pipefail\nmkdir -p {merged_dir}\n" + "\n".join(merges) + "\n")
        print(f"{len(merges) // 2} samples span several lanes: run  bash {script}  first.",
              file=sys.stderr)
    if unparsed:
        print("Not recognised (add these rows by hand):\n  " + "\n  ".join(unparsed), file=sys.stderr)


if __name__ == "__main__":
    main()
