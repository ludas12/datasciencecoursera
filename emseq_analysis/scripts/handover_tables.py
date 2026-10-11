#!/usr/bin/env python3
"""Build the QC-review handover tables from an analysed OUTDIR (standard library only).

usage: handover_tables.py OUTDIR HANDOVER_DIR SAMPLE_INFO CLIP_R1 CLIP_R2
Writes into HANDOVER_DIR:
  alignment_qc.csv, spikein_controls.csv, conversion_proxies.csv,
  fragment_lengths_<label>.tsv (0-600 bp), methylation_1Mb_bins.csv,
  chr20_cpg/<label>_chr20_CpG.bedGraph.gz
"""
import csv
import gzip
import json
import re
import sys
from pathlib import Path

outdir, hdir, sample_info = Path(sys.argv[1]), Path(sys.argv[2]), Path(sys.argv[3])
default_clip = (int(sys.argv[4]), int(sys.argv[5]))
BIN = 1_000_000
MAIN = re.compile(r"^chr([0-9]+|X|Y)$")


def rows_tsv(p):
    if not p.exists():
        return []
    with open(p) as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def text(p):
    return p.read_text() if p.exists() else ""


samples = [r for r in rows_tsv(sample_info) if r.get("include", "yes") == "yes"]
qc = {r["sample"]: r for r in rows_tsv(outdir / "summary" / "qc_summary.tsv")}
fai = [l.split("\t") for l in text(outdir / "reference" / "genome_with_controls.fa.fai").splitlines()]
main_len = sum(int(f[1]) for f in fai if MAIN.match(f[0]))


def clip_of(s):
    t = text(outdir / "samples" / s / "qc" / f"{s}.clip5p.txt").split()
    return (int(t[0]), int(t[1])) if len(t) == 2 else default_clip


def g(r, k):
    v = r.get(k, "NA") if r else "NA"
    return "NA" if v in ("", None) else v


# ---- C. alignment and duplication ---------------------------------------------------
align, spike, conv = [], [], []
for r in samples:
    s, lab = r["sample"], r["label"]
    q = outdir / "samples" / s / "qc"
    fp = json.loads(text(q / f"{s}.fastp.json") or "{}").get("summary", {})
    flag = text(q / f"{s}.flagstat.txt")

    def fcount(label):
        m = re.search(rf"^(\d+) \+ \d+ {label}\b", flag, re.M)
        return int(m.group(1)) if m else None
    primary, pmapped, pdup = fcount("primary"), fcount("primary mapped"), fcount("primary duplicates")
    idx = {f[0]: int(f[2]) for f in (l.split("\t") for l in text(q / f"{s}.idxstats.txt").splitlines()) if len(f) > 3}
    mapped_total = sum(idx.values()) or None
    st = text(q / f"{s}.samtools_stats.txt")
    bases = re.search(r"^SN\tbases mapped \(cigar\):\t(\d+)", st, re.M)
    hist = {}
    for m in re.finditer(r"^IS\t(\d+)\t(\d+)", st, re.M):
        hist[int(m.group(1))] = int(m.group(2))
    tot = sum(hist.values())
    med = None
    if tot:
        acc = 0
        for k in sorted(hist):
            acc += hist[k]
            if acc >= tot / 2:
                med = k
                break
    c1, c2 = clip_of(s)
    dup = float(g(qc.get(s), "pct_duplicates")) if g(qc.get(s), "pct_duplicates") != "NA" else None
    align.append({
        "library": lab, "sample": s,
        "raw_read_pairs": fp.get("before_filtering", {}).get("total_reads", 0) // 2 or "NA",
        "post_trim_pairs": fp.get("after_filtering", {}).get("total_reads", 0) // 2 or "NA",
        "mapping_rate_pct": g(qc.get(s), "pct_mapped"),
        "duplicate_rate_pct": g(qc.get(s), "pct_duplicates"),
        "optical_duplicate_pct": g(qc.get(s), "pct_optical_dup"),
        "pcr_duplicate_pct": g(qc.get(s), "pct_pcr_dup"),
        "unique_dedup_pairs": (pmapped - pdup) // 2 if pmapped is not None and pdup is not None else "NA",
        "mean_genome_coverage_approx": round(int(bases.group(1)) * (1 - (dup or 0) / 100) / main_len, 3)
        if bases and main_len else "NA",
        "median_insert_tlen_bp": med if med is not None else "NA",
        "read_start_clip_bp": f"{c1}+{c2}",
        "median_fragment_corrected_bp": med + c1 + c2 if med is not None else "NA",
        "pct_chrM": round(100 * idx.get("chrM", 0) / mapped_total, 4) if mapped_total else "NA",
        "pct_unmapped": round(100 - 100 * pmapped / primary, 3) if primary and pmapped is not None else "NA",
    })
    # fragment length table: step 09 output if present, else samtools TLEN (+ clip correction)
    with open(hdir / f"fragment_lengths_{lab}.tsv", "w") as fh:
        fh.write("length_bp\tcount\tsource\n")
        for L in range(0, 601):
            src = L - c1 - c2
            fh.write(f"{L}\t{hist.get(src, 0)}\tsamtools_TLEN+{c1 + c2}bp\n")

    q_ = qc.get(s, {})
    spike.append({
        "library": lab, "lambda_reads_mapped": idx.get("phage_lambda", "NA"),
        "puc19_reads_mapped": idx.get("plasmid_puc19c", "NA"),
        "lambda_CpG_calls": g(q_, "lambda_CpG_calls"), "lambda_pct_meth_CpG": g(q_, "lambda_pct_meth_CpG"),
        "lambda_pct_meth_nonCpG": g(q_, "lambda_pct_meth_nonCpG"),
        "lambda_conversion_pct_all_C": g(q_, "conversion_efficiency_pct"),
        "puc19_CpG_calls": g(q_, "puc19_CpG_calls"), "puc19_pct_meth_CpG": g(q_, "puc19_pct_meth_CpG")})

    chrm = text(hdir / "work" / f"{s}.chrM.txt").strip()
    conv.append({
        "library": lab, "global_mCpG_pct": g(q_, "genome_pct_meth_CpG"),
        "mCHH_pct_chr22": g(q_, "genome_pct_meth_CHH"), "CHH_calls_chr22": g(q_, "CHH_calls"),
        "mCHG_pct": "not computed", "chrM_mCpG_pct": chrm or "not computed",
        "LINE1_mCpG_pct": "not computed",
        "imprinted_promoters_mCpG_pct": g(q_, "imprinted_pct_meth"),
        "imprinted_genes_used": g(q_, "imprinted_genes_used"),
        "CpGs_ge1x": g(q_, "CpGs_covered"), "CpGs_ge5x": g(q_, "CpGs_cov5"), "CpGs_ge10x": g(q_, "CpGs_cov10"),
        "mean_calls_per_covered_CpG": g(q_, "mean_CpG_depth")})


def write(name, rows):
    if rows:
        with open(hdir / name, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=list(rows[0]))
            w.writeheader()
            w.writerows(rows)


write("alignment_qc.csv", align)
write("spikein_controls.csv", spike)
write("conversion_proxies.csv", conv)

# ---- F. 1 Mb binned methylation matrix + chr20 per-CpG ----------------------------------
bins = {}
labels = []
(hdir / "chr20_cpg").mkdir(exist_ok=True)
for r in samples:
    s, lab = r["sample"], r["label"]
    bg = outdir / "samples" / s / "meth" / f"{s}_CpG.bedGraph.gz"
    if not bg.exists():
        continue
    labels.append(lab)
    with gzip.open(bg, "rt") as fh, gzip.open(hdir / "chr20_cpg" / f"{lab}_chr20_CpG.bedGraph.gz", "wt") as c20:
        for line in fh:
            if line.startswith("track"):
                c20.write(line)
                continue
            f = line.split("\t", 6)
            if not MAIN.match(f[0]):
                continue
            if f[0] == "chr20":
                c20.write(line)
            key = (f[0], int(f[1]) // BIN)
            d = bins.setdefault(key, {})
            m, u = d.get(lab, (0, 0))
            d[lab] = (m + int(f[4]), u + int(f[5]))


def ckey(k):
    c = k[0][3:]
    return (int(c) if c.isdigit() else 100 + ord(c[0]), k[1])


with open(hdir / "methylation_1Mb_bins.csv", "w", newline="") as fh:
    w = csv.writer(fh)
    w.writerow(["chrom", "start", "end"] + [f"{l}_beta" for l in labels] + [f"{l}_calls" for l in labels])
    for k in sorted(bins, key=ckey):
        d = bins[k]
        betas, calls = [], []
        for l in labels:
            m, u = d.get(l, (0, 0))
            betas.append(f"{m / (m + u):.4f}" if m + u else "NA")
            calls.append(m + u)
        w.writerow([k[0], k[1] * BIN, (k[1] + 1) * BIN] + betas + calls)
print(f"handover tables written for {len(samples)} libraries", file=sys.stderr)
