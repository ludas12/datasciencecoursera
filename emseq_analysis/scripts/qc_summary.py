#!/usr/bin/env python3
"""Collect per-sample QC metrics into one table and flag problem samples.

usage: qc_summary.py OUTDIR CONFIG_KEY=VALUE ...
Reads OUTDIR/samples/<sample>/{qc,meth,mbias} and writes
OUTDIR/summary/qc_summary.tsv.
"""
import gzip
import json
import re
import sys
from pathlib import Path


def read_text(path):
    return path.read_text() if path.exists() else ""


def fastp_metrics(qc, s):
    p = qc / f"{s}.fastp.json"
    if not p.exists():
        return {}
    j = json.loads(p.read_text())["summary"]
    before, after = j["before_filtering"], j["after_filtering"]
    return {
        "raw_reads": before["total_reads"],
        "reads_after_trim": after["total_reads"],
        "pct_reads_kept": 100 * after["total_reads"] / max(before["total_reads"], 1),
        "pct_q30_raw": 100 * before.get("q30_rate", 0),
        "gc_pct_trimmed": 100 * after.get("gc_content", 0),
    }


def flagstat_metrics(qc, s):
    txt = read_text(qc / f"{s}.flagstat.txt")
    def count(label):
        m = re.search(rf"^(\d+) \+ \d+ {label}\b", txt, re.M)
        return int(m.group(1)) if m else None
    total = count("primary") or count("in total")
    mapped = count("primary mapped") or count("mapped")
    dups = count("primary duplicates") or count("duplicates")
    out = {}
    if total:
        out["pct_mapped"] = 100 * mapped / total if mapped is not None else None
    if mapped:
        out["pct_duplicates"] = 100 * dups / mapped if dups is not None else None
    return out


def bismark_metrics(qc, s):
    out = {}
    rep = "".join(read_text(p) for p in qc.glob(f"{s}_*E_report.txt"))
    m = re.search(r"Mapping efficiency:\s*([\d.]+)%", rep)
    if m:
        out["pct_mapped"] = float(m.group(1))
    dd = "".join(read_text(p) for p in qc.glob("*deduplication_report.txt"))
    m = re.search(r"duplicated alignments removed:\s*\d+\s*\(([\d.]+)%\)", dd)
    if m:
        out["pct_duplicates"] = float(m.group(1))
    return out


def bedgraph_counts(path):
    """Sum methylated / unmethylated calls in a MethylDackel bedGraph."""
    meth = unmeth = 0
    if not path.exists():
        return None
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt") as fh:
        for line in fh:
            if line.startswith("track"):
                continue
            f = line.split("\t")
            meth += int(f[4]); unmeth += int(f[5])
    return meth, unmeth


def pct(counts):
    if not counts or sum(counts) == 0:
        return None
    return 100 * counts[0] / sum(counts)


def control_metrics(meth, s, unmeth_ctrl, meth_ctrl):
    out = {}
    if unmeth_ctrl:
        pre = f"{s}.control_{unmeth_ctrl}"
        cpg = bedgraph_counts(meth / f"{pre}_CpG.bedGraph")
        chg = bedgraph_counts(meth / f"{pre}_CHG.bedGraph")
        chh = bedgraph_counts(meth / f"{pre}_CHH.bedGraph")
        out["lambda_pct_meth_CpG"] = pct(cpg)
        non_cpg = None
        if chg or chh:
            non_cpg = tuple(sum(x) for x in zip(*(c for c in (chg, chh) if c)))
        out["lambda_pct_meth_nonCpG"] = pct(non_cpg)
        allc = [c for c in (cpg, chg, chh) if c]
        if allc:
            tot = tuple(sum(x) for x in zip(*allc))
            p = pct(tot)
            out["conversion_efficiency_pct"] = None if p is None else 100 - p
        out["lambda_CpG_calls"] = sum(cpg) if cpg else 0
    if meth_ctrl:
        cpg = bedgraph_counts(meth / f"{s}.control_{meth_ctrl}_CpG.bedGraph")
        out["puc19_pct_meth_CpG"] = pct(cpg)
        out["puc19_CpG_calls"] = sum(cpg) if cpg else 0
    return out


def genome_cpg_metrics(meth, s, controls):
    path = meth / f"{s}_CpG.bedGraph.gz"
    if not path.exists():
        path = meth / f"{s}_CpG.bedGraph"
    if not path.exists():
        return {}
    m = u = n = n5 = n10 = 0
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt") as fh:
        for line in fh:
            if line.startswith("track"):
                continue
            f = line.split("\t")
            if f[0] in controls:
                continue
            a, b = int(f[4]), int(f[5])
            d = a + b
            m += a; u += b; n += 1
            n5 += d >= 5; n10 += d >= 10
    return {
        "genome_pct_meth_CpG": 100 * m / (m + u) if m + u else None,
        "CpGs_covered": n,
        "CpGs_cov5": n5,
        "CpGs_cov10": n10,
        "mean_CpG_depth": (m + u) / n if n else None,
    }


def mbias_suggestion(mb, s):
    txt = read_text(mb / f"{s}.mbias_suggestion.txt")
    m = re.search(r"Suggested inclusion options:\s*(.*)", txt)
    return m.group(1).strip() if m else ""


def fmt(v):
    if v is None:
        return "NA"
    if isinstance(v, float):
        return f"{v:.2f}"
    return str(v)


def main():
    outdir = Path(sys.argv[1])
    cfg = dict(a.split("=", 1) for a in sys.argv[2:])
    unmeth_ctrl, meth_ctrl = cfg.get("CONTROL_UNMETH", ""), cfg.get("CONTROL_METH", "")
    controls = {c for c in (unmeth_ctrl, meth_ctrl) if c}
    controls |= {"phage_lambda", "plasmid_puc19c", "phage_T4", "phage_Xp12"}
    thr = {k: float(cfg[k]) for k in ("MAX_LAMBDA_CPG_METH", "MIN_PUC19_CPG_METH",
                                     "MIN_MAPPING_RATE", "MAX_DUP_RATE")}

    sheet = outdir / "samplesheet.tsv"
    samples = [l.split("\t")[0] for l in sheet.read_text().splitlines()[1:] if l.strip()]
    cols = ["sample", "status", "raw_reads", "reads_after_trim", "pct_reads_kept", "pct_q30_raw",
            "gc_pct_trimmed", "pct_mapped", "pct_duplicates", "conversion_efficiency_pct",
            "lambda_pct_meth_CpG", "lambda_pct_meth_nonCpG", "lambda_CpG_calls",
            "puc19_pct_meth_CpG", "puc19_CpG_calls", "genome_pct_meth_CpG", "CpGs_covered",
            "CpGs_cov5", "CpGs_cov10", "mean_CpG_depth", "flags", "mbias_suggestion"]
    rows = []
    for s in samples:
        d = outdir / "samples" / s
        r = {"sample": s}
        if not (d / ".done.extract").exists():
            r["status"] = "INCOMPLETE"
        r.update(fastp_metrics(d / "qc", s))
        r.update(flagstat_metrics(d / "qc", s))
        r.update(bismark_metrics(d / "qc", s))
        r.update(control_metrics(d / "meth", s, unmeth_ctrl, meth_ctrl))
        r.update(genome_cpg_metrics(d / "meth", s, controls))
        r["mbias_suggestion"] = mbias_suggestion(d / "mbias", s)

        flags = []
        def check(key, bad, msg):
            v = r.get(key)
            if v is not None and bad(v):
                flags.append(msg.format(v))
        check("lambda_pct_meth_CpG", lambda v: v > thr["MAX_LAMBDA_CPG_METH"],
              "lambda CpG meth {:.2f}% (incomplete conversion)")
        check("puc19_pct_meth_CpG", lambda v: v < thr["MIN_PUC19_CPG_METH"],
              "pUC19 CpG meth {:.1f}% (over-conversion / protection failure)")
        check("pct_mapped", lambda v: v < thr["MIN_MAPPING_RATE"], "mapping {:.1f}%")
        check("pct_duplicates", lambda v: v > thr["MAX_DUP_RATE"], "duplicates {:.1f}%")
        if unmeth_ctrl and not r.get("lambda_CpG_calls"):
            flags.append("no lambda reads (spike-in missing?)")
        if meth_ctrl and not r.get("puc19_CpG_calls"):
            flags.append("no pUC19 reads (spike-in missing?)")
        r.setdefault("status", "WARN" if flags else "PASS")
        r["flags"] = "; ".join(flags)
        rows.append(r)

    out = outdir / "summary" / "qc_summary.tsv"
    out.parent.mkdir(parents=True, exist_ok=True)
    with open(out, "w") as fh:
        fh.write("\t".join(cols) + "\n")
        for r in rows:
            fh.write("\t".join(fmt(r.get(c)) for c in cols) + "\n")

    show = [("sample", "sample"), ("status", "status"), ("raw_reads", "raw reads"),
            ("pct_mapped", "%map"), ("pct_duplicates", "%dup"),
            ("conversion_efficiency_pct", "conv%"), ("lambda_pct_meth_CpG", "lambda CpG%"),
            ("puc19_pct_meth_CpG", "pUC19 CpG%"), ("genome_pct_meth_CpG", "genome CpG%"),
            ("CpGs_cov10", "CpGs>=10x")]
    table = [[h for _, h in show]] + [[fmt(r.get(k)) for k, _ in show] for r in rows]
    widths = [max(len(row[i]) for row in table) for i in range(len(show))]
    for row in table:
        print("  ".join(v.ljust(w) for v, w in zip(row, widths)))
    for r in rows:
        if r["flags"]:
            print(f"  ! {r['sample']}: {r['flags']}")
    print(f"\nFull table: {out}")


if __name__ == "__main__":
    main()
