"""extract keeps forward-reverse pairs whether or not the aligner flagged them as proper.

bwa sets the proper-pair flag from each library's insert-size spread, so long plasma fragments
often lack it. Here two long pairs lack the flag: they are kept by default and dropped with
--proper-only. Reverse-forward pairs and pairs split across chromosomes are never kept.
"""
import argparse
import gzip
import os
import sys
import tempfile

import numpy as np
import pandas as pd
import pysam

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "fragmentomics"))
import cfdna_frag as cf  # noqa: E402

READ = 50


def pair(header, name, chrom, start, length, proper, rf=False, mate_chrom=None):
    """Two reads of one fragment [start, start + length); rf=True flips the orientation."""
    tid = header.get_tid(chrom)
    mtid = header.get_tid(mate_chrom) if mate_chrom else tid
    r1, r2 = pysam.AlignedSegment(header), pysam.AlignedSegment(header)
    p1, p2 = start, start + length - READ
    for r, pos, mpos, rev, read1 in ((r1, p1, p2, rf, True), (r2, p2, p1, not rf, False)):
        r.query_name, r.query_sequence = name, "A" * READ
        r.query_qualities = pysam.qualitystring_to_array("I" * READ)
        r.reference_id = tid if r is r1 else mtid
        r.next_reference_id = mtid if r is r1 else tid
        r.reference_start, r.next_reference_start = pos, mpos
        r.cigarstring, r.mapping_quality = f"{READ}M", 60
        r.is_paired, r.is_proper_pair = True, proper
        r.is_read1, r.is_read2 = read1, not read1
        r.is_reverse, r.mate_is_reverse = rev, not rev
        r.template_length = (length if r is r1 else -length) if r.reference_id == r.next_reference_id else 0
        r.set_tag("MQ", 60)
    return [r1, r2]


def main():
    d = tempfile.mkdtemp()
    rng = np.random.default_rng(0)
    seq = "".join(rng.choice(list("ACGT"), 20_000))
    with open(f"{d}/ref.fa", "w") as fh:
        fh.write(">chr1\n" + seq + "\n>chr2\n" + seq + "\n")
    pysam.faidx(f"{d}/ref.fa")
    header = pysam.AlignmentHeader.from_dict({"HD": {"VN": "1.6", "SO": "coordinate"},
                                              "SQ": [{"SN": "chr1", "LN": 20_000}, {"SN": "chr2", "LN": 20_000}]})
    reads = (pair(header, "mono", "chr1", 1000, 167, True)
             + pair(header, "short", "chr1", 3000, 50, True)          # both reads start together
             + pair(header, "long1", "chr1", 5000, 520, False)
             + pair(header, "long2", "chr1", 8000, 880, False)
             + pair(header, "rf", "chr1", 11000, 300, False, rf=True)
             + pair(header, "split", "chr1", 14000, 200, False, mate_chrom="chr2"))
    reads.sort(key=lambda r: (r.reference_id, r.reference_start))
    with pysam.AlignmentFile(f"{d}/t.unsorted.bam", "wb", header=header) as out:
        for r in reads:
            out.write(r)
    pysam.sort("-o", f"{d}/t.bam", f"{d}/t.unsorted.bam")
    pysam.index(f"{d}/t.bam")

    def run(proper_only):
        args = argparse.Namespace(bam=f"{d}/t.bam", reference=f"{d}/ref.fa", out=f"{d}/f.tsv.gz",
                                  blacklist=None, min_mapq=30, min_len=20, chroms=r"^chr[12]$",
                                  proper_only=proper_only)
        cf.extract(args)
        with gzip.open(args.out, "rt") as fh:
            return sorted(pd.read_csv(fh, sep="\t").length)

    default, proper = run(False), run(True)
    print("default:", default, " proper-only:", proper)
    assert default == [50, 167, 520, 880], default
    assert proper == [50, 167], proper
    print("OK")


if __name__ == "__main__":
    main()
