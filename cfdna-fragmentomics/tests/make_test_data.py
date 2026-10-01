"""Simulate a small reference + paired-end BAM with cfDNA-like fragment sizes, for smoke tests."""
import sys
import numpy as np
import pysam

def main(outdir, n_frag=60000, short_shift=0, seed=1, name="sim"):
    chroms = {"chr1": 3_000_000, "chr2": 2_000_000}
    ref = f"{outdir}/ref.fa"
    ref_rng = np.random.default_rng(0)  # fixed: every sample shares one reference
    seqs = {}
    for c, L in chroms.items():
        gc = np.repeat(ref_rng.uniform(0.35, 0.6, L // 100000), 100000)
        s = np.where(ref_rng.random(L) < gc, ref_rng.choice(list("GC"), L), ref_rng.choice(list("AT"), L))
        seqs[c] = "".join(s)
    with open(ref, "w") as fh:
        for c, s in seqs.items():
            fh.write(f">{c}\n" + "\n".join(s[i:i+60] for i in range(0, len(s), 60)) + "\n")
    pysam.faidx(ref)
    rng = np.random.default_rng(seed)
    header = {"HD": {"VN": "1.6", "SO": "coordinate"},
              "SQ": [{"SN": c, "LN": L} for c, L in chroms.items()]}
    reads = []
    for i in range(n_frag):
        c = "chr1" if rng.random() < 0.6 else "chr2"
        mono = rng.random() < 0.9
        L = int(rng.normal(167 - short_shift, 15) if mono else rng.normal(334, 25))
        L = max(60, L)
        st = int(rng.integers(0, chroms[c] - L))
        rl = min(100, L)
        tid = list(chroms).index(c)
        for k, (pos, rev) in enumerate([(st, False), (st + L - rl, True)]):
            a = pysam.AlignedSegment()
            a.query_name = f"f{i}"; a.reference_id = tid; a.reference_start = pos
            a.next_reference_id = tid; a.next_reference_start = st if k else st + L - rl
            a.query_sequence = seqs[c][pos:pos + rl]; a.query_qualities = [30] * rl
            a.cigarstring = f"{rl}M"; a.mapping_quality = 60
            a.flag = 1 | 2 | (16 if rev else 32) | (64 if k == 0 else 128)
            a.template_length = L if k == 0 else -L
            a.set_tag("MQ", 60)
            reads.append(a)
    reads.sort(key=lambda a: (a.reference_id, a.reference_start))
    bam = f"{outdir}/{name}.bam"
    with pysam.AlignmentFile(bam, "wb", header=header) as out:
        for a in reads:
            out.write(a)
    pysam.index(bam)
    with open(f"{outdir}/blacklist.bed", "w") as fh:
        fh.write("chr1\t1000000\t1100000\n")

if __name__ == "__main__":
    main(sys.argv[1], short_shift=int(sys.argv[2]) if len(sys.argv) > 2 else 0,
         name=sys.argv[3] if len(sys.argv) > 3 else "sim", seed=int(sys.argv[4]) if len(sys.argv) > 4 else 1)
