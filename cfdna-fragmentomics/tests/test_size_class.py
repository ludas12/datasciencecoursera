"""Size-class end motifs: same-process long fragments sit at the noise floor; random ends do not."""
import itertools
import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "fragmentomics"))
import cfdna_frag as cf  # noqa: E402

rng = np.random.default_rng(1)
motifs = np.array(["".join(m) for m in itertools.product("ACGT", repeat=4)])
pref = rng.dirichlet(np.full(256, 3.0))          # a skewed, nuclease-like end profile


def frags(n, lo, hi, p):
    m5 = rng.choice(motifs, n, p=p)
    m3 = rng.choice(motifs, n, p=p)
    return pd.DataFrame({"length": rng.integers(lo, hi + 1, n), "motif_5p": m5, "motif_3p": m3})


same = pd.concat([frags(200_000, 150, 220, pref), frags(20_000, 401, 600, pref)])
out = cf.size_class_motifs(same, n_boot=100)
r = out.loc["401_600"]
assert r.js_p_vs_null > 0.05, r                       # indistinguishable from counting noise
assert abs(r.js_excess) < r.js_null_95

rand = pd.concat([frags(200_000, 150, 220, pref), frags(20_000, 401, 600, np.full(256, 1 / 256))])
r2 = cf.size_class_motifs(rand, n_boot=100).loc["401_600"]
assert r2.js_p_vs_null < 0.05 and r2.js_excess > 10 * r.js_null_median, r2
print("size-class motifs: same-process", round(r.js_vs_mono, 5), "noise", round(r.js_null_median, 5),
      "| random ends", round(r2.js_vs_mono, 5), "OK")
