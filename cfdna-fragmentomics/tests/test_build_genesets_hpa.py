"""Builder test on a mock HPA table (format as in proteinatlas.tsv: 'tissue: nTPM;...' cells)."""
import os
import subprocess
import sys
import tempfile
import zipfile

import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
TOOL = os.path.join(HERE, "..", "scripts", "build_genesets_hpa.py")


def main():
    rows = []
    def add(prefix, n, **kw):
        for i in range(n):
            r = {"Gene": f"{prefix}{i}", "Ensembl": "ENSG0", "RNA tissue specificity": "Low tissue specificity",
                 "RNA tissue distribution": "Detected in many", "RNA tissue specific nTPM": "",
                 "RNA single cell type specificity": "Low cell type specificity",
                 "RNA single cell type specific nTPM": "",
                 "RNA blood cell specificity": "Low immune cell specificity", "RNA blood cell specific nTPM": ""}
            r.update(kw)
            rows.append(r)
    add("AT2_", 30, **{"RNA single cell type specificity": "Cell type enriched",
                       "RNA single cell type specific nTPM": "alveolar cells type 2: 900.1"})
    add("CLUB_", 25, **{"RNA single cell type specificity": "Group enriched",
                        "RNA single cell type specific nTPM": "club cells: 300;ciliated cells: 200"})
    # lung single-cell gene that is ALSO neutrophil-specific -> excluded from lung_epithelium
    add("LUNGNEU_", 5, **{"RNA single cell type specificity": "Cell type enriched",
                          "RNA single cell type specific nTPM": "alveolar cells type 1: 50",
                          "RNA blood cell specificity": "Immune cell enriched",
                          "RNA blood cell specific nTPM": "neutrophil: 400"})
    add("LUNGT_", 25, **{"RNA tissue specificity": "Tissue enriched", "RNA tissue specific nTPM": "lung: 500"})
    add("NEU_", 30, **{"RNA blood cell specificity": "Immune cell enriched", "RNA blood cell specific nTPM": "neutrophil: 800"})
    add("MONO_", 25, **{"RNA blood cell specificity": "Group enriched",
                        "RNA blood cell specific nTPM": "classical monocyte: 300;non-classical monocyte: 200"})
    add("TCELL_", 25, **{"RNA blood cell specificity": "Immune cell enriched", "RNA blood cell specific nTPM": "naive CD4 T-cell: 100"})
    add("LIV_", 25, **{"RNA tissue specificity": "Tissue enriched", "RNA tissue specific nTPM": "liver: 2000"})
    add("HK_", 60, **{"RNA tissue distribution": "Detected in all"})
    add("TES_", 40, **{"RNA tissue specificity": "Tissue enriched", "RNA tissue specific nTPM": "testis: 50"})
    add("TESGRP_", 10, **{"RNA tissue specificity": "Group enriched", "RNA tissue specific nTPM": "testis: 50;lung: 30"})
    add("ENH_", 10, **{"RNA tissue specificity": "Tissue enhanced", "RNA tissue specific nTPM": "liver: 20"})

    d = tempfile.mkdtemp()
    tsv = os.path.join(d, "proteinatlas.tsv")
    pd.DataFrame(rows).to_csv(tsv, sep="\t", index=False)
    with zipfile.ZipFile(tsv + ".zip", "w") as z:
        z.write(tsv, "proteinatlas.tsv")
    out = os.path.join(d, "sets.tsv")
    subprocess.run([sys.executable, TOOL, "--hpa", tsv + ".zip", "-o", out, "--n-housekeeping", "50"], check=True)
    s = pd.read_csv(out, sep="\t", comment="#")
    got = s.groupby("set").gene.apply(set)
    pref = lambda st: {g.split("_")[0] for g in got[st]}  # noqa: E731
    assert pref("lung_epithelium") == {"AT2", "CLUB"}, pref("lung_epithelium")   # LUNGNEU excluded
    assert pref("lung_tissue") == {"LUNGT", "TESGRP"}   # group-enriched lung+testis is lung-enriched
    assert pref("neutrophil") == {"NEU", "LUNGNEU"}
    assert pref("monocyte_macrophage") == {"MONO"}
    assert pref("lymphocyte") == {"TCELL"}
    assert pref("liver") == {"LIV"}                       # 'Tissue enhanced' not included
    assert pref("inactive_control") == {"TES"}            # group-enriched testis excluded
    assert len(got["housekeeping"]) == 50 and pref("housekeeping") == {"HK"}
    assert len(got["megakaryocyte_erythroid"]) > 5
    print("OK")


if __name__ == "__main__":
    main()
