# Capture panel draft

`build_capture_panel.py` drafts a targeted methylation capture panel from the human cell-type methylation atlas of Loyfer et al. 2023 (*Nature* 613:355), using the UXM U250 marker set. The panel is aimed at measuring lung-derived DNA in plasma and BAL from lung transplant recipients.

| Tier | Cell types | Markers per type | How markers are chosen |
|---|---|---|---|
| lung | alveolar and bronchial epithelium | 250 | Ranked by the contrast with the highest unmethylated fraction in any plasma/BAL background cell type (erythroid progenitors excluded; their value is reported as `eryth_U`) |
| background | neutrophils, monocytes/macrophages, T, B, NK, megakaryocytes, erythroid, endothelium, liver | 100 | Largest methylation difference |
| stromal | fibroblasts, smooth and skeletal muscle, cardiomyocytes, adipocytes | 40 | Largest methylation difference |
| other | all remaining atlas cell types | 25 | Largest methylation difference |

The stromal and "other" tiers are there so that deconvolution doesn't force unrelated signal into the main compartments. Markers need at least 4 CpGs, which is what UXM counts. Capture targets are padded ±60 bp and merged.

```bash
# atlas files come with UXM (cloned by downstream/setup_downstream.sh)
U=$OUTDIR/tools/UXM_deconv/supplemental
python3 panel/build_capture_panel.py $U/Markers.U250.hg38.tsv $U/Atlas.U250.l4.hg38.full.tsv panel_out/
```

Outputs:
- `panel_markers.tsv`: selected markers with statistics
- `panel_targets.bed`: hg38 capture targets for the probe vendor
- `panel_summary.tsv`: per cell type, including how many markers have ≤ 2 % / ≤ 5 % background

This is a draft for discussion with a probe vendor and a bioinformatician. Before ordering:
- Add probes for the conversion controls (lambda, pUC19).
- Consider adding SNP sites if donor/recipient genotyping is wanted.
- Check probe design in repeats with the vendor.
