# cf-EpiTracing — Plasma Tn5 Dose Pilot

**H. Ludewick · Version 1.0 · 9 November 2026**
Standalone bench worksheet (capture steps folded in from Bench Protocol v1.3). Tn5 dose titration to optimise library complexity. H3K4me3 + IgG, one donor.

---

## 1. Objective and rationale

- **Why:** the in-house Tn5 is confirmed active, correctly loaded and at the right concentration (Library Complexity Plan, B.10). Complexity is limited by tagmentation **dose**, not a broken enzyme.
- **Direction:** on naked DNA the enzyme under-cuts at low enzyme:DNA ratio; on the tiny cfDNA input (~1–5 ng / 200 µl plasma) 100 nM is a **high** ratio and may over-tagment (short, unmappable fragments). **Bracket both sides of 100 nM.**
- **Goal:** find the Tn5 dose that maximises unique mapped fragments without losing enrichment, then adopt it in the bench protocol.

## 2. Design

| Element | Choice |
|---|---|
| Plasma | One donor; ≥ 2 ml pooled from the same draw, mixed, split into 200 µl aliquots |
| Mark | H3K4me3 (most robust), single per dose — no replicate (conjugate-limited); IgG at the top dose |
| Doses (final, in 10 µl) | 12.5, 25, 50, 100, 200 nM |
| S2 spike-in | 1 ng S2 chromatin per 200 µl plasma (target dm6 5–20%) |
| Libraries | 6 total (5 H3K4me3 + 1 IgG) — limited by conjugated antibody |
| Sequencing | ~5–10 M read pairs each; compare at a matched subsample (e.g. 3 M pairs) |

> **Note:** singlicate is a deliberate trade-off for a dose-finding pilot — 5 single doses define the complexity-vs-dose trend better than fewer doses in duplicate, and one technical failure still leaves a 4-point curve. IgG sits at the top dose (200 nM): background rises with dose, so a clean 200 nM IgG means every lower dose is clean. Add replicates and more IgG doses in the confirmatory run.

### Reaction map (assign a unique dual index to each — see CSV)

| # | Library ID | Mark | Tn5 dose | Index (fill in) |
|---|---|---|---|---|
| 1 | PILOT_K4_12.5 | H3K4me3 | 12.5 nM | |
| 2 | PILOT_K4_25 | H3K4me3 | 25 nM | |
| 3 | PILOT_K4_50 | H3K4me3 | 50 nM | |
| 4 | PILOT_K4_100 | H3K4me3 | 100 nM | |
| 5 | PILOT_K4_200 | H3K4me3 | 200 nM | |
| 6 | PILOT_IgG_200 | IgG | 200 nM | |

## 2a. 96-well plate layout and equipment

All 6 reactions run in **one row** of a 96-well plate (300 µl wells). cf-EpiTracing is CUT&Tag-type chemistry, so your colleagues' plate magnet + thermomixer workflow transfers directly.

| Well | Library | Tn5 dose |
|---|---|---|
| A1 | K4_12.5 | 12.5 nM |
| A2 | K4_25 | 25 nM |
| A3 | K4_50 | 50 nM |
| A4 | K4_100 | 100 nM |
| A5 | K4_200 | 200 nM |
| A6 | IgG_200 | 200 nM |

- **Equipment:** low-bind 96-well plate, 300 µl wells; adhesive seal film; 96-well plate magnet; **Eppendorf ThermoMixer C with the ThermoTop heated lid** and the SmartBlock for your plate; thermocycler for PCR only. With the ThermoTop fitted, **every step up to PCR runs on the ThermoMixer C** — no evaporation of the small volumes.
- **Magnet check:** confirm the magnet pellets the 2.8 µm M-270 beads cleanly in your wells (test with blank beads). Allow 3–5 min for the first separation out of viscous plasma.
- **Capture (5.1):** 300 µl plate on the ThermoMixer C at **5 °C, overnight**, ThermoTop on, sealed. Mix to keep M-270 suspended: ~1,000–1,500 rpm, or intermittent (1,500 rpm 10 s every 2 min) — check beads don't settle.
- **Washes (5.2):** on the plate magnet, **200 µl per wash**; use **4× wash buffer + 1× 10 mM Tris** (one extra wash to offset the smaller volume). Resuspend off-magnet with a brief ThermoMixer pulse (~1,500 rpm) or by pipetting.
- **Tagmentation / lysis / quench (5.3–5.5):** same plate on the **ThermoMixer C with the ThermoTop** (37 / 55 / 37 °C), sealed, gentle mix (~500–800 rpm) or static. The ThermoTop prevents evaporation/condensation of the 10–17 µl volumes.
- **PCR (5.6):** transfer each **17 µl quenched bead slurry** to a 96-well PCR plate, add PCR mix, cycle in the thermocycler. (If your 300 µl plate is thermocycler-rated, add PCR mix and cycle in place — no transfer.)
- **Sealing:** use a **peelable adhesive plate film** (optically clear or pierceable foil) for the capture and enzymatic steps — peel and apply a **fresh** film each time you add a reagent (lysis, then quench); press flat with a sealing roller so it survives capture mixing and seats evenly under the ThermoTop. For **PCR** use a **thermocycling-rated seal** (PCR foil or clear PCR film rated ~100 °C, or heat-seal). Avoid Parafilm and non-temperature-rated storage films. Seal first, then fit the ThermoTop.
- **First run:** keep **1–2 reactions in tube format in parallel** as a format control.

## 3. Buffers and stocks (make fresh where noted)

| Buffer | Composition | Notes |
|---|---|---|
| Wash buffer (cold) | 50 mM Tris-HCl pH 7.4, 150 mM NaCl, 0.5% Triton X-100, 2 mM EDTA, 10 mM Na-butyrate | Keep on ice |
| Capture additives | cOmplete protease inhibitor (1×), 10 mM Na-butyrate | Add as a master mix to all wells |
| 2× tagmentation salts | 20 mM MgCl₂, 20 mM TAPS-NaOH pH 8.5 | 1:1 with 2× Tn5 → 10 mM MgCl₂, 10 mM TAPS, dose |
| Tn5 dilution buffer | 10 mM Tris-HCl pH 7.5, 50 mM NaCl, 0.1 mM EDTA, 0.1% Triton X-100, 50% glycerol | **Add DTT to 1 mM FRESH** (1 µl 1 M per 1 ml); discard at day end |
| Lysis buffer | 60 mM Tris-HCl pH 8.0, 60 mM EDTA, 0.25% SDS, **0.6 mg/ml proteinase K** | PK = 3 µl of 20 mg/ml per 100 µl; add PK last |
| Quench buffer | 60 mM MgCl₂, 0.36% Triton X-100, **4 mM PMSF** | **Make fresh**, PMSF last (hydrolyses in ~30 min at pH 8) |
| PMSF stock | **100 mM in ethanol or isopropanol** (17.4 mg/ml) | **TOXIC** — gloves, fume hood. Do NOT use a 10 mM stock (too much ethanol → PCR inhibition) |

**Lysis buffer (100 µl):** 6 µl 1 M Tris pH 8.0 + 12 µl 0.5 M EDTA + 2.5 µl 10% SDS + 3 µl 20 mg/ml proteinase K + water to 100 µl. Add 2 µl per reaction.

**Quench buffer (per 5 µl reaction; scale ×N, add PMSF last):** 0.3 µl 1 M MgCl₂ + 0.18 µl 10% Triton X-100 + 0.2 µl 100 mM PMSF + water to 5 µl. Add 5 µl per reaction.

## 4. Tn5 working stocks (2× the final dose)

Each reaction = 5 µl of 2× Tn5 working + 5 µl 2× tagmentation salts = 10 µl at the final dose. Make a 5 µM intermediate, then the four working stocks. Loaded Tn5 assumed ~37.5 µM (use your best-performing loaded stock); keep on ice, dilution buffer with FRESH 1 mM DTT.

| Step | Recipe | Result |
|---|---|---|
| 5 µM intermediate | 4 µl loaded Tn5 (37.5 µM) + 26 µl dilution buffer | ~5 µM |
| 25 nM working (→ 12.5 nM) | 50 µl of the 50 nM working + 50 µl dilution buffer | 25 nM |
| 50 nM working (→ 25 nM) | 1 µl of 5 µM + 99 µl dilution buffer | 50 nM |
| 100 nM working (→ 50 nM) | 2 µl of 5 µM + 98 µl | 100 nM |
| 200 nM working (→ 100 nM) | 4 µl of 5 µM + 96 µl | 200 nM |
| 400 nM working (→ 200 nM) | 8 µl of 5 µM + 92 µl | 400 nM |

## 5. Procedure

### 5.1 Capture — Day 1 (overnight)
1. Vortex the coupled bead stock (1 mg/ml). For each reaction take **50 µl (0.05 mg, ~1 µg antibody)**; place on magnet, remove storage buffer.
2. To each **200 µl plasma** aliquot add the beads, **1× cOmplete protease inhibitor**, **10 mM Na-butyrate**, and **1 ng S2 spike-in** (1 µl of 1 ng/µl). Prepare additives + S2 as one master mix for all wells.
3. Rotate **overnight at 4 °C**.

### 5.2 Washes — Day 2
4. Wash **3× with cold wash buffer** (allow 3–5 min on the magnet for the first separation — plasma is viscous).
5. Wash **1× with 10 mM Tris-HCl pH 7.4**. Remove all liquid. Keep beads on ice; proceed immediately.
   > No wash after tagmentation — specificity comes entirely from these pre-tagmentation washes.

### 5.3 Tagmentation (dose-specific)
6. For each reaction make 10 µl tagmentation mix: **5 µl 2× tagmentation salts + 5 µl of the matching 2× Tn5 working stock** (final 10 mM MgCl₂, 10 mM TAPS pH 8.5, dose nM Tn5).
7. Resuspend the washed beads in the 10 µl mix for that library's dose.
8. Incubate **37 °C, 30 min**.

### 5.4 Stop and digest
9. Add **2 µl lysis buffer** (0.6 mg/ml PK). Incubate **55 °C, 15 min**.

### 5.5 Quench (prepare fresh with 100 mM PMSF stock)
10. Add **5 µl quench buffer**. Incubate **37 °C, 15 min**. (Inactivates PK so it cannot degrade the PCR polymerase.)

### 5.6 PCR enrichment
11. Add PCR mix directly to the bead reaction: KAPA HiFi (GC), dNTPs, assigned i7 + i5 index primers (~50 µl total). **21 cycles for all libraries** — do not adjust per sample.

### 5.7 Clean-up (double 1.2× SPRI)
12. Remove the M-270 beads: magnet 2 min, transfer the ~50 µl supernatant to a fresh tube.
13. **Round 1:** 60 µl SPRI (1.2×), 5 min RT, magnet, 2× 80% ethanol wash, air-dry 2–3 min, elute in 40 µl 10 mM Tris pH 8.0.
14. **Round 2:** 48 µl SPRI (1.2×) on the 40 µl, repeat; final elute in 22 µl, recover 20 µl. Store −20 °C.

### 5.8 QC
15. TapeStation **D5000** (or HS D5000 if < 1 ng/µl); record yield (nM, 170–1,000 bp region) and fragment-size profile per library.

## 6. Read-outs (pipeline as for run 20261026)
- Library yield after 21 cycles (TapeStation nM) and fragment-size profile.
- Unique human fragments (MAPQ ≥ 30, duplicates removed) and duplication % at a matched subsample.
- TSS/flank, FRiP, % reads at TSS ± 1 kb; IgG background and peak count.
- Sequenced fragment-size mode (was ~55–65 bp; watch whether lower doses shift it up toward ~150 bp).
- dm6 (S2) fraction per library (target 5–20%).

## 7. Decision
- Pick the dose with the **most unique fragments** at matched depth, provided TSS/flank and FRiP fall by ≤ ~10% vs 100 nM and IgG stays clean.
- **Complexity peaks at a lower dose (25–50 nM)** with the fragment mode shifting up toward ~150 bp → 100 nM was over-tagmenting the cfDNA; adopt the lower dose.
- **Complexity rises with dose to 200 nM** with enrichment held → adopt the higher dose.
- **No dose improves unique fragments by ≥ 1.5×** → the limit is capture/recovery, not Tn5 dose; move to antibody/bead amount and wash losses.
- Adopt the chosen dose in Bench Protocol v1.3 and record the batch-specific Tn5 working dilution.

## 8. Safety and handling
- **PMSF is highly toxic** (irreversible serine inhibitor) — gloves, fume hood. 100 mM ethanol/isopropanol stock; make quench buffer fresh daily.
- **DTT:** add to dilution buffer fresh (1 mM) on the day; discard spiked buffer at day end.
- Never vortex Tn5; mix by gentle pipetting; keep on ice/cold block.
- Load Tn5 fresh from the premix in single-use aliquots; the loaded transposome (~37% glycerol) is perishable.
