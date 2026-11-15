# cf-EpiTracing — 96-well plate workflow (BioRender figure brief)

**H. Ludewick · Draft v0.1 · 15 November 2026**
Schematic of the plate-based cf-EpiTracing method (96-well + ThermoMixer C). **Draft — method not yet validated by sequencing; Tn5 dose to be finalised from the pilot.** Keep it descriptive (steps + hardware), no performance claims yet.

---

## Layout
Horizontal left-to-right flow, 8 numbered steps joined by arrows, all shown happening in **one row of a 96-well plate**. A hardware strip runs along the bottom. Use one consistent colour each for: beads (grey/brown), DNA (blue), antibody (green Y), Tn5 (orange), adapters (red/purple).

## Steps (left → right)

| # | Panel | Icon(s) | Label |
|---|---|---|---|
| 1 | Plasma input | plasma/blood tube → 96-well plate, one well highlighted | **200 µl plasma per well** (one donor, pooled) |
| 2 | Capture + S2 spike-in | Dynabead M-270 with a green antibody (Y) binding a nucleosome (histone octamer + wrapped DNA); small *Drosophila* icon for S2 | **Anti-H3K4me3 (or IgG) on Dynabeads M-270 + 1 ng S2 spike-in** · 5 °C overnight, ThermoMixer C |
| 3 | Wash | plate on ring magnet, beads pelleted, supernatant arrow out | **Magnetic wash** · 4× wash buffer + 1× Tris (plate magnet) |
| 4 | Tagmentation | loaded Tn5 transposome (orange Tn5 dimer + two red/purple ME adapters) inserting into the bead-bound chromatin DNA | **Tn5 tagmentation** (barcoded) · 37 °C 30 min, ThermoTop · *optimised dose (titration)* |
| 5 | Lysis / stop | SDS + proteinase K scissors releasing DNA from protein | **Lysis** (SDS / proteinase K) · 55 °C 15 min |
| 6 | Quench | PMSF icon / stop symbol | **Quench** (PMSF) · 37 °C 15 min · inactivates proteinase K |
| 7 | PCR | index primers (i5/i7) on fragment → amplified library; arrow to a PCR plate | **21-cycle dual-index PCR** (KAPA HiFi) · transfer to PCR plate, thermocycler |
| 8 | Clean-up → sequencing | SPRI beads, then a flow-cell / sequencer icon | **Double 1.2× SPRI → pool → 2×50 bp, ~8–10 M pairs/library** |

## Hardware strip (bottom row of icons, labelled)
- **ThermoMixer C + ThermoTop** — temperature (5–55 °C) + mixing + no-evaporation lid; runs steps 2–6
- **96-well 300 µl low-bind plate** — capture → quench
- **Ring-shaped plate magnet** — washes
- **VACUSAFE** — bulk liquid removal (pipette the final wash)
- **Thermocycler** — PCR only (PCR-96M2-HS-C plate)

## "What's new vs the tube method" callout box (optional side panel)
- 96-well plate format (scalable)
- Single-instrument capture→quench on the ThermoMixer C + ThermoTop
- Ring magnet + VACUSAFE liquid handling
- Same chemistry as Chen et al. (antibody-bead capture → free Tn5 tagmentation → direct PCR)

## BioRender assembly notes
- **Icon search terms:** "blood collection tube", "magnetic bead", "antibody", "nucleosome", "Tn5 / transposase" (or build from a dimer + DNA), "96 well plate", "heat block / thermomixer", "magnet", "PCR tube", "DNA fragment with adapters", "Illumina flow cell / sequencer".
- Number the steps 1–8; connect with arrows; keep the plate motif under every step so the "all in one plate" message reads at a glance.
- Footer: "Draft v0.1 — workflow schematic; performance not yet shown. Based on Chen et al., Nature 2026 (cf-EpiTracing)."
- Leave the Tn5 dose as a variable/placeholder until the pilot reads out, then replace with the chosen value.

## When to finalise
Add a results panel (e.g. complexity-vs-dose, tube-vs-plate equivalence) only **after** the pilot sequencing confirms the format works and the dose is set — then promote from "draft" to a method figure.
