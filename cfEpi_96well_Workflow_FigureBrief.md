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

---

## BioRender asset map (templates + icon packs)

Editing templates/icons in the BioRender editor is free (no AI credits). Open → duplicate → adapt.

### Icon packs (toolbox — grab individual icons)
- **Heating Instruments** → ThermoMixer C / heat block — https://app.biorender.com/biorender-templates/details/t-649608880c342739f6a86036?source=mcp
- **Magnets in Chemistry** → plate / ring magnet — https://app.biorender.com/biorender-templates/details/t-640e70ab386eba8e0f74d815?source=mcp
- **Tubes** → plasma tube, PCR tube — https://app.biorender.com/biorender-templates/details/t-633f362403394e304035dc31?source=mcp
- **Cell Culture & Plates (Editable Icons)** → 96-well plate — https://app.biorender.com/biorender-templates/details/t-608c4592da41a600a3a682c8?source=mcp
- **Lab and Machines** → thermocycler, vacuum aspirator, sequencer — https://app.biorender.com/biorender-templates/details/t-5fb581cc00cf4e00b0cefc84?source=mcp
- **Primer** → i5/i7 index primers — https://app.biorender.com/biorender-templates/details/t-63a38430aebcba548dcc924f?source=mcp
- **Epigenetics** → nucleosome / histone marks — https://app.biorender.com/biorender-templates/details/t-6064c3233d2ca1009f4b0f1f?source=mcp

### Element → source
| Figure element | Best source | Template ID |
|---|---|---|
| Panel/step scaffold | Blank Panels, 4 Step | 6175caa057e93400aeaf2637 |
| Panel/step scaffold (alt) | qPCR 4 Steps | 63a26af456a41b7e85601c1c |
| 1. Plasma input | Tubes pack + Liquid Biopsy | 5e1f35b66bca230087062b15 |
| 2. Bead capture + S2 | CUT&RUN Procedure | 5f872753c2ef9200ab24e4c3 |
| 2. (bead IP detail) | Immunoprecipitation With Magnetic Beads | 67acf95b4937b4a7c2099a5c |
| 3. Magnetic wash | Magnets pack + plate icons | 640e70ab386eba8e0f74d815 |
| 4. Tn5 tagmentation (base) | ATAC Sequencing | 6009a62628bf69021fdcdc04 |
| 4. (transposase mechanism) | Mechanism of Transposable Mobilization | 5fff153e7fd8e600a4dac80f |
| 5/6. Lysis / Quench | Tubes pack + STOP/scissors | 633f362403394e304035dc31 |
| 7. Dual-index PCR | Polymerase Chain Reaction (PCR) | 5cf004269cf1ff003377d8cd |
| 8. SPRI → cleanup | DNA Extraction | 5e68048a10234378502a316b |
| 8. Flow cell / sequencing | Next Generation Sequencing (Illumina) | 5ef134a11c72b100ad8d13ac |
| Sequencing workflow (alt) | Next Generation Sequencing Workflow | 62ccc36459dc10d0850c42df |
| ChIP framing (alt base) | ChIP sequencing | 5f20a148563a0600ad9d9f87 |
| 96-well plate motif | 96-well Plate with Callout | 6175c62132211700aefc962c |

### Recommended build
1. Open **ATAC Sequencing** (6009a62628bf69021fdcdc04) as the base — keeps the Tn5 tagmentation visuals.
2. Add antibody-bead capture from **CUT&RUN Procedure** (5f872753c2ef9200ab24e4c3) as steps 2–3.
3. Lay the 8 steps on the **Blank 4-Step panel** scaffold (6175caa057e93400aeaf2637), duplicated for 8, over the **96-well plate** motif.
4. Pull nucleosome/histone icons from the **Epigenetics** pack; equipment icons from **Heating Instruments / Magnets / Lab and Machines**; flow cell from **NGS (Illumina)**.
5. Relabel the 8 steps from the brief above; keep the Draft / dose-TBD framing.

### AI option
To generate the custom figure in one shot, top up 5 BioRender AI credits (https://app.biorender.com/account/ai-credits), then it is a one-step generation from the prepared prompt.
