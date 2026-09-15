# Intrinsic surface complexation constants for Li⁺ and Co²⁺ on dolomite

Determination of intrinsic stability constants for lithium and cobalt binding on dolomite in
NaCl brine, by constant capacitance surface complexation modelling. All code is Python; no
external software is used, and no data is read from a file.

## The workflow this follows

| Source | What it supplies |
|---|---|
| **Belova et al. (2014)**, *Nickel adsorption on chalk and calcite*, J. Contam. Hydrol. 170, 1–9 | **The workflow.** Capacity and saturation screening before any fitting, FITEQL inversion, electrostatic and non-electrostatic constants side by side, the exchange constant as the transferable quantity |
| **Pokrovsky, Schott & Thomas (1999)**, GCA 63, 3133–3143 | **The surface model.** Three dolomite sites in 1:1:2 stoichiometry, their densities, the eleven fixed constants, α |
| **Goldberg (1985)**, SSSAJ 49, 851–856 | **The formalism.** Assumptions, σ–ψ, site and charge balances, intrinsic vs conditional |
| **Westall (1982)**, FITEQL 2.0 | **The objective function**, WSOS/DF |

Belova et al. (2014) is the paper to check this against.

## Files

| File | What it is |
|---|---|
| **`LogK_dolomite_SCM.ipynb`** | **Start here.** The complete analysis, executed, in 11 numbered parts. Generated from `logK_dolomite_SCM.py`. |
| `logK_dolomite_SCM.py` | The same thing as a cell-marked Python file. This is the one to edit; the notebook is regenerated from it. |
| `Paper_LiCo_dolomite_carbonate_site_SCM.md` / `.docx` | Manuscript: every reaction, equation, parameter, result and reference. |
| `Goldberg1985_Procedure_Dolomite.ipynb` / `.md` / `.docx` | Side explainer: the Goldberg (1985) constant capacitance procedure transcribed step for step onto this system. Method only, no new results. |
| `figures/` | Four figures, 300 dpi PNG and vector PDF. |
| `superseded/` | Earlier attempts, kept only as a record. Do not use. |

## The headline result

**Neither constant is determined by this data set**, and the two fail for different reasons.

Γ is the difference of two concentrations each measured to ±3 %. At the best of twelve points
the uptake is 6.1 mg/L against 3.3 mg/L of noise — a signal-to-noise of 1.9, where a constant
needs about 3. Relative error on Γ runs from 54 % to 2968 %, median 220 %. **Zero points** clear
a 30 % bar, for either metal.

Two further checks confirm it. A Monte Carlo test recovers a *known* constant with bias within
±0.06 whenever uptake exceeds ~2 %, so the analysis code is not the limitation. And the two
indicative values land within 0.04 log units of each other — a monovalent alkali and a divalent
transition metal binding the same site with the same strength is not chemically credible, which
is what an inversion pinned by noise looks like.

The capacity test does settle one thing from the earlier study: removal of 70 mg/L Co or 8.5 % of
138 mg/L Li exceeds the carbonate-site density by 1.7 to 2.4 times, so those recoveries are
dominated by precipitation, not adsorption.

## What would determine the constants

Part 10 computes this rather than asserting it. **Surface area per litre is the controlling
variable for cobalt** — dose barely matters — and **dose is the controlling variable for lithium**,
which is currently biased low by 0.45 log units because the surface holds only 3 % of the lithium
present.

Target: **300 m² L⁻¹** (this study has 46), Co 20 mg/L, Li 14 mg/L, a 6 to 8 point pH edge at
fixed pH, an isotherm at one pH, solid-free blanks, alkalinity at every sampling, the solid
pre-equilibrated for a week, and a desorption step at the end. That reaches ±0.02 to ±0.04.

## Run it

```bash
pip install numpy pandas scipy matplotlib
python logK_dolomite_SCM.py          # or open LogK_dolomite_SCM.ipynb and run all
```
