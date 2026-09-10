# The constant capacitance procedure of Goldberg (1985), applied to Li⁺ and Co²⁺ on dolomite

**A step-for-step transcription of the method, with this system's sites, reactions and parameters**

Marwa Elshebli, Javier Vilcáez, James Smay — Boone Pickens School of Geology, Oklahoma State University

---

## Purpose of this document

Goldberg (1985) is the clearest published statement of how the constant capacitance model is set
up, solved and fitted. This document follows that paper in its own order, with its own equation
numbering, substituting the dolomite surface for goethite and Li⁺ and Co²⁺ for phosphate,
selenite and silicate. Goldberg's equation numbers appear in square brackets throughout, so every
line can be checked against the source. Nothing here is a summary: every reaction, every mass
balance, every constant expression and every parameter is written out.

The reason for choosing that paper as the template goes beyond clarity. Its structure is this
project's structure. Goldberg fits surface complexation constants on **single-adsorbate** systems
with the acid–base constants and the capacitance held fixed, and then uses those constants, with
**no parameter fitted at all**, to predict a **competitive** two-adsorbate system. That is exactly
the path from the single-ion Li and Co batches to produced water.

**Reference.** Goldberg, S. (1985). Chemical modeling of anion competition on goethite using the
constant capacitance model. *Soil Science Society of America Journal, 49*(4), 851–856.

---

## 1. The four assumptions

Goldberg states four assumptions in his Data and Methods section. Three carry over unchanged. One
does not, and the place where the mimic breaks is more informative than the places where it holds.

| | Goldberg (1985), goethite | This system, dolomite |
|---|---|---|
| (i) | adsorption proceeds by a **ligand exchange** mechanism | adsorption proceeds by **proton exchange** at the carbonate site. The adsorbates here are cations, not anions, so the exchanged ion is H⁺ rather than OH⁻. The formalism is unchanged; only the reaction differs |
| (ii) | all surface complexes formed are **inner-sphere** | assumed here too, and it is precisely this assumption that Zachara et al. (1991) put in doubt for strongly hydrated divalent cations |
| (iii) | **no complexes are formed with ions in the background electrolyte** | **violated, and it cannot be waived.** The background is 0.684 M NaCl, in which CoCl⁺ carries about 20 % of the dissolved cobalt and NaCl ion pairing lowers the effective ionic strength. Goldberg could neglect this at 0.1 M; here an aqueous speciation stage must precede the surface calculation |
| (iv) | a **linear** relation holds between net surface charge and surface potential | assumed here too, Eq. [1] |

---

## 2. Charge and potential [Goldberg Eq. 1]

Goldberg writes the constant capacitance relation on a volume-of-suspension basis,

> σ = (C · S · a / F) · ψ  [1]

where σ is the net surface charge in mol_c m⁻³, C is the capacitance density in F m⁻², S is the
specific surface area in m² kg⁻¹, a is the concentration of solid in suspension in kg m⁻³, F is
the Faraday constant in C mol⁻¹, and ψ is the surface potential in V.

This is the same statement as the area-based form σ_area = C·ψ used elsewhere in this project,
multiplied by S·a/F to place charge on a per-volume basis. Both conventions are reported side by
side in the notebook, because the confusion between them, and between either of them and a
total-ion mass balance, is the origin of the recurring impression that the computed surface charge
is too large.

For this system C = √I / α with α = 0.004 (Pokrovsky et al., 1999, §3.3), giving C = 207 F m⁻² at
I = 0.684 M, S = 760 m² kg⁻¹, a = 60 kg m⁻³, and therefore C·S·a/F = 97.8 mol_c m⁻³ V⁻¹.

---

## 3. Surface functional groups and their acid–base reactions [Goldberg Eqs. 2–3]

Goethite carries a single amphoteric group, FeOH, and Goldberg needs two reactions for it.
Dolomite carries **three** primary hydration sites in 1 : 1 : 2 stoichiometry (Pokrovsky et al.,
1999), and only one of them binds metal cations. Goldberg's pair therefore becomes five reactions:

> \>CaOH⁰ + H⁺ ⇌ >CaOH₂⁺  [2a]

> \>MgOH⁰ + H⁺ ⇌ >MgOH₂⁺  [2b]

> \>CaOH⁰ ⇌ >CaO⁻ + H⁺  [3a]

> \>MgOH⁰ ⇌ >MgO⁻ + H⁺  [3b]

> \>CO₃H⁰ ⇌ >CO₃⁻ + H⁺  [3c]

Reactions [2a] to [3b] are Goldberg's [2] and [3] written once for each metal site. Reaction [3c]
has no counterpart in the goethite system: it is the deprotonation of the carbonate site, and the
deprotonated >CO₃⁻ is what binds the cation.

---

## 4. Adsorbate surface complexation reactions [Goldberg Eqs. 4–7]

Goldberg's Eqs. [4] to [7] exchange an anion for a surface hydroxyl. The cation analogue exchanges
a proton at the carbonate site. The two reactions this work sets out to quantify are

> \>CO₃H⁰ + Li⁺ ⇌ >CO₃Li⁰ + H⁺  [4]

> \>CO₃H⁰ + Co²⁺ ⇌ >CO₃Co⁺ + H⁺  [5]

Two further reactions of identical form are **not** adjustable, because the dolomite supplies
their cations by dissolving. They are the competitors for the same site, and Goldberg has no
equivalent because goethite does not dissolve:

> \>CO₃H⁰ + Ca²⁺ ⇌ >CO₃Ca⁺ + H⁺  [6]

> \>CO₃H⁰ + Mg²⁺ ⇌ >CO₃Mg⁺ + H⁺  [7]

The two hydroxyl sites form carbonate and bicarbonate complexes. These bind no metal, but two of
them carry charge and therefore enter the charge balance:

> \>MeOH⁰ + CO₃²⁻ + 2H⁺ ⇌ >MeHCO₃⁰ + H₂O  (Me = Ca, Mg)  [8]

> \>MeOH⁰ + CO₃²⁻ + H⁺ ⇌ >MeCO₃⁻ + H₂O  (Me = Ca, Mg)  [9]

---

## 5. Intrinsic conditional equilibrium constants [Goldberg Eqs. 8–13]

Goldberg writes one mass-action quotient for each reaction, in surface concentrations. Writing {}
for aqueous activity and [] for surface concentration in mol L⁻¹ of suspension, the quotients for
[4] to [7] are

> K_Li = [>CO₃Li⁰]·{H⁺} / ( [>CO₃H⁰]·{Li⁺} )  [10]

> K_Co = [>CO₃Co⁺]·{H⁺} / ( [>CO₃H⁰]·{Co²⁺} )  [11]

> K_Ca = [>CO₃Ca⁺]·{H⁺} / ( [>CO₃H⁰]·{Ca²⁺} )  [12]

> K_Mg = [>CO₃Mg⁺]·{H⁺} / ( [>CO₃H⁰]·{Mg²⁺} )  [13]

and correspondingly for [2a] to [3c], [8] and [9]. Each of these is an **intrinsic** constant only
after the electrostatic term of Eq. [16] has been applied.

---

## 6. Mass balance on the surface functional groups [Goldberg Eq. 14]

Goldberg needs one site balance. Dolomite needs three, one for each site:

> [>CO₃]_T = [>CO₃H⁰] + [>CO₃⁻] + [>CO₃Li⁰] + [>CO₃Co⁺] + [>CO₃Ca⁺] + [>CO₃Mg⁺]  [14a]

> [>CaOH]_T = [>CaOH⁰] + [>CaOH₂⁺] + [>CaO⁻] + [>CaHCO₃⁰] + [>CaCO₃⁻]  [14b]

> [>MgOH]_T = [>MgOH⁰] + [>MgOH₂⁺] + [>MgO⁻] + [>MgHCO₃⁰] + [>MgCO₃⁻]  [14c]

The totals follow from the crystallographic site densities of 14 µmol m⁻² for the carbonate site
and 7 µmol m⁻² for each metal site, multiplied by the surface area per litre of suspension.

---

## 7. Charge balance [Goldberg Eq. 15]

Goldberg's charge balance sums the charged surface species and nothing else. The same holds here,
and what the equation **excludes** matters as much as what it contains:

> σ = [>CaOH₂⁺] + [>MgOH₂⁺] + [>CO₃Co⁺] + [>CO₃Ca⁺] + [>CO₃Mg⁺]
> − [>CO₃⁻] − [>CaO⁻] − [>MgO⁻] − [>CaCO₃⁻] − [>MgCO₃⁻]  [15]

The species >CO₃Li⁰, >CaHCO₃⁰ and >MgHCO₃⁰ are neutral and contribute nothing. Dissolved Ca²⁺ and
Mg²⁺ released by dolomite dissolution appear nowhere in Eq. [15]. A surface charge computed from a
total-ion mass balance on the solution is a different quantity, it carries the dissolution signal,
and at acidic pH it exceeds the entire site capacity many times over. It is not the σ of Eq. [15]
and the two must never be substituted for one another.

---

## 8. Intrinsic versus conditional constants [Goldberg Eqs. 16–17]

Goldberg gives the conversion for the protonation constant, K₊(int) = ᶜK₊ exp(Fψ/RT). The general
form is

> K(int) = ᶜK · exp( ΔZ · F · ψ / R · T )  [16]

> ᶜK = the same mass-action quotient in concentrations, without the exponential  [17]

where ΔZ is the change in the charge of the surface species in the reaction as written, R is the
molar gas constant and T the absolute temperature. The value of ΔZ decides which constants feel
the electrostatics at all:

| Reaction | surface charge, before → after | ΔZ | electrostatic factor |
|---|---|---|---|
| [4] lithium | 0 → 0 | **0** | **none** |
| [5] cobalt | 0 → +1 | +1 | exp(+Fψ/RT) |
| [6], [7] calcium, magnesium | 0 → +1 | +1 | exp(+Fψ/RT) |
| [2a], [2b] protonation | 0 → +1 | +1 | exp(+Fψ/RT) |
| [3a], [3b], [3c] deprotonation | 0 → −1 | −1 | exp(−Fψ/RT) |
| [8] neutral bicarbonate complex | 0 → 0 | 0 | none |
| [9] carbonate complex | 0 → −1 | −1 | exp(−Fψ/RT) |

Because ΔZ = 0 for lithium, the lithium constant is electrostatically blind: its constant
capacitance and non-electrostatic values are identical to within rounding. This is a property of
the reaction stoichiometry and must not be read as evidence that the lithium constant is well
constrained.

---

## 9. Model parameters, in the layout of Goldberg's Table 1

| Parameter | Value | Units | Source |
|---|---|---|---|
| Capacitance density, C | 207 | F m⁻² | C = √I/α, α = 0.004 (Pokrovsky et al., 1999, §3.3) |
| Specific surface area, S | 760 | m² kg⁻¹ | BET, this study (0.76 m² g⁻¹) |
| Solid concentration, a | 60 | kg m⁻³ | 6 g per 100 mL |
| Surface area per litre, S·a | 45.6 | m² L⁻¹ | product of the two above |
| Ionic strength, I | 0.684 | mol L⁻¹ | 40 g L⁻¹ NaCl |
| C·S·a/F, the coefficient in Eq. [1] | 97.8 | mol_c m⁻³ V⁻¹ | computed |
| Maximum adsorption density [>CO₃]_T | 0.638 | mol m⁻³ | 14 µmol m⁻² × S·a |
| [>CaOH]_T | 0.319 | mol m⁻³ | 7 µmol m⁻² × S·a |
| [>MgOH]_T | 0.319 | mol m⁻³ | 7 µmol m⁻² × S·a |

**Fixed intrinsic constants**, all from Pokrovsky et al. (1999), Table 3, dolomite columns, held
constant throughout exactly as Goldberg holds his protonation constants and capacitance fixed:

| Reaction | log K(int) |
|---|---|
| [3c] >CO₃H⁰ ⇌ >CO₃⁻ + H⁺ | −4.8 ± 0.2 |
| [6] >CO₃H⁰ + Ca²⁺ ⇌ >CO₃Ca⁺ + H⁺ | −1.8 ± 0.2 |
| [7] >CO₃H⁰ + Mg²⁺ ⇌ >CO₃Mg⁺ + H⁺ | −2.0 ± 0.2 |
| [2a] >CaOH⁰ + H⁺ ⇌ >CaOH₂⁺ | +11.5 ± 0.2 |
| [2b] >MgOH⁰ + H⁺ ⇌ >MgOH₂⁺ | +10.6 ± 0.2 |
| [3a] >CaOH⁰ ⇌ >CaO⁻ + H⁺ | −12 ± 2 |
| [3b] >MgOH⁰ ⇌ >MgO⁻ + H⁺ | −12 ± 2 |
| [8] >CaOH⁰ + CO₃²⁻ + 2H⁺ ⇌ >CaHCO₃⁰ + H₂O | −4.0 ± 0.5 |
| [8] >MgOH⁰ + CO₃²⁻ + 2H⁺ ⇌ >MgHCO₃⁰ + H₂O | −3.5 ± 0.5 |
| [9] >CaOH⁰ + CO₃²⁻ + H⁺ ⇌ >CaCO₃⁻ + H₂O | +16.6 ± 0.2 |
| [9] >MgOH⁰ + CO₃²⁻ + H⁺ ⇌ >MgCO₃⁻ + H₂O | +15.4 ± 0.2 |

**Adjustable constants.** Goldberg fits three for phosphate and two each for selenite and silicate.
This work fits **one per metal**, reactions [4] and [5], because the data cannot support more.

---

## 10. Numerical solution and parameter estimation

Goldberg used FITEQL (Westall, 1982) to fit the surface complexation constants and MICROQL
(Westall, 1979) to run the competitive prediction. No software is used in this project, so both
are reimplemented in Python:

1. **Equilibrium solution.** Equations [14a] to [14c] and [15], together with [1] and the mass-action
   expressions [10] to [13], are solved simultaneously at each data point. The surface potential
   enters as its own unknown through the Boltzmann factor of Eq. [16], and the loop between σ and ψ
   is closed by Eq. [1].
2. **Fitting.** With the acid–base constants of §9 and the capacitance held fixed, the adsorbate
   constants of [4] and [5] are adjusted to minimise the FITEQL objective,
   WSOS/DF = Σ (Y_i/s_i)² / (N_obs − N_param), where Y_i are the mass-balance residuals and s_i
   their estimated errors. A value near unity means the model reproduces the data to within
   measurement error.
3. **Prediction.** With every constant fixed and nothing refitted, the same solver is run forward.
   This is Goldberg's competitive step, and it is the basis on which single-ion constants can be
   carried to produced water.

---

## 11. Where this mimic is exact, and where it is not

**Exact.** Assumptions (i), (ii) and (iv); Eqs. [1], [14], [15], [16] and [17]; acid–base constants
and capacitance held fixed with only the adsorbate constants adjustable; a single-adsorbate fit
followed by a competitive prediction in which no parameter is refitted.

**Not exact**, each for a stated reason:

1. **Background electrolyte.** Assumption (iii) fails at 0.684 M NaCl, so an aqueous speciation
   stage is inserted ahead of the surface calculation. Goldberg has no equivalent stage.
2. **The solid dissolves.** Goethite does not; dolomite does, and it supplies the competing cations
   of reactions [6] and [7]. Their surface fractions therefore vary from point to point with the
   measured calcium and magnesium, rather than being fixed.
3. **Capacitance.** Goldberg fixes C = 1.06 F m⁻², the γ-Al₂O₃ value of Westall and Hohl (1980).
   Carbonates use the ionic-strength form, which gives a far stiffer double layer and surface
   potentials of a few millivolts.
4. **Constraint.** Goldberg fits adsorption envelopes of many points spanning pH 3 to 12. This data
   set provides five usable lithium points and two usable cobalt points spanning 0.08 pH units. The
   procedure is identical; what it can deliver is not, and establishing that difference is the
   result of the present work rather than a shortcoming of the method.

---

## References

Goldberg, S. (1985). Chemical modeling of anion competition on goethite using the constant
capacitance model. *Soil Science Society of America Journal, 49*(4), 851–856.
https://doi.org/10.2136/sssaj1985.03615995004900040013x

Goldberg, S., & Sposito, G. (1984). A chemical model of phosphate adsorption by soils: I.
Reference oxide minerals. *Soil Science Society of America Journal, 48*(4), 772–778.

Pokrovsky, O. S., Schott, J., & Thomas, F. (1999). Dolomite surface speciation and reactivity in
aquatic systems. *Geochimica et Cosmochimica Acta, 63*(19–20), 3133–3143.

Westall, J. C. (1979). *MICROQL: A chemical equilibrium program in BASIC*. Swiss Federal Institute
for Water Resources and Water Pollution Control, Dübendorf.

Westall, J. C. (1982). *FITEQL: A computer program for determination of chemical equilibrium
constants from experimental data* (Report 82-01). Department of Chemistry, Oregon State University.

Westall, J. C., & Hohl, H. (1980). A comparison of electrostatic models for the oxide/solution
interface. *Advances in Colloid and Interface Science, 12*(4), 265–294.

Zachara, J. M., Cowan, C. E., & Resch, C. T. (1991). Sorption of divalent metals on calcite.
*Geochimica et Cosmochimica Acta, 55*(6), 1549–1562.
