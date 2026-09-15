# Constant capacitance model constants for lithium and cobalt on dolomite, fitted simultaneously to an adsorption envelope and an isotherm

**Marwa Elshebli**¹, Javier Vilcáez¹, James Smay²

¹ Boone Pickens School of Geology, Oklahoma State University, Stillwater, OK 74078, USA
² School of Materials, Mechatronics & Manufacturing Engineering, Oklahoma State University, Tulsa, OK 74106, USA

*All numbers in the Results were produced by the accompanying notebook, `CCM_FITEQL_LiCo_dolomite.ipynb`. Every datum is embedded in that notebook and nothing is read from an external file.*

---

## Abstract

Dolomite recovers Li⁺ and Co²⁺ from petroleum produced water, and the quantitative split between sorption and carbonate mineral formation was left open by Elshebli et al. (2025). Here that split is resolved and intrinsic surface complexation constants are reported for both metals, using the constant capacitance model under the parameter selection discipline of Hayes et al. (1991), the simultaneous isotherm and envelope fitting protocol of Goldberg (2004), and the per point reporting of Ioannou and Dimirkou (1997). Two experiments are combined for the first time: the single ion batches, in which dolomite dissolution carries the suspension from pH 2.25 to 8.26 at fixed total metal and which therefore constitute an adsorption envelope, and the produced water batches, which probe a different total metal against 4186 mg L⁻¹ of competing calcium and therefore constitute an isotherm point. Strontium serves as an internal standard for the produced water: because it does not sorb on carbonate surfaces, strontianite equilibrium fixes the carbonate activity at log a(CO₃²⁻) = −5.85 ± 0.01 across every sampling and both batches, which in turn fixes the solubility floor of every other metal and separates the precipitated from the sorbed fraction without the alkalinity measurement the experiments did not record. Three parameter choices are shown to matter very differently. The capacitance is immaterial anywhere inside the bracket that Hayes Eq. 22 allows, 0.12 to 1.92 F m⁻²: the fitted constants move by 0.008 log units for cobalt and 0.005 for lithium across that whole range. The site total is the dominant sensitivity, moving the cobalt constant by 1.54 log units across three published site densities, so every constant here is quoted with the site total it was fitted at. The activity convention is worth 0.25 log units for cobalt and 0.16 for lithium at this ionic strength. At C = 1.06 F m⁻² and N_t = 0.638 mmol L⁻¹, the fitted constants are log K(>CO₃Co⁺) = −1.60 with WSOS/DF = 0.39 over eleven points, and log K(>CO₃Li⁰) = −1.69 with WSOS/DF = 0.14 over eight points. The cobalt constant is supported independently: inverting one constant from each envelope point and averaging in the manner of Ioannou and Dimirkou gives −1.46 ± 0.27 across pH 7.40 to 8.09, which agrees with the simultaneous fit to 0.14 log units, and the value falls 0.20 log units from the calcium constant Pokrovsky et al. (1999a) measured for the same site on the same mineral, the order that the Irving and Williams series predicts and that nothing in the fit constrained. The lithium constant is weaker and is reported as a lower bound, because lithium uptake reaches 0.91 of a monolayer in the single ion batches and 1.48 times the site inventory in produced water, which no surface complexation constant can reproduce; the excess is surface zabuyelite, consistent with the XRD of Elshebli et al. (2025) and with the bulk solution being four orders of magnitude undersaturated with respect to that phase. For cobalt in produced water the surface model predicts less than one per cent coverage because calcium holds the exchange sites, which is the same conclusion the solubility analysis reaches by an independent route: cobalt recovery from produced water is carbonate mineral formation, not sorption.

---

## 1. Introduction

Surface complexation constants are the transferable form of a sorption measurement. A distribution coefficient describes one experiment; an intrinsic constant, paired with a site density and a capacitance, can be carried into a speciation or transport calculation for a different water. That is why Elshebli et al. (2025), having established that dolomite recovers lithium and cobalt from produced water and that the two metals behave differently, left the quantitative mechanism split as the open question. This work answers it and reports the constants.

The obstacle is that neither of the two available data sets can carry a constant on its own. The single ion batches have the right mechanism, with uptake well inside the monolayer bound, but the uptake is small: at ±3 per cent analytical precision on each concentration, the best cobalt point reaches a signal to noise ratio of only 1.9, because the measured quantity is a difference between two large and nearly equal numbers. The produced water batches have twenty times the signal but most of the removal is carbonate precipitation rather than sorption, and separating the two requires a carbonate activity that was never measured.

Three published studies supply the way through, and the method here is assembled from them rather than invented.

Hayes, Redden, Ela and Leckie (1991) evaluated how surface complexation model parameters should be estimated from data, using FITEQL and titration data for three oxides. Three of their conclusions govern the present work. First, the constant capacitance model is the appropriate choice above about 0.01 M ionic strength, where the diffuse layer model overpredicts the interfacial potential; both experiments here run far above that, at 0.68 and 1.0 M. Second, the capacitance cannot be measured and must be treated as a fitting parameter bounded by theory, and their Eq. 22 brackets it at 0.1 to 2.0 F m⁻² for any physically possible dielectric constant and ion approach distance. Third, and most consequential, the fitted constants move directly with the assumed total site concentration, so a constant reported without the site total it was fitted at cannot be transferred, because whoever transfers it will use a different one.

Goldberg (2004) supplies the fitting protocol. Her advance over previous constant capacitance applications was to fit adsorption isotherms and adsorption envelopes simultaneously with a single constant set, rather than fitting one form and testing against the other. Her reason applies directly here: an isotherm constrains the magnitude of the affinity and an envelope constrains its pH dependence, and sparse data in either form leaves one of the two unconstrained.

Ioannou and Dimirkou (1997) supply the reporting discipline. Modelling phosphate adsorption on hematite and kaolinite, they took the surface protonation constants from a published compilation rather than optimising them, held the capacitance fixed at 1.06 F m⁻², fitted one complexation constant at each of six pH values, and reported the mean and standard deviation across pH. The spread across pH is the honest measure of whether a constant is intrinsic, because a constant that drifts systematically with pH is not one.

---

## 2. Materials and methods

### 2.1 Experiments

Dolomite at 60 g L⁻¹, six grams per 100 mL, BET specific surface area 0.76 m² g⁻¹, giving 45.6 m² per litre of suspension. Two experiments at 25 °C, both sampled at days 0, 2, 4 and 6, both analysed by ICP-OES with a relative precision of ±3 per cent on each concentration.

The **envelope arm** is the single ion series: cobalt at 80.5 mg L⁻¹ or lithium at 138.4 mg L⁻¹ in 40 g L⁻¹ NaCl, ionic strength 0.684 M, in batches started at nominal pH 2 and pH 6. Those labels describe the starting condition only. Dolomite dissolution buffers the suspension, so the measured pH runs from 2.25 to 8.26 over the six days, and the measured value at each sampling is what enters the model. The series is therefore an adsorption envelope acquired by drift rather than by titration, and it is treated as one.

The **isotherm arm** is the produced water series: six metals together, each near 80 mg L⁻¹, in synthetic produced water of ionic strength 1.0 M containing 4186 mg L⁻¹ calcium and 1006 mg L⁻¹ magnesium. It probes a different total metal in a different matrix and therefore samples a different point on the same isotherm.

### 2.2 Choice of model

The constant capacitance model is adopted on the ground Hayes et al. (1991) give. Their analysis restricts the diffuse layer model to low ionic strength, because it significantly overpredicts the diffuse layer potential above roughly 0.1 M, and restricts the constant capacitance model to high ionic strength, because its charge and potential relation, σ = Cψ, carries no ionic strength dependence and is justified only where the interfacial potential is small. Both experiments here sit well inside the range where the constant capacitance model is the indicated choice. The triple layer model was not used: it carries seven adjustable parameters including two electrolyte binding constants, and Hayes et al. note that FITEQL will not converge when protolysis and electrolyte binding constants are declared adjustable together. Nineteen usable measurements cannot resolve seven parameters.

### 2.3 Capacitance

Hayes et al. Eq. 22 gives the integral capacitance as C₁ = ε_r ε₀ / d. Taking their range of 6 to 50 for the dielectric constant in the double layer region and 2.3 to 4.3 Å for the distance of closest approach of a hydrated univalent ion, the capacitance cannot lie outside 0.12 to 1.92 F m⁻². A value of 1.06 F m⁻² is used here, the optimum Westall and Hohl (1980) found for aluminium oxide and the value both Goldberg (2004) and Ioannou and Dimirkou (1997) adopt. Two independent checks support it. It lies inside the Hayes bracket, and the diffuse layer capacitance at this ionic strength, ε_r ε₀ κ, is 1.89 F m⁻², which also lies inside the bracket. Hayes et al. note that when C₁ equals ε_r ε₀ κ the diffuse layer model reduces exactly to the constant capacitance model, so the two routes to a capacitance for this brine agree on where it belongs.

This is a change from earlier work on this system. The carbonate literature convention C = √I / α with α = 0.004, which Pokrovsky et al. (1999a) fitted for dolomite, gives 207 F m⁻² at this ionic strength, which is 107 times the top of the Hayes bracket and corresponds to no physically possible combination of dielectric constant and approach distance. Section 3.1 quantifies what that costs.

### 2.4 Site totals

Three independently sourced site densities are carried, and none is privileged in advance:

* **Goldberg (2004)**, 2.31 sites nm⁻², the value Davis and Kent (1990) recommend for natural materials, converted by her Eq. 10 to 0.175 mmol L⁻¹ here.
* **Zachara et al. (1991)**, 3.46 × 10⁻⁶ mol g⁻¹, measured directly on calcite by ⁴⁵Ca isotopic exchange, giving 0.208 mmol L⁻¹.
* **Pokrovsky et al. (1999a)**, 14 µmol m⁻² for the carbonate site on dolomite, giving 0.638 mmol L⁻¹.

That the first two agree to 16 per cent is a real cross check, because a generic soil site density and an isotopically measured carbonate exchange capacity have no reason to agree that closely by chance. Pokrovsky's value is about three times higher because it counts every surface carbonate group rather than only those that exchange.

### 2.5 Surface and aqueous reactions

The surface reactions are Pokrovsky et al. (1999a) Table 3, dolomite columns, taken as published and held fixed (Table 1). Following Ioannou and Dimirkou, who took their protonation constants from a published compilation of measured values rather than optimising them, only the two metal reactions are adjustable:

>CO₃H⁰ + Co²⁺ ⇌ >CO₃Co⁺ + H⁺ and >CO₃H⁰ + Li⁺ ⇌ >CO₃Li⁰ + H⁺.

Mass action takes the constant capacitance form of Goldberg Eqs. 5 to 7, in which a surface species of charge z carries the factor exp(−zFψ/RT); the site balance is her Eq. 8, the charge balance her Eq. 9, and the closure σ = Cψ is Hayes Eq. 14. Surface charge is assembled from surface species only, in the sense of Hayes Table Ib, with σ₀ = B(Σ z_i [surface species]) and B = F/(C_s·SA); dissolved calcium and magnesium belong to dolomite dissolution and carry no surface charge. Aqueous complexes for each metal are carried with NIST and MINTEQ formation constants.

Goldberg records that in the FITEQL constant capacitance model activities are set equal to concentrations and no activity corrections are performed. That convention was written for soil suspensions near 0.1 M. These experiments run at 0.68 and 1.0 M, so Davies coefficients are applied as the base case and the FITEQL convention is carried as an explicit sensitivity case rather than either being assumed silently.

Every surface species and the resulting charge are computed in logarithmic space, with site fractions formed by a normalised exponential, so that a strongly polarised surface where one species outweighs the reference form by many orders produces neither an overflow nor an indeterminate ratio. Because every positively charged surface species carries the Boltzmann factor and every negatively charged one carries its reciprocal, σ falls monotonically with ψ while Cψ rises, so the root of σ(ψ) = Cψ is unique and the interval ±1.5 V brackets it at any capacitance in the Hayes range.

### 2.6 Separating precipitation from sorption in the isotherm arm

The produced water removes metal by two routes and the surface model describes only one. Zachara et al. (1991) measured log ᶜK_ex(Sr) = −2.04 on calcite and describe strontium and barium as effectively non sorbing on carbonate surfaces. Strontium is therefore the one metal in the suite whose dissolved concentration is set by precipitation alone, and strontianite equilibrium fixes the carbonate activity as a(CO₃²⁻) = K_sp(SrCO₃)/a(Sr²⁺). With the carbonate activity known, every other metal has a solubility floor,

[Me]_floor = K_sp(MeCO₃) / (a(CO₃²⁻) γ₂),

with the square root taken for Li₂CO₃. Three outcomes are possible at each sampling and each means something different. Where the measured concentration lies above the floor, precipitation has not run to completion and no sorption is required. Where it lies below, more metal has left solution than its carbonate can account for and the excess is the sorbed quantity. Where the floor exceeds the initial concentration, the carbonate cannot form at all and the whole removal is a surface process. Only the sorbed quantity is passed to the surface model.

### 2.7 Parameter estimation

Two estimates are made independently.

The **per point inversion** follows Ioannou and Dimirkou. At each measurement the sorbed quantity is known, so the site balance and the closure σ = Cψ determine the surface potential, and the mass action expression then yields log K directly with no optimisation. Because the test asks whether one constant holds across pH within a single experiment, it is applied to the envelope arm alone, where the matrix is fixed and only pH moves. The mean and standard deviation across pH are reported, as in their Table 3.

The **simultaneous fit** follows Goldberg. A single constant per metal is optimised against the envelope and isotherm arms together, minimising the FITEQL objective

WSOS/DF = Σ (Y_i/s_i)² / (n_obs − n_param),  s_i = √(σ_abs² + (σ_rel C_i)²),

with σ_abs = 0.05 mg L⁻¹ and σ_rel = 0.03. Herbelin and Westall take WSOS/DF between about 0.1 and 20 as an acceptable fit. One adjustable parameter per metal is used, within the limit of about two that both Hayes et al. and Ioannou and Dimirkou identify as the point beyond which FITEQL convergence becomes unreliable.

A point whose measured uptake exceeds the total site inventory is excluded, with the reason recorded. No value of log K can place more metal on a monolayer than the monolayer holds, so such a point does not constrain the constant; it only drives the optimiser to its bound.

### 2.8 Sensitivity protocol

Hayes et al. recommend that a unique parameter set be chosen by mapping the trade offs and selecting from the region where the constants are stable, rather than by adopting one set and reporting a single number. Three sweeps are run: the capacitance from 0.1 to 207 F m⁻², the three published site densities, and the two activity conventions. Each is reported in full.

---

## 3. Results

### 3.1 The surface model before any metal constant is asked of it

At C = 1.06 F m⁻² the surface potential reaches 108 mV at pH 4 and −31 mV at pH 9, and the Boltzmann factor spans 0.015 to 3.28, so the electrostatic term is worth up to 1.83 log units in the units a constant is reported in. At C = 207 F m⁻² the potential never exceeds 8.7 mV and the Boltzmann factor spans only 0.80 to 1.40, a correction of at most 0.15 log units. The high capacitance convention therefore suppresses the electrostatic term by a factor of twelve and leaves what is in substance a mass action model carrying the label of an electrostatic one.

Two checks precede any fitting. The model crosses zero charge near pH 8, against the pH 7.5 to 8.5 point of zero charge reported for dolomite. Its charge reaches 0.0012 mmol m⁻² over the pH 4 to 9 window of these experiments, about an order of magnitude below the 0.01 to 0.02 mmol m⁻² Pokrovsky et al. report at the extremes of titrations that ran to pH 3 and pH 11. Both are the right sign and the right order, which is as much as a transferred constant set can be asked to deliver.

### 3.2 The carbonate activity and the precipitation correction

Strontianite equilibrium gives log a(CO₃²⁻) = −5.85 ± 0.01 across all five post zero samplings and both batches. A solubility controlled system should give a steady value and it does, which is the first check that using strontium this way is sound.

With that carbonate activity, lithium is the clean case. Zabuyelite would require 4.1 × 10⁵ mg L⁻¹ of lithium before it could nucleate in the bulk solution, against the 77 mg L⁻¹ present, so lithium cannot precipitate from this water and its entire removal enters the surface model uncorrected. Cobalt is the opposite case: sphaerocobaltite carries most of the removal, and only the fraction that drives the solution below the solubility floor, 0.066 mmol L⁻¹ at day 6 of the pH 6 batch, is passed to the surface model.

### 3.3 Per point constants

Table 2 gives one constant inverted from each envelope measurement. Over pH 7.40 to 8.09 and a threefold range of uptake, cobalt gives −1.46 ± 0.27, and over pH 7.65 to 8.26 lithium gives −2.23 ± 0.67, with no systematic trend in either.

The batches started at pH 2 are excluded from these averages and the reason is quantitative rather than discretionary. They spent most of the run below pH 6, where the model itself predicts almost no uptake, so their measured differences carry a mean uptake to noise ratio of 0.30 against 0.89 for the batches started at pH 6. Their scatter in log K, 1.89 for cobalt and 1.17 for lithium, is a property of the ICP rather than of the dolomite surface.

These spreads are far wider than the ±0.03 log units Ioannou and Dimirkou obtained across six pH values. The difference is not method but signal: their phosphate uptake was most of the phosphate added, while ours is a few per cent of the metal added.

### 3.4 Simultaneous fit

Fitting the envelope and isotherm arms together gives

* log K(>CO₃Co⁺) = **−1.60**, WSOS/DF = 0.39, over eleven points, six envelope and five isotherm;
* log K(>CO₃Li⁰) = **−1.69**, WSOS/DF = 0.14, over eight points, six envelope and two isotherm, with three lithium isotherm points excluded for coverages of 1.20 to 1.48 monolayers.

Both sit inside the accepted WSOS/DF band. Values below one indicate that the residuals are smaller than the analytical error assigned to them, which with ±3 per cent on a difference of two large concentrations is what should happen: the data are consistent with the model but do not pin it tightly.

The two estimates are independent in method, in the data they use and in their weighting, and for cobalt they agree to 0.14 log units, well inside the per point scatter of 0.27. That agreement is the strongest internal evidence in this work that the cobalt constant is real. For lithium they differ by 0.54 log units, still inside the per point scatter of 0.67, but lithium sits at far higher coverage and its margin is correspondingly thinner.

One result is visible in the point by point fit and is worth stating plainly. In produced water the model places almost no cobalt on the surface, predicting under one per cent coverage against 18 per cent in the single ion batch at a comparable pH, because 4186 mg L⁻¹ of calcium holds the exchange sites. The surface model and the solubility analysis of §3.2 reach the same conclusion by independent routes: cobalt recovery from produced water is carbonate mineral formation and not sorption.

### 3.5 Sensitivity

**Capacitance.** Inside the Hayes bracket the cobalt constant moves by 0.008 log units and the lithium constant by 0.005. The capacitance is the one parameter in this model that can be chosen freely without consequence, which reproduces the insensitivity Hayes et al. reported for the protolysis constants and is the reason they treat it as a fitting parameter bounded by theory rather than as a quantity to be measured. Pushed past the bracket it does begin to matter: at 207 F m⁻² cobalt shifts by −0.37 and lithium by −0.51 log units.

**Site total.** This is the dominant sensitivity in the whole model. Across the three published site densities the cobalt constant spans 1.54 log units, reading −0.06 at the Goldberg density, −0.23 at the Zachara density and −1.60 at the Pokrovsky density. That is roughly thirty times what the capacitance is worth, and it is precisely the point Hayes et al. made.

The lithium column decides which site density to adopt. At the two lower densities the measured lithium uptake already exceeds a monolayer, at 2.61 and 2.20 times the inventory respectively, so no constant exists at any value: the model cannot hold that much lithium on a surface that small. Only the Pokrovsky carbonate site density leaves both metals inside the range a surface reaction can describe, and it is also the only one of the three measured on dolomite rather than transferred from calcite or from a generic soil. Those two reasons, one from this data set and one from the literature, are why it is the base case, and every constant reported here is quoted with it attached.

**Activity convention.** Davies coefficients give −1.60 and −1.69; the FITEQL convention of setting activities equal to concentrations gives −1.85 and −1.54. The convention is worth 0.25 log units for cobalt and 0.16 for lithium at this ionic strength, which is why it is stated rather than left implicit.

---

## 4. Discussion

### 4.1 The cobalt constant against the literature

Pokrovsky et al. (1999a) measured the same reaction on the same site of the same mineral for the two metals that build it: log K = −1.8 for calcium and −2.0 for magnesium. The cobalt constant fitted here, −1.60, falls 0.20 log units from calcium and 0.40 from magnesium. A divalent transition metal binding a surface carbonate group slightly more strongly than the alkaline earths that build the mineral is the order the Irving and Williams series predicts, and nothing in the fit constrained it to land there. That is the main external check available on this result, and it is the check the earlier treatment of this system could not offer, because a capacitance two orders outside the admissible range had removed the electrostatic term that makes the comparison meaningful.

### 4.2 Why lithium is a lower bound

Lithium uptake reaches 0.91 of a monolayer in the single ion batches and 1.48 times the site inventory in produced water. The second of those is not a surface complexation result at all, since no constant can reproduce it, and it is the strongest evidence in either data set that lithium is removed by more than the exchange sites. The XRD of Elshebli et al. (2025) identifies zabuyelite in the reacted solid, and §3.2 shows the bulk solution to be four orders of magnitude undersaturated with respect to it, so that phase must be forming at the surface, where the carbonate activity released by dolomite dissolution is far above the bulk value. The constant of −1.69 therefore describes the surface complexation component of lithium uptake and should be read as a lower bound on the total surface affinity.

This refines rather than contradicts the published conclusion. Lithium recovery is a surface process, as Elshebli et al. concluded, but it is two surface processes rather than one, and only the smaller of them is surface complexation.

### 4.3 What limits the constants, and what would tighten them

The limit throughout is that uptake is measured as a difference between two large concentrations. At ±3 per cent precision the best single ion point reaches a signal to noise ratio of 1.9, and that single number explains every wide interval in this work. Two changes fix it without any new chemistry. Raising the solid loading scales the site inventory and therefore the uptake directly. Lowering the initial metal concentration raises the fraction removed. Both move the measurement away from the difference of two nearly equal numbers, which is the only reason the constants reported here carry the uncertainty they do.

A third change would help the lithium constant specifically. Because lithium already sits near monolayer coverage at 138 mg L⁻¹, a lithium series run at a tenth of that concentration would fall in the range where the isotherm is steep and a constant is best determined, and would separate surface complexation from surface zabuyelite formation by starting well below the coverage at which the second appears.

---

## 5. Conclusions

1. The constant capacitance model, parameterised under the Hayes et al. (1991) discipline and fitted under the Goldberg (2004) protocol, yields intrinsic constants for both metals on the dolomite carbonate site: log K(>CO₃Co⁺) = −1.60 and log K(>CO₃Li⁰) = −1.69, at C = 1.06 F m⁻² and N_t = 0.638 mmol L⁻¹.

2. The cobalt constant is supported by two independent estimates that agree to 0.14 log units, and it falls 0.20 log units from the calcium constant measured for the same site, in the order the Irving and Williams series predicts.

3. The lithium constant is a lower bound. Lithium uptake reaches 1.48 times the site inventory in produced water, which no surface complexation constant can reproduce, and the excess is surface zabuyelite formed where the carbonate activity exceeds the bulk value.

4. Cobalt recovery from produced water is carbonate mineral formation rather than sorption. The surface model predicts under one per cent coverage there because calcium holds the exchange sites, and the solubility analysis using strontium as an internal standard reaches the same conclusion independently. This answers the question left open by Elshebli et al. (2025).

5. The site total is the dominant model sensitivity, worth 1.54 log units for cobalt against 0.008 for the capacitance anywhere inside its physically admissible bracket. Surface complexation constants for carbonate minerals should not be reported, and cannot be transferred, without the site total they were fitted at.

---

## Tables

**Table 1.** Surface reactions and intrinsic constants used as fixed parameters, from Pokrovsky et al. (1999a) Table 3, dolomite columns. Only the two metal reactions are adjustable.

| Surface reaction | log K | Source |
|---|---|---|
| >CO₃H⁰ = >CO₃⁻ + H⁺ | −4.8 | Pokrovsky et al. (1999a) |
| >CO₃H⁰ + Ca²⁺ = >CO₃Ca⁺ + H⁺ | −1.8 | Pokrovsky et al. (1999a) |
| >CO₃H⁰ + Mg²⁺ = >CO₃Mg⁺ + H⁺ | −2.0 | Pokrovsky et al. (1999a) |
| >CO₃H⁰ + Co²⁺ = >CO₃Co⁺ + H⁺ | fitted | this work |
| >CO₃H⁰ + Li⁺ = >CO₃Li⁰ + H⁺ | fitted | this work |
| >CaOH⁰ + H⁺ = >CaOH₂⁺ | +11.5 | Pokrovsky et al. (1999a) |
| >CaOH⁰ = >CaO⁻ + H⁺ | −12.0 | Pokrovsky et al. (1999a) |
| >CaOH⁰ + CO₃²⁻ + 2H⁺ = >CaHCO₃ + H₂O | −4.0 | Pokrovsky et al. (1999a) |
| >CaOH⁰ + CO₃²⁻ + H⁺ = >CaCO₃⁻ + H₂O | +16.6 | Pokrovsky et al. (1999a) |
| >MgOH⁰ + H⁺ = >MgOH₂⁺ | +10.6 | Pokrovsky et al. (1999a) |
| >MgOH⁰ = >MgO⁻ + H⁺ | −12.0 | Pokrovsky et al. (1999a) |
| >MgOH⁰ + CO₃²⁻ + 2H⁺ = >MgHCO₃ + H₂O | −3.5 | Pokrovsky et al. (1999a) |
| >MgOH⁰ + CO₃²⁻ + H⁺ = >MgCO₃⁻ + H₂O | +15.4 | Pokrovsky et al. (1999a) |

**Table 2.** Per point constants from the envelope arm, in the form of Ioannou and Dimirkou (1997) Table 3. Batches started at pH 2 are listed but excluded from the averages, for the reason given in §3.3.

| Metal | Batch | pH | Uptake (mg L⁻¹) | Uptake / noise | Coverage N | log K |
|---|---|---|---|---|---|---|
| Co | pH 6 | 7.40 | 2.61 | 0.78 | 0.069 | −1.51 |
| Co | pH 6 | 7.48 | 4.10 | 1.23 | 0.109 | −1.17 |
| Co | pH 6 | 8.09 | 6.12 | 1.86 | 0.163 | −1.70 |
| Co | pH 2 | 4.12 | 0.85 | 0.25 | 0.023 | +1.26 |
| Co | pH 2 | 5.66 | 0.72 | 0.21 | 0.019 | −0.54 |
| Co | pH 2 | 7.11 | 0.12 | 0.03 | 0.003 | −2.52 |
| Li | pH 6 | 7.65 | 3.77 | 0.65 | 0.852 | −1.51 |
| Li | pH 6 | 7.92 | 1.71 | 0.29 | 0.386 | −2.83 |
| Li | pH 6 | 8.26 | 3.17 | 0.55 | 0.715 | −2.34 |
| Li | pH 2 | 4.26 | 2.17 | 0.37 | 0.489 | +0.80 |
| Li | pH 2 | 5.73 | 2.47 | 0.43 | 0.558 | −0.46 |
| Li | pH 2 | 7.05 | 2.84 | 0.49 | 0.641 | −1.54 |
| **Co mean, pH 6 series** | | **7.40 to 8.09** | | | | **−1.46 ± 0.27** |
| **Li mean, pH 6 series** | | **7.65 to 8.26** | | | | **−2.23 ± 0.67** |

**Table 3.** Reported constants. 25 °C, C = 1.06 F m⁻², N_t = 0.638 mmol L⁻¹ (8.4 sites nm⁻², Pokrovsky et al. 1999a carbonate site), Davies activity coefficients.

| Reaction | log K | Per point SD | WSOS/DF | n | Status |
|---|---|---|---|---|---|
| >CO₃H⁰ + Co²⁺ = >CO₃Co⁺ + H⁺ | −1.60 | 0.27 | 0.39 | 11 | determined |
| >CO₃H⁰ + Li⁺ = >CO₃Li⁰ + H⁺ | −1.69 | 0.67 | 0.14 | 8 | lower bound |

**Table 4.** Sensitivity of the fitted constants, following the protocol of Hayes et al. (1991).

| Parameter varied | Range | Δ log K, Co | Δ log K, Li |
|---|---|---|---|
| Capacitance, inside the Hayes Eq. 22 bracket | 0.12 to 1.92 F m⁻² | 0.008 | 0.005 |
| Capacitance, extended to the carbonate convention | 0.1 to 207 F m⁻² | 0.38 | 0.51 |
| Site total, three published densities | 0.175 to 0.638 mmol L⁻¹ | 1.54 | not defined below 0.638 |
| Activity convention, Davies against FITEQL | — | 0.25 | 0.16 |

---

## Figure captions

**Figure 1.** Adsorption envelopes for cobalt and lithium on dolomite. Symbols are measured uptake with ±3 per cent analytical error bars, filled for batches started at pH 6 and open for batches started at pH 2. Lines are the constant capacitance model at the fitted constants of Table 3. The cobalt edge rises through pH 7 to 8 and turns over above pH 8.7 as carbonate complexation removes free Co²⁺ from solution, the same shape Goldberg (2004) reports for boron adsorption envelopes.

**Figure 2.** Predicted against measured uptake for both arms, with the 1:1 line. Horizontal bars are the ±3 per cent analytical error on the uptake, which is a difference of two concentrations. The scatter is of the size that error implies, which is the honest statement of what these data can constrain.

**Figure 3.** Sensitivity of the fitted constants, following Hayes et al. (1991). Left: capacitance, with the Hayes Eq. 22 bracket shaded and the adopted value marked; the constants are flat across the whole admissible range and drift only when the capacitance is pushed orders of magnitude beyond it. Right: site total, at the three published site densities; the cobalt constant spans 1.54 log units and lithium has no solution at the two lower densities because the measured uptake exceeds a monolayer there.

---

## Data and code availability

Every datum, reaction, constant and equation is embedded in `CCM_FITEQL_LiCo_dolomite.ipynb`, which reproduces every number, table and figure in this manuscript from a clean kernel with no external file dependency.

---

## References

Davis, J. A., & Kent, D. B. (1990). Surface complexation modeling in aqueous geochemistry. *Reviews in Mineralogy, 23*, 177–260.

Elshebli, M., Vilcáez, J., & Smay, J. (2025). Lithium and cobalt recovery from petroleum produced water using dolomite: Impact of pH, temperature, and surface area. *Science of the Total Environment, 1003*, 180743. https://doi.org/10.1016/j.scitotenv.2025.180743

Goldberg, S. (2004). Modeling boron adsorption isotherms and envelopes using the constant capacitance model. *Vadose Zone Journal, 3*(2), 676–680. https://doi.org/10.2136/vzj2004.0676

Hayes, K. F., Redden, G., Ela, W., & Leckie, J. O. (1991). Surface complexation models: An evaluation of model parameter estimation using FITEQL and oxide mineral titration data. *Journal of Colloid and Interface Science, 142*(2), 448–469. https://doi.org/10.1016/0021-9797(91)90075-J

Herbelin, A. L., & Westall, J. C. (1996). *FITEQL: A computer program for determination of chemical equilibrium constants from experimental data* (Report 96-01, Version 3.2). Department of Chemistry, Oregon State University.

Ioannou, A., & Dimirkou, A. (1997). Phosphate adsorption on hematite, kaolinite, and kaolinite–hematite (k–h) systems as described by a constant capacitance model. *Journal of Colloid and Interface Science, 192*(1), 119–128. https://doi.org/10.1006/jcis.1997.4970

Pokrovsky, O. S., Schott, J., & Thomas, F. (1999a). Dolomite surface speciation and reactivity in aquatic systems. *Geochimica et Cosmochimica Acta, 63*(19–20), 3133–3143. https://doi.org/10.1016/S0016-7037(99)00240-9

Smith, R. M., & Martell, A. E. (2004). *NIST critically selected stability constants of metal complexes database* (Version 8.0, NIST Standard Reference Database 46). National Institute of Standards and Technology.

Stumm, W., & Morgan, J. J. (1996). *Aquatic chemistry: Chemical equilibria and rates in natural waters* (3rd ed.). Wiley.

Westall, J., & Hohl, H. (1980). A comparison of electrostatic models for the oxide/solution interface. *Advances in Colloid and Interface Science, 12*(4), 265–294. https://doi.org/10.1016/0001-8686(80)80012-1

Zachara, J. M., Cowan, C. E., & Resch, C. T. (1991). Sorption of divalent metals on calcite. *Geochimica et Cosmochimica Acta, 55*(6), 1549–1562. https://doi.org/10.1016/0016-7037(91)90126-Q
