# Constant capacitance model constants for lithium and cobalt on dolomite, fitted simultaneously to an adsorption envelope and an isotherm under the Pokrovsky parameterisation

**Marwa Elshebli**¹, Javier Vilcáez¹, James Smay²

¹ Boone Pickens School of Geology, Oklahoma State University, Stillwater, OK 74078, USA
² School of Materials, Mechatronics & Manufacturing Engineering, Oklahoma State University, Tulsa, OK 74106, USA

*All numbers in the Results were produced by the accompanying notebook, `CCM_FITEQL_LiCo_dolomite.ipynb`. Every datum is embedded in that notebook and nothing is read from an external file.*

---

## Abstract

Dolomite recovers Li⁺ and Co²⁺ from petroleum produced water, and the quantitative split between sorption and carbonate mineral formation was left open by Elshebli et al. (2025). Here that split is resolved and intrinsic surface complexation constants are reported for both metals. The chemical model is that of Pokrovsky et al. (1999a) for dolomite, adopted whole: their eleven surface constants, their site densities of 14, 7 and 7 µmol m⁻² for the carbonate, calcium and magnesium sites, and their capacitance relation C = √I/α with α = 0.004, which gives 207 F m⁻² for the single ion brine and 250 F m⁻² for the produced water. Adopting the set whole is the point. Their constants were fitted under their capacitance, and transferring the constants under a capacitance taken from the oxide literature makes the set internally inconsistent: at 1.06 F m⁻², the value appropriate to oxides, the 0.02 mmol m⁻² surface charge Pokrovsky et al. measured on dolomite would require a surface potential of 1820 mV. Under their own relation the model reproduces their measured charge, reaching 0.019 mmol m⁻² at the ends of the pH 4 to 9 window, and crosses zero charge near pH 8 against the pH 7.5 to 8.5 point of zero charge reported for dolomite. The fitting protocol is that of Goldberg (2004), in which adsorption isotherms and envelopes are fitted simultaneously with one constant set in FITEQL, and the reporting follows Ioannou and Dimirkou (1997), who inverted one constant per pH and averaged across pH. Two experiments are combined for the first time: the single ion batches, in which dolomite dissolution carries the suspension from pH 2.25 to 8.26 at fixed total metal and which therefore constitute an adsorption envelope, and the produced water batches, which probe a different total metal against 4186 mg L⁻¹ of competing calcium and therefore constitute an isotherm point. Strontium serves as an internal standard for the produced water: because it does not sorb on carbonate surfaces, strontianite equilibrium fixes the carbonate activity at log a(CO₃²⁻) = −5.85 ± 0.01 across every sampling and both batches, which in turn fixes the solubility floor of every other metal and separates the precipitated from the sorbed fraction without the alkalinity measurement the experiments did not record. The fitted constants are log K(>CO₃Co⁺) = −1.98 with WSOS/DF = 0.56 over eleven points, and log K(>CO₃Li⁰) = −2.24 with WSOS/DF = 0.16 over eight points. The cobalt constant is supported independently: inverting one constant from each envelope point and averaging gives −1.93 ± 0.24 across pH 7.40 to 8.09, which agrees with the simultaneous fit to 0.05 log units. It is also supported externally. Cobalt lands 0.02 log units from the magnesium constant Pokrovsky et al. fitted for the same site on the same mineral, and the Shannon radius of Co²⁺, 0.745 Å, is within 3 per cent of that of Mg²⁺, 0.72 Å, and 26 per cent from that of Ca²⁺; an affinity that lands on magnesium's rather than calcium's is what the crystal chemistry predicts, and nothing in the fit was constrained to produce it. The lithium constant is reported as a lower bound, because lithium uptake reaches 0.77 of a monolayer in the single ion batches and 1.48 times the site inventory in produced water, which no surface complexation constant can reproduce; the excess is surface zabuyelite, consistent with the XRD of Elshebli et al. (2025) and with the bulk solution being four orders of magnitude undersaturated with respect to that phase. Sensitivity analysis after Hayes et al. (1991) shows the site total to be the dominant parameter, worth 0.77 log units for cobalt across three published site densities, five times what the capacitance is worth across the carbonate range, so every constant is quoted with its site total attached. For cobalt in produced water the surface model predicts less than one per cent coverage because calcium holds the exchange sites, which is the same conclusion the solubility analysis reaches by an independent route: cobalt recovery from produced water is carbonate mineral formation, not sorption.

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

Pokrovsky et al. (1999a) fitted every dolomite constant used here under the carbonate convention C = √I/α with α = 0.004. Their constants and their capacitance are one parameter set. Transferring the constants under a capacitance taken from somewhere else would make the set internally inconsistent, so the whole parameterisation is adopted, and because the relation is a function of ionic strength the two experiments take different values: 207 F m⁻² at I = 0.684 M and 250 F m⁻² at I = 1.0 M.

That is far above the 0.12 to 1.92 F m⁻² that Hayes et al. Eq. 22, C₁ = ε_r ε₀/d, allows for ε_r between 6 and 50 and d between 2.3 and 4.3 Å. The difference is a property of the mineral rather than a discrepancy. Hayes et al. derived their bracket for oxides, whose site densities are a few sites per square nanometre. Dolomite carries 8.4 carbonate sites per square nanometre, and Pokrovsky et al. measured its surface charge at 0.01 to 0.02 mmol m⁻², which is 0.96 to 1.93 C m⁻². A charge that large at the tens of millivolts a carbonate surface actually reaches implies C = σ/ψ of order 50 to 200 F m⁻², so the high capacitance follows directly from their own measurement. The check is symmetric and decisive: at the oxide value of 1.06 F m⁻², a charge of 0.02 mmol m⁻² would require a surface potential of 1820 mV, which no interface sustains. Carbonate surfaces do not sit inside a bracket derived for oxides.

The capacitance is nevertheless carried through the sensitivity protocol of §2.8, so what the choice is worth is reported rather than asserted.

### 2.4 Site totals

Three independently sourced site densities are carried, and none is privileged in advance:

* **Goldberg (2004)**, 2.31 sites nm⁻², the value Davis and Kent (1990) recommend for natural materials, converted by her Eq. 10 to 0.175 mmol L⁻¹ here.
* **Zachara et al. (1991)**, 3.46 × 10⁻⁶ mol g⁻¹, measured directly on calcite by ⁴⁵Ca isotopic exchange, giving 0.208 mmol L⁻¹.
* **Pokrovsky et al. (1999a)**, 14 µmol m⁻² for the carbonate site on dolomite, giving 0.638 mmol L⁻¹.

That the first two agree to 16 per cent is a real cross check, because a generic soil site density and an isotopically measured carbonate exchange capacity have no reason to agree that closely by chance. Pokrovsky's value is about three times higher because it counts every surface carbonate group rather than only those that exchange. It is the base case here for two reasons given in §3.5: it is the only one of the three measured on dolomite itself, and it is the only one that leaves lithium inside a monolayer.

### 2.5 Surface and aqueous reactions

The surface reactions are Pokrovsky et al. (1999a) Table 3, dolomite columns, taken as published and held fixed (Table 1). Two of them are corrected here: reaction 5, >MeOH⁰ + CO₃²⁻ + 2H⁺ ⇌ >MeHCO₃⁰ + H₂O, was carried in earlier drafts of this work as −4.0 for calcium and −3.5 for magnesium, and the published values reproduced in Table 1 of Ebrahimi and Vilcáez (2018) are +24.0 and +23.5. Thermodynamics alone identifies the error: reaction 5 is reaction 6 with one further proton on the left, so its constant must be the larger of the two. The correction moves the fitted constants by 0.01 log units for cobalt and 0.03 for lithium, because these two species sit on the calcium and magnesium sites and reach the metal only through the surface potential, which is small at Pokrovsky's capacitance. Following Ioannou and Dimirkou, who took their protonation constants from a published compilation of measured values rather than optimising them, only the two metal reactions are adjustable:

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

Under Pokrovsky's capacitance the surface potential runs from +5.6 mV at pH 4 to −8.7 mV at pH 9, and the Boltzmann factor spans 0.80 to 1.40, a correction of up to 0.15 log units in the units a constant is reported in. The electrostatic term on a carbonate surface is genuinely small, which is a consequence of the high site density rather than a deficiency of the model: the same charge that forces the capacitance to be large forces the potential to be small. At the oxide value of 1.06 F m⁻² the same model would place the potential at +108 mV at pH 4, which no measurement of a carbonate surface supports.

Two checks precede any fitting. The model crosses zero charge near pH 8, against the pH 7.5 to 8.5 point of zero charge reported for dolomite. Its charge reaches 0.019 mmol m⁻² at the ends of the pH 4 to 9 window, which falls inside the 0.01 to 0.02 mmol m⁻² Pokrovsky et al. measured by titration. Holding their capacitance together with their constants reproduces their own surface charge; separating the two does not. That is the internal consistency check that decides the parameterisation.

### 3.2 The carbonate activity and the precipitation correction

Strontianite equilibrium gives log a(CO₃²⁻) = −5.85 ± 0.01 across all five post zero samplings and both batches. A solubility controlled system should give a steady value and it does, which is the first check that using strontium this way is sound.

With that carbonate activity, lithium is the clean case. Zabuyelite would require 4.1 × 10⁵ mg L⁻¹ of lithium before it could nucleate in the bulk solution, against the 77 mg L⁻¹ present, so lithium cannot precipitate from this water and its entire removal enters the surface model uncorrected. Cobalt is the opposite case: sphaerocobaltite carries most of the removal, and only the fraction that drives the solution below the solubility floor, 0.066 mmol L⁻¹ at day 6 of the pH 6 batch, is passed to the surface model.

### 3.3 Per point constants

Table 2 gives one constant inverted from each envelope measurement. Over pH 7.40 to 8.09 and a threefold range of uptake, cobalt gives −1.93 ± 0.24, and over pH 7.65 to 8.26 lithium gives −2.42 ± 0.48, with no systematic trend in either.

The batches started at pH 2 are excluded from these averages and the reason is quantitative rather than discretionary. They spent most of the run below pH 6, where the model itself predicts almost no uptake, so their measured differences carry a mean uptake to noise ratio of 0.30 against 0.89 for the batches started at pH 6. Their scatter in log K, 1.89 for cobalt and 1.17 for lithium, is a property of the ICP rather than of the dolomite surface.

These spreads are far wider than the ±0.03 log units Ioannou and Dimirkou obtained across six pH values. The difference is not method but signal: their phosphate uptake was most of the phosphate added, while ours is a few per cent of the metal added.

### 3.4 Simultaneous fit

Fitting the envelope and isotherm arms together gives

* log K(>CO₃Co⁺) = **−1.98**, WSOS/DF = 0.56, over eleven points, six envelope and five isotherm;
* log K(>CO₃Li⁰) = **−2.24**, WSOS/DF = 0.16, over eight points, six envelope and two isotherm, with three lithium isotherm points excluded for coverages of 1.20 to 1.48 monolayers.

Both sit inside the accepted WSOS/DF band. Values below one indicate that the residuals are smaller than the analytical error assigned to them, which with ±3 per cent on a difference of two large concentrations is what should happen: the data are consistent with the model but do not pin it tightly.

The two estimates are independent in method, in the data they use and in their weighting, and for cobalt they agree to 0.05 log units, well inside the per point scatter of 0.24. That agreement is the strongest internal evidence in this work that the cobalt constant is real. For lithium they differ by 0.18 log units, still inside the per point scatter of 0.48, but lithium sits at far higher coverage and its margin is correspondingly thinner.

One result is visible in the point by point fit and is worth stating plainly. In produced water the model places almost no cobalt on the surface, predicting under one per cent coverage against 10 per cent in the single ion batch at a comparable pH, because 4186 mg L⁻¹ of calcium holds the exchange sites. The surface model and the solubility analysis of §3.2 reach the same conclusion by independent routes: cobalt recovery from produced water is carbonate mineral formation and not sorption.

### 3.5 Sensitivity

**Capacitance.** Across the carbonate range, 50 F m⁻² and above, the cobalt constant moves by 0.14 log units and the lithium constant by 0.31, so the exact value of α is not critical once the capacitance is in the range a carbonate surface implies. Hayes et al. reported the same kind of plateau for oxides above 1.2 F m⁻²; in both cases the constant stops responding to the capacitance once the potential has become small. Moving all the way down to the oxide value of 1.06 F m⁻² shifts cobalt by +0.37 and lithium by +0.52 log units, which is the price of pairing a constant set fitted under one capacitance with a capacitance taken from a different mineral class.

**Site total.** This is the dominant sensitivity in the whole model. Across the three published site densities the cobalt constant spans 0.77 log units, reading −1.21 at the Goldberg density, −1.32 at the Zachara density and −1.98 at the Pokrovsky density. That is five times what the capacitance is worth across the carbonate range, and it is precisely the point Hayes et al. made.

The lithium column decides which site density to adopt. At the two lower densities the measured lithium uptake already exceeds a monolayer, at 2.61 and 2.20 times the inventory respectively, so no constant exists at any value: the model cannot hold that much lithium on a surface that small. Only the Pokrovsky carbonate site density leaves both metals inside the range a surface reaction can describe, and it is also the only one of the three measured on dolomite rather than transferred from calcite or from a generic soil. Those two reasons, one from this data set and one from the literature, are why it is the base case, and every constant reported here is quoted with it attached.

**Activity convention.** Davies coefficients give −1.98 and −2.24; the FITEQL convention of setting activities equal to concentrations gives −2.30 and −2.26. The convention is worth 0.32 log units for cobalt and 0.02 for lithium at this ionic strength, which is why it is stated rather than left implicit.

---

## 4. Discussion

### 4.1 The cobalt constant against the literature

Pokrovsky et al. (1999a) measured the same reaction on the same site of the same mineral for the two metals that build it: log K = −1.8 for calcium and −2.0 for magnesium. The cobalt constant fitted here, −1.98, falls 0.02 log units from magnesium and 0.18 from calcium.

The comparison that matters is with ionic size, because substitution into a surface carbonate site is a size question before it is anything else. On the Shannon (1976) six coordinate scale, Co²⁺ has a radius of 0.745 Å, within 3 per cent of the 0.72 Å of Mg²⁺ and 26 per cent from the 1.00 Å of Ca²⁺. An affinity that lands on magnesium's rather than calcium's is therefore what the crystal chemistry predicts. Sphaerocobaltite and magnesite are isostructural for the same reason, and cobalt substitutes for magnesium in natural carbonates far more readily than for calcium.

The strength of this check lies in its independence. Pokrovsky's magnesium constant was fitted by different workers, from different data, on the same site of the same mineral, and it enters the present model as a fixed parameter in the competition term rather than as anything the metal fit can adjust. Nothing in the optimisation was constrained to approach it, and it is approached to 0.02 log units.

Lithium is also close to magnesium in size, at 0.76 Å, but carries one charge rather than two, and its constant comes out 0.24 log units weaker, which is the direction a halved charge on the same site predicts.

### 4.2 Why lithium is a lower bound

Lithium uptake reaches 0.77 of a monolayer in the single ion batches and 1.48 times the site inventory in produced water. The second of those is not a surface complexation result at all, since no constant can reproduce it, and it is the strongest evidence in either data set that lithium is removed by more than the exchange sites. The XRD of Elshebli et al. (2025) identifies zabuyelite in the reacted solid, and §3.2 shows the bulk solution to be four orders of magnitude undersaturated with respect to it, so that phase must be forming at the surface, where the carbonate activity released by dolomite dissolution is far above the bulk value. The constant of −2.24 therefore describes the surface complexation component of lithium uptake and should be read as a lower bound on the total surface affinity.

This refines rather than contradicts the published conclusion. Lithium recovery is a surface process, as Elshebli et al. concluded, but it is two surface processes rather than one, and only the smaller of them is surface complexation.

### 4.3 What limits the constants, and what would tighten them

The limit throughout is that uptake is measured as a difference between two large concentrations. At ±3 per cent precision the best single ion point reaches a signal to noise ratio of 1.9, and that single number explains every wide interval in this work. Two changes fix it without any new chemistry. Raising the solid loading scales the site inventory and therefore the uptake directly. Lowering the initial metal concentration raises the fraction removed. Both move the measurement away from the difference of two nearly equal numbers, which is the only reason the constants reported here carry the uncertainty they do.

A third change would help the lithium constant specifically. Because lithium already sits near monolayer coverage at 138 mg L⁻¹, a lithium series run at a tenth of that concentration would fall in the range where the isotherm is steep and a constant is best determined, and would separate surface complexation from surface zabuyelite formation by starting well below the coverage at which the second appears.

---

## 5. Conclusions

1. The Pokrovsky et al. (1999a) dolomite parameterisation, adopted whole and fitted under the Goldberg (2004) protocol, yields intrinsic constants for both metals on the carbonate site: log K(>CO₃Co⁺) = −1.98 and log K(>CO₃Li⁰) = −2.24, at C = √I/0.004 and N_t = 0.638 mmol L⁻¹.

2. The surface constants and the capacitance of a published model must be transferred together. Pairing Pokrovsky's constants with the capacitance appropriate to oxides would shift the cobalt constant by 0.37 log units and would imply a surface potential of 1820 mV for the surface charge those same authors measured.

3. The cobalt constant is supported by two independent estimates that agree to 0.05 log units, and it falls 0.02 log units from the magnesium constant measured for the same site, which is where the near equality of the Co²⁺ and Mg²⁺ ionic radii places it.

4. The lithium constant is a lower bound. Lithium uptake reaches 1.48 times the site inventory in produced water, which no surface complexation constant can reproduce, and the excess is surface zabuyelite formed where the carbonate activity exceeds the bulk value.

5. Cobalt recovery from produced water is carbonate mineral formation rather than sorption. The surface model predicts under one per cent coverage there because calcium holds the exchange sites, and the solubility analysis using strontium as an internal standard reaches the same conclusion independently. This answers the question left open by Elshebli et al. (2025).

6. The site total is the dominant model sensitivity, worth 0.77 log units for cobalt against 0.14 for the capacitance across the carbonate range. Surface complexation constants for carbonate minerals should not be reported, and cannot be transferred, without the site total they were fitted at.

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
| >CaOH⁰ + CO₃²⁻ + 2H⁺ = >CaHCO₃ + H₂O | +24.0 | Pokrovsky et al. (1999a) |
| >CaOH⁰ + CO₃²⁻ + H⁺ = >CaCO₃⁻ + H₂O | +16.6 | Pokrovsky et al. (1999a) |
| >MgOH⁰ + H⁺ = >MgOH₂⁺ | +10.6 | Pokrovsky et al. (1999a) |
| >MgOH⁰ = >MgO⁻ + H⁺ | −12.0 | Pokrovsky et al. (1999a) |
| >MgOH⁰ + CO₃²⁻ + 2H⁺ = >MgHCO₃ + H₂O | +23.5 | Pokrovsky et al. (1999a) |
| >MgOH⁰ + CO₃²⁻ + H⁺ = >MgCO₃⁻ + H₂O | +15.4 | Pokrovsky et al. (1999a) |

**Table 2.** Per point constants from the envelope arm, in the form of Ioannou and Dimirkou (1997) Table 3. Batches started at pH 2 are listed but excluded from the averages, for the reason given in §3.3.

| Metal | Batch | pH | Uptake (mg L⁻¹) | Uptake / noise | Coverage N | log K |
|---|---|---|---|---|---|---|
| Co | pH 6 | 7.40 | 2.61 | 0.78 | 0.069 | −1.93 |
| Co | pH 6 | 7.48 | 4.10 | 1.23 | 0.109 | −1.65 |
| Co | pH 6 | 8.09 | 6.12 | 1.86 | 0.163 | −2.18 |
| Co | pH 2 | 4.12 | 0.85 | 0.25 | 0.023 | +1.26 |
| Co | pH 2 | 5.66 | 0.72 | 0.21 | 0.019 | −0.54 |
| Co | pH 2 | 7.11 | 0.12 | 0.03 | 0.003 | −2.52 |
| Li | pH 6 | 7.65 | 3.77 | 0.65 | 0.852 | −1.73 |
| Li | pH 6 | 7.92 | 1.71 | 0.29 | 0.386 | −2.97 |
| Li | pH 6 | 8.26 | 3.17 | 0.55 | 0.715 | −2.56 |
| Li | pH 2 | 4.26 | 2.17 | 0.37 | 0.489 | +0.80 |
| Li | pH 2 | 5.73 | 2.47 | 0.43 | 0.558 | −0.46 |
| Li | pH 2 | 7.05 | 2.84 | 0.49 | 0.641 | −1.54 |
| **Co mean, pH 6 series** | | **7.40 to 8.09** | | | | **−1.93 ± 0.24** |
| **Li mean, pH 6 series** | | **7.65 to 8.26** | | | | **−2.42 ± 0.48** |

**Table 3.** Reported constants. 25 °C, Pokrovsky et al. (1999a) parameterisation throughout: C = √I/0.004, giving 207 F m⁻² for the envelope arm and 250 F m⁻² for the isotherm arm, N_t = 0.638 mmol L⁻¹ (8.4 sites nm⁻², carbonate site), Davies activity coefficients.

| Reaction | log K | Per point SD | WSOS/DF | n | Status |
|---|---|---|---|---|---|
| >CO₃H⁰ + Co²⁺ = >CO₃Co⁺ + H⁺ | −1.98 | 0.24 | 0.56 | 11 | determined |
| >CO₃H⁰ + Li⁺ = >CO₃Li⁰ + H⁺ | −2.24 | 0.48 | 0.16 | 8 | lower bound |

For comparison, on the same site of the same mineral, Pokrovsky et al. (1999a) give log K = −1.8 for Ca²⁺ and −2.0 for Mg²⁺.

**Table 4.** Sensitivity of the fitted constants, following the protocol of Hayes et al. (1991).

| Parameter varied | Range | Δ log K, Co | Δ log K, Li |
|---|---|---|---|
| Capacitance, across the carbonate range | 50 to 400 F m⁻² | 0.14 | 0.31 |
| Capacitance, down to the oxide value | 1.06 against 207 F m⁻² | 0.37 | 0.52 |
| Site total, three published densities | 0.175 to 0.638 mmol L⁻¹ | 0.77 | not defined below 0.638 |
| Activity convention, Davies against FITEQL | — | 0.32 | 0.02 |

---

## Figure captions

**Figure 1.** Adsorption envelopes for cobalt and lithium on dolomite. Symbols are measured uptake with ±3 per cent analytical error bars, filled for batches started at pH 6 and open for batches started at pH 2. Lines are the constant capacitance model at the fitted constants of Table 3, under the Pokrovsky parameterisation. The cobalt edge rises through pH 7 to 8 and turns over above pH 8.7 as carbonate complexation removes free Co²⁺ from solution, the same shape Goldberg (2004) reports for boron adsorption envelopes.

**Figure 2.** Predicted against measured uptake for both arms, with the 1:1 line. Horizontal bars are the ±3 per cent analytical error on the uptake, which is a difference of two concentrations. The scatter is of the size that error implies, which is the honest statement of what these data can constrain.

**Figure 3.** Sensitivity of the fitted constants, following Hayes et al. (1991). Left: capacitance, with the Hayes Eq. 22 bracket for oxides shaded and Pokrovsky's value for this brine marked; the constants plateau at each end and the whole span between the oxide and carbonate regimes is 0.37 log units for cobalt. Right: site total, at the three published site densities; the cobalt constant spans 0.77 log units and lithium has no solution at the two lower densities because the measured uptake exceeds a monolayer there.

---

## Data and code availability

Every datum, reaction, constant and equation is embedded in `CCM_FITEQL_LiCo_dolomite.ipynb`, which reproduces every number, table and figure in this manuscript from a clean kernel with no external file dependency.

---

## References

Davis, J. A., & Kent, D. B. (1990). Surface complexation modeling in aqueous geochemistry. *Reviews in Mineralogy, 23*, 177–260.

Elshebli, M., Vilcáez, J., & Smay, J. (2025). Lithium and cobalt recovery from petroleum produced water using dolomite: Impact of pH, temperature, and surface area. *Science of the Total Environment, 1003*, 180743. https://doi.org/10.1016/j.scitotenv.2025.180743

Ebrahimi, P., & Vilcáez, J. (2018). Effect of brine salinity and guar gum on the transport of barium through dolomite rocks: Implications for unconventional oil and gas wastewater disposal. *Journal of Environmental Management, 214*, 370–378. https://doi.org/10.1016/j.jenvman.2018.03.008

Goldberg, S. (2004). Modeling boron adsorption isotherms and envelopes using the constant capacitance model. *Vadose Zone Journal, 3*(2), 676–680. https://doi.org/10.2136/vzj2004.0676

Hayes, K. F., Redden, G., Ela, W., & Leckie, J. O. (1991). Surface complexation models: An evaluation of model parameter estimation using FITEQL and oxide mineral titration data. *Journal of Colloid and Interface Science, 142*(2), 448–469. https://doi.org/10.1016/0021-9797(91)90075-J

Herbelin, A. L., & Westall, J. C. (1996). *FITEQL: A computer program for determination of chemical equilibrium constants from experimental data* (Report 96-01, Version 3.2). Department of Chemistry, Oregon State University.

Ioannou, A., & Dimirkou, A. (1997). Phosphate adsorption on hematite, kaolinite, and kaolinite–hematite (k–h) systems as described by a constant capacitance model. *Journal of Colloid and Interface Science, 192*(1), 119–128. https://doi.org/10.1006/jcis.1997.4970

Pokrovsky, O. S., Schott, J., & Thomas, F. (1999a). Dolomite surface speciation and reactivity in aquatic systems. *Geochimica et Cosmochimica Acta, 63*(19–20), 3133–3143. https://doi.org/10.1016/S0016-7037(99)00240-9

Shannon, R. D. (1976). Revised effective ionic radii and systematic studies of interatomic distances in halides and chalcogenides. *Acta Crystallographica A, 32*(5), 751–767. https://doi.org/10.1107/S0567739476001551

Smith, R. M., & Martell, A. E. (2004). *NIST critically selected stability constants of metal complexes database* (Version 8.0, NIST Standard Reference Database 46). National Institute of Standards and Technology.

Stumm, W., & Morgan, J. J. (1996). *Aquatic chemistry: Chemical equilibria and rates in natural waters* (3rd ed.). Wiley.

Westall, J., & Hohl, H. (1980). A comparison of electrostatic models for the oxide/solution interface. *Advances in Colloid and Interface Science, 12*(4), 265–294. https://doi.org/10.1016/0001-8686(80)80012-1

Zachara, J. M., Cowan, C. E., & Resch, C. T. (1991). Sorption of divalent metals on calcite. *Geochimica et Cosmochimica Acta, 55*(6), 1549–1562. https://doi.org/10.1016/0016-7037(91)90126-Q
