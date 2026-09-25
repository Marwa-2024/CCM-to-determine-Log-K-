# Intrinsic surface complexation constants for lithium and cobalt on dolomite: single ion in sodium chloride against produced water

**Marwa Elshebli**¹, Javier Vilcáez¹, James Smay²

¹ Boone Pickens School of Geology, Oklahoma State University, Stillwater, OK 74078, USA
² School of Materials, Mechatronics & Manufacturing Engineering, Oklahoma State University, Tulsa, OK 74106, USA

*Every number in the Results was produced by the accompanying notebook, `LogK_LiCo_singleion_and_PW.ipynb`, in which all data are embedded and nothing is read from an external file.*

---

## Abstract

Dolomite recovers Li⁺ and Co²⁺ from petroleum produced water, and the quantitative split between surface complexation and carbonate mineral formation was left open by Elshebli et al. (2025). Here intrinsic surface complexation constants are determined for both metals, in two experiments: the metal alone in 40 g L⁻¹ NaCl, and the same metals in produced water carrying 4186 mg L⁻¹ calcium, 1006 mg L⁻¹ magnesium and lead, cadmium, barium and strontium at about 80 mg L⁻¹ each. The chemical model is that of Pokrovsky, Schott and Thomas (1999) for dolomite, adopted whole: their three site model with its published constants held fixed, their site densities of 14 µmol m⁻² for the carbonate site and 7 µmol m⁻² for the calcium and magnesium sites, and their capacitance relation C = √I/α with α = 0.004. A single reaction is added for each metal, >CO₃H⁰ + Meᶻ⁺ ⇌ >CO₃Me⁽ᶻ⁻¹⁾⁺ + H⁺, and its constant is the only adjustable quantity. One entry in the source table needs care and is settled here by three independent checks: reaction 5, >MeOH⁰ + CO₃²⁻ + 2H⁺ ⇌ >MeHCO₃⁰ + H₂O, is positive (+24.0 for calcium, +23.5 for magnesium), because the alternative reading would make the protonation of the surface carbonate complex unfavourable by twenty orders of magnitude, would place >MeHCO₃⁰ at 10⁻²⁵ relative abundance despite Pokrovsky et al. naming it one of five species the model requires, and would give a scattered rather than a consistent protonation constant across the four mineral columns of that table. Aqueous chloro complexation is carried explicitly for every metal, because at this salinity it sets how much of the dissolved metal is present as the free ion the surface binds. Dolomite dissolution is not simulated: pH, calcium and magnesium were measured at every sampling and are used directly. Saturation is screened before any constant is formed. In the single ion system neither metal can precipitate at five of six samplings for cobalt and all six for lithium; the exception is cobalt at day 6 of the batch started at pH 6, where the pH has reached 8.09 and the sphaerocobaltite solubility floor has fallen to 5.7 mg L⁻¹ against 74.4 mg L⁻¹ measured, so the largest uptake in the cobalt series is the one point that must be dropped. In produced water cobalt is supersaturated with sphaerocobaltite at five of six samplings and three of the six lithium points exceed a monolayer, at coverages up to 1.48. The single ion experiment therefore yields both constants and produced water yields neither cleanly. Inverting the constant at each sampling and averaging gives log K_int = −2.39 ± 0.32 for lithium over all six samplings across pH 4.26 to 8.26, and −2.34 ± 0.61 for cobalt over five samplings across pH 4.12 to 7.48, with no trend in either. Both fall close to the constants Pokrovsky et al. measured for the same reaction on the same site: cobalt is 0.34 log units from magnesium and lithium 0.39, and the Shannon radius of Co²⁺, 0.745 Å, is within 3 per cent of that of Mg²⁺ and 26 per cent from that of Ca²⁺, so an affinity close to magnesium's is what the crystal chemistry predicts. Those constants enter the model only as the fixed competition term, which the fit cannot adjust. In produced water the same inversion returns values 1.26 and 2.99 log units higher, and the cause is identified rather than assumed: calcium and magnesium occupy 96 per cent of the carbonate sites there against 6 per cent in the single ion system, so removal that is not surface complexation is forced into the constant as an inflated affinity. The practical consequence is that the greater recovery produced water achieves is not stronger sorption. Sorption is weaker there. The additional recovery is carbonate mineral formation for cobalt and uptake beyond a monolayer for lithium, and the two respond to different process levers.

---

## 1. Introduction

Surface complexation constants are the transferable form of a sorption measurement. A percentage removed describes one experiment; an intrinsic constant, quoted with the site density and capacitance it was obtained at, can be carried into a speciation or transport calculation for a different water. Elshebli et al. (2025) established that dolomite recovers lithium and cobalt from produced water and that the two metals behave differently, and left the quantitative mechanism split open. This work determines the constants and answers that question.

The obstacle is that produced water is not a clean sorption experiment. It contains enough calcium to saturate the surface and, over six days, dolomite dissolution releases enough carbonate for several metals to precipitate. Both processes remove metal from solution and only one of them is surface complexation. A constant inverted from total removal in such a system is not a sorption constant, whatever it is called.

The procedure used here is the one Ebrahimi and Vilcáez (2018) set out for obtaining an intrinsic constant on dolomite in a sodium chloride brine: adopt the Pokrovsky et al. (1999) model with its constants fixed, add one reaction for the metal on the carbonate site, carry the aqueous chloro complexes explicitly, account for dolomite dissolution, and estimate the single unknown by fitting the measured sorption profile. Their procedure is followed; none of their results is used, and every number reported below comes from the experiments of this study.

---

## 2. Materials and methods

### 2.1 Experiments

Dolomite at 60 g L⁻¹, six grams per 100 mL, BET specific surface area 0.76 m² g⁻¹, giving 45.6 m² per litre of suspension. All experiments at 25 °C, sampled at days 0, 2, 4 and 6, analysed by ICP-OES with a relative precision of ±3 per cent on each concentration.

* **Single ion.** Cobalt at 80.5 mg L⁻¹ or lithium at 138.4 mg L⁻¹ alone in 40 g L⁻¹ NaCl, ionic strength 0.684 M, in batches started at nominal pH 2 and pH 6.
* **Produced water.** Lead, cadmium, cobalt, barium, strontium and lithium together at about 80 mg L⁻¹ each, in a brine of ionic strength 1.0 M carrying 4186 mg L⁻¹ calcium and 1006 mg L⁻¹ magnesium, in the same two starting conditions.

The batch labels describe the starting condition only. Dolomite dissolution buffers the suspension and the measured pH runs from 2.25 to 8.26 over the six days; the measured value at each sampling is what enters the model.

### 2.2 Chemical model

Table 1 gives the surface reactions, from Pokrovsky et al. (1999) Table 3, dolomite columns, held fixed. One reaction for each metal is added and its constant is the only adjustable quantity in the work:

>CO₃H⁰ + Co²⁺ ⇌ >CO₃Co⁺ + H⁺ and >CO₃H⁰ + Li⁺ ⇌ >CO₃Li⁰ + H⁺,

with K_int = [>CO₃Me⁺][H⁺] / ([>CO₃H⁰][Me²⁺]).

Site densities are Pokrovsky's own: 14 µmol m⁻² for the carbonate site and 7 µmol m⁻² for the calcium and magnesium sites, a 1:1:2 stoichiometry corresponding to eight sites per square nanometre, giving 0.638 mmol L⁻¹ of carbonate sites here. Surface charge is assembled from surface species only, as his Eq. 5; the apparent and intrinsic constants are related by his Eq. 2, K_s = K_int exp(−zFψ₀/RT); and the constant capacitance closure is his Eq. 3, ψ₀ = σ/C, with C = √I/α and α = 0.004 C mol^½ V⁻¹ m³, the value he fitted to reproduce the high surface charge of a carbonate and its weak response to ionic strength. At the ionic strengths here that gives 207 and 250 F m⁻².

### 2.3 Reading reaction 5

Text extracted from the source PDF renders the leading minus sign of a negative entry as a digit, so reaction 5 of Table 3 is ambiguous between +24.0 and −4.0 for the calcium site and +23.5 and −3.5 for magnesium, with the same ambiguity that rows 1 to 3 carry where the values are certainly negative. Three independent checks settle it, and all point the same way.

1. Reaction 5 minus reaction 6 is the protonation >MeCO₃⁻ + H⁺ ⇌ >MeHCO₃⁰, which must be positive and of order 10, since aqueous HCO₃⁻ formation is 10.33. The positive reading gives +6.4 to +8.1 across all four mineral columns of that table, a tight and sensible spread. The negative reading gives −16.8 to −20.6, scattered and thermodynamically impossible.
2. Pokrovsky et al. name the five species their model requires: >CO₃⁻, >CO₃Me⁺, >MeOH₂⁺, >MeHCO₃⁰ and >MeCO₃⁻. At the conditions of these experiments the two readings put >MeHCO₃⁰ at 10⁺³·¹ and 10⁻²⁴·⁹ relative to >MeOH⁰. A species at 10⁻²⁵ cannot be one of five a model is built on.
3. Under the negative reading >MeHCO₃⁰ disappears from the surface entirely and the calcium and magnesium sites are speciated wrongly at every pH in these experiments.

The positive values are used throughout. Earlier drafts of this work carried the negative reading, and the correction is recorded here because it changes the speciation of the metal sites completely.

### 2.4 Aqueous speciation

Each metal carries its chloride, carbonate and hydrolysis complexes with formation constants from the NIST and MINTEQ compilations. The constant is written on the free ion, so the free ion fraction is computed at every point from the measured pH, the carbonate activity of §2.5 and the chloride of the brine. This is the step Ebrahimi and Vilcáez identified as the dominant control at brine salinity, and it matters here: cobalt loses about half of its free ion pool to chloride and carbonate complexing, while lithium barely complexes at all.

### 2.5 Dolomite dissolution and the saturation screen

Ebrahimi and Vilcáez solved dolomite dissolution with a transition state rate law because they did not measure the pH and the dissolved calcium and magnesium through the run. These experiments measured all three at every sampling, so the measured values are used directly and no rate law is required. That removes two kinetic parameters from the problem and replaces them with data.

What six days of dissolution also does, which 300 minutes would not, is release enough carbonate for metals to precipitate. Saturation is therefore screened at every sampling before any constant is formed. For the single ion system the carbonate activity is taken from open system equilibrium with atmospheric pCO₂. For produced water it is fixed by strontium: Zachara et al. (1991) measured log ᶜK_ex(Sr) = −2.04 on calcite and describe strontium as effectively non sorbing, so strontianite equilibrium alone sets its dissolved concentration, giving log a(CO₃²⁻) = −5.85 ± 0.01 across every sampling and both batches. Where the measured concentration lies above the solubility floor of the metal's own carbonate, precipitation is incomplete and no sorption is required; where it lies below, only the excess is sorbed; where the floor exceeds the initial concentration, the carbonate cannot form and the whole removal is a surface process.

### 2.6 The dissolution rate law as a check on the carbonate activity

The rate law is not used to generate the pH, but the term in braces within it is not kinetic: it is the dolomite saturation state, and it can be evaluated from the measured calcium, magnesium and pH together with the carbonate activity §2.5 assumes. It therefore tests the one quantity in this work that had to be inferred rather than measured, and it is run for that reason.

Two points of care. The published form carries [H⁺] to the first power, which does not balance the reaction it is written for, CaMg(CO₃)₂ + 2H⁺ ⇌ Ca²⁺ + Mg²⁺ + 2HCO₃⁻; the squared form is used here. The equilibrium constant is built from the dolomite solubility product and the bicarbonate constant, which is equivalent to the tabulated value: the source gives log K_eq = +2.525, and that is exactly log K_sp(dolomite-ord) = −18.13 plus twice 10.325. An earlier draft of this work read that constant as negative and reported it as an error; it is not, and the misreading was ours.

### 2.7 Estimating the constant

The sorbed quantity is known at each sampling from §2.5, so the site balances and the closure σ = Cψ fix the surface potential and the mass action expression returns the constant directly, with no optimisation. A constant is formed only where the sorbed quantity is positive and below one monolayer. Constants are reported as a mean and standard deviation over the usable samplings of each system, because a constant that drifts systematically with pH is not an intrinsic constant and the spread across pH is the evidence either way.

---

## 3. Results

### 3.1 The model before any metal constant is asked of it

At Pokrovsky's own conditions, 0.01 M NaCl with 10⁻³ M calcium and magnesium in an open system, the model crosses zero surface charge at pH 8.5 against the pH 8.0 ± 0.1 he measured for both the point of zero charge and the isoelectric point. The residual of half a pH unit reflects the aqueous model and the site density rather than the constants, but it places the transferred set correctly on the one quantity he determined independently of his own fit, which is the check worth making before using it.

### 3.2 The saturation state of dolomite

Every single ion sampling is undersaturated with dolomite, from −18.6 to −1.5, with implied bicarbonate from below a micromolar to 1.17 mmol L⁻¹. Both are what a suspension dissolving dolomite into a dilute brine should give, so the open system carbonate activity used for that experiment passes the test.

Produced water sits steadily at about +2.2, the expected sign for a brine that dissolves dolomite and cannot reprecipitate it at any useful rate. The implied bicarbonate, however, fails in the acid start batch, reaching 1907 mmol L⁻¹ at pH 4.29. No solution holds that much bicarbonate. Strontianite equilibrium fixes a(CO₃²⁻) almost independently of pH, so the bicarbonate it implies rises as the pH falls, and below about pH 7 the assumption becomes inadmissible: strontium is not solubility controlled there, its removal in that batch being 1.9 to 2.6 per cent, within the analytical noise.

No reported number changes. The single ion constants, which are the recommended values, never use strontium. In produced water lithium needs no precipitation correction at any carbonate activity, so a(CO₃²⁻) reaches its constant only through the minor LiCO₃⁻ complex, and a tenfold change in it moves the lithium constant by less than 0.001 log units. The single usable produced water cobalt point sits at pH 7.86, where the implied bicarbonate is 0.54 mmol L⁻¹ and physically reasonable. The check tightens the account of the analysis without moving a result, and it bounds where the strontium reference may be used.

### 3.3 The saturation screen

Lithium is usable at every sampling of both systems. Zabuyelite would require lithium at 1.9 × 10⁵ mg L⁻¹ at its most soluble point in this set, against the 77 to 138 mg L⁻¹ present, so it cannot form and no correction is applied anywhere.

Cobalt is the metal the carbonate chemistry interferes with, and it does so differently in the two systems. In the single ion experiment cobalt is undersaturated at five of six samplings; the exception is day 6 of the batch started at pH 6, where the pH has risen to 8.09 and the sphaerocobaltite floor has fallen to 5.7 mg L⁻¹ against the 74.4 mg L⁻¹ measured. That is the largest uptake in the whole cobalt series, and it is the one point that has to be dropped. In produced water cobalt is supersaturated at five of six samplings.

### 3.4 Competition for the site

Before any trace metal arrives, calcium and magnesium occupy 96 per cent of the carbonate sites in produced water and 6 per cent in the single ion system. This is the largest single difference between the two experiments, and it means a constant fitted in produced water must be large simply to account for any uptake at all.

### 3.5 The constants

**Table 2** gives the constants. The single ion experiment delivers both:

* log K_int(>CO₃Li⁰) = **−2.39 ± 0.32**, all six samplings, pH 4.26 to 8.26
* log K_int(>CO₃Co⁺) = **−2.34 ± 0.61**, five of six samplings, pH 4.12 to 7.48

Both hold across four pH units and two batches with no trend. Produced water delivers neither cleanly: three of six lithium points exceed a monolayer and five of six cobalt points are supersaturated, and where a value can be formed at all it comes out 1.26 log units higher for lithium and 2.99 higher for cobalt.

### 3.6 Forward check

The constants were obtained by inversion, so the fair test is to put them back into the model and predict the sorption profiles they were never shown (Figure 1). For the single ion system the root mean square difference between measured and modelled removal is 1.3 percentage points against removals that reach 7.6 per cent. For produced water it is 25 points, which is the same conclusion as §3.4 reached from coverage and saturation, arrived at independently: the produced water profiles are not surface complexation profiles and no constant reproduces them.

### 3.7 Against the constants measured for the same site

Pokrovsky et al. measured the same reaction on the same site of the same mineral for the two metals that build it: −1.8 for calcium and −2.0 for magnesium. Cobalt comes out 0.34 log units from magnesium and 0.54 from calcium, and lithium 0.39 from magnesium.

Ordered by ionic radius on the Shannon (1976) six coordinate scale, the four constants read Mg²⁺ 0.72 Å at −2.0, Co²⁺ 0.745 Å at −2.34, Li⁺ 0.76 Å at −2.39 and Ca²⁺ 1.00 Å at −1.8. Cobalt sits within 3 per cent of magnesium in size and 26 per cent from calcium, so an affinity close to magnesium's is what the crystal chemistry predicts; sphaerocobaltite and magnesite are isostructural for the same reason, and cobalt substitutes for magnesium in natural carbonates far more readily than for calcium. Lithium is similar in size but carries one charge rather than two and comes out weaker, the direction a halved charge on the same site predicts.

Neither comparison was built into the fit. The calcium and magnesium constants enter the model only as the fixed competition term, which the metal constant cannot adjust, and both metals were free to come out anywhere over fourteen log units.

### 3.8 Sensitivity

Raising the carbonate site density from Pokrovsky's 14 µmol m⁻² to the 20 µmol m⁻² of Brady et al. (1999), a 43 per cent change, moves lithium by 0.36 log units and cobalt by 0.20. A factor of four in the capacitance is worth 0.03 log units for lithium and 0.09 for cobalt: on a carbonate surface the capacitance is high and the potential correspondingly small, so the electrostatic term does little work, which is why α can be treated as a single empirical parameter. The constants are therefore quoted with the site density they were obtained at, and the capacitance is not critical.

---

## 4. Discussion

The two systems give different constants for the same metal on the same mineral, and the difference is not a measurement problem. It is the signature of removal that is not surface complexation.

Three things separate them, and all three are visible in the data rather than assumed. Coverage: every usable single ion point sits inside a monolayer, while three of the six produced water lithium points reach up to 1.48 monolayers, which no surface complexation constant can reproduce at any value. Precipitation: cobalt is undersaturated at five of six single ion samplings and one of six in produced water. Competition: calcium and magnesium hold 96 per cent of the sites in produced water against 6 per cent in the single ion system, so any removal not actually due to sorption is forced into the constant as an inflated affinity.

The practical reading is the useful one. Dolomite recovers about four times more lithium and twelve times more cobalt from produced water than from a simple brine, and the constants show that the additional recovery is not stronger sorption. Sorption is weaker there, because calcium occupies the sites. The additional recovery is carbonate mineral formation for cobalt and, for lithium, uptake beyond what a monolayer can hold. That distinction matters for process design, because the two respond to different levers: mineral formation to carbonate supply and residence time, surface complexation to surface area and to keeping calcium out of the way.

The uncertainties on the reported constants are set by one thing, that uptake is measured as a difference between two large concentrations. At ±3 per cent precision the best single ion point reaches a signal to noise of 1.9. Raising the solid loading scales the site inventory and therefore the uptake; lowering the initial metal concentration raises the fraction removed. Measuring alkalinity at each sampling would remove the one quantity in this analysis that had to be inferred rather than measured.

---

## 5. Conclusions

1. Intrinsic surface complexation constants on the dolomite carbonate site, at the Pokrovsky et al. (1999) site density of 14 µmol m⁻² and capacitance C = √I/0.004: log K_int(>CO₃Li⁰) = −2.39 ± 0.32 and log K_int(>CO₃Co⁺) = −2.34 ± 0.61, both from the single ion experiment.

2. Both fall close to the constants measured for the same reaction on the same site, 0.39 and 0.34 log units from magnesium, which is where the near equality of the Li⁺, Co²⁺ and Mg²⁺ ionic radii places them, and nothing in the fit constrained them to.

3. Produced water does not yield transferable constants. Three of six lithium points exceed a monolayer, five of six cobalt points are supersaturated with sphaerocobaltite, and calcium and magnesium hold 96 per cent of the sites.

4. The greater recovery achieved from produced water is not stronger sorption. It is carbonate mineral formation for cobalt and uptake beyond a monolayer for lithium. This answers the question left open by Elshebli et al. (2025).

5. Reaction 5 of Pokrovsky et al. Table 3 is positive, +24.0 and +23.5. Under the negative reading >MeHCO₃⁰ vanishes from the surface and the metal sites are speciated wrongly at every pH.

6. Setting the carbonate by equilibrium with atmospheric pCO₂, as the source method specifies, is preferable to the strontianite reference used in an earlier draft of this work: the two agree to 0.07 log units at pH 7.86, but the strontium route is valid only above about pH 7. Below that it implies bicarbonate concentrations no solution can hold, because strontianite equilibrium fixes a(CO₃²⁻) almost independently of pH. None of the reported constants depends on it, but the bound should be stated by anyone reusing the method.

---

## Tables

**Table 1.** Surface reactions and intrinsic constants used as fixed parameters, from Pokrovsky et al. (1999) Table 3, dolomite columns. Only the two metal reactions are adjustable.

| # | Reaction on the dolomite surface | log K, Ca site | log K, Mg site | Uncertainty |
|---|---|---|---|---|
| 1 | >CO₃H⁰ = >CO₃⁻ + H⁺ | −4.8 | −4.8 | ±0.2 |
| 2 | >CO₃H⁰ + Me²⁺ = >CO₃Me⁺ + H⁺ | −1.8 | −2.0 | ±0.2 |
| 3 | >MeOH⁰ = >MeO⁻ + H⁺ | −12 | −12 | ±2 |
| 4 | >MeOH⁰ + H⁺ = >MeOH₂⁺ | +11.5 | +10.6 | ±0.2 |
| 5 | >MeOH⁰ + CO₃²⁻ + 2H⁺ = >MeHCO₃⁰ + H₂O | +24.0 | +23.5 | ±0.5 |
| 6 | >MeOH⁰ + CO₃²⁻ + H⁺ = >MeCO₃⁻ + H₂O | +16.6 | +15.4 | ±0.2 |
| — | >CO₃H⁰ + Li⁺ = >CO₃Li⁰ + H⁺ | fitted | — | this work |
| — | >CO₃H⁰ + Co²⁺ = >CO₃Co⁺ + H⁺ | fitted | — | this work |

**Table 2.** Constants obtained. 25 °C, carbonate site density 14 µmol m⁻² (0.638 mmol L⁻¹ here), C = √I/0.004, Davies activity coefficients.

| Reaction | System | log K_int | SD | Usable samplings | Status |
|---|---|---|---|---|---|
| >CO₃H⁰ + Li⁺ = >CO₃Li⁰ + H⁺ | single ion, 40 g L⁻¹ NaCl | **−2.39** | 0.32 | 6 of 6 | determined |
| >CO₃H⁰ + Co²⁺ = >CO₃Co⁺ + H⁺ | single ion, 40 g L⁻¹ NaCl | **−2.34** | 0.61 | 5 of 6 | determined |
| >CO₃H⁰ + Li⁺ = >CO₃Li⁰ + H⁺ | produced water | −1.13 | 0.63 | 3 of 6 | upper bound, not transferable |
| >CO₃H⁰ + Co²⁺ = >CO₃Co⁺ + H⁺ | produced water | +0.65 | — | 1 of 6 | upper bound, not transferable |

For comparison, on the same site of the same mineral, Pokrovsky et al. (1999) give −1.8 for Ca²⁺ and −2.0 for Mg²⁺.

**Table 3.** Why produced water does not yield a transferable constant.

| Criterion | Single ion | Produced water |
|---|---|---|
| Carbonate sites held by Ca and Mg before the trace metal | 6 % | 96 % |
| Cobalt samplings undersaturated with sphaerocobaltite | 5 of 6 | 1 of 6 |
| Lithium samplings inside one monolayer | 6 of 6 | 3 of 6 |
| Maximum lithium coverage | 0.85 monolayer | 1.48 monolayers |

**Table 4.** Sensitivity of the single ion constants.

| Parameter varied | Range | Δ log K, Li | Δ log K, Co |
|---|---|---|---|
| Carbonate site density | 14 to 20 µmol m⁻² | 0.36 | 0.20 |
| Capacitance, α halved to doubled | 103 to 414 F m⁻² | 0.03 | 0.09 |

---

## Figure captions

**Figure 1.** The fitted constant put back into the model and compared with the sorption profiles it was not shown, for both metals in both systems. Symbols are measured removal, filled for the batches started at pH 6 and open for those started at pH 2; lines are the model at the constant given in each panel title. The curve is not a kinetic trajectory: it is the equilibrium sorption the model predicts at the pH, calcium and magnesium measured at that moment, interpolated across the six days, so it does not start at zero and a measured point lying on it is a batch that has reached equilibrium. Grey rings mark the samplings excluded before fitting, either because the metal's carbonate had become supersaturated or because the uptake exceeded one monolayer. Panels A and B, the single ion system, track the measured profiles to a root mean square error of 1.3 percentage points; panels C and D, produced water, miss by 25 points, because most of the removal there is not surface complexation and no value of the constant can reproduce it. Panel B shows the largest cobalt uptake of the single ion series lying well above the model and ringed: that is the sampling at which sphaerocobaltite becomes supersaturated.

**Figure 2.** The constant inverted at each sampling, against the measured pH, on one shared scale. Filled symbols are the single ion experiment and open symbols produced water; the line and shaded band are the single ion mean and standard deviation; the dashed and dotted lines are the constants Pokrovsky et al. (1999) measured for calcium and magnesium on the same site. The single ion determinations scatter about the band across four pH units with no trend, which is what an intrinsic constant should do; the produced water points sit well above it.

**Figure 3.** Removal against time for both metals in both systems and both starting conditions, as measured. Produced water removes far more of each metal, which §3.4 and §3.5 show is not stronger sorption.

**Figure 4.** Why the two systems differ. Left: occupancy of the carbonate site before the trace metal arrives, showing calcium and magnesium holding 96 per cent of it in produced water against 6 per cent in the single ion system. Right: the constants obtained in each system, against the calcium and magnesium values for the same site.

---

## Data and code availability

Every datum, reaction, constant and equation is embedded in `LogK_LiCo_singleion_and_PW.ipynb`, which reproduces every number, table and figure in this manuscript from a clean kernel with no external file dependency.

---

## References

Brady, P. V., Papenguth, H. W., & Kelly, J. W. (1999). Metal sorption to dolomite surfaces. *Applied Geochemistry, 14*(5), 569–579. https://doi.org/10.1016/S0883-2927(98)00085-7

Ebrahimi, P., & Vilcáez, J. (2018). Effect of brine salinity and guar gum on the transport of barium through dolomite rocks: Implications for unconventional oil and gas wastewater disposal. *Journal of Environmental Management, 214*, 370–378. https://doi.org/10.1016/j.jenvman.2018.03.008

Elshebli, M., Vilcáez, J., & Smay, J. (2025). Lithium and cobalt recovery from petroleum produced water using dolomite: Impact of pH, temperature, and surface area. *Science of the Total Environment, 1003*, 180743. https://doi.org/10.1016/j.scitotenv.2025.180743

Ioannou, A., & Dimirkou, A. (1997). Phosphate adsorption on hematite, kaolinite, and kaolinite–hematite (k–h) systems as described by a constant capacitance model. *Journal of Colloid and Interface Science, 192*(1), 119–128. https://doi.org/10.1006/jcis.1997.4970

Pokrovsky, O. S., Schott, J., & Thomas, F. (1999). Dolomite surface speciation and reactivity in aquatic systems. *Geochimica et Cosmochimica Acta, 63*(19–20), 3133–3143. https://doi.org/10.1016/S0016-7037(99)00240-9

Shannon, R. D. (1976). Revised effective ionic radii and systematic studies of interatomic distances in halides and chalcogenides. *Acta Crystallographica A, 32*(5), 751–767. https://doi.org/10.1107/S0567739476001551

Smith, R. M., & Martell, A. E. (2004). *NIST critically selected stability constants of metal complexes database* (Version 8.0, NIST Standard Reference Database 46). National Institute of Standards and Technology.

Zachara, J. M., Cowan, C. E., & Resch, C. T. (1991). Sorption of divalent metals on calcite. *Geochimica et Cosmochimica Acta, 55*(6), 1549–1562. https://doi.org/10.1016/0016-7037(91)90126-Q
