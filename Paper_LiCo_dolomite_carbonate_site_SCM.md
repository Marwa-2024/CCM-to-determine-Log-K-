# Surface complexation of lithium and cobalt on dolomite in a sodium chloride brine: a carbonate-site model, capacity constraints, and the separation of sorption from mineralization

**Marwa Elshebli**¹, Javier Vilcáez¹, James Smay²

¹ Boone Pickens School of Geology, Oklahoma State University, Stillwater, OK 74078, USA
² School of Materials, Mechatronics & Manufacturing Engineering, Oklahoma State University, Tulsa, OK 74106, USA

*Draft prepared from the single-ion batch data and the model workflow in the accompanying notebook (`LogK_LiCo_dolomite_carbonate_site_STEPWISE.ipynb`). All numbers in the Results were produced by that notebook; nothing is read from external files.*

---

## Abstract

Dolomite has been shown to recover Li⁺ and Co²⁺ from petroleum produced water, with Li⁺ removal governed mainly by sorption and Co²⁺ removal by coupled sorption and carbonate mineralization (Elshebli et al., 2025). Here the single-ion batch data (dolomite 60 g L⁻¹, 40 g L⁻¹ NaCl, 25 °C, six days) are interpreted with the three-site constant capacitance surface complexation model of Pokrovsky et al. (1999), in which metal cations bind the carbonate site, >CO₃H⁰ + Meᶻ⁺ ⇌ >CO₃Me⁽ᶻ⁻¹⁾⁺ + H⁺. A capacity test applied before any fitting shows that the high-recovery runs reported previously (≈70 mg L⁻¹ Co removed on 0.84 m² g⁻¹ dolomite) correspond to 1.7 to 2.4 times the crystallographic carbonate-site density and therefore cannot be adsorption, consistent with the zabuyelite and sphaerocobaltite detected by XRD. The present low-uptake data set lies inside the monolayer bound (Co 16 %, Li 72 % of carbonate sites), but the cobalt "pH 6" equilibrium point is supersaturated with sphaerocobaltite (SI = +0.93). Aqueous speciation at the measured equilibrium pH (7.05 to 8.26, not the nominal 2 and 6) shows 67 to 69 % of dissolved cobalt as free Co²⁺ and 27 % as CoCl⁺, and 90 % of lithium as free Li⁺. Solving the carbonate-site mass balance with the measured Ca and Mg as competitors, and closing the constant capacitance electrostatic loop (C = 138 F m⁻²), gives point-by-point constants that are flat in pH for lithium and a fitted log K(>CO₃Li⁰) = −2.42 ± 0.30 (electrostatic) and −2.40 (non-electrostatic), conditional on the borrowed site parameters. For cobalt, a fitted value would be a lumped parameter; instead the constant is taken from the aqueous–surface analogy of Van Cappellen et al. (1993) and Pokrovsky and Schott (2002), log K(>CO₃Co⁺) = −0.79, and used to predict the sorbable fraction. Sorption cannot account for more than 41 % of the cobalt removed in the best previously reported run even if every carbonate site were occupied, which quantifies the coupled sorption–mineralization mechanism proposed earlier. The experimental design required for a genuinely transferable constant is specified.

**Keywords:** dolomite; surface complexation; constant capacitance model; lithium; cobalt; produced water; sphaerocobaltite; zabuyelite

---

## 1. Introduction

Lithium and cobalt are strategic metals whose recovery from saline produced water using low-cost carbonate minerals has attracted interest (Elshebli et al., 2025). Interpreting such uptake quantitatively requires a surface complexation model (SCM) that accounts for aqueous speciation, competition for surface sites, and the electrostatics of the mineral–water interface (Dzombak & Morel, 1990; Stumm & Morgan, 1996). For carbonate minerals the appropriate framework is the hydration-site model developed for calcite by Van Cappellen et al. (1993) and extended to magnesite and dolomite by Pokrovsky et al. (1999a, 1999b), in which three primary surface groups, >CO₃H⁰, >CaOH⁰ and >MgOH⁰, undergo protonation, deprotonation and complexation with the lattice ions and dissolved carbonate. In this model metal cations adsorb by forming >CO₃Me complexes on the deprotonated carbonate site, whereas anions and ligands bind the metal hydroxyl sites (Pokrovsky et al., 1999a). Belova et al. (2014) applied the same approach to nickel adsorption on calcite and chalk, determining the >CO₃Ni⁺ constant with FITEQL (Herbelin & Westall, 1999) and demonstrating the importance of preventing surface precipitation during the experiment.

An intrinsic stability constant describes a monolayer reaction at an interface. Where a large share of the metal has left solution as a discrete carbonate phase, a constant fitted to that removal is not an interfacial property but a lumped parameter that absorbs the solubility product, the carbonate supply from dolomite dissolution and the nucleation behaviour of the precipitate. The earlier study established that cobalt removal by dolomite is coupled sorption and mineralization and left the quantitative split as an open question (Elshebli et al., 2025). The purpose of this work is therefore twofold: to specify the carbonate-site SCM for lithium and cobalt on dolomite in brine, with every reaction, equation and parameter stated and sourced, and to determine what the existing data can and cannot constrain, using a capacity test and a saturation test before any fitting, and using the model in prediction mode to quantify the sorption–mineralization split.

---

## 2. Materials and methods

### 2.1 Dolomite and batch experiments

Dolomite from the Arbuckle Group (98 % dolomite by XRD) was ground and sieved; the single-ion batches used 6 g of dolomite per 100 mL of solution (60 g L⁻¹) in 40 g L⁻¹ NaCl at 25 °C, shaken at 150 rpm, with a nominal initial pH of 2 or 6 (Elshebli et al., 2025). Lithium (138.4 mg L⁻¹ from LiCl) or cobalt (80.5 mg L⁻¹ from CoCl₂·6H₂O) was added as the only trace metal. Aliquots were taken at 2, 4 and 6 days and analysed by ICP-OES for Li, Co, Ca and Mg; pH was recorded at each sampling. The specific surface area of the powder used here is 0.76 m² g⁻¹ (mean particle diameter 2.8 µm, density 2850 kg m⁻³), giving a surface concentration

S = 60 g L⁻¹ × 0.76 m² g⁻¹ = 45.6 m² L⁻¹.   (1)

The value of 0.046 to 0.049 m² L⁻¹ quoted in the methods of Elshebli et al. (2025) is inconsistent with 0.84 m² g⁻¹ × 60 g L⁻¹ ≈ 50 m² L⁻¹ by a factor of 10³ and is most likely m² mL⁻¹; every surface normalization in an SCM divides by this quantity, so the corrected value is used throughout.

### 2.2 Experimental data

Table 1 gives the embedded data set. The measured pH at day 6 (7.05 to 8.26) differs from the nominal initial pH because dolomite dissolution buffers the suspension; the measured value is the one that enters the model.

**Table 1.** Single-ion batch data (dolomite 60 g L⁻¹, 40 g L⁻¹ NaCl, 25 °C).

| Metal | Batch | Day | Ca (mg L⁻¹) | Mg (mg L⁻¹) | Metal (mg L⁻¹) | pH (measured) | Removal (%) |
|---|---|---|---|---|---|---|---|
| Co | pH 6 | 0 | 0.10 | 0.10 | 80.513 | 6.32 | 0 |
| Co | pH 6 | 2 | 0.10 | 0.10 | 77.905 | 7.40 | 3.24 |
| Co | pH 6 | 4 | 0.10 | 8.14 | 76.411 | 7.48 | 5.09 |
| Co | pH 6 | 6 | 0.13 | 12.94 | 74.393 | 8.09 | 7.60 |
| Co | pH 2 | 0 | 0.10 | 0.10 | 80.513 | 2.35 | 0 |
| Co | pH 2 | 2 | 57.48 | 17.46 | 79.661 | 4.12 | 1.06 |
| Co | pH 2 | 4 | 24.18 | 45.59 | 79.795 | 5.66 | 0.89 |
| Co | pH 2 | 6 | 21.66 | 51.38 | 80.398 | 7.11 | 0.14 |
| Li | pH 6 | 0 | 0.66 | 0.10 | 138.430 | 6.42 | 0 |
| Li | pH 6 | 2 | 37.63 | 0.10 | 134.657 | 7.65 | 2.73 |
| Li | pH 6 | 4 | 54.82 | 0.10 | 136.722 | 7.92 | 1.23 |
| Li | pH 6 | 6 | 49.89 | 0.87 | 135.262 | 8.26 | 2.29 |
| Li | pH 2 | 0 | 0.66 | 0.10 | 138.430 | 2.25 | 0 |
| Li | pH 2 | 2 | 39.19 | 0.10 | 136.263 | 4.26 | 1.57 |
| Li | pH 2 | 4 | 79.45 | 21.32 | 135.956 | 5.73 | 1.79 |
| Li | pH 2 | 6 | 87.74 | 20.55 | 135.592 | 7.05 | 2.05 |

### 2.3 Aqueous speciation

Free-ion activities were computed at each data point from the measured pH, the measured total Ca and Mg, the NaCl medium and the measured total metal, at fixed atmospheric pCO₂ = 10⁻³·⁵ atm (open system, as in Pokrovsky et al., 1999a). Activity coefficients follow the Davies equation (Davies, 1962), valid to I ≈ 0.7 M,

log γᵢ = −A zᵢ² [ √I / (1 + √I) − 0.3 I ],  A = 0.509 at 25 °C.   (2)

The carbonate system is fixed by pH and pCO₂:

CO₂(g) ⇌ CO₂(aq), log K_H = −1.469;  CO₂(aq) + H₂O ⇌ HCO₃⁻ + H⁺, log K = −6.345;  CO₃²⁻ + H⁺ ⇌ HCO₃⁻, log K = 10.329   (3)

(Wolery, 1992). Aqueous complexes were written as formation reactions from the free ions with the constants in Table 2 (Smith & Martell, 2004; Wolery, 1992). For each cation X with analytical total T_X the mass balance

T_X = [X] + Σᵢ [Xᵢ-complex]   (4)

was solved by fixed-point iteration together with the chloride balance and the ionic strength. Saturation indices were computed as SI = log(IAP/K_sp) with log K_sp = −8.48 for calcite (Plummer & Busenberg, 1982), −17.09 for ordered dolomite (Parkhurst & Appelo, 2013), −9.98 for sphaerocobaltite CoCO₃ (Smith & Martell, 2004) and ≈ −2.5 for zabuyelite Li₂CO₃.

**Table 2.** Aqueous complexation reactions (formation from free ions, 25 °C, I = 0).

| Reaction | log β | Reaction | log β |
|---|---|---|---|
| Ca²⁺ + Cl⁻ ⇌ CaCl⁺ | −0.70 | Co²⁺ + Cl⁻ ⇌ CoCl⁺ | 0.30 |
| Ca²⁺ + CO₃²⁻ ⇌ CaCO₃⁰ | 3.22 | Co²⁺ + 2Cl⁻ ⇌ CoCl₂⁰ | −0.20 |
| Ca²⁺ + H⁺ + CO₃²⁻ ⇌ CaHCO₃⁺ | 11.43 | Co²⁺ + CO₃²⁻ ⇌ CoCO₃⁰ | 4.23 |
| Ca²⁺ + H₂O ⇌ CaOH⁺ + H⁺ | −12.70 | Co²⁺ + H⁺ + CO₃²⁻ ⇌ CoHCO₃⁺ | 12.20 |
| Mg²⁺ + Cl⁻ ⇌ MgCl⁺ | −0.14 | Co²⁺ + H₂O ⇌ CoOH⁺ + H⁺ | −9.65 |
| Mg²⁺ + CO₃²⁻ ⇌ MgCO₃⁰ | 2.98 | Co²⁺ + 2H₂O ⇌ Co(OH)₂⁰ + 2H⁺ | −18.80 |
| Mg²⁺ + H⁺ + CO₃²⁻ ⇌ MgHCO₃⁺ | 11.40 | Li⁺ + Cl⁻ ⇌ LiCl⁰ | −0.50 |
| Mg²⁺ + H₂O ⇌ MgOH⁺ + H⁺ | −11.79 | Li⁺ + H₂O ⇌ LiOH⁰ + H⁺ | −13.64 |
| Na⁺ + Cl⁻ ⇌ NaCl⁰ | −0.78 | Li⁺ + CO₃²⁻ ⇌ LiCO₃⁻ | 0.90 |
| Na⁺ + CO₃²⁻ ⇌ NaCO₃⁻ | 1.27 | | |
| Na⁺ + H⁺ + CO₃²⁻ ⇌ NaHCO₃⁰ | 10.48 | | |
| Na⁺ + H₂O ⇌ NaOH⁰ + H⁺ | −14.18 | | |

### 2.4 Surface complexation model

The dolomite surface is described by three primary hydration sites in the ratio 1 : 1 : 2, with site densities of 7 µmol m⁻² for >CaOH⁰, 7 µmol m⁻² for >MgOH⁰ and 14 µmol m⁻² for >CO₃H⁰ (Pokrovsky et al., 1999a). The total carbonate-site concentration in suspension is

S_T = 14 × 10⁻⁶ mol m⁻² × 45.6 m² L⁻¹ = 6.38 × 10⁻⁴ mol L⁻¹.   (5)

The surface reactions and their intrinsic constants (Table 3) are those of Pokrovsky et al. (1999a, Table 3). Following the analogy between surface and solution complexation (Schindler & Stumm, 1987; Van Cappellen et al., 1993), metal cations bind the deprotonated carbonate site,

>CO₃H⁰ + Li⁺ ⇌ >CO₃Li⁰ + H⁺,  K_Li  (ΔZ = 0)   (6a)
>CO₃H⁰ + Co²⁺ ⇌ >CO₃Co⁺ + H⁺,  K_Co  (ΔZ = +1)   (6b)

These are the two unknowns. Lithium forms a neutral surface species and therefore carries no Boltzmann term; cobalt forms a +1 species and does.

**Table 3.** Surface reactions and intrinsic constants (25 °C, I = 0) used as fixed parameters.

| Reaction | log K_int | Source |
|---|---|---|
| >CO₃H⁰ ⇌ >CO₃⁻ + H⁺ | −4.8 | Pokrovsky et al. (1999a) |
| >CO₃H⁰ + Ca²⁺ ⇌ >CO₃Ca⁺ + H⁺ | −1.8 | Pokrovsky et al. (1999a) |
| >CO₃H⁰ + Mg²⁺ ⇌ >CO₃Mg⁺ + H⁺ | −2.0 | Pokrovsky et al. (1999a) |
| >CaOH⁰ + H⁺ ⇌ >CaOH₂⁺ | 11.5 | Pokrovsky et al. (1999a) |
| >CaOH⁰ ⇌ >CaO⁻ + H⁺ | −12.0 | Pokrovsky et al. (1999a) |
| >CaOH⁰ + CO₃²⁻ + 2H⁺ ⇌ >CaHCO₃⁰ + H₂O | −4.0 | Pokrovsky et al. (1999a) |
| >CaOH⁰ + CO₃²⁻ + H⁺ ⇌ >CaCO₃⁻ + H₂O | 16.6 | Pokrovsky et al. (1999a) |
| >MgOH⁰ + H⁺ ⇌ >MgOH₂⁺ | 10.6 | Pokrovsky et al. (1999a) |
| >MgOH⁰ ⇌ >MgO⁻ + H⁺ | −12.0 | Pokrovsky et al. (1999a) |
| >MgOH⁰ + CO₃²⁻ + 2H⁺ ⇌ >MgHCO₃⁰ + H₂O | −3.5 | Pokrovsky et al. (1999a) |
| >MgOH⁰ + CO₃²⁻ + H⁺ ⇌ >MgCO₃⁻ + H₂O | 15.4 | Pokrovsky et al. (1999a) |
| Site densities (Ca : Mg : CO₃) | 7 : 7 : 14 µmol m⁻² | Pokrovsky et al. (1999a) |
| Capacitance C = √I / α, α = 0.006 | 138 F m⁻² at I = 0.684 M | Pokrovsky & Schott (2002); Belova et al. (2014) |

The apparent and intrinsic constants are related through the Boltzmann factor (Pokrovsky et al., 1999a, Eq. 2; Gustafsson, 2014, Eq. 4.4),

K_int = K_app · exp(ΔZ F ψ₀ / RT),   (7)

where ΔZ is the change in surface charge in the reaction, ψ₀ the surface potential, F the Faraday constant, R the gas constant and T the temperature. In concentration form a surface species of charge z carries the factor exp(−z F ψ₀/RT); for example

[>CO₃Ca⁺] = K_Ca [>CO₃H⁰] {Ca²⁺}/{H⁺} · exp(−F ψ₀/RT).   (8)

The carbonate-site mass balance,

S_T = [>CO₃H⁰] + [>CO₃⁻] + [>CO₃Ca⁺] + [>CO₃Mg⁺] + [>CO₃Li⁰] + [>CO₃Co⁺],   (9)

collapses to one equation in [>CO₃H⁰] because every other term is proportional to it. Calcium and magnesium are treated as competitors whose concentrations follow the measured solution chemistry at each point. The two metal hydroxyl sites are solved likewise because they carry surface charge. The net surface charge is computed from the charged surface species only (Charlet et al., 1990; Turner & Fein, 2006, Eq. 2.2.3),

σ = (F / S) Σᵢ zᵢ [iₛᵤᵣf] = (F/S) { [>CaOH₂⁺] + [>MgOH₂⁺] + [>CO₃Ca⁺] + [>CO₃Mg⁺] + [>CO₃Co⁺] − [>CO₃⁻] − [>CaO⁻] − [>MgO⁻] − [>CaCO₃⁻] − [>MgCO₃⁻] },   (10)

and related to the potential by the constant capacitance model (Pokrovsky et al., 1999a, Eqs. 3–4; Gustafsson, 2014, Eq. 4.11),

ψ₀ = σ / C,  C = √I / α.   (11)

Equations 9 to 11 are implicit and were solved by iteration on ψ₀ until it changed by less than 0.1 mV. A non-electrostatic model (ψ₀ = 0) was run in parallel. Note that Eq. 10 never contains dissolved Ca²⁺ or Mg²⁺: the mass-balance surface charge of Charlet et al. (1990) applied to a dissolving mineral carries the dissolution flux and, in the present titration data, reaches 1.7 mmol m⁻², sixty times the site density, so it cannot be used as a surface charge.

### 2.5 Surface coverage and its uncertainty

The measured coverage at each point is

Γ = (C₀ − C_eq) V / (m · SSA)  (mol m⁻²),   (12)

with V = 0.1 L, m = 6 g and SSA = 0.76 m² g⁻¹. Taking the ICP-OES precision as ±3 % on both C₀ and C_eq, the uncertainty of the difference is σ_Γ = √[(0.03 C₀)² + (0.03 C_eq)²] V/(m·SSA), which is large relative to Γ when removal is small; such points were down-weighted, not discarded. No solid-free blanks were available, so wall losses could not be subtracted.

### 2.6 Capacity and saturation criteria

Two tests precede any fitting. The capacity test compares Γ with the crystallographic site density: Γ/N_s > 1 is not reachable by adsorption. The saturation test requires SI < 0 for the metal carbonate at the point being fitted (Belova et al., 2014).

### 2.7 Parameter estimation

Two procedures were used. Point by point, with Γ known, Eq. 6 is solved algebraically at each datum,

K_Co = ( [>CO₃Co⁺] {H⁺} ) / ( [>CO₃H⁰] {Co²⁺} ) · exp(F ψ₀ / RT),   (13)

with [>CO₃Co⁺] = Γ·S and [>CO₃H⁰] from Eq. 9, giving one constant per point without an optimizer. Globally, the forward model predicts Γ for a trial log K and the weighted sum of squares over degrees of freedom,

WSOS/DF = Σᵢ [ (Γ_pred,ᵢ − Γ_meas,ᵢ) / σ_Γ,ᵢ ]² / (N − p),   (14)

was minimized (Westall, 1982; Herbelin & Westall, 1999), as in the FITEQL treatment of Belova et al. (2014) with a dummy adsorbed component. The diagnostic of Step 6 is the plot of point-by-point log K against pH and against Γ: a constant independent of both is a transferable constant, drift with pH signals a wrong surface species or a missing coulombic term, and drift with loading signals surface precipitation.

### 2.8 Prediction mode and the mechanism split

For cobalt the constant was taken from the linear free-energy analogy between surface and aqueous carbonate complexes (Van Cappellen et al., 1993; Pokrovsky & Schott, 2002),

log K(>CO₃Co⁺) ≈ log K(>CO₃Ca⁺) + [log β(CoCO₃⁰) − log β(CaCO₃⁰)] = −1.8 + (4.23 − 3.22) = −0.79,   (15)

and the model was run forward with all constants fixed to predict the sorbable cobalt at each time point. The difference between measured removal and predicted sorption, read together with SI(CoCO₃), is the mineralized fraction. Lithium, whose removal the earlier study attributed to sorption, was treated as the single adjustable parameter.

All calculations were performed in Python 3 (NumPy, SciPy, pandas, Matplotlib); the notebook reproduces every table and figure.

---

## 3. Results

### 3.1 Capacity test

**Table 4.** Surface coverage relative to the carbonate-site density.

| Case | Removed (mg L⁻¹) | Area (m² L⁻¹) | Γ (µmol m⁻²) | Γ / 14 µmol m⁻² | Γ / 8.22 µmol m⁻² |
|---|---|---|---|---|---|
| This data, Co, pH 6 batch, day 6 | 6.12 | 45.6 | 2.28 | 0.16 | 0.28 |
| This data, Li, pH 6 batch, day 6 | 3.17 | 45.6 | 10.0 | 0.72 | 1.22 |
| Elshebli et al. (2025), Co, fine dolomite, pH 6, 25 °C | ≈70 | 50.4 | 23.6 | 1.68 | 2.87 |
| Elshebli et al. (2025), Li, fine dolomite, 8.5 % of 138 mg L⁻¹ | 11.8 | 50.4 | 33.6 | 2.40 | 4.09 |

The high-recovery runs of the earlier study exceed the carbonate-site density by a factor of 1.7 to 2.4 (2.9 to 4.1 relative to the calcite value of Belova et al., 2014) even if every site were occupied by the metal alone. They therefore cannot be adsorption, in agreement with the 8 % zabuyelite and 4 % sphaerocobaltite found by XRD. The present low-uptake data set lies within the bound for cobalt (16 %) and at the margin for lithium (72 %).

### 3.2 Aqueous speciation and saturation

At the measured equilibrium pH, 67 to 69 % of dissolved cobalt is free Co²⁺, 27 % is CoCl⁺ and 3 % CoCl₂⁰; carbonate and hydroxide complexes are below 1.5 %. Lithium is 90 % free Li⁺ and 10 % LiCl⁰. Ionic strength is 0.65 to 0.66 M. The saturation index of sphaerocobaltite in the cobalt pH 6 batch rises from −0.42 (day 2) to −0.27 (day 4) to **+0.93 at day 6**; the pH 2 batch stays undersaturated (−6.96 to −0.98). Calcite approaches saturation in the lithium pH 6 batch at day 6 (SI = −0.08); dolomite remains undersaturated at all points (SI ≤ −1.6); Li₂CO₃ is far undersaturated (SI ≤ −6.4).

### 3.3 Coverage and site occupancy

Cobalt coverage at day 6 is 2.28 µmol m⁻² in the pH 6 batch (relative error 54 %) and 0.04 µmol m⁻² in the pH 2 batch (relative error > 1000 %, i.e. within noise). Lithium coverage is 9 to 10 µmol m⁻² in both batches (relative error ≈ 200 % at ±3 % ICP precision). At the equilibrium points the carbonate site is dominated by >CO₃⁻ (74 to 77 % in the cobalt batches, 21 to 25 % in the lithium batches), with >CO₃Ca⁺ at 0.1 to 8.7 %, >CO₃Mg⁺ at 0.1 to 18 %, and the metal complex at 16 % (Co, pH 6) and 64 to 72 % (Li). The surface potential is −0.1 to +8 mV and the Boltzmann factor 0.73 to 1.0.

### 3.4 Stability constants

**Table 5.** Point-by-point log K on the carbonate site (constant capacitance model).

| Metal | Batch | Day | pH | Γ (µmol m⁻²) | log K (CCM) | log K (NEM) | Note |
|---|---|---|---|---|---|---|---|
| Co | pH 6 | 2 | 7.40 | 0.97 | −2.34 | −2.37 | pre-equilibrium |
| Co | pH 6 | 4 | 7.48 | 1.53 | −2.05 | −2.13 | pre-equilibrium |
| Co | pH 6 | 6 | 8.09 | 2.28 | −1.89 | −1.89 | SI(CoCO₃) = +0.93 |
| Co | pH 2 | 2 | 4.12 | 0.32 | −1.94 | −2.10 | within noise |
| Co | pH 2 | 4 | 5.66 | 0.27 | −2.63 | −2.76 | within noise |
| Co | pH 2 | 6 | 7.11 | 0.04 | −3.47 | −3.59 | within noise |
| Li | pH 6 | 2 | 7.65 | 11.9 | −1.97 | −2.05 | |
| Li | pH 6 | 4 | 7.92 | 5.4 | −2.95 | −2.97 | |
| Li | pH 6 | 6 | 8.26 | 10.0 | −2.37 | −2.38 | equilibrium |
| Li | pH 2 | 2 | 4.26 | 6.9 | −2.23 | −2.26 | |
| Li | pH 2 | 4 | 5.73 | 7.8 | −2.51 | −2.55 | |
| Li | pH 2 | 6 | 7.05 | 9.0 | −2.37 | −2.41 | equilibrium |

**Table 6.** Reported constants (25 °C, I = 0.68 M NaCl).

| Reaction | Electrostatic (CCM) | Non-electrostatic | Basis |
|---|---|---|---|
| >CO₃H⁰ + Li⁺ ⇌ >CO₃Li⁰ + H⁺ | −2.42 ± 0.30 (all six points); −2.37 (day 6) | −2.40 | fitted, conditional on Table 3 site parameters |
| >CO₃H⁰ + Co²⁺ ⇌ >CO₃Co⁺ + H⁺ | −0.79 | −0.79 | predicted from Eq. 15 |
| >CO₃H⁰ + Co²⁺ ⇌ >CO₃Co⁺ + H⁺ (apparent, this data) | −2.06 ± 0.79 | −2.13 | lumped; not an intrinsic constant |

The lithium point-by-point values are flat in pH (slope 0.00 per pH unit across the equilibrium points, spread 0.30 log units across all six points), and the electrostatic and non-electrostatic values agree to 0.02 log units, as expected at C = 138 F m⁻². The cobalt values drift with pH (slope +1.6 per pH unit), driven by the pH 2 batch whose uptake is within analytical noise, and the only cobalt point with a measurable coverage is supersaturated with sphaerocobaltite. A fitted cobalt constant is therefore reported only as an apparent value.

### 3.5 Prediction and mechanism split

With the constants of Table 3 and log K(>CO₃Co⁺) = −0.79, the forward model predicts that sorption alone would remove 21 to 27 % of the added cobalt at the measured equilibrium chemistry, whereas 0.1 to 7.6 % was removed. In this coarse, low-surface-area system the observed removal is therefore *below* what the literature analogy allows for sorption, and no mineralization is required to explain it; the apparent constant of −2.06 is 1.3 log units weaker than the analogy, comparable to the >CO₃Ca⁺ and >CO₃Mg⁺ values of −1.8 and −2.0. For the high-recovery runs of Elshebli et al. (2025), the absolute sorption ceiling at 50.4 m² L⁻¹ is 14 µmol m⁻² × 50.4 m² L⁻¹ = 0.71 mmol L⁻¹, against 1.19 mmol L⁻¹ of cobalt removed: at least 41 % of that removal is mineralization even if every carbonate site held cobalt, and a much larger share once protons, calcium and magnesium occupy their share of sites.

---

## 4. Discussion

The capacity test settles the question left open in the earlier study. Removal of 70 mg L⁻¹ of cobalt or 8.5 % of 138 mg L⁻¹ of lithium on 50 m² L⁻¹ of dolomite is one to three monolayers beyond what any surface complexation model can hold, so those recoveries are dominated by precipitation of sphaerocobaltite and zabuyelite, as the XRD already indicated. A constant fitted to such data would change with pH, temperature, loading and duration and would not be an interfacial property. The present single-ion data set is different: its uptake is small enough to sit inside the monolayer bound, and for lithium the point-by-point constant is independent of pH and of the electrostatic model, which is the behaviour of a transferable constant. Its uncertainty is nevertheless large because lithium uptake (2 to 3 %) is within about one standard deviation of the ICP measurement, and the value is conditional on the borrowed site density and competitor constants; Belova et al. (2014) reported ±0.45 to ±0.66 log units for nickel on calcite under far better controlled conditions.

For cobalt the analogy of Eq. 15 predicts a constant about one log unit stronger than that of calcium, consistent with the higher stability of CoCO₃⁰ over CaCO₃⁰. The observed uptake on the coarse dolomite is weaker than this prediction. Three explanations are open and cannot be separated with the present design: a reactive area smaller than the BET area (Belova et al., 2014, found that assuming 70 % of the BET area raised the fitted constant), a slower approach to equilibrium on the coarse solid, or competition by calcium and magnesium at the concentrations released by dissolution (up to 88 and 51 mg L⁻¹). The saturation history shows that the cobalt pH 6 batch crossed sphaerocobaltite saturation between day 4 and day 6, so the day 6 point is at the onset of precipitation and the earlier points, still undersaturated, give the cleaner estimate of sorption.

Three features of the experimental design limit what can be fitted irrespective of precipitation: two initial pH values rather than an edge, so the proton stoichiometry is not constrained; pH that drifts during the run, so each point is not an equilibrium at a defined pH; and a single initial concentration, so there is no isotherm and site density cannot be separated from binding constant even in principle. A genuinely transferable constant requires pre-equilibrating the dolomite with the brine until the pH is stable, filtering and using the solution without storage, adding metal at concentrations kept undersaturated with respect to CoCO₃ and Li₂CO₃ (micrograms per litre for cobalt), holding pH at six or more values across the uptake edge, varying the initial concentration at three or four levels at each pH, limiting contact to 24 h to stay in the fast adsorption regime, and running solid-free blanks at every condition (Belova et al., 2014).

---

## 5. Conclusions

1. The carbonate-site surface complexation model of Pokrovsky et al. (1999a), with metal cations on >CO₃⁻ and anions on >MeOH sites, is the appropriate framework for lithium and cobalt on dolomite; all reactions, equations and parameters are stated and sourced above.
2. A capacity test must precede any fitting. The high-recovery runs of Elshebli et al. (2025) exceed the carbonate-site density by 1.7 to 2.4 times and are dominated by mineralization; at least 41 % of the best cobalt recovery is precipitation even at the sorption ceiling.
3. For the present low-uptake data set, lithium behaves as a sorbing ion with a pH-independent constant, log K(>CO₃Li⁰) = −2.42 ± 0.30 (electrostatic) and −2.40 (non-electrostatic), conditional on the borrowed site parameters.
4. Cobalt cannot be fitted reliably from these data because the measurable point is supersaturated with sphaerocobaltite; its constant is taken from the aqueous–surface analogy, log K(>CO₃Co⁺) = −0.79, and the apparent value from the data (−2.06) is reported only as a lumped parameter.
5. The published surface concentration of 0.046 to 0.049 m² L⁻¹ should be corrected to ≈ 50 m² L⁻¹ (m² mL⁻¹ in the original).
6. The experimental design for a transferable constant is specified and is a matter of weeks rather than months.

---

## Data and code availability

The complete workflow, with every reaction, constant and datum embedded, is in `LogK_LiCo_dolomite_carbonate_site_STEPWISE.ipynb` (Python 3) in the project repository. A Wolfram Mathematica notebook implementing the same model is provided alongside.

---

## References

Belova, D. A., Lakshtanov, L. Z., Carneiro, J. F., & Stipp, S. L. S. (2014). Nickel adsorption on chalk and calcite. *Journal of Contaminant Hydrology, 170*, 1–9. https://doi.org/10.1016/j.jconhyd.2014.09.007

Charlet, L., Wersin, P., & Stumm, W. (1990). Surface charge of MnCO₃ and FeCO₃. *Geochimica et Cosmochimica Acta, 54*(8), 2329–2336. https://doi.org/10.1016/0016-7037(90)90059-T

Davies, C. W. (1962). *Ion association*. Butterworths.

Dzombak, D. A., & Morel, F. M. M. (1990). *Surface complexation modeling: Hydrous ferric oxide*. Wiley.

Elshebli, M., Vilcáez, J., & Smay, J. (2025). Lithium and cobalt recovery from petroleum produced water using dolomite: Impact of pH, temperature, and surface area. *Science of the Total Environment, 1003*, 180743. https://doi.org/10.1016/j.scitotenv.2025.180743

Gustafsson, J. P. (2014). *Visual MINTEQ 3.1 user guide*. KTH Royal Institute of Technology.

Herbelin, A. L., & Westall, J. C. (1999). *FITEQL 4.0: A computer program for determination of chemical equilibrium constants from experimental data* (Report 99-01). Department of Chemistry, Oregon State University.

Parkhurst, D. L., & Appelo, C. A. J. (2013). *Description of input and examples for PHREEQC version 3: A computer program for speciation, batch-reaction, one-dimensional transport, and inverse geochemical calculations* (Techniques and Methods 6-A43). U.S. Geological Survey. https://doi.org/10.3133/tm6A43

Plummer, L. N., & Busenberg, E. (1982). The solubilities of calcite, aragonite and vaterite in CO₂-H₂O solutions between 0 and 90 °C, and an evaluation of the aqueous model for the system CaCO₃-CO₂-H₂O. *Geochimica et Cosmochimica Acta, 46*(6), 1011–1040. https://doi.org/10.1016/0016-7037(82)90056-4

Pokrovsky, O. S., & Schott, J. (2002). Surface chemistry and dissolution kinetics of divalent metal carbonates. *Environmental Science & Technology, 36*(3), 426–432. https://doi.org/10.1021/es010925u

Pokrovsky, O. S., Schott, J., & Thomas, F. (1999a). Dolomite surface speciation and reactivity in aquatic systems. *Geochimica et Cosmochimica Acta, 63*(19–20), 3133–3143. https://doi.org/10.1016/S0016-7037(99)00240-9

Pokrovsky, O. S., Schott, J., & Thomas, F. (1999b). Processes at the magnesium-bearing carbonates/solution interface. I. A surface speciation model for magnesite. *Geochimica et Cosmochimica Acta, 63*(6), 863–880. https://doi.org/10.1016/S0016-7037(99)00008-3

Schindler, P. W., & Stumm, W. (1987). The surface chemistry of oxides, hydroxides, and oxide minerals. In W. Stumm (Ed.), *Aquatic surface chemistry* (pp. 83–110). Wiley.

Smith, R. M., & Martell, A. E. (2004). *NIST critically selected stability constants of metal complexes database* (Version 8.0, NIST Standard Reference Database 46). National Institute of Standards and Technology.

Stumm, W., & Morgan, J. J. (1996). *Aquatic chemistry: Chemical equilibria and rates in natural waters* (3rd ed.). Wiley.

Turner, B. F., & Fein, J. B. (2006). Protofit: A program for determining surface protonation constants from titration data. *Computers & Geosciences, 32*(9), 1344–1356. https://doi.org/10.1016/j.cageo.2005.12.005

Van Cappellen, P., Charlet, L., Stumm, W., & Wersin, P. (1993). A surface complexation model of the carbonate mineral-aqueous solution interface. *Geochimica et Cosmochimica Acta, 57*(15), 3505–3518. https://doi.org/10.1016/0016-7037(93)90135-J

Westall, J. C. (1982). *FITEQL: A computer program for determination of chemical equilibrium constants from experimental data* (Report 82-01). Department of Chemistry, Oregon State University.

Wolery, T. J. (1992). *EQ3/6, a software package for geochemical modeling of aqueous systems: Package overview and installation guide* (UCRL-MA-110662 PT I). Lawrence Livermore National Laboratory.
