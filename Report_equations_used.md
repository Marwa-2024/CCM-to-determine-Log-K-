# Determination of intrinsic surface complexation constants for Li⁺ and Co²⁺ on dolomite
## Report: every equation used, where it comes from, and where it enters the calculation

**Marwa Elshebli**¹, Javier Vilcáez¹, James Smay²

¹ Boone Pickens School of Geology, Oklahoma State University, Stillwater, OK 74078, USA
² School of Materials, Mechatronics & Manufacturing Engineering, Oklahoma State University, Tulsa, OK 74106, USA

*Every equation below is implemented in `LogK_LiCo_TST_full_method.ipynb`, which reproduces every number in this report from a clean kernel. Nothing is quoted that is not executed.*

---

## 1. Scope

This report lists the complete set of equations used to obtain intrinsic surface complexation constants for lithium and cobalt on dolomite from the batch experiments of this study, in two systems: the metal alone in 40 g L⁻¹ NaCl, and the same metals in produced water. The equations are grouped by what they do. For each, the report gives the form as implemented, its source, the parameter values used, and the place in the calculation where it acts.

Two departures from the source formulations are documented in §8, with the evidence for each.

---

## 2. Surface complexation reactions and their mass action expressions

### 2.1 The reaction set

Six reactions on each of the three hydration sites of dolomite, from Pokrovsky, Schott and Thomas (1999) Table 3, dolomite columns, held fixed throughout. Me denotes Ca or Mg.

| # | Reaction | log K_int, Ca site | log K_int, Mg site |
|---|---|---|---|
| E1 | >CO₃H⁰ ⇌ >CO₃⁻ + H⁺ | −4.8 ± 0.2 | −4.8 ± 0.2 |
| E2 | >CO₃H⁰ + Me²⁺ ⇌ >CO₃Me⁺ + H⁺ | −1.8 ± 0.2 | −2.0 ± 0.2 |
| E3 | >MeOH⁰ ⇌ >MeO⁻ + H⁺ | −12 ± 2 | −12 ± 2 |
| E4 | >MeOH⁰ + H⁺ ⇌ >MeOH₂⁺ | 11.5 ± 0.2 | 10.6 ± 0.2 |
| E5 | >MeOH⁰ + CO₃²⁻ + 2H⁺ ⇌ >MeHCO₃⁰ + H₂O | 24.0 ± 0.5 | 23.5 ± 0.5 |
| E6 | >MeOH⁰ + CO₃²⁻ + H⁺ ⇌ >MeCO₃⁻ + H₂O | 16.6 ± 0.2 | 15.4 ± 0.2 |

### 2.2 The two fitted reactions

One reaction is added per metal. Its constant is the only quantity fitted in the surface model.

**E7**  >CO₃H⁰ + Li⁺ ⇌ >CO₃Li⁰ + H⁺

**E8**  >CO₃H⁰ + Co²⁺ ⇌ >CO₃Co⁺ + H⁺

### 2.3 The constant

**E9**   K_int = ([>CO₃Me⁺][H⁺]) / ([>CO₃H⁰][Me²⁺])

This is the quantity reported in §9. In the code the inverted form is used, since the sorbed quantity is known at each sampling:

**E10**   log K_int = log[>CO₃Me] + log a(H⁺) − log[>CO₃H⁰] − log a(Me^z⁺) − (z − 1)·log b

where b is the Boltzmann factor of E13. The last term is the electrostatic correction; setting it to zero recovers E9 exactly as written in the source.

### 2.4 Mass action for every surface species

Each species is computed relative to the neutral reference form of its site. For a species formed by reaction with stoichiometry aⱼ on component j and carrying charge z:

**E11**   [species] / [reference form] = K_int · Π_j (a_j)^{ν_j} · b^z

Implemented in logarithmic space so that a strongly polarised surface cannot overflow the exponentials:

**E12**   log([species]/[reference]) = log K_int + Σ_j ν_j log a_j + z·log b

---

## 3. Electrostatics: the constant capacitance model

**E13**  Boltzmann factor, from Pokrovsky et al. (1999) Eq. 2:

  K_s = K_int · exp(−zFψ₀/RT),  so  b ≡ exp(−Fψ₀/RT)

**E14**  Constant capacitance closure, their Eq. 3:

  ψ₀ = σ / C

**E15**  Capacitance as a function of ionic strength, their Eq. 4:

  C = √I / α,  with α = 0.004 C mol^½ V⁻¹ m³

  Giving C = 207 F m⁻² at I = 0.684 M and 250 F m⁻² at I = 1.00 M.

**E16**  Net surface charge from surface species only, their Eq. 5:

  σ = (F / S_A) · Σ_k z_k [k]

  where S_A = 45.6 m² L⁻¹ is the surface per litre of suspension and the sum runs over
  >CO₃⁻, >CO₃Ca⁺, >CO₃Mg⁺, >CO₃Me^(z−1)⁺, >CaOH₂⁺, >CaO⁻, >CaCO₃⁻, >MgOH₂⁺, >MgO⁻, >MgCO₃⁻.
  Dissolved Ca²⁺ and Mg²⁺ do not appear: they belong to dolomite dissolution, not to the surface.

E14 and E16 are solved together for ψ₀ at every evaluation. Because every positively charged surface species carries b and every negatively charged one carries b⁻¹, σ falls monotonically with ψ₀ while Cψ₀ rises, so the root is unique and is bracketed by ψ₀ ∈ [−1.5, +1.5] V at any capacitance in range.

---

## 4. Site mass balances

**E17**   T(>CO₃H) = [>CO₃H⁰] + [>CO₃⁻] + [>CO₃Ca⁺] + [>CO₃Mg⁺] + [>CO₃Me^(z−1)⁺]

**E18**   T(>CaOH) = [>CaOH⁰] + [>CaOH₂⁺] + [>CaO⁻] + [>CaHCO₃⁰] + [>CaCO₃⁻]

**E19**   T(>MgOH) = [>MgOH⁰] + [>MgOH₂⁺] + [>MgO⁻] + [>MgHCO₃⁰] + [>MgCO₃⁻]

**E20**  Total sites from the site density, with N_s = 2 × 10⁻⁵ mol m⁻² on each site (Table 2.2, from Brady et al. 1999):

  T_site = N_s · SSA · a_solid = 2 × 10⁻⁵ × 0.76 × 60 = 0.912 mmol L⁻¹ on each of the three sites

A second value, N_s = 14 µmol m⁻² on the carbonate site and 7 µmol m⁻² on the metal sites (Pokrovsky et al. 1999, §3.3, their 1:1:2 stoichiometry, eight sites per nm²), gives 0.638 mmol L⁻¹ and is carried as the sensitivity case of §9.3.

---

## 5. Aqueous speciation

**E21**  Davies activity coefficient, applied to every dissolved species:

  log γ_i = −A z_i² [ √I / (1 + √I) − 0.3 I ],  A = 0.509

  Giving γ(1±) = 0.748 and γ(2±) = 0.313 at I = 0.684 M, and 0.791 and 0.392 at I = 1.00 M.

**E22**  Water dissociation (Table 2.3): OH⁻ + H⁺ ⇌ H₂O, log K = 13.991

**E23**  Carbonate from the gas phase at fixed pCO₂ (Table 2.3): CO₂(g) + H₂O ⇌ H⁺ + HCO₃⁻, log K = −7.809

  a(HCO₃⁻) = 10^(−7.809) · pCO₂ / a(H⁺),  with pCO₂ = 3.15 × 10⁻⁴ bar (Table 2.2)

**E24**  Carbonate speciation (Table 2.3): CO₃²⁻ + H⁺ ⇌ HCO₃⁻, log K = 10.325

  a(CO₃²⁻) = a(HCO₃⁻) / (10^10.325 · a(H⁺))

**E25**  Aqueous complexation, general form. Every complex is written as formation from the free ions, log β, and its abundance relative to the free ion is

  [complex] / [free ion] = 10^(log β) · Π (a_Cl)^n (a_H)^m (a_CO₃)^p

The complexes carried are listed in Table R1. Calcium, magnesium and carbonate constants are exactly those of Table 2.3, converted from the dissociation form in which that table is written. Lithium and cobalt do not appear in that table and are added from the same database family, marked as such.

**E26**  Free ion fraction, which is what the surface reaction sees:

  f_free = 1 / (1 + Σ_i [complex_i]/[free ion])

  a(Me^z⁺) = C_eq · f_free · γ_z

At pH 7.86 in produced water this gives 43 per cent of dissolved cobalt as free Co²⁺ and 85 per cent of dissolved lithium as free Li⁺; the remainder is chloride and carbonate complexes.

---

## 6. pH from charge balance

Table 2.2 sets chloride by charge balance at the imposed initial pH, then holds it. The same closure is used here.

**E27**  Charge balance on the solution:

  Σ_i z_i c_i = 0

  expanded as

  [H⁺] − [OH⁻] − [HCO₃⁻] − 2[CO₃²⁻] + [Na⁺] − [Cl⁻] + Σ_metals ( z·[free] + Σ_complexes z_c·[complex] ) = 0

**E28**  Initial closure: at t = 0 the measured pH is imposed and E27 is solved for [Cl⁻].

**E29**  Thereafter [Cl⁻] is held and E27 is solved for pH at every step of the integration.

This is the mechanism to which the method attributes the pH rise, so it is solved rather than assumed. In the single ion batches the balancing chloride comes out at 686 to 709 mmol L⁻¹, and within each metal pair the acid batch needs 6 to 8 mmol L⁻¹ more than the one started near pH 6, which is the hydrochloric acid used to set the starting point. The direction and the size are both right.

---

## 7. Dolomite dissolution

**E30**  The dissolution reaction (source Eq. 5):

  CaMg(CO₃)₂ + 2H⁺ ⇌ Ca²⁺ + Mg²⁺ + 2HCO₃⁻,  log K_eq = +2.525 (Table 2.2)

**E31**  Transition state rate law (source Eq. 6, after Lasaga 1984), as implemented:

  −R_dolomite = A · k_m · { 1 − ( [HCO₃⁻]²[Ca²⁺][Mg²⁺] ) / ( [H⁺]² K_eq ) }

  with A = 45.6 m² L⁻¹ and log k_m = −9 mol m⁻² s⁻¹ (Table 2.2), neither fitted.

**E32**  Integration, with congruent stoichiometry from E30:

  d[Ca]_T/dt = d[Mg]_T/dt = −R_dolomite

  integrated over six days by an implicit stiff solver, with pH from E29 at every step.

**E33**  Saturation state, the quotient inside the braces of E31:

  Ω = ( a(HCO₃⁻)² a(Ca²⁺) a(Mg²⁺) ) / ( a(H⁺)² K_eq )

  The rate goes to zero at Ω = 1, which is why the pH plateaus. The plateau is therefore thermodynamic and its value does not depend on k_m.

---

## 8. The two departures from the printed equations, and the evidence

### 8.1 The [H⁺] prefactor of Eq. 6

The source prints the rate as −R = A·k_m[H⁺]{1 − …}. With the tabulated k_m = 10⁻⁹ mol m⁻² s⁻¹ that gives an areal rate of order 10⁻¹⁵ mol m⁻² s⁻¹ near neutral pH. Three independent facts identify the multiplier as the error rather than the constant:

1. It dissolves nothing. Run as printed over six days, the predicted pH does not move from its starting value in any batch (Table R2).
2. It contradicts the source's own result, in which the pH rises and levels within about an hour.
3. Without it, k_m is itself the areal rate, and 10⁻⁹ mol m⁻² s⁻¹ is exactly the 10⁻¹⁰ to 10⁻⁹ mol m⁻² s⁻¹ measured for dolomite near neutral pH.

Removing the multiplier and changing nothing else, E31 predicts that the batches started near pH 6 reach 90 per cent of their pH rise in 1.8 hours and plateau at pH 8.29, against 8.17 measured at day 6, with total dissolved Ca + Mg of 0.79 mmol L⁻¹ against 0.91 measured. No parameter is adjusted anywhere in the dissolution model.

### 8.2 The exponent inside the affinity quotient

The source prints [H⁺] to the first power inside the quotient. Reaction E30 carries 2H⁺, so the quotient with the first power is not the saturation state of the reaction K_eq belongs to. The squared form of E33 is used.

### 8.3 A verification, not a departure: log K_eq = +2.525

The tabulated value is confirmed rather than replaced. Combined with the same table's K(HCO₃⁻) = 10.325 it corresponds to

  log K_eq = log K_sp(dolomite) + 2 log K(HCO₃⁻) = −18.13 + 2(10.325) = +2.52

which is the EQ3/6 value for ordered dolomite. An earlier draft of this work read the sign as negative and reported the constant as an error; that was our misreading and is corrected here.

### 8.4 A verification, not a departure: reaction E5 is positive

Reaction E5 is +24.0 and +23.5, read directly from Table 2.1. Three checks agree: E5 minus E6 is the protonation >MeCO₃⁻ + H⁺ ⇌ >MeHCO₃⁰, which must be positive and of order 10, and the positive reading gives +6.4 to +8.1 across all four mineral columns against −16.8 to −20.6 for the negative reading; Pokrovsky et al. name >MeHCO₃⁰ as one of five species their model requires, which a species at 10⁻²⁵ relative abundance cannot be; and under the negative reading >MeHCO₃⁰ vanishes from the surface at every pH of these experiments.

---

## 9. Fitting, screening and goodness of fit

### 9.1 Saturation screen, applied before any constant is formed

**E34**  Solubility floor of the metal's own carbonate:

  [Me]_floor = K_sp(MeCO₃) / ( a(CO₃²⁻) γ₂ )   and   [Li]_floor = √( K_sp(Li₂CO₃) / a(CO₃²⁻) ) / γ₁

  with log K_sp = −9.98 for sphaerocobaltite and −2.50 for zabuyelite.

**E35**  Assignment of the removal at each sampling:

  if [Me]_floor > C₀:  sorbed = C₀ − C_eq   (the carbonate cannot form; all removal is surface)
  if C_eq < [Me]_floor:  sorbed = [Me]_floor − C_eq   (only the excess below the floor is surface)
  otherwise:  sorbed = 0   (precipitation accounts for the removal)

**E36**  Coverage screen. A constant is formed only where

  0 < N = [>CO₃Me] / T(>CO₃H) < 1

### 9.2 The two fitting routes

**Route 1, point by point.** E10 is inverted at each usable sampling and the results averaged, with the standard deviation across pH reported as the measure of whether the constant is intrinsic.

**Route 2, whole profile.** The model is run forward through the simulated dissolution and log K_int is chosen to minimise

**E37**   RMSE = √( (1/n) Σ_j ( S_model(t_j) − S_measured(t_j) )² )

  where S is the percentage removed. The same statistic is reported for the pH.

**E38**  Sorption predicted at each moment, solved self consistently with E9 to E20:

  S(t) = 100 · ( T_Me − C_dissolved(t) ) / T_Me

### 9.3 Sensitivity

**E39**  Site density, E20 evaluated at 14 and 20 µmol m⁻²: moves lithium by 0.36 and cobalt by 0.20 log units.

**E40**  Capacitance, E15 with α halved and doubled: moves lithium by 0.03 and cobalt by 0.09 log units.

---

## 10. Results obtained from these equations

**Table R3** gives the constants. The recommended values are from the single ion experiment, at the site density of E20 and the capacitance of E15.

Both routes agree. Route 1 gives Li −2.69 ± 0.25 and Co −2.26 ± 0.61; Route 2 gives Li −2.83 ± 0.12 and Co −2.34 ± 0.41. The constants do not depend on which is used.

For context, E2 gives −1.8 for calcium and −2.0 for magnesium on the same site of the same mineral. Cobalt falls 0.26 log units from magnesium and lithium 0.69. The Shannon (1976) six coordinate radii are Mg²⁺ 0.72 Å, Co²⁺ 0.745 Å, Li⁺ 0.76 Å and Ca²⁺ 1.00 Å, so an affinity close to magnesium's is what the crystal chemistry predicts. Those two constants enter the model only through E2 as the fixed competition term, which the fit cannot adjust.

---

## Tables

**Table R1.** Aqueous complexes carried in E25, as formation constants from the free ions.

| Metal | Complexes | Source |
|---|---|---|
| Ca²⁺ | CaOH⁺, CaCO₃⁰, CaHCO₃⁺, CaCl⁺ | Table 2.3 |
| Mg²⁺ | MgOH⁺, MgCO₃⁰, MgHCO₃⁺, MgCl⁺ | Table 2.3 |
| Li⁺ | LiCl⁰, LiCO₃⁻, LiOH⁰ | same database family, added here |
| Co²⁺ | CoCl⁺, CoCl₂⁰, CoCO₃⁰, CoHCO₃⁺, CoOH⁺, Co(OH)₂⁰ | same database family, added here |

**Table R2.** Eq. 6 (E31) run both ways over six days, at the tabulated k_m, nothing fitted.

| Batch | pH at start | pH at day 6, with the [H⁺] prefactor | pH at day 6, without | pH measured at day 6 |
|---|---|---|---|---|
| Li, pH 6 start | 6.42 | 6.42 | 8.30 | 8.26 |
| Li, pH 2 start | 2.25 | 2.28 | 7.91 | 7.05 |
| Co, pH 6 start | 6.32 | 6.32 | 8.28 | 8.09 |
| Co, pH 2 start | 2.35 | 2.38 | 7.95 | 7.11 |

**Table R3.** Constants from E9 and E10, at T_site = 0.912 mmol L⁻¹ (E20) and C = √I/0.004 (E15), Davies coefficients (E21).

| Reaction | System | log K_int | SD | Usable samplings | Status |
|---|---|---|---|---|---|
| E7, >CO₃H⁰ + Li⁺ = >CO₃Li⁰ + H⁺ | single ion, 40 g L⁻¹ NaCl | **−2.69** | 0.25 | 6 of 6 | determined |
| E8, >CO₃H⁰ + Co²⁺ = >CO₃Co⁺ + H⁺ | single ion, 40 g L⁻¹ NaCl | **−2.26** | 0.61 | 5 of 6 | determined |
| E7 | produced water | −1.11 | 0.78 | 5 of 6 | upper bound, not transferable |
| E8 | produced water | +0.77 | 1.38 | 5 of 6 | upper bound, not transferable |

**Table R4.** Whole profile fits by E37, single ion experiment, with the dissolution predicted by E31 at the tabulated k_m.

| Batch | log k_m | pH RMSE | log K_int | Sorption RMSE |
|---|---|---|---|---|
| Li, pH 6 start | −9, tabulated | 0.43 | −2.91 | 0.63 % |
| Li, pH 2 start | −9, tabulated | 2.51 | −2.74 | 0.20 % |
| Co, pH 6 start | −9, tabulated | 0.70 | −2.05 | 1.79 % |
| Co, pH 2 start | −9, tabulated | 2.62 | −2.63 | 0.40 % |

---

## 11. Equation index

| Tag | Equation | Source |
|---|---|---|
| E1 to E6 | Surface reactions on the three sites | Pokrovsky et al. (1999) Table 3 |
| E7, E8 | The two fitted metal reactions | this work, in the form of source Eq. 3 |
| E9 | The intrinsic constant | source Eq. 4 |
| E10 | Inverted form of E9, with the electrostatic term | this work |
| E11, E12 | Mass action for a surface species | Pokrovsky et al. (1999) |
| E13 | Boltzmann factor | Pokrovsky et al. (1999) Eq. 2 |
| E14 | Constant capacitance closure | Pokrovsky et al. (1999) Eq. 3 |
| E15 | Capacitance against ionic strength | Pokrovsky et al. (1999) Eq. 4 |
| E16 | Net surface charge | Pokrovsky et al. (1999) Eq. 5 |
| E17 to E19 | Site mass balances | Pokrovsky et al. (1999) |
| E20 | Total sites from site density | Table 2.2, Brady et al. (1999) |
| E21 | Davies activity coefficient | Davies (1962) |
| E22 to E24 | Water and carbonate equilibria | Table 2.3 |
| E25, E26 | Aqueous complexation and free ion fraction | Table 2.3 |
| E27 to E29 | Charge balance and its closure | Table 2.2 |
| E30 | Dolomite dissolution reaction | source Eq. 5 |
| E31 | Transition state rate law | source Eq. 6, after Lasaga (1984) |
| E32 | Integration of E31 | this work |
| E33 | Saturation state | the affinity term of E31 |
| E34 to E36 | Solubility floor, assignment, coverage screen | this work |
| E37 | Goodness of fit | this work |
| E38 | Predicted sorption | this work |
| E39, E40 | Sensitivity to E20 and E15 | after Hayes et al. (1991) |

---

## References

Brady, P. V., Papenguth, H. W., & Kelly, J. W. (1999). Metal sorption to dolomite surfaces. *Applied Geochemistry, 14*(5), 569–579.

Davies, C. W. (1962). *Ion association*. Butterworths.

Elshebli, M., Vilcáez, J., & Smay, J. (2025). Lithium and cobalt recovery from petroleum produced water using dolomite: Impact of pH, temperature, and surface area. *Science of the Total Environment, 1003*, 180743.

Hayes, K. F., Redden, G., Ela, W., & Leckie, J. O. (1991). Surface complexation models: An evaluation of model parameter estimation using FITEQL and oxide mineral titration data. *Journal of Colloid and Interface Science, 142*(2), 448–469.

Lasaga, A. C. (1984). Chemical kinetics of water–rock interactions. *Journal of Geophysical Research, 89*(B6), 4009–4025.

Pokrovsky, O. S., Schott, J., & Thomas, F. (1999). Dolomite surface speciation and reactivity in aquatic systems. *Geochimica et Cosmochimica Acta, 63*(19–20), 3133–3143.

Shannon, R. D. (1976). Revised effective ionic radii and systematic studies of interatomic distances in halides and chalcogenides. *Acta Crystallographica A, 32*(5), 751–767.

Smith, R. M., & Martell, A. E. (2004). *NIST critically selected stability constants of metal complexes database* (Version 8.0). National Institute of Standards and Technology.

Steefel, C. I., & Molins, S. (2016). *CrunchFlow: Software for modeling multicomponent reactive flow and transport*. Lawrence Berkeley National Laboratory.
