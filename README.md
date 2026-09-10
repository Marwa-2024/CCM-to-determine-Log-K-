# CCM to determine log K — Li and Co on dolomite

Two self-contained Python scripts that determine the intrinsic surface complexation
constants (log K) for lithium and cobalt binding on dolomite, using the constant
capacitance three-site model of Pokrovsky, Schott & Thomas (1999).

## Start here: the stepwise notebook and the paper

* **`LogK_LiCo_dolomite_carbonate_site_STEPWISE.ipynb`** — the complete procedure in one
  executed notebook (outputs embedded): capacity check → aqueous speciation and saturation
  indices at the *measured* equilibrium pH → surface coverage with propagated error → carbonate
  site balance with Ca/Mg competition → constant capacitance electrostatic loop (and the
  non-electrostatic case) → point-by-point log K → global WSOS/DF fit → the flat-line
  diagnostic → kinetics, adsorption edge, isotherm, K_d and surface speciation figures →
  prediction from literature constants and the sorption/mineralisation split. Source:
  `logK_carbonate_site_stepwise.py` (cell-marked; the notebook is generated from it).
* **`Paper_LiCo_dolomite_carbonate_site_SCM.md` / `.docx`** — publication-style write-up with
  every reaction, equation, parameter table, data table, result and APA reference.

Headline: the high-recovery runs in Elshebli et al. (2025) exceed the carbonate-site density
by 1.7 to 2.4 times, so they are mineralisation, not adsorption. In the low-uptake single-ion
data set, lithium gives a pH-independent log K(>CO3Li0) = −2.42 ± 0.30; cobalt cannot be
fitted (its measurable point is supersaturated with sphaerocobaltite) and is predicted from
the aqueous–surface analogy, log K(>CO3Co+) = −0.79.

## Notebooks (all data embedded, no external files)

* **`FITEQL_LiCo_dolomite_Colab.ipynb`** — Google Colab notebook. Open it at
  colab.research.google.com (File, Upload notebook) and choose Runtime, Run all. Every
  reaction, equation, and datum is written out in the first cells; the rest is the
  FITEQL-style solver in Python.
* **`FITEQL_LiCo_dolomite_Mathematica.nb`** — Wolfram Mathematica notebook with the same
  content and data, implemented in the Wolfram Language. Open it in Mathematica and choose
  Evaluation, Evaluate Notebook.

Both reproduce log K(Co) = −1.89, log K(Li) = +0.16, WSOS/DF = 0.11.

## The binding site (important)

Per Pokrovsky (1999), metal **cations** adsorb on the **carbonate** site,
`>CO3H0 + Me(z+) = >CO3Me(z-1) + H+` (the deprotonated `>CO3-` binds the cation). The
`>CaOH`/`>MgOH` hydroxyl sites take anions and ligands, not metal cations. The fitted
cobalt constant, **−1.89**, sits right next to Pokrovsky's own `>CO3H0 + Ca2+ = >CO3Ca+`
value of **−1.8**, which is the expected result for a divalent cation on the carbonate
site and a good check on the model.

## Which script to run

* **`fiteql_logK_Li_Co_dolomite.py`** — the correct, recommended method. Pure-Python
  reproduction of FITEQL / Visual MINTEQ: a component/species tableau, a Newton-Raphson
  solve of the full coupled equilibrium at each pH (surface potential as its own Boltzmann
  component), and a fit of the metal constants by minimising WSOS/DF (Westall 1982;
  Herbelin & Westall 1999; Visual MINTEQ Eqs 4.4, 4.6, 4.11). Metal cations bind the
  carbonate site. Run this one.
* **`determine_logK_Li_Co_dolomite.py`** — a first attempt that (a) fit by a direct match
  to uptake and (b) placed the metal on the `>CaOH`/`>MgOH` sites. **Superseded** by the
  script above; kept only as a record of the earlier approach.

## The FITEQL method, in three steps

1. **Tableau.** Every species, aqueous and surface, is written as a product of a few
   components: `C_i = K_i · Π_j X_j^(a_ij)`. The surface potential enters as its own
   component `P = exp(−Fψ/RT)`, and each surface species carries `P` to its charge.
2. **Inner loop (equilibrium).** For fixed log K, the mole-balance residuals
   `Y_j = Σ_i a_ij C_i − T_j = 0` are solved together by Newton-Raphson with the
   analytical Jacobian `Z_jk = Σ_i a_ij a_ik C_i`, plus the capacitance term on the
   electrostatic component.
3. **Outer loop (fit).** The unknown constants are adjusted to minimise
   `WSOS/DF = Σ(Y_i/s_i)² / (N_obs − N_param)`; a value near 1 means the model matches
   the data to within its measurement error.

---

Below refers to the earlier `determine_logK_Li_Co_dolomite.py` pipeline.

The workflow supports the study *"Lithium and cobalt recovery from petroleum produced
water using dolomite"* (Elshebli, Vilcáez & Smay). This stage covers the single-ion Li
and Co experiments; the code is structured so the same framework extends to the full
produced-water cation suite later.

## Run it

```bash
pip install numpy pandas scipy matplotlib
python determine_logK_Li_Co_dolomite.py
```

Everything is embedded in the one file — no notebook or Excel file is required. If the
original Excel workbook sits alongside the script it is used only for a cross-check.

## What the script does

It walks through the workflow in order, printing a short plain explanation before each
block. Step 0 writes out every equation the model solves (also saved to
`MODEL_EQUATIONS.txt`) so the model can be checked by hand.

0. The complete equation set: aqueous reactions, surface reactions, mass and site
   balances, the constant capacitance relations, the Boltzmann correction, and the two
   separate surface-charge definitions (predicted from surface species vs measured from
   the solution mass balance).
1. Organise and validate the measured datasets (σ₀, ζ, Li and Co adsorption).
2. Define the three-site surface complexation model and the four metal reactions.
3. Compute aqueous speciation and free-ion activities at each pH.
4. Fix the background (proton) constants from Pokrovsky 1999 and calibrate the EDL
   capacitance against the measured surface charge.
5. Build the least-squares objective.
6. Provide the forward prediction function (self-consistent CCM solver).
7. Fit the metal constants with `scipy.optimize.least_squares`.
8. Validate: R², residuals, and figures.
9. Sensitivity analysis (±0.5 log K on each constant).
10. Compare to literature and interpret.
11. Confidence intervals from the fit covariance.

## Outputs

* `MODEL_EQUATIONS.txt` — every equation the code solves, written out precisely.
* `results_logK_summary.txt` — the final constants, fit quality, and interpretation.
* `surface_charge_fit.png`, `adsorption_fit.png`, `surface_speciation.png`,
  `sensitivity.png` — validation figures.

## Predicted vs measured surface charge

The predicted surface charge is built only from the surface species the model solves; it
reads nothing from the spreadsheet. The spreadsheet σ is the separate measured quantity
(the Charlet mass balance, your Eq 9), which also carries the effect of dolomite
dissolution — at acidic pH it reaches roughly sixty times the total site capacity, so it
cannot be a surface charge there. The two are compared on the point of zero charge and on
magnitude near it, never fed into one another.

## Honest scope

The single-ion data constrains binding at two pH points per metal, so one constant per
metal is fitted (the Ca-site value) and the Mg-site constant is tied to it with a fixed
offset taken from Pokrovsky's Ca vs Mg protonation difference. The model represents
surface sorption only; cobalt's carbonate mineralisation pathway is discussed but not
fitted. The capacitance uses Pokrovsky's published α = 0.004 (high-capacitance regime,
surface potential near zero); the fitted constants shift if α is changed. The script
reports what the data can and cannot constrain rather than implying a precision the
measurements do not support.
