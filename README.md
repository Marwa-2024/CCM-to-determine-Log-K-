# CCM to determine log K — Li and Co on dolomite

A single, self-contained Python pipeline that determines the intrinsic surface
complexation constants (log K) for lithium and cobalt binding on dolomite, using the
constant capacitance three-site model of Pokrovsky, Schott & Thomas (1999).

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
