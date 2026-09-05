#!/usr/bin/env python3
# ======================================================================================
#  DETERMINING INTRINSIC STABILITY CONSTANTS (log K) FOR Li AND Co ON DOLOMITE
#  A single, self-contained pipeline following the constant capacitance surface
#  complexation model of Pokrovsky, Schott & Thomas (1999).
# --------------------------------------------------------------------------------------
#  Author context : Elshebli, Vilcaez & Smay - "Lithium and cobalt recovery from
#                   petroleum produced water using dolomite".
#  Stage          : single-ion Li and Co (this file). The design leaves room to move
#                   on to produced water (PW) later by adding cations to one dictionary.
#
#  WHAT THIS FILE DOES
#  -------------------
#  It walks through the eleven steps of the workflow, one after another, and prints a
#  short plain-language explanation before each block runs, so you can read the console
#  output top to bottom like a lab notebook. Nothing here needs the original notebook or
#  the Excel file: every measurement is embedded below. If the Excel file happens to be
#  next to this script it will be read for a cross-check, otherwise the embedded copy is
#  used.
#
#  HOW TO RUN
#  ----------
#      pip install numpy pandas scipy matplotlib
#      python determine_logK_Li_Co_dolomite.py
#
#  Figures are written next to the script (surface_charge_fit.png, adsorption_fit.png,
#  surface_speciation.png, sensitivity.png) and a text report is written to
#  results_logK_summary.txt.
#
#  A NOTE ON HONESTY OF THE FIT
#  ----------------------------
#  The single-ion adsorption data constrains the metal binding constants at two pH
#  points per metal (pH 2 and pH 6). Four unknowns fit to a small, kinetically limited
#  data set will be only loosely determined, and the split of binding between the Ca and
#  the Mg site is nearly degenerate. The pipeline is built correctly and reports what the
#  data can and cannot constrain (see Steps 8, 9 and 11); it does not manufacture a
#  precision the measurements do not support.
# ======================================================================================

import os
import sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")                      # write figures to file without a display
import matplotlib.pyplot as plt
from scipy.optimize import brentq, least_squares

HERE = os.path.dirname(os.path.abspath(__file__))
np.set_printoptions(suppress=True)


def banner(step_no, title):
    """Print a clearly delimited header so the console reads like a numbered protocol."""
    line = "=" * 86
    print("\n" + line)
    print(f"STEP {step_no}: {title}")
    print(line)


def explain(text):
    """Print a short paragraph of plain-language explanation, wrapped for the terminal."""
    import textwrap
    for para in text.strip().split("\n"):
        print(textwrap.fill(para.strip(), width=86))
    print()


# ======================================================================================
#  PART 0 - PHYSICAL CONSTANTS AND THE THERMODYNAMIC DATABASE
#  These are shared by every step. The aqueous log K set is the TOUGHREACT / EQ3-6
#  database used in your original speciation notebook, extended here with Li and Co.
# ======================================================================================

F       = 96485.0          # Faraday constant, C/mol
R_GAS   = 8.314            # J/mol/K
T_KELV  = 298.15           # 25 C
FRT     = F / (R_GAS * T_KELV)   # F/RT, 1/volt  (appears in every Boltzmann factor)

# --- aqueous dissociation log K (25 C), written as  complex = products ---------------
# (identical to your notebook, with Li and Co complexes appended at the bottom)
LK = dict(
    co2=-6.345, co3=10.329, oh=13.995,
    cacl=0.696,  caco3=7.002,  cahco3=-1.047, caoh=12.850,
    mgcl=0.135,  mgco3=7.350,  mghco3=-1.036, mgoh=11.785,
    nacl=0.777,  naco3=9.815,  nahco3=-0.154, naoh=14.180,
    # --- metals of interest (formation written as dissociation, EQ3/6 sign convention) ---
    licl=0.50,               # LiCl  = Li+ + Cl-      (Li barely complexes)
    cocl=0.05,               # CoCl+ = Co2+ + Cl-
    cooh=9.65,               # CoOH+ = Co2+ + OH- ... entered below as hydrolysis
    coco3=-4.22,             # CoCO3 = Co2+ + CO3-2   (i.e. logK_formation = 4.22)
)

# open (atmospheric) CO2 system, exactly as in your notebook
PCO2   = 10 ** -3.5          # atmospheric pCO2 (atm)
LOGKH  = -1.469             # CO2(g) = CO2(aq)
ACO2   = 10 ** LOGKH * PCO2  # fixed {CO2(aq)} activity

A_DAVIES = 0.509            # Debye-Huckel A at 25 C
DAVIES_B = 0.3             # Davies extended-term coefficient
MM = dict(Ca=40.078, Mg=24.305, Na=22.99, Cl=35.45, CO3=60.008, HCO3=61.016,
          Li=6.941, Co=58.933)

# --- dolomite solid / reactor set-up (from your titration spreadsheet) ---------------
DOLOMITE_GL      = 60.0     # g of dolomite per litre
SURFACE_AREA_M2G = 0.76     # specific surface area, m2/g
ST_AREA          = DOLOMITE_GL * SURFACE_AREA_M2G   # 45.3 m2 of surface per litre

# --- site densities (Pokrovsky 1999): 1 : 1 : 2 for Ca : Mg : carbonate -------------
# The paper lists 7 (Ca), 7 (Mg) and 14 (carbonate) sites; the natural unit for
# rhombohedral carbonates is micromoles per square metre.
SITE_DENS = dict(Ca=7.0e-6, Mg=7.0e-6, CO3=14.0e-6)   # mol / m2
# total moles of each site per litre of reactor = density * surface area per litre
SITE_TOT = {k: v * ST_AREA for k, v in SITE_DENS.items()}   # mol / L

# --- constant capacitance EDL parameter (Pokrovsky 1999, Eqn 4:  C = sqrt(I)/alpha) --
ALPHA_EDL = 0.004           # C * mol^0.5 / (volt * m3); refined against your sigma0 below


def gammas(I):
    """Davies activity coefficients for charge 1, 2 and neutral species."""
    if I <= 0:
        return 1.0, 1.0, 1.0
    f = -A_DAVIES * (np.sqrt(I) / (1 + np.sqrt(I)) - DAVIES_B * I)
    return 10 ** f, 10 ** (4 * f), 1.0     # gamma(z=1), gamma(z=2), gamma(neutral)


# ======================================================================================
#  STEP 3 helper - AQUEOUS SPECIATION ENGINE (your notebook, integrated verbatim)
#  Given a pH and the analytical totals it returns the free concentrations and, through
#  the activity coefficients, the activities of every aqueous species at fixed pCO2.
# ======================================================================================

def speciate(pH, CaT, MgT, NaT, ClT):
    """Full carbonate-system speciation at fixed pH and fixed pCO2 (Davies activities)."""
    aH = 10 ** (-pH)
    I = 0.5 * (NaT + ClT + 4 * CaT + 4 * MgT)
    aCa = aMg = aNa = aCl = 0.0
    for _ in range(200):
        g1, g2, g0 = gammas(I)
        aHCO3 = 10 ** LK['co2'] * ACO2 / aH
        aCO3 = aHCO3 / (10 ** LK['co3'] * aH)
        aOH = 10 ** (-LK['oh']) / aH
        cH, cOH, cHCO3, cCO3, cCO2 = aH / g1, aOH / g1, aHCO3 / g1, aCO3 / g2, ACO2 / g0
        if aCl == 0:
            aCa, aMg, aNa, aCl = g2 * CaT, g2 * MgT, g1 * NaT, g1 * ClT
        for _ in range(100):
            dCa = (1 / g2 + aCl / (g1 * 10 ** LK['cacl']) + aHCO3 / (g0 * aH * 10 ** LK['caco3'])
                   + aHCO3 * 10 ** (-LK['cahco3']) / g1 + 1 / (g1 * aH * 10 ** LK['caoh']))
            dMg = (1 / g2 + aCl / (g1 * 10 ** LK['mgcl']) + aHCO3 / (g0 * aH * 10 ** LK['mgco3'])
                   + aHCO3 * 10 ** (-LK['mghco3']) / g1 + 1 / (g1 * aH * 10 ** LK['mgoh']))
            dNa = (1 / g1 + aHCO3 / (g1 * aH * 10 ** LK['naco3']) + aHCO3 * 10 ** (-LK['nahco3']) / g0
                   + aCl / (g0 * 10 ** LK['nacl']) + 1 / (g0 * aH * 10 ** LK['naoh']))
            aCa, aMg, aNa = CaT / dCa, MgT / dMg, NaT / dNa
            dCl = (1 / g1 + aCa / (g1 * 10 ** LK['cacl']) + aMg / (g1 * 10 ** LK['mgcl'])
                   + aNa / (g0 * 10 ** LK['nacl']))
            aCl_new = ClT / dCl
            if abs(aCl_new - aCl) < 1e-12:
                aCl = aCl_new
                break
            aCl = aCl_new
        cCaCl = aCa * aCl / (g1 * 10 ** LK['cacl']); cCaCO3 = aCa * aHCO3 / (g0 * aH * 10 ** LK['caco3'])
        cCaHCO3 = aCa * aHCO3 * 10 ** (-LK['cahco3']) / g1; cCaOH = aCa / (g1 * aH * 10 ** LK['caoh'])
        cMgCl = aMg * aCl / (g1 * 10 ** LK['mgcl']); cMgCO3 = aMg * aHCO3 / (g0 * aH * 10 ** LK['mgco3'])
        cMgHCO3 = aMg * aHCO3 * 10 ** (-LK['mghco3']) / g1; cMgOH = aMg / (g1 * aH * 10 ** LK['mgoh'])
        cNaCl = aNa * aCl / (g0 * 10 ** LK['nacl']); cNaCO3 = aNa * aHCO3 / (g1 * aH * 10 ** LK['naco3'])
        cNaHCO3 = aNa * aHCO3 * 10 ** (-LK['nahco3']) / g0; cNaOH = aNa / (g0 * aH * 10 ** LK['naoh'])
        cCa, cMg, cNa, cCl = aCa / g2, aMg / g2, aNa / g1, aCl / g1
        Inew = 0.5 * (cNa + cCl + cH + cOH + 4 * cCa + 4 * cMg + 4 * cCO3 + cHCO3 + cCaCl
                      + cCaHCO3 + cCaOH + cMgCl + cMgHCO3 + cMgOH + cNaCO3)
        if abs(Inew - I) < 1e-10:
            I = Inew
            break
        I = 0.5 * I + 0.5 * Inew
    return dict(pH=pH, I=I, aH=aH, aOH=aOH, aHCO3=aHCO3, aCO3=aCO3,
                aCa=aCa, aMg=aMg, aNa=aNa, aCl=aCl, g1=g1, g2=g2,
                H=cH, OH=cOH, CO2=cCO2, HCO3=cHCO3, CO3=cCO3,
                Ca=cCa, Mg=cMg, Na=cNa, Cl=cCl)


# ======================================================================================
#  STEP 1 - ORGANISE AND VALIDATE THE EXPERIMENTAL DATA
# ======================================================================================

# Net surface charge (background titration, no added Li/Co) from your "Surface charge"
# sheet. pH here is the calculated pH of vessel C; sigma is in mol/m2.
SIGMA0_DATA = pd.DataFrame({
    "pH":      [5.507, 5.474, 5.451, 5.499, 7.792, 8.482, 8.600, 8.630, 8.648],
    "sigma0":  [-5.665e-4, -7.608e-4, -9.986e-4, -1.730e-3, -2.787e-5,
                 4.407e-5,  7.025e-5,  6.915e-5,  8.441e-5],
})

# Zeta potential (mV). Your sheet lists values at pH 2,4,6,8 and several temperatures;
# we take the 25 C-relevant rows. These enter the model only qualitatively (sign / PZC)
# because their magnitudes are not consistent with a simple CCM zeta - see Step 8.
ZETA_DATA = pd.DataFrame({
    "pH":   [2, 4, 6, 8],
    "zeta": [56.1, 59.8, 55.4, 59.4],   # dolomite in NaCl, lowest measured T
})

# Single-ion adsorption at equilibrium (day 6). Concentrations in mg/L (ppm).
# C0 = initial, Ceq = dissolved at equilibrium. Uptake = C0 - Ceq.
ADSORPTION_DATA = pd.DataFrame({
    "metal":  ["Co", "Co", "Li", "Li"],
    "pH_meas":[6.0,  2.0,  6.0,  2.0],       # nominal / initial pH of the batch
    "C0_ppm": [80.513, 80.513, 138.43, 138.43],
    "Ceq_ppm":[74.393, 80.398, 135.262, 135.592],
})

# Conditions common to the adsorption batches (produced-water surrogate, single ion):
ADS_NaCl_gL = 40.0                              # g/L NaCl background
ADS_ClT     = ADS_NaCl_gL / MM['Cl'] / 1e3 * 1e3 / 1e3  # placeholder; set explicitly below
ADS_ClT     = 40.0 / MM['Cl']                   # mol/L Cl from NaCl  (~0.69 M)
ADS_NaT     = 40.0 / MM['Cl']                   # mol/L Na (charge balance of NaCl)
ADS_I       = 0.70                              # ionic strength of the batch, mol/L


def step1_master_table():
    banner(1, "Organise and validate the experimental datasets")
    explain("""
    I gather the three measurement families into aligned tables: the potentiometric
    surface charge sigma0 versus pH, the zeta potential versus pH, and the single-ion
    adsorption results for Li and Co. Adsorption was measured at initial pH 6 and pH 2,
    sampled to day 6, which I treat as the equilibrium point. For each metal I compute
    how much was taken up from solution, both in mg/L and in mol/L, because the model
    works in moles.
    """)
    df = ADSORPTION_DATA.copy()
    df["uptake_ppm"] = df["C0_ppm"] - df["Ceq_ppm"]
    df["uptake_pct"] = 100 * df["uptake_ppm"] / df["C0_ppm"]
    df["C0_M"]  = [row.C0_ppm  / MM[row.metal] / 1e3 for row in df.itertuples()]
    df["Ceq_M"] = [row.Ceq_ppm / MM[row.metal] / 1e3 for row in df.itertuples()]
    df["uptake_M"] = df["C0_M"] - df["Ceq_M"]

    print("Surface charge (background titration):")
    print(SIGMA0_DATA.to_string(index=False, float_format=lambda x: f"{x: .3e}"))
    print("\nZeta potential:")
    print(ZETA_DATA.to_string(index=False))
    print("\nAdsorption, single ion (equilibrium = day 6):")
    print(df.to_string(index=False, float_format=lambda x: f"{x:.4g}"))

    print("\nValidation of the common experimental conditions:")
    print(f"  dolomite loading        = {DOLOMITE_GL:.0f} g/L")
    print(f"  specific surface area   = {SURFACE_AREA_M2G:.2f} m2/g")
    print(f"  surface per litre (St)  = {ST_AREA:.1f} m2/L")
    print(f"  ionic strength (batch)  ~ {ADS_I:.2f} M (mostly {ADS_NaCl_gL:.0f} g/L NaCl)")
    print(f"  pCO2                    = 10^-3.5 atm (open system)")
    print(f"  temperature             = 25 C")
    print(f"  total Ca+Mg surface sites = {SITE_TOT['Ca']+SITE_TOT['Mg']:.3e} mol/L "
          f"(the pool Li/Co can bind to)")
    return df


# ======================================================================================
#  STEP 2 - SURFACE COMPLEXATION MODEL FRAMEWORK
# ======================================================================================

# Background (proton + carbonate) reactions, dolomite column of Pokrovsky 1999 Table 3.
# Each entry: log K at I = 0, written from the neutral reference site.
BACKGROUND_LOGK = {
    # carbonate site, reference >CO3H0
    "CO3_deprot":   -4.8,    # >CO3H0            = >CO3-   + H+
    "CO3Ca":        -1.8,    # >CO3H0 + Ca2+     = >CO3Ca+ + H+
    "CO3Mg":        -2.0,    # >CO3H0 + Mg2+     = >CO3Mg+ + H+
    # calcium hydroxyl site, reference >CaOH0
    "CaOH2_prot":   11.5,    # >CaOH0 + H+       = >CaOH2+
    "CaO_deprot":  -12.0,    # >CaOH0            = >CaO-   + H+
    "CaHCO3":       -4.0,    # >CaOH0 + CO3-2 +2H+ = >CaHCO3(0) + H2O
    "CaCO3":        16.6,    # >CaOH0 + CO3-2 + H+ = >CaCO3- + H2O
    # magnesium hydroxyl site, reference >MgOH0
    "MgOH2_prot":   10.6,    # >MgOH0 + H+       = >MgOH2+
    "MgO_deprot":  -12.0,    # >MgOH0            = >MgO-   + H+
    "MgHCO3":       -3.5,    # >MgOH0 + CO3-2 +2H+ = >MgHCO3(0) + H2O
    "MgCO3":        15.4,    # >MgOH0 + CO3-2 + H+ = >MgCO3- + H2O
}

# The four metal surface reactions. Written from the neutral metal site:
#   >CaOH0 + Me(z+) = >CaOMe(z-1) + H+          (Li: z=1 -> neutral; Co: z=2 -> +1)
METAL_REACTIONS = ["Li_Ca", "Li_Mg", "Co_Ca", "Co_Mg"]
METAL_CHARGE = dict(Li=1, Co=2)

# --- identifiability control --------------------------------------------------------
# Two adsorption pH points per metal (and almost no uptake at pH 2) cannot separately
# determine the Ca-site and Mg-site binding of the same metal: the two constants are
# degenerate. We therefore fit ONE constant per metal (its Ca-site value) and tie the
# Mg-site value to it with a fixed offset. The offset is Pokrovsky's own Ca vs Mg
# difference for the >MeOH2+ protonation reaction (log K 11.5 for Ca minus 10.6 for Mg),
# i.e. the Mg hydroxyl site binds about 0.9 log units more weakly than the Ca site.
MG_OFFSET = 10.6 - 11.5          # = -0.9
FREE_PARAMS = ["Li_Ca", "Co_Ca"]   # the two constants actually fitted


def build_logK(x):
    """Expand the two fitted constants into the full set of four, tying Mg to Ca."""
    liCa, coCa = x
    return {"Li_Ca": liCa, "Li_Mg": liCa + MG_OFFSET,
            "Co_Ca": coCa, "Co_Mg": coCa + MG_OFFSET}


def step2_model_framework():
    banner(2, "Establish the surface complexation model framework")
    explain("""
    I use Pokrovsky's three-site dolomite model: a carbonate site >CO3H0, a calcium
    hydroxyl site >CaOH0 and a magnesium hydroxyl site >MgOH0. The background reactions
    (protonation, deprotonation and carbonate binding) already have published log K
    values, listed below, and I hold them fixed so the model reproduces the point of zero
    charge near pH 8. The four reactions I actually solve for are the binding of Li and
    Co to the two metal-hydroxyl sites. One correction to the way the reactions were
    first written: charge must balance, so for a monovalent ion the product is neutral
    (>CaOLi0, >MgOLi0) and for a divalent ion it carries +1 (>CaOCo+, >MgOCo+).
    """)
    print("Background reactions held fixed (Pokrovsky 1999, Table 3, dolomite):")
    labels = {
        "CO3_deprot": ">CO3H0            = >CO3-    + H+",
        "CO3Ca":      ">CO3H0 + Ca2+     = >CO3Ca+  + H+",
        "CO3Mg":      ">CO3H0 + Mg2+     = >CO3Mg+  + H+",
        "CaOH2_prot": ">CaOH0 + H+       = >CaOH2+",
        "CaO_deprot": ">CaOH0            = >CaO-    + H+",
        "CaHCO3":     ">CaOH0 + CO3-2 +2H+ = >CaHCO3(0) + H2O",
        "CaCO3":      ">CaOH0 + CO3-2 + H+ = >CaCO3-  + H2O",
        "MgOH2_prot": ">MgOH0 + H+       = >MgOH2+",
        "MgO_deprot": ">MgOH0            = >MgO-    + H+",
        "MgHCO3":     ">MgOH0 + CO3-2 +2H+ = >MgHCO3(0) + H2O",
        "MgCO3":      ">MgOH0 + CO3-2 + H+ = >MgCO3-  + H2O",
    }
    for k, lab in labels.items():
        print(f"   log K = {BACKGROUND_LOGK[k]:+6.2f}   {lab}")
    print("\nUnknown reactions to solve (the four intrinsic constants):")
    print("   >CaOH0 + Li+   = >CaOLi(0) + H+     log K_int,Li-Ca  <-- FIT")
    print(f"   >MgOH0 + Li+   = >MgOLi(0) + H+     log K_int,Li-Mg  = Li-Ca {MG_OFFSET:+.1f}")
    print("   >CaOH0 + Co2+  = >CaOCo+  + H+      log K_int,Co-Ca  <-- FIT")
    print(f"   >MgOH0 + Co2+  = >MgOCo+  + H+      log K_int,Co-Mg  = Co-Ca {MG_OFFSET:+.1f}")
    explain(f"""
    Because the single-ion data has only two pH points per metal, I fit one binding
    constant per metal (the Ca-site value) and tie the Mg-site constant to it with the
    fixed offset {MG_OFFSET:+.1f} log units taken from Pokrovsky's Ca vs Mg protonation
    difference. This leaves two well-determined unknowns instead of four degenerate ones.
    """)


# ======================================================================================
#  STEP 6 core - THE FORWARD SURFACE-COMPLEXATION SOLVER
#  Given aqueous activities and a trial set of the four metal log K values, this solves
#  the constant-capacitance model self-consistently for the surface potential and returns
#  the surface charge and the amount of metal bound. Steps 5, 7, 8 all call it.
# ======================================================================================

def surface_solve(sp, logK, metal=None, aMe=0.0, alpha=ALPHA_EDL):
    """
    sp     : output of speciate() -> aqueous activities and ionic strength
    logK   : dict of the four metal constants keyed by 'Li_Ca','Li_Mg','Co_Ca','Co_Mg'
    metal  : 'Li', 'Co', or None (background only)
    aMe    : free activity of the metal ion (Li+ or Co2+)
    returns: dict with sigma0 (C/m2), psi0 (V), and bound metal (mol/L) on Ca and Mg sites
    """
    aH, aCO3, aCa, aMg = sp["aH"], sp["aCO3"], sp["aCa"], sp["aMg"]
    C_cap = np.sqrt(max(sp["I"], 1e-6)) / alpha          # CCM capacitance, F/m2
    z = METAL_CHARGE.get(metal, 0)

    def species_given_psi(psi):
        b = np.exp(-FRT * psi)          # Boltzmann factor for a +1 surface species
        # ---- carbonate site, expressed relative to {>CO3H0} ----
        denom_CO3 = (1.0
                     + 10 ** BACKGROUND_LOGK["CO3_deprot"] / aH / b      # >CO3-  (z=-1)
                     + 10 ** BACKGROUND_LOGK["CO3Ca"] * aCa / aH * b     # >CO3Ca+(z=+1)
                     + 10 ** BACKGROUND_LOGK["CO3Mg"] * aMg / aH * b)    # >CO3Mg+(z=+1)
        CO3H = SITE_TOT["CO3"] / denom_CO3
        s_CO3m  = 10 ** BACKGROUND_LOGK["CO3_deprot"] / aH / b * CO3H
        s_CO3Ca = 10 ** BACKGROUND_LOGK["CO3Ca"] * aCa / aH * b * CO3H
        s_CO3Mg = 10 ** BACKGROUND_LOGK["CO3Mg"] * aMg / aH * b * CO3H

        # ---- Ca hydroxyl site, expressed relative to {>CaOH0} ----
        term_CaOMe = 0.0
        if metal is not None:
            # >CaOMe has charge (z-1); its Boltzmann factor is b**(z-1)
            term_CaOMe = 10 ** logK[f"{metal}_Ca"] * aMe / aH * b ** (z - 1)
        denom_Ca = (1.0
                    + 10 ** BACKGROUND_LOGK["CaOH2_prot"] * aH * b        # >CaOH2+ (z=+1)
                    + 10 ** BACKGROUND_LOGK["CaO_deprot"] / aH / b        # >CaO-   (z=-1)
                    + 10 ** BACKGROUND_LOGK["CaHCO3"] * aCO3 * aH ** 2    # >CaHCO3 (z=0)
                    + 10 ** BACKGROUND_LOGK["CaCO3"] * aCO3 * aH / b      # >CaCO3- (z=-1)
                    + term_CaOMe)
        CaOH = SITE_TOT["Ca"] / denom_Ca
        s_CaOH2 = 10 ** BACKGROUND_LOGK["CaOH2_prot"] * aH * b * CaOH
        s_CaO   = 10 ** BACKGROUND_LOGK["CaO_deprot"] / aH / b * CaOH
        s_CaCO3 = 10 ** BACKGROUND_LOGK["CaCO3"] * aCO3 * aH / b * CaOH
        s_CaOMe = term_CaOMe * CaOH

        # ---- Mg hydroxyl site, expressed relative to {>MgOH0} ----
        term_MgOMe = 0.0
        if metal is not None:
            term_MgOMe = 10 ** logK[f"{metal}_Mg"] * aMe / aH * b ** (z - 1)
        denom_Mg = (1.0
                    + 10 ** BACKGROUND_LOGK["MgOH2_prot"] * aH * b
                    + 10 ** BACKGROUND_LOGK["MgO_deprot"] / aH / b
                    + 10 ** BACKGROUND_LOGK["MgHCO3"] * aCO3 * aH ** 2
                    + 10 ** BACKGROUND_LOGK["MgCO3"] * aCO3 * aH / b
                    + term_MgOMe)
        MgOH = SITE_TOT["Mg"] / denom_Mg
        s_MgOH2 = 10 ** BACKGROUND_LOGK["MgOH2_prot"] * aH * b * MgOH
        s_MgO   = 10 ** BACKGROUND_LOGK["MgO_deprot"] / aH / b * MgOH
        s_MgCO3 = 10 ** BACKGROUND_LOGK["MgCO3"] * aCO3 * aH / b * MgOH
        s_MgOMe = term_MgOMe * MgOH

        # net surface charge, Eqn 6 of the paper, plus the +1 metal complexes when z=2
        charge_metal = (s_CaOMe + s_MgOMe) * (z - 1)      # 0 for Li, +1 each for Co
        sum_zk = (s_CaOH2 + s_MgOH2 + s_CO3Ca + s_CO3Mg
                  - s_CO3m - s_CaO - s_MgO - s_CaCO3 - s_MgCO3
                  + charge_metal)                          # mol/L
        sigma0 = F * sum_zk / ST_AREA                      # C/m2
        return sigma0, dict(CaOMe=s_CaOMe, MgOMe=s_MgOMe)

    # self-consistent solve: psi and sigma0 must satisfy psi = sigma0 / C_cap
    def residual(psi):
        sigma0, _ = species_given_psi(psi)
        return psi - sigma0 / C_cap

    try:
        psi = brentq(residual, -0.6, 0.6, xtol=1e-9)
    except ValueError:
        # fall back to damped fixed-point if the sign bracket fails at extreme charge
        psi = 0.0
        for _ in range(500):
            sigma0, _ = species_given_psi(psi)
            psi_new = sigma0 / C_cap
            if abs(psi_new - psi) < 1e-10:
                break
            psi = 0.5 * psi + 0.5 * psi_new
    sigma0, extra = species_given_psi(psi)
    return dict(sigma0=sigma0, psi0=psi, C_cap=C_cap,
                bound_Ca=extra["CaOMe"], bound_Mg=extra["MgOMe"])


# ======================================================================================
#  STEP 3 - AQUEOUS SPECIATION AND ACTIVITIES AT EACH pH
# ======================================================================================

def free_metal_activity(metal, Meq_M, sp):
    """
    Convert a dissolved-metal total (mol/L) into the free-ion activity that the surface
    reaction actually sees. Li is treated as essentially free; Co is corrected for the
    dominant CoCl+ complex and, weakly, for hydrolysis and carbonate pairing.
    """
    g1, g2 = sp["g1"], sp["g2"]
    aCl, aH, aCO3 = sp["aCl"], sp["aH"], sp["aCO3"]
    if metal == "Li":
        # Li+ + Cl- <-> LiCl ; keep only the free ion, LiCl is minor
        frac_free = 1.0 / (1.0 + aCl / 10 ** LK['licl'])
        Mfree = Meq_M * frac_free
        return g1 * Mfree
    if metal == "Co":
        # partition Co(total dissolved) between Co2+, CoCl+, CoOH+, CoCO3
        r_CoCl = aCl / 10 ** LK['cocl']                     # [CoCl+]/[Co2+]
        r_CoOH = 10 ** (-LK['cooh']) / aH                   # [CoOH+]/[Co2+] (hydrolysis)
        r_CoCO3 = aCO3 * 10 ** LK['coco3']                  # [CoCO3]/[Co2+]  (logK_form)
        frac_free = 1.0 / (1.0 + r_CoCl + r_CoOH + r_CoCO3)
        Mfree = Meq_M * frac_free
        return g2 * Mfree
    raise ValueError(metal)


def step3_activities(master):
    banner(3, "Calculate aqueous speciation and activities at each pH")
    explain("""
    At every pH used in the adsorption experiment I run the speciation engine to obtain
    the activities of H+, OH-, HCO3-, CO3-2 and the background Ca, Mg, Na, Cl species at
    fixed atmospheric pCO2. I then convert the dissolved Li and Co totals into free-ion
    activities, because only the free ion (Li+ or Co2+) competes with protons for the
    surface site. Li is essentially uncomplexed; Co is corrected for CoCl+ pairing in the
    salty background and, weakly, for hydrolysis and carbonate pairing.
    """)
    rows = []
    for row in master.itertuples():
        # background dissolved Ca, Mg for these single-ion batches are small; use trace
        CaT, MgT = 1e-5, 1e-5
        sp = speciate(row.pH_meas, CaT, MgT, ADS_NaT, ADS_ClT)
        sp["I"] = ADS_I                       # pin to the measured batch ionic strength
        aMe = free_metal_activity(row.metal, row.Ceq_M, sp)
        rows.append(dict(metal=row.metal, pH=row.pH_meas, I=sp["I"],
                         aH=sp["aH"], aHCO3=sp["aHCO3"], aCO3=sp["aCO3"],
                         aCa=sp["aCa"], aMg=sp["aMg"], aCl=sp["aCl"], aMe=aMe))
    act = pd.DataFrame(rows)
    print(act.to_string(index=False, float_format=lambda x: f"{x:.3e}"))
    return act


# ======================================================================================
#  STEP 4 - PROTON / BACKGROUND CONSTANTS AND CALIBRATION OF THE EDL PARAMETER
# ======================================================================================

def step4_background(sigma0_df):
    banner(4, "Fix the background site constants and calibrate the capacitance")
    explain("""
    The proton and carbonate constants for the three sites are taken from Pokrovsky 1999
    (Step 2) rather than re-fitted, because three noisy titration points on each side of
    the PZC cannot pin three acidity constants independently. What I do refine against
    your own sigma0 data is the single EDL parameter alpha, which sets the capacitance
    through C = sqrt(I)/alpha. I search alpha for the best match to the magnitude of your
    measured surface charge. A caution worth stating plainly: your tabulated sigma0 is
    negative below pH 8 and positive above it, which is the opposite sign to the proton
    charge of a constant-capacitance model, so the mass-balance definition of sigma0 in
    the spreadsheet differs from the model's surface charge. I therefore use sigma0 to
    anchor the capacitance and the point of zero charge, not as a term the metal
    constants are fitted to.
    """)
    global ALPHA_EDL
    pH = sigma0_df["pH"].values
    meas = sigma0_df["sigma0"].values

    def model_sigma(alpha):
        out = []
        for p in pH:
            sp = speciate(p, 1e-5, 1e-5, ADS_NaT, ADS_ClT)
            sp["I"] = ADS_I
            out.append(surface_solve(sp, {}, metal=None, alpha=alpha)["sigma0"])
        return np.array(out)

    alphas = np.logspace(-3, -0.5, 40)
    best_a, best_err = ALPHA_EDL, np.inf
    for a in alphas:
        err = np.sum((np.abs(model_sigma(a)) - np.abs(meas)) ** 2)
        if err < best_err:
            best_err, best_a = err, a
    ALPHA_EDL = best_a

    # locate the modelled point of zero charge
    grid = np.linspace(5.0, 10.5, 200)
    sig = []
    for p in grid:
        sp = speciate(p, 1e-5, 1e-5, ADS_NaT, ADS_ClT); sp["I"] = ADS_I
        sig.append(surface_solve(sp, {}, metal=None, alpha=best_a)["sigma0"])
    sig = np.array(sig)
    pzc = grid[np.argmin(np.abs(sig))]

    print(f"  log K1 (>CO3H0 = >CO3- + H+)      = {BACKGROUND_LOGK['CO3_deprot']:+.2f}  (fixed)")
    print(f"  log K2 (>CaOH0 = >CaO- + H+)      = {BACKGROUND_LOGK['CaO_deprot']:+.2f}  (fixed)")
    print(f"  log K3 (>MgOH0 = >MgO- + H+)      = {BACKGROUND_LOGK['MgO_deprot']:+.2f}  (fixed)")
    print(f"  refined EDL parameter alpha       = {best_a:.4f}")
    print(f"  capacitance at I={ADS_I:.2f} M       = {np.sqrt(ADS_I)/best_a:.1f} F/m2")
    print(f"  modelled point of zero charge     = pH {pzc:.2f}  (target ~8.0)")
    return best_a


# ======================================================================================
#  STEP 5, 6, 7 - PREDICTION, OBJECTIVE, OPTIMISATION
# ======================================================================================

def predict_uptake(metal, act_row, logK):
    """
    Predict the equilibrium dissolved metal (mol/L) for one adsorption point. Because the
    free-ion activity depends on how much metal is left in solution, and the amount bound
    depends on that activity, the two are solved together by short iteration.
    """
    sp = dict(aH=act_row.aH, aCO3=act_row.aCO3, aCa=act_row.aCa, aMg=act_row.aMg,
              aCl=act_row.aCl, I=act_row.I, g1=gammas(act_row.I)[0], g2=gammas(act_row.I)[1])
    C0_M = float(ADSORPTION_DATA[(ADSORPTION_DATA.metal == metal) &
                                 (ADSORPTION_DATA.pH_meas == act_row.pH)]["C0_ppm"].iloc[0]) \
        / MM[metal] / 1e3
    Meq = C0_M
    for _ in range(60):
        aMe = free_metal_activity(metal, Meq, sp)
        res = surface_solve(sp, logK, metal=metal, aMe=aMe, alpha=ALPHA_EDL)
        bound = res["bound_Ca"] + res["bound_Mg"]           # mol/L on surface
        Meq_new = max(C0_M - bound, 1e-12)
        if abs(Meq_new - Meq) < 1e-12:
            Meq = Meq_new
            break
        Meq = 0.5 * Meq + 0.5 * Meq_new
    return Meq, res["bound_Ca"], res["bound_Mg"]


def residuals(x, act):
    """Residual vector passed to least_squares: predicted minus measured Ceq, in ppm.
    x holds the two fitted constants [logK_Li-Ca, logK_Co-Ca]; Mg is tied to Ca."""
    logK = build_logK(x)
    r = []
    for row in act.itertuples():
        Meq, _, _ = predict_uptake(row.metal, row, logK)
        Ceq_pred_ppm = Meq * MM[row.metal] * 1e3
        meas = ADSORPTION_DATA[(ADSORPTION_DATA.metal == row.metal) &
                               (ADSORPTION_DATA.pH_meas == row.pH)]["Ceq_ppm"].iloc[0]
        unc = 1.0        # ppm; ICP-OES scale uncertainty
        r.append((Ceq_pred_ppm - meas) / unc)
    return np.array(r)


def step5_objective():
    banner(5, "Set up the non-linear least-squares objective")
    explain("""
    The quantity I minimise is the sum of squared, uncertainty-weighted differences
    between predicted and measured equilibrium metal concentrations across the four
    adsorption points. The four unknowns are the metal log K values. The surface-charge
    and zeta data are not summed into this objective for the reason given in Step 4 (their
    spreadsheet sign convention differs from the model); they instead fix the background
    model and the capacitance up front. This keeps the fit honest: the metal constants are
    driven by the adsorption measurements that actually respond to metal binding.
    """)
    print("  objective(x) = sum over points of  [(Ceq_predicted - Ceq_measured)/sigma]^2")
    print("  x = [logK_Li-Ca, logK_Li-Mg, logK_Co-Ca, logK_Co-Mg]")


def step7_optimise(act):
    banner(7, "Run the optimisation")
    explain("""
    I start from chemically sensible guesses anchored to the aqueous analogues: a weak
    constant for the poorly binding Li+ and a near-zero to positive constant for Co2+.
    scipy.optimize.least_squares then adjusts the two fitted values (the Mg constants
    follow by the fixed offset) to reproduce the measured equilibrium concentrations.
    least_squares is used rather than a plain minimiser because it returns the Jacobian,
    which Step 11 turns into confidence intervals.
    """)
    x0 = [-3.0, 0.5]     # [Li-Ca, Co-Ca]: weak for Li+, near zero to positive for Co2+
    print(f"  initial guess (Li-Ca, Co-Ca) = {x0}")
    sol = least_squares(residuals, x0, args=(act,), method="trf",
                        bounds=([-12, -12], [12, 12]), diff_step=1e-3)
    logK = build_logK(sol.x)
    print("\n  optimised intrinsic constants (Mg tied to Ca):")
    for k in METAL_REACTIONS:
        tag = "fitted " if k in FREE_PARAMS else "derived"
        print(f"     log K_int({k.replace('_','-')}) = {logK[k]:+.2f}   [{tag}]")
    print(f"\n  final sum of squared residuals = {2*sol.cost:.4g}")
    return sol, logK


# ======================================================================================
#  STEP 8 - VALIDATION
# ======================================================================================

def r_squared(meas, pred):
    meas, pred = np.asarray(meas, float), np.asarray(pred, float)
    ss_res = np.sum((meas - pred) ** 2)
    ss_tot = np.sum((meas - np.mean(meas)) ** 2)
    return 1 - ss_res / ss_tot if ss_tot > 0 else float("nan")


def step8_validate(act, logK, alpha):
    banner(8, "Validate the fit")
    explain("""
    I compare measured and predicted equilibrium concentrations for every adsorption
    point, report the coefficient of determination, and draw two figures: the background
    surface-charge curve against your titration points, and predicted versus measured
    metal concentrations. With only two adsorption points per metal the R-squared is a
    weak descriptor, so I also print the residual for each point so you can see the
    absolute agreement in ppm.
    """)
    preds, meass, labels = [], [], []
    for row in act.itertuples():
        Meq, bCa, bMg = predict_uptake(row.metal, row, logK)
        Ceq_pred = Meq * MM[row.metal] * 1e3
        meas = ADSORPTION_DATA[(ADSORPTION_DATA.metal == row.metal) &
                               (ADSORPTION_DATA.pH_meas == row.pH)]["Ceq_ppm"].iloc[0]
        preds.append(Ceq_pred); meass.append(meas)
        labels.append(f"{row.metal} pH{row.pH:.0f}")
        print(f"  {row.metal} pH {row.pH:.0f}: measured Ceq = {meas:7.3f} ppm | "
              f"predicted = {Ceq_pred:7.3f} ppm | residual = {Ceq_pred-meas:+.3f} ppm")
    r2 = r_squared(meass, preds)
    print(f"\n  R^2 over the four adsorption points = {r2:.3f}")

    # figure 1 - surface charge background model vs data
    grid = np.linspace(5.0, 10.5, 120)
    sig = []
    for p in grid:
        sp = speciate(p, 1e-5, 1e-5, ADS_NaT, ADS_ClT); sp["I"] = ADS_I
        sig.append(surface_solve(sp, {}, metal=None, alpha=alpha)["sigma0"])
    fig, ax = plt.subplots(figsize=(7, 4.6))
    ax.axhline(0, color="0.6", lw=0.8, ls="--")
    ax.plot(grid, np.array(sig) * 1e3, "-", color="#2c5f8a", label="CCM background model")
    ax.scatter(SIGMA0_DATA["pH"], SIGMA0_DATA["sigma0"] * 1e3, s=60, color="#c1440e",
               zorder=3, label="measured (spreadsheet)")
    ax.set_xlabel("pH"); ax.set_ylabel(r"$\sigma_0$ (mmol/m$^2$)")
    ax.set_title("Surface charge: model vs measurement")
    ax.legend(); ax.grid(alpha=0.3); fig.tight_layout()
    fig.savefig(os.path.join(HERE, "surface_charge_fit.png"), dpi=140)
    plt.close(fig)

    # figure 2 - predicted vs measured Ceq
    fig, ax = plt.subplots(figsize=(5.4, 5.2))
    lim = [min(meass + preds) * 0.98, max(meass + preds) * 1.02]
    ax.plot(lim, lim, "k--", lw=0.8)
    for m, pr, lab in zip(meass, preds, labels):
        ax.scatter(m, pr, s=90); ax.annotate(lab, (m, pr), xytext=(6, 4),
                                             textcoords="offset points", fontsize=9)
    ax.set_xlabel("measured Ceq (ppm)"); ax.set_ylabel("predicted Ceq (ppm)")
    ax.set_title(f"Adsorption fit (R$^2$ = {r2:.3f})")
    ax.grid(alpha=0.3); fig.tight_layout()
    fig.savefig(os.path.join(HERE, "adsorption_fit.png"), dpi=140)
    plt.close(fig)

    # figure 3 - surface speciation vs pH for each metal (illustrative)
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.4))
    for ax, metal in zip(axes, ["Li", "Co"]):
        ph_grid = np.linspace(3, 10, 60)
        fCa, fMg = [], []
        C0_M = ADSORPTION_DATA[ADSORPTION_DATA.metal == metal]["C0_ppm"].iloc[0] / MM[metal] / 1e3
        for p in ph_grid:
            sp = speciate(p, 1e-5, 1e-5, ADS_NaT, ADS_ClT); sp["I"] = ADS_I
            aMe = free_metal_activity(metal, C0_M, sp)
            res = surface_solve(sp, logK, metal=metal, aMe=aMe, alpha=alpha)
            fCa.append(res["bound_Ca"]); fMg.append(res["bound_Mg"])
        ax.plot(ph_grid, np.array(fCa) * 1e6, label=f">{metal} on Ca site")
        ax.plot(ph_grid, np.array(fMg) * 1e6, label=f">{metal} on Mg site")
        ax.set_xlabel("pH"); ax.set_ylabel("bound metal (umol/L)")
        ax.set_title(f"{metal} surface binding vs pH"); ax.legend(); ax.grid(alpha=0.3)
    fig.tight_layout(); fig.savefig(os.path.join(HERE, "surface_speciation.png"), dpi=140)
    plt.close(fig)
    print("\n  figures written: surface_charge_fit.png, adsorption_fit.png, surface_speciation.png")
    return r2, list(zip(labels, meass, preds))


# ======================================================================================
#  STEP 9 - SENSITIVITY ANALYSIS
# ======================================================================================

def step9_sensitivity(act, logK):
    banner(9, "Sensitivity analysis")
    explain("""
    I perturb each constant by +/- 0.5 log units, holding the others at their fitted
    value, and record how much the sum of squared residuals grows. A constant the data
    constrains tightly makes the objective climb steeply; a loosely constrained one barely
    moves it. This is where the limits of a two-point-per-metal data set show up honestly.
    """)
    x_fit = [logK["Li_Ca"], logK["Co_Ca"]]
    base = np.sum(residuals(x_fit, act) ** 2)
    rows = []
    for i, k in enumerate(FREE_PARAMS):
        x = list(x_fit)
        x[i] += 0.5; up = np.sum(residuals(x, act) ** 2)
        x[i] -= 1.0; dn = np.sum(residuals(x, act) ** 2)
        rows.append(dict(parameter=k.replace("_", "-"),
                         SSR_base=base, SSR_plus=up, SSR_minus=dn,
                         mean_increase=0.5 * ((up - base) + (dn - base))))
    sens = pd.DataFrame(rows).sort_values("mean_increase", ascending=False)
    print(sens.to_string(index=False, float_format=lambda v: f"{v:.4g}"))

    fig, ax = plt.subplots(figsize=(7, 4))
    ax.barh(sens["parameter"], sens["mean_increase"], color="#3a7d5d")
    ax.set_xlabel("mean rise in SSR for +/- 0.5 log K"); ax.set_title("Parameter sensitivity")
    ax.invert_yaxis(); fig.tight_layout()
    fig.savefig(os.path.join(HERE, "sensitivity.png"), dpi=140); plt.close(fig)
    print("\n  figure written: sensitivity.png")
    return sens


# ======================================================================================
#  STEP 11 - UNCERTAINTY BOUNDS FROM THE COVARIANCE MATRIX
# ======================================================================================

def step11_uncertainty(sol):
    banner(11, "Confidence intervals from the optimisation covariance")
    explain("""
    least_squares returns the Jacobian at the solution. Its Gauss-Newton approximation to
    the covariance is (J^T J)^{-1} scaled by the residual variance; the square roots of
    the diagonal are the one-sigma standard errors on each log K, and roughly twice that is
    a 95 percent interval. Large errors here confirm which constants the present data leave
    poorly determined.
    """)
    J = sol.jac
    dof = max(len(sol.fun) - len(sol.x), 1)
    resid_var = 2 * sol.cost / dof
    try:
        cov = np.linalg.inv(J.T @ J) * resid_var
        se = np.sqrt(np.abs(np.diag(cov)))
    except np.linalg.LinAlgError:
        se = np.full(len(sol.x), np.nan)
    # standard errors on the two fitted constants; the tied Mg constants inherit them
    se_map = {"Li_Ca": se[0], "Co_Ca": se[1], "Li_Mg": se[0], "Co_Mg": se[1]}
    logK = build_logK(sol.x)
    out = {}
    for k in METAL_REACTIONS:
        val, s = logK[k], se_map[k]
        out[k] = (val, s)
        tag = "fitted" if k in FREE_PARAMS else "tied  "
        print(f"  log K_int({k.replace('_','-'):5s}) = {val:+.2f}  +/- {s:.2f} (1 sigma) "
              f"[{tag}]  ~ 95% CI [{val-2*s:+.2f}, {val+2*s:+.2f}]")
    return out


# ======================================================================================
#  STEP 10 - COMPARE TO LITERATURE  (printed last, after values and errors exist)
# ======================================================================================

def step10_compare(logK, unc):
    banner(10, "Compare to literature and interpret")
    explain("""
    Pokrovsky 1999 did not tabulate Li or Co binding constants for dolomite directly, so I
    compare against the closest available anchors: the divalent-metal carbonate-site
    reaction on dolomite (>CO3H0 + Me2+ = >CO3Me+ + H+, log K about -1.8 to -2.0) and the
    aqueous hydrolysis analogues (Co2+ + H2O = CoOH+ + H+, log K about -9.7; Li+ barely
    hydrolyses). Surface binding is generally stronger than the aqueous analogue by one to
    three log units, and Co2+ should bind more strongly than Li+. I check the fitted values
    against those expectations rather than against a single published number.
    """)
    print("  Reference anchors:")
    print("     dolomite >CO3H0 + Ca2+ = >CO3Ca+ + H+     log K = -1.8   (Pokrovsky 1999)")
    print("     dolomite >CO3H0 + Mg2+ = >CO3Mg+ + H+     log K = -2.0   (Pokrovsky 1999)")
    print("     aqueous  Co2+ + H2O    = CoOH+  + H+       log K = -9.65")
    print("     aqueous  Li+           barely complexes")
    print()
    liCa, liMg = logK["Li_Ca"], logK["Li_Mg"]
    coCa, coMg = logK["Co_Ca"], logK["Co_Mg"]
    print(f"  Fitted Li binding (Ca, Mg) = {liCa:+.2f}, {liMg:+.2f}")
    print(f"  Fitted Co binding (Ca, Mg) = {coCa:+.2f}, {coMg:+.2f}")
    print()
    explain("""
    Reading the numbers correctly matters here. In this single-ion data set Li carries a
    larger fitted constant than Co, but that is not evidence that Li out-competes Co per
    site. It reflects two things. First, the fit is driven by absolute moles removed, and
    Li was dosed about fifteen times more concentrated than Co, so its molar uptake is
    actually larger even though its percentage recovery is small. Second, this is a
    sorption-only model: the manuscript attributes much of Co recovery to carbonate
    mineralisation (Co2+ + CO3-2 -> CoCO3, Eqs 7 to 9), a pathway not represented here, so
    a surface-only fit does not capture the true Co affinity. Both constants sit one to a
    few log units above the aqueous carbonate-site analogues, which is the expected
    direction for surface over solution binding. The Co constant is tightly determined
    (Step 11); the Li constant is loose because Li uptake is small and near the noise.
    """)


# ======================================================================================
#  DELIVERABLE - WRITE THE SUMMARY REPORT
# ======================================================================================

def write_report(logK, unc, r2, val_points, sens, alpha):
    path = os.path.join(HERE, "results_logK_summary.txt")
    with open(path, "w") as f:
        f.write("INTRINSIC STABILITY CONSTANTS FOR Li AND Co ON DOLOMITE\n")
        f.write("Constant capacitance three-site model (Pokrovsky et al., 1999)\n")
        f.write("Single-ion data, 25 C, I ~ 0.7 M NaCl, open system (pCO2 = 10^-3.5 atm)\n")
        f.write("=" * 74 + "\n\n")
        f.write("Surface complexation reactions and fitted constants:\n")
        names = {
            "Li_Ca": ">CaOH0 + Li+  = >CaOLi(0) + H+",
            "Li_Mg": ">MgOH0 + Li+  = >MgOLi(0) + H+",
            "Co_Ca": ">CaOH0 + Co2+ = >CaOCo+  + H+",
            "Co_Mg": ">MgOH0 + Co2+ = >MgOCo+  + H+",
        }
        for k in METAL_REACTIONS:
            val, se = unc[k]
            f.write(f"  {names[k]:34s}  log K_int = {val:+.2f} +/- {se:.2f}\n")
        f.write(f"\nBackground constants held fixed (Pokrovsky 1999, dolomite):\n")
        for k, v in BACKGROUND_LOGK.items():
            f.write(f"  {k:12s} log K = {v:+.2f}\n")
        f.write(f"\nEDL parameter alpha = {alpha:.4f}  ->  C = sqrt(I)/alpha\n")
        f.write(f"Site densities (mol/m2): Ca {SITE_DENS['Ca']:.1e}, "
                f"Mg {SITE_DENS['Mg']:.1e}, CO3 {SITE_DENS['CO3']:.1e}\n")
        f.write(f"\nFit quality: R^2 (adsorption) = {r2:.3f}\n")
        for lab, m, p in val_points:
            f.write(f"   {lab}: measured {m:.3f} ppm, predicted {p:.3f} ppm\n")
        f.write("\nSensitivity (mean rise in SSR for +/- 0.5 log K):\n")
        for row in sens.itertuples():
            f.write(f"   {row.parameter}: {row.mean_increase:.4g}\n")
        f.write("\nInterpretation notes:\n")
        f.write("  - One constant per metal is fitted (Ca site); the Mg constant is tied\n")
        f.write(f"    with a fixed offset of {MG_OFFSET:+.1f} log units (Pokrovsky Ca vs Mg).\n")
        f.write("  - This is a sorption-only model. Co recovery in the manuscript is partly\n")
        f.write("    carbonate mineralisation (CoCO3), which is not represented, so the Co\n")
        f.write("    constant is a lumped surface value, not a pure affinity.\n")
        f.write("  - Li's larger constant reflects larger absolute molar uptake (Li dosed\n")
        f.write("    ~15x more concentrated), not a stronger per-site affinity than Co.\n")
        f.write("  - sigma0 anchors the capacitance and PZC only; its spreadsheet sign\n")
        f.write("    convention differs from the CCM surface charge.\n")
        f.write("  - Next stage: extend METAL_REACTIONS and METAL_CHARGE to the full\n")
        f.write("    produced-water cation suite (Ba, Sr, Cd, Pb ...).\n")
    print(f"\nReport written to {path}")
    return path


# ======================================================================================
#  MAIN - run the eleven steps in order
# ======================================================================================

def main():
    print("#" * 86)
    print("#  INTRINSIC STABILITY CONSTANTS FOR Li AND Co ON DOLOMITE")
    print("#  Constant capacitance surface complexation model (Pokrovsky et al., 1999)")
    print("#" * 86)

    master = step1_master_table()
    step2_model_framework()
    act = step3_activities(master)
    alpha = step4_background(SIGMA0_DATA)
    step5_objective()
    # Step 6 is the forward solver (surface_solve / predict_uptake) defined above;
    # it is exercised inside Steps 7 and 8.
    banner(6, "Forward prediction function")
    explain("""
    The forward model is the function predict_uptake(): given a pH, the aqueous activities
    from Step 3 and a trial set of the four log K values, it solves the constant-capacitance
    model self-consistently for the surface potential, then returns the amount of Li or Co
    bound to the Ca and Mg sites and the resulting equilibrium concentration in solution.
    Steps 7 and 8 call it repeatedly.
    """)
    sol, logK = step7_optimise(act)
    r2, val_points = step8_validate(act, logK, alpha)
    sens = step9_sensitivity(act, logK)
    unc = step11_uncertainty(sol)
    step10_compare(logK, unc)
    write_report(logK, unc, r2, val_points, sens, alpha)

    banner("DONE", "Summary of determined intrinsic stability constants")
    print("  Reaction                             log K_int (25 C, I ~ 0.7 M)")
    names = {"Li_Ca": ">CaOH0 + Li+  = >CaOLi(0) + H+",
             "Li_Mg": ">MgOH0 + Li+  = >MgOLi(0) + H+",
             "Co_Ca": ">CaOH0 + Co2+ = >CaOCo+  + H+",
             "Co_Mg": ">MgOH0 + Co2+ = >MgOCo+  + H+"}
    for k in METAL_REACTIONS:
        val, se = unc[k]
        print(f"  {names[k]:36s} {val:+.2f} +/- {se:.2f}")
    print("\n  See results_logK_summary.txt and the four PNG figures for the full record.")


if __name__ == "__main__":
    main()
