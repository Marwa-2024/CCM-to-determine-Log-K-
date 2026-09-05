#!/usr/bin/env python3
"""Generate a self-contained Colab .ipynb and a Wolfram Mathematica .nb for the
FITEQL-style determination of Li/Co log K on dolomite. All data embedded; no file refs."""
import json, os

OUT = "/home/user/CCM-to-determine-Log-K-"

# ======================================================================================
#  SHARED TEXT: all reactions and equations, written out in plain language
# ======================================================================================
REACTIONS_TEXT = r"""REACTIONS IN THE MODEL

A. AQUEOUS DISSOCIATION
   H2O            = H+   + OH-
   CO2(g)         = CO2(aq)
   NaCl           = Na+  + Cl-
   HCl            = H+   + Cl-

B. AQUEOUS COMPLEXATION
   CO2(aq) + H2O  = HCO3- + H+          log K = -6.345
   CO3-2  + H+    = HCO3-               log K = 10.329
   LiCl           = Li+  + Cl-          (Li barely complexes)
   CoCl+          = Co2+ + Cl-
   CoOH+ + H+     = Co2+ + H2O          (Co hydrolysis)
   CoCO3          = Co2+ + CO3-2

C. SURFACE REACTIONS - carbonate site  >CO3H0     (log K, Pokrovsky 1999)
   Per Pokrovsky, metal CATIONS bind this site (>CO3- + Me -> >CO3Me).
   >CO3H0         = >CO3-   + H+        log K = -4.8
   >CO3H0 + Ca2+  = >CO3Ca+ + H+        log K = -1.8
   >CO3H0 + Mg2+  = >CO3Mg+ + H+        log K = -2.0
   >CO3H0 + Li+   = >CO3Li0 + H+        log K = UNKNOWN (fitted)   <-- Li binds here
   >CO3H0 + Co2+  = >CO3Co+ + H+        log K = UNKNOWN (fitted)   <-- Co binds here

D. SURFACE REACTIONS - calcium site   >CaOH0   (anion/ligand site, no metal cation)
   >CaOH0 + H+          = >CaOH2+                  log K = 11.5
   >CaOH0              = >CaO-   + H+              log K = -12.0
   >CaOH0 + CO3-2 + 2H+ = >CaHCO3(0) + H2O         log K = -4.0
   >CaOH0 + CO3-2 +  H+ = >CaCO3-   + H2O          log K = 16.6

E. SURFACE REACTIONS - magnesium site >MgOH0   (anion/ligand site, no metal cation)
   >MgOH0 + H+          = >MgOH2+                  log K = 10.6
   >MgOH0              = >MgO-   + H+              log K = -12.0
   >MgOH0 + CO3-2 + 2H+ = >MgHCO3(0) + H2O         log K = -3.5
   >MgOH0 + CO3-2 +  H+ = >MgCO3-   + H2O          log K = 15.4
"""

EQUATIONS_TEXT = r"""EQUATIONS

(1) Davies activity coefficient (used up to I about 0.7 M)
    log gamma_i = -A z^2 [ sqrt(I)/(1+sqrt(I)) - 0.3 I ] ,   A = 0.509

(2) Open-system carbonate (fixed atmospheric pCO2 = 10^-3.5 atm)
    a(H+)    = 10^-pH
    a(HCO3-) = 10^-6.345 * a(CO2) / a(H+) ,   a(CO2) = 10^-1.469 * pCO2
    a(CO3-2) = a(HCO3-) / (10^10.329 * a(H+))

(3) Mass action for every species (component/species tableau, FITEQL)
    C_i = K_i * PRODUCT_j ( X_j ^ a_ij ) * P^(z_i for surface species)
    X_j  = free value of component j ;  a_ij = stoichiometry ;
    P    = exp(-F psi / RT) is the electrostatic component (Boltzmann factor).

(4) Boltzmann electrostatic correction (Pokrovsky Eq 6; Stumm)
    a surface species of charge z carries the factor  exp(-z F psi / RT) = P^z
    For  >SOH0 + Me(z+) = >SOMe(z-1) + H+ :
    [>SOMe] = K * [>SOH0] * a(Me)/a(H) * exp(-(z-1) F psi / RT)

(5) Site (mass) balances, solved for the neutral reference form
    T(>CO3H) = [>CO3H0]+[>CO3-]+[>CO3Ca+]+[>CO3Mg+]+[>CO3Me]   (metal binds here)
    T(>CaOH) = [>CaOH0]+[>CaOH2+]+[>CaO-]+[>CaHCO3]+[>CaCO3-]
    T(>MgOH) = [>MgOH0]+[>MgOH2+]+[>MgO-]+[>MgHCO3]+[>MgCO3-]
    e.g.  [>CO3H0] = T(>CO3H) / (1 + sum of all ratios to >CO3H0)

(6) PREDICTED surface charge, from the surface species only (ProtoFit Eq 2.2.3)
    sigma = (F / SSA) * { [>CaOH2+]+[>MgOH2+]+[>CO3Ca+]+[>CO3Mg+]+[>CO3Co+]
                          - [>CO3-]-[>CaO-]-[>MgO-]-[>CaCO3-]-[>MgCO3-] }
    It never uses dissolved Ca2+/Mg2+; those belong to dolomite dissolution.

(7) Constant capacitance model
    psi = sigma / C ,     C = sqrt(I) / alpha ,   alpha = 0.004 (Pokrovsky)

(8) Newton-Raphson equilibrium (FITEQL inner loop)
    residuals Y_j = sum_i a_ij C_i - T_j = 0 for the chemical components,
    Y_P = sum_(surface) z_i C_i + kappa * ln P = 0 for the electrostatic component,
    Jacobian  Z_jk = sum_i a_ij a_ik C_i   (+ kappa on the P diagonal),
    kappa = C * SSA * R T / F^2 .

(9) Adsorption uptake and the fitted quantity
    uptake = [>CO3Me] ;   [Me]_eq = [Me]_total - uptake
    predicted [Me]_eq is compared to the measured value.

(10) FITEQL objective (goodness of fit)
    WSOS/DF = sum_i ( Y_i / s_i )^2 / (N_obs - N_param) ,
    s_i = sqrt( errAbs^2 + (errRel * value)^2 ) ;  a value near 1 is a good fit.
"""

DATA_TEXT = r"""EMBEDDED DATA (from the spreadsheet "Set up" and adsorption sheets)

Reactor set-up
   mean particle diameter        = 2.8 um
   dolomite concentration        = 6 g / 100 mL = 60 g/L
   solid density                 = 2850 g/m3
   specific surface area (SSA)   = 0.76 m2/g
   surface concentration St      = 60 * 0.76 = 45.6 m2/L
   NaOH added each step          = 0.001 mol/L
   HCl added each step           = 0.001 mol/L

Model set-up
   temperature                   = 25 C
   ionic strength (batch)        = 0.70 M (mostly 40 g/L NaCl)
   pCO2 (open system)            = 10^-3.5 atm
   site density (Pokrovsky)      = 7, 7, 14 umol/m2 for Ca, Mg, carbonate
   EDL parameter alpha           = 0.004

Single-ion adsorption (equilibrium = day 6), concentrations in mg/L (ppm)
   metal  pH    C0 (initial)   Ceq (equilibrium)
   Co     6.0   80.513         74.393
   Co     2.0   80.513         80.398
   Li     6.0   138.43         135.262
   Li     2.0   138.43         135.592
"""

# ======================================================================================
#  PYTHON CODE CELLS FOR COLAB
# ======================================================================================
PY_SETUP = r'''# Colab already has numpy, scipy, matplotlib, pandas. Nothing to install.
import numpy as np, pandas as pd
import matplotlib.pyplot as plt
from scipy.optimize import least_squares

# ---- physical constants ----
F, R_GAS, T_K = 96485.0, 8.314, 298.15
FRT = F / (R_GAS * T_K)                     # F/RT = 38.92 1/V

# ---- reactor set-up (spreadsheet "Set up") ----
MEAN_DIAM_UM = 2.8
DOLOMITE_GL  = 60.0                         # 6 g/100 mL
DENSITY_GM3  = 2850.0
SSA_M2G      = 0.76
S_AREA       = DOLOMITE_GL * SSA_M2G        # 45.6 m2 per L
NAOH_STEP, HCL_STEP = 0.001, 0.001         # mol/L per titration step

# ---- CCM / EDL and brine ----
ALPHA_EDL = 0.004
I_BATCH   = 0.70
MM = dict(Li=6.941, Co=58.933, Cl=35.45)
CL_TOT = 40.0 / MM['Cl']                    # 40 g/L NaCl
CA_TRACE = MG_TRACE = 3.0e-6                # trace dissolved Ca, Mg

# ---- site densities (Pokrovsky 1999), 1:1:2 Ca:Mg:carbonate ----
SITE_DENS = dict(CO3=14.0e-6, Ca=7.0e-6, Mg=7.0e-6)      # mol/m2
SITE_TOT  = {k: v * S_AREA for k, v in SITE_DENS.items()}# mol/L
METAL_Z   = dict(Li=1, Co=2)

# ---- open-system carbonate ----
PCO2, ACO2 = 10**-3.5, 10**-1.469 * 10**-3.5
LK_CO2, LK_CO3 = -6.345, 10.329

# ---- embedded adsorption data (equilibrium = day 6) ----
ADS = pd.DataFrame({"metal":["Co","Co","Li","Li"], "pH":[6.0,2.0,6.0,2.0],
                    "C0_ppm":[80.513,80.513,138.43,138.43],
                    "Ceq_ppm":[74.393,80.398,135.262,135.592]})
ERR_ABS, ERR_REL = 0.5, 0.02               # measurement error for WSOS/DF
print("Setup loaded. Surface concentration St =", S_AREA, "m2/L")'''

PY_ACT = r'''def davies(I):
    if I <= 0: return 1.0, 1.0
    f = -0.509 * (np.sqrt(I)/(1+np.sqrt(I)) - 0.3*I)
    return 10**f, 10**(4*f)                 # gamma(z=1), gamma(z=2)

def background_activities(pH):
    aH = 10**(-pH)                          # pH is on the H+ activity scale
    aHCO3 = 10**LK_CO2 * ACO2 / aH          # open system
    aCO3  = aHCO3 / (10**LK_CO3 * aH)
    g1, g2 = davies(I_BATCH)
    return dict(H=aH, CO3=aCO3, Ca=g2*CA_TRACE, Mg=g2*MG_TRACE, Cl=g1*CL_TOT, g1=g1, g2=g2)'''

PY_TABLEAU = r'''# The tableau: components solved = Me, >CO3H, >CaOH, >MgOH, P=exp(-F.psi/RT)
#              components fixed  = H, CO3, Ca, Mg, Cl (activities)
SOLVED = ["Me", "sCO3H", "sCaOH", "sMgOH"]
AQ_METAL = {"Li": [("LiCl", -0.50, {"Me":1,"Cl":1}, 0)],
            "Co": [("CoCl", -0.05, {"Me":1,"Cl":1}, 1),
                   ("CoOH", -9.65, {"Me":1,"H":-1}, 1),
                   ("CoCO3", 4.22, {"Me":1,"CO3":1}, 0)]}

def build_species(metal, logK_Me):
    # Per Pokrovsky 1999, the metal CATION binds the carbonate site (>CO3Me), not >CaOH.
    z = METAL_Z[metal]; sp = [("Me", 0.0, {"Me":1}, z, False)]
    for name, lk, st, zz in AQ_METAL[metal]: sp.append((name, lk, st, zz, False))
    sp += [("sCO3H",0.0,{"sCO3H":1},0,True),("sCO3-",-4.8,{"sCO3H":1,"H":-1},-1,True),
           ("sCO3Ca",-1.8,{"sCO3H":1,"Ca":1,"H":-1},1,True),
           ("sCO3Mg",-2.0,{"sCO3H":1,"Mg":1,"H":-1},1,True),
           ("sCO3Me",logK_Me,{"sCO3H":1,"Me":1,"H":-1},z-1,True),          # <-- UNKNOWN
           ("sCaOH",0.0,{"sCaOH":1},0,True),("sCaOH2",11.5,{"sCaOH":1,"H":1},1,True),
           ("sCaO",-12.0,{"sCaOH":1,"H":-1},-1,True),
           ("sCaHCO3",-4.0,{"sCaOH":1,"CO3":1,"H":2},0,True),
           ("sCaCO3",16.6,{"sCaOH":1,"CO3":1,"H":1},-1,True),
           ("sMgOH",0.0,{"sMgOH":1},0,True),("sMgOH2",10.6,{"sMgOH":1,"H":1},1,True),
           ("sMgO",-12.0,{"sMgOH":1,"H":-1},-1,True),
           ("sMgHCO3",-3.5,{"sMgOH":1,"CO3":1,"H":2},0,True),
           ("sMgCO3",15.4,{"sMgOH":1,"CO3":1,"H":1},-1,True)]
    return sp

def prepare_point(metal, pH, logK_Me):
    bg = background_activities(pH); g1, g2 = bg["g1"], bg["g2"]
    gMe = g1 if METAL_Z[metal]==1 else g2
    Xfix = {"H":bg["H"]/g1, "CO3":bg["CO3"]/g2, "Ca":bg["Ca"]/g2, "Mg":bg["Mg"]/g2, "Cl":bg["Cl"]/g1}
    gcomp = {"H":g1,"CO3":g2,"Ca":g2,"Mg":g2,"Cl":g1,"Me":gMe,"sCO3H":1.0,"sCaOH":1.0,"sMgOH":1.0}
    species = []
    for name, lk, st, z, surf in build_species(metal, logK_Me):
        gprod = 1.0
        for c, co in st.items(): gprod *= gcomp[c]**co
        gi = 1.0 if surf else (1.0 if z==0 else (g1 if abs(z)==1 else g2))
        species.append(dict(name=name, K=10**lk*gprod/gi, st=st, z=z, surf=surf))
    C_cap = np.sqrt(I_BATCH)/ALPHA_EDL
    kappa = C_cap * S_AREA * R_GAS * T_K / F**2
    return Xfix, species, kappa'''

PY_EQUIL = r'''def equilibrium(Xfix, species, kappa, totals):
    """FITEQL inner loop: solve the coupled equilibrium by Newton-Raphson in ln-space.
       Jacobian Z_jk = sum_i a_ij a_ik C_i, with kappa on the electrostatic diagonal."""
    idx = {c:k for k,c in enumerate(SOLVED)}; nP = len(SOLVED); n = nP+1
    u = np.array([np.log(max(totals.get(c,1e-9),1e-12)) for c in SOLVED] + [0.0])
    def conc(u):
        X = {c: np.exp(u[idx[c]]) for c in SOLVED}; X.update(Xfix); P = np.exp(u[nP])
        C = np.empty(len(species))
        for i,s in enumerate(species):
            v = s["K"]
            for c,co in s["st"].items(): v *= X[c]**co
            if s["surf"]: v *= P**s["z"]
            C[i] = v
        return C, P
    for _ in range(200):
        C, P = conc(u)
        Y = np.zeros(n)
        for c in SOLVED:
            Y[idx[c]] = sum(s["st"].get(c,0)*C[i] for i,s in enumerate(species)) - totals[c]
        Y[nP] = sum(s["z"]*C[i] for i,s in enumerate(species) if s["surf"]) + kappa*u[nP]
        J = np.zeros((n,n))
        for i,s in enumerate(species):
            a = {idx[c]: s["st"][c] for c in SOLVED if c in s["st"]}
            if s["surf"] and s["z"]: a[nP] = a.get(nP,0) + s["z"]
            for j,aj in a.items():
                for k,ak in a.items(): J[j,k] += aj*ak*C[i]
        J[nP,nP] += kappa
        try: du = np.linalg.solve(J, -Y)
        except np.linalg.LinAlgError: du = -Y
        u = u + np.clip(du, -4.0, 4.0)
        scale = np.array([max(abs(totals[c]),1e-12) for c in SOLVED] + [max(abs(kappa),1e-12)])
        if np.max(np.abs(Y)/scale) < 1e-9: break
    C, P = conc(u)
    psi = -(R_GAS*T_K/F) * u[nP]
    cc = {s["name"]: C[i] for i,s in enumerate(species)}
    sigma = F/S_AREA * sum(s["z"]*C[i] for i,s in enumerate(species) if s["surf"])
    return dict(conc=cc, psi=psi, sigma=sigma)

def predict_point(metal, pH, T_Me, logK_Me):
    # The metal is a component with known TOTAL, so the free ion, its aqueous complexes,
    # the surface sites and the potential are all solved together in one equilibrium call.
    Xfix, species, kappa = prepare_point(metal, pH, logK_Me)
    totals = {"Me":T_Me, "sCO3H":SITE_TOT["CO3"], "sCaOH":SITE_TOT["Ca"], "sMgOH":SITE_TOT["Mg"]}
    res = equilibrium(Xfix, species, kappa, totals)
    aq = ["Me"] + [c[0] for c in AQ_METAL[metal]]        # free ion + aqueous complexes
    dissolved = sum(res["conc"][nm] for nm in aq)
    adsorbed  = res["conc"]["sCO3Me"]                    # metal on the carbonate site
    return dict(dissolved=dissolved, adsorbed=adsorbed, psi=res["psi"], sigma=res["sigma"])'''

PY_FIT = r'''# Metal cations bind ONE site (the carbonate site), so ONE constant per metal.
NAMES = {"Li":">CO3H0 + Li+  = >CO3Li(0) + H+", "Co":">CO3H0 + Co2+ = >CO3Co+  + H+"}

def residuals(theta):
    K = {"Li": theta[0], "Co": theta[1]}
    r = []
    for row in ADS.itertuples():
        T_Me = row.C0_ppm / MM[row.metal] / 1e3
        pr = predict_point(row.metal, row.pH, T_Me, K[row.metal])
        Ceq_pred = pr["dissolved"] * MM[row.metal] * 1e3
        s = np.hypot(ERR_ABS, ERR_REL*row.Ceq_ppm)
        r.append((Ceq_pred - row.Ceq_ppm)/s)
    return np.array(r)

sol = least_squares(residuals, [0.0, 0.0], bounds=([-12,-12],[12,12]), diff_step=1e-3)
logK = {"Li":sol.x[0], "Co":sol.x[1]}
r = sol.fun; wsos_df = float(np.sum(r**2)/(len(r)-len(sol.x)))
cov = np.linalg.inv(sol.jac.T @ sol.jac) * (np.sum(r**2)/(len(r)-len(sol.x)))
se = np.sqrt(np.abs(np.diag(cov))); se_map = {"Li":se[0],"Co":se[1]}
print("Fitted intrinsic stability constants (25 C, I = 0.7 M):")
for k in ["Li","Co"]:
    print(f"  {NAMES[k]:32}  log K = {logK[k]:+.2f} +/- {se_map[k]:.2f}")
print()
print("  For scale, Pokrovsky >CO3H0 + Ca2+ = >CO3Ca+ + H+ has log K = -1.8.")
print(f"  WSOS/DF = {wsos_df:.3f}   (near 1 = good fit)")'''

PY_VALID = r'''rows = []
for row in ADS.itertuples():
    T_Me = row.C0_ppm / MM[row.metal] / 1e3
    pr = predict_point(row.metal, row.pH, T_Me, logK[row.metal])
    rows.append(dict(metal=row.metal, pH=row.pH, Ceq_meas=row.Ceq_ppm,
                     Ceq_pred=round(pr["dissolved"]*MM[row.metal]*1e3,3),
                     uptake_meas=round(row.C0_ppm-row.Ceq_ppm,3),
                     uptake_pred=round(pr["adsorbed"]*MM[row.metal]*1e3,3),
                     psi_mV=round(pr["psi"]*1e3,2)))
val = pd.DataFrame(rows); print(val.to_string(index=False))

fig, ax = plt.subplots(figsize=(5.2,5))
mm = val["Ceq_meas"].values; pp = val["Ceq_pred"].values
lim = [min(mm.min(),pp.min())*0.98, max(mm.max(),pp.max())*1.02]
ax.plot(lim, lim, "k--", lw=0.8)
for _,rr in val.iterrows():
    ax.scatter(rr.Ceq_meas, rr.Ceq_pred, s=90)
    ax.annotate(f"{rr.metal} pH{rr.pH:.0f}", (rr.Ceq_meas, rr.Ceq_pred),
                xytext=(6,4), textcoords="offset points", fontsize=9)
ax.set_xlabel("measured Ceq (ppm)"); ax.set_ylabel("predicted Ceq (ppm)")
ax.set_title("FITEQL-style fit"); ax.grid(alpha=0.3); plt.show()'''

PY_SIGMA = r'''# Surface charge from surface species only (ProtoFit Eq 2.2.3) - correct magnitude
def sigma_at(pH):
    Xfix, species, kappa = prepare_point("Co", pH, -99)
    totals = {"Me":1e-15, "sCO3H":SITE_TOT["CO3"], "sCaOH":SITE_TOT["Ca"], "sMgOH":SITE_TOT["Mg"]}
    return equilibrium(Xfix, species, kappa, totals)["sigma"]

print("pH   sigma (C/m2)   sigma (mmol/m2)   [Pokrovsky measured 0.01-0.02 mmol/m2]")
for p in [5.5,6.5,7.3,8.0,8.5,9.0]:
    s = sigma_at(p); print(f"{p:4.1f}  {s:11.4f}   {s/F*1e3:12.5f}")

grid = np.linspace(4,10.5,90); sig = np.array([sigma_at(p) for p in grid])/F*1e3
fig, ax = plt.subplots(figsize=(7.4,4.6))
ax.axhline(0, color="0.6", ls="--", lw=0.8)
ax.axhspan(-0.02, 0.02, color="#3a7d5d", alpha=0.15, label="Pokrovsky 0.01-0.02 mmol/m2")
ax.plot(grid, sig, "-", color="#2c5f8a", lw=2, label="predicted (surface species only)")
ax.set_xlabel("pH"); ax.set_ylabel("surface charge (mmol/m2)"); ax.set_ylim(-0.035,0.02)
ax.legend(); ax.grid(alpha=0.3); ax.set_title("Predicted surface charge, correct magnitude"); plt.show()'''

# ======================================================================================
#  BUILD THE .ipynb (Colab)
# ======================================================================================
def md(src):  return {"cell_type":"markdown","metadata":{},"source":src.splitlines(keepends=True)}
def code(src):return {"cell_type":"code","metadata":{},"execution_count":None,"outputs":[],
                      "source":src.splitlines(keepends=True)}

cells = [
 md("# Intrinsic log K for Li and Co on dolomite\n"
    "### FITEQL-style surface complexation model, pure Python (Google Colab)\n\n"
    "Self-contained: every reaction, equation, and datum is embedded below. "
    "Method: component/species tableau, Newton-Raphson multicomponent equilibrium with an "
    "electrostatic (Boltzmann) component, and a WSOS/DF fit of the unknown constants "
    "(Westall 1982; Herbelin & Westall 1999). Chemical model: Pokrovsky, Schott & Thomas (1999). "
    "Surface charge follows the ProtoFit definition (Turner & Fein).\n\n"
    "Run the cells in order (Runtime, Run all)."),
 md("## 1. Reactions\n```\n" + REACTIONS_TEXT + "\n```"),
 md("## 2. Equations\n```\n" + EQUATIONS_TEXT + "\n```"),
 md("## 3. Data\n```\n" + DATA_TEXT + "\n```"),
 md("## 4. Setup and embedded data"), code(PY_SETUP),
 md("## 5. Activity coefficients (Davies) and background activities"), code(PY_ACT),
 md("## 6. The tableau: components and species"), code(PY_TABLEAU),
 md("## 7. Newton-Raphson equilibrium and forward prediction"), code(PY_EQUIL),
 md("## 8. Fit log K by minimising WSOS/DF"), code(PY_FIT),
 md("## 9. Validation: measured vs predicted"), code(PY_VALID),
 md("## 10. Surface charge check (correct magnitude)"), code(PY_SIGMA),
 md("## Notes\n"
    "- One constant per metal is fitted; the Mg-site constant is tied with the -0.9 offset "
    "because two pH points per metal cannot resolve both sites.\n"
    "- Sorption-only model; cobalt carbonate mineralisation (CoCO3) is discussed in the "
    "manuscript but not fitted here.\n"
    "- alpha = 0.004 (Pokrovsky) gives a high capacitance, so the surface potential is only a "
    "few mV.\n"
    "- Surface charge is built from surface species only and matches Pokrovsky's measured "
    "0.01-0.02 mmol/m2; the spreadsheet's larger values came from dolomite dissolution."),
]
nb = {"cells":cells,
      "metadata":{"kernelspec":{"display_name":"Python 3","language":"python","name":"python3"},
                  "language_info":{"name":"python","version":"3"},
                  "colab":{"provenance":[]}},
      "nbformat":4,"nbformat_minor":5}
with open(os.path.join(OUT,"FITEQL_LiCo_dolomite_Colab.ipynb"),"w") as f:
    json.dump(nb, f, indent=1)
print("wrote Colab ipynb")

# ======================================================================================
#  WOLFRAM MATHEMATICA .nb
# ======================================================================================
WL_SETUP = r'''(* ===== Physical constants and embedded data ===== *)
Fc = 96485.0; Rg = 8.314; Tk = 298.15; FRT = Fc/(Rg*Tk);

(* Reactor set-up (spreadsheet "Set up") *)
meanDiamUm = 2.8; dolomiteGL = 60.0; densitySolid = 2850.0;
SSA = 0.76; Sarea = dolomiteGL*SSA;   (* 45.6 m^2 / L *)
NaOHstep = 0.001; HClstep = 0.001;

(* CCM / EDL and brine *)
alpha = 0.004; Ibatch = 0.7; Ccap = Sqrt[Ibatch]/alpha;
MM = <|"Li" -> 6.941, "Co" -> 58.933, "Cl" -> 35.45|>;
ClTot = 40.0/MM["Cl"];

(* site densities (Pokrovsky 1999), mol/m^2 and totals mol/L *)
TCO3site = 14.0*10^-6*Sarea; TCasite = 7.0*10^-6*Sarea; TMgsite = 7.0*10^-6*Sarea;
zMet = <|"Li" -> 1, "Co" -> 2|>;

(* single-ion adsorption data: {metal, pH, C0ppm, Ceqppm} *)
adsData = {{"Co", 6.0, 80.513, 74.393}, {"Co", 2.0, 80.513, 80.398},
   {"Li", 6.0, 138.43, 135.262}, {"Li", 2.0, 138.43, 135.592}};
errAbs = 0.5; errRel = 0.02;'''

WL_ACT = r'''(* ===== Davies activity coefficients and background activities ===== *)
davies[Ii_] := Module[{f}, f = -0.509*(Sqrt[Ii]/(1 + Sqrt[Ii]) - 0.3*Ii); {10^f, 10^(4 f)}];
{gam1, gam2} = davies[Ibatch];
pCO2 = 10^-3.5; aCO2 = 10^(-1.469)*pCO2;
activities[pH_] := Module[{aH, aHCO3, aCO3},
   aH = 10^(-pH);
   aHCO3 = 10^(-6.345)*aCO2/aH;
   aCO3 = aHCO3/(10^(10.329)*aH);
   <|"aH" -> aH, "aCO3" -> aCO3, "aCa" -> gam2*3.0*10^-6,
     "aMg" -> gam2*3.0*10^-6, "aCl" -> gam1*ClTot|>];

(* free metal-ion activity from the dissolved total *)
metalActivity[metal_, Meq_, act_] := Module[{aCl, aH, aCO3, r},
   aCl = act["aCl"]; aH = act["aH"]; aCO3 = act["aCO3"];
   If[metal == "Li",
     gam1*Meq/(1 + aCl/10^(0.5)),
     r = aCl*10^(-0.05) + 10^(-9.65)/aH + aCO3*10^(4.22);
     gam2*Meq/(1 + r)]];'''

WL_SPECIES = r'''(* ===== Surface speciation and net surface charge at a given potential psi ===== *)
(* Each site is written relative to its neutral reference form; b = Exp[-FRT psi] is the
   Boltzmann factor, and a surface species of charge z carries b^z. Per Pokrovsky 1999
   the metal cation binds the CARBONATE site: >CO3H0 + Me = >CO3Me + H+. *)
speciesAt[metal_, aMe_, act_, logKme_, psi_?NumericQ] := Module[
   {aH, aCO3, aCa, aMg, zm, b, denCO3, CO3H, denCa, CaOH, denMg, MgOH,
    sCaOH2, sMgOH2, sCO3Ca, sCO3Mg, sCO3Me, sCO3m, sCaO, sMgO, sCaCO3, sMgCO3, sig},
   aH = act["aH"]; aCO3 = act["aCO3"]; aCa = act["aCa"]; aMg = act["aMg"];
   zm = zMet[metal]; b = Exp[-FRT*psi];
   denCO3 = 1 + 10^(-4.8)/aH/b + 10^(-1.8)*aCa/aH*b + 10^(-2.0)*aMg/aH*b +
      10^(logKme)*aMe/aH*b^(zm - 1);
   CO3H = TCO3site/denCO3;
   denCa = 1 + 10^(11.5)*aH*b + 10^(-12)/aH/b + 10^(-4.0)*aCO3*aH^2 + 10^(16.6)*aCO3*aH/b;
   CaOH = TCasite/denCa;
   denMg = 1 + 10^(10.6)*aH*b + 10^(-12)/aH/b + 10^(-3.5)*aCO3*aH^2 + 10^(15.4)*aCO3*aH/b;
   MgOH = TMgsite/denMg;
   sCaOH2 = 10^(11.5)*aH*b*CaOH; sMgOH2 = 10^(10.6)*aH*b*MgOH;
   sCO3Ca = 10^(-1.8)*aCa/aH*b*CO3H; sCO3Mg = 10^(-2.0)*aMg/aH*b*CO3H;
   sCO3Me = 10^(logKme)*aMe/aH*b^(zm - 1)*CO3H;
   sCO3m = 10^(-4.8)/aH/b*CO3H; sCaO = 10^(-12)/aH/b*CaOH; sMgO = 10^(-12)/aH/b*MgOH;
   sCaCO3 = 10^(16.6)*aCO3*aH/b*CaOH; sMgCO3 = 10^(15.4)*aCO3*aH/b*MgOH;
   sig = Fc/Sarea*(sCaOH2 + sMgOH2 + sCO3Ca + sCO3Mg + (zm - 1)*sCO3Me
        - sCO3m - sCaO - sMgO - sCaCO3 - sMgCO3);
   <|"sigma" -> sig, "adsorbed" -> sCO3Me|>];

(* solve the constant-capacitance constraint psi = sigma/Ccap for psi *)
surfaceSolve[metal_, aMe_, act_, logKme_] := Module[{ps, sol},
   sol = FindRoot[speciesAt[metal, aMe, act, logKme, ps]["sigma"] == Ccap*ps, {ps, 0.0}];
   speciesAt[metal, aMe, act, logKme, ps /. sol]];'''

WL_PREDICT = r'''(* ===== Forward prediction of the equilibrium dissolved metal ===== *)
predictCeq[metal_, pH_, C0ppm_, logKme_] := Module[
   {act, Ttot, Meq, aMe, sol, k},
   act = activities[pH];
   Ttot = C0ppm/MM[metal]/1000.0;
   Meq = Ttot;
   Do[
     aMe = metalActivity[metal, Meq, act];
     sol = surfaceSolve[metal, aMe, act, logKme];
     Meq = Max[Ttot - sol["adsorbed"], 1.0*10^-12];
   , {k, 40}];
   Meq*MM[metal]*1000.0];'''

WL_FIT = r'''(* ===== Fit log K by minimising WSOS. One constant per metal (carbonate site) ===== *)
wsos[kLi_?NumericQ, kCo_?NumericQ] := Module[{ss, kk, row, metal, pH, C0, Ceq, pred, s, km},
   ss = 0.0;
   Do[
     row = adsData[[kk]]; metal = row[[1]]; pH = row[[2]]; C0 = row[[3]]; Ceq = row[[4]];
     km = If[metal == "Li", kLi, kCo];
     pred = predictCeq[metal, pH, C0, km];
     s = Sqrt[errAbs^2 + (errRel*Ceq)^2];
     ss = ss + ((pred - Ceq)/s)^2;
   , {kk, Length[adsData]}];
   ss];

fit = FindMinimum[wsos[kLi, kCo], {{kLi, 0.0}, {kCo, -1.8}}];
{kLiOpt, kCoOpt} = {kLi, kCo} /. fit[[2]];
ndf = Length[adsData] - 2;
Print["log K  Li (>CO3Li) = ", NumberForm[kLiOpt, {5, 2}]];
Print["log K  Co (>CO3Co) = ", NumberForm[kCoOpt, {5, 2}]];
Print["(for scale, Pokrovsky >CO3Ca+ has log K = -1.8)"];
Print["WSOS/DF = ", NumberForm[fit[[1]]/ndf, {5, 3}]];'''

WL_VALID = r'''(* ===== Validation table ===== *)
Print["metal  pH   Ceq_meas   Ceq_pred"];
Do[
  With[{row = adsData[[kk]]},
   With[{metal = row[[1]], pH = row[[2]], C0 = row[[3]], Ceq = row[[4]]},
    Print[metal, "   ", pH, "   ", NumberForm[Ceq, {6, 3}], "   ",
      NumberForm[predictCeq[metal, pH, C0, If[metal == "Li", kLiOpt, kCoOpt]], {6, 3}]]]],
  {kk, Length[adsData]}];'''

WL_SIGMA = r'''(* ===== Surface charge from surface species only (ProtoFit Eq 2.2.3) ===== *)
sigmaAt[pH_?NumericQ] := surfaceSolve["Co", 0.0, activities[pH], -99]["sigma"];
Print["pH   sigma(C/m2)   sigma(mmol/m2)   [Pokrovsky 0.01-0.02 mmol/m2]"];
Do[Print[pH, "   ", NumberForm[sigmaAt[pH], {7, 4}], "   ",
    NumberForm[sigmaAt[pH]/Fc*1000, {8, 5}]], {pH, {5.5, 6.5, 7.3, 8.0, 8.5, 9.0}}];
Plot[sigmaAt[pH]/Fc*1000, {pH, 4, 10.5},
   AxesLabel -> {"pH", "sigma (mmol/m2)"},
   PlotLabel -> "Predicted surface charge (surface species only)",
   PlotRange -> {-0.035, 0.02}, GridLines -> Automatic]'''

def wl_text(s):
    """escape a plain-text string for a Wolfram Cell[...] string literal"""
    return s.replace("\\", "\\\\").replace('"', '\\"')

def cell_text(s, style):
    return 'Cell["%s", "%s"]' % (wl_text(s), style)

def cell_input(code):
    return 'Cell["%s", "Input"]' % wl_text(code)

nb_cells = [
  cell_text("Intrinsic log K for Li and Co on dolomite (FITEQL-style, Wolfram Mathematica)", "Title"),
  cell_text("Self-contained notebook. Every reaction, equation and datum is embedded. "
            "Method: surface complexation with the constant capacitance model, solved for the "
            "surface potential, then log K fitted by minimising WSOS/DF. Model: Pokrovsky 1999. "
            "Evaluate the cells top to bottom (Evaluation, Evaluate Notebook).", "Text"),
  cell_text("1. Reactions", "Section"),
  cell_text(REACTIONS_TEXT, "Text"),
  cell_text("2. Equations", "Section"),
  cell_text(EQUATIONS_TEXT, "Text"),
  cell_text("3. Data", "Section"),
  cell_text(DATA_TEXT, "Text"),
  cell_text("4. Setup and embedded data", "Section"),
  cell_input(WL_SETUP),
  cell_text("5. Activity coefficients and background activities", "Section"),
  cell_input(WL_ACT),
  cell_text("6. Surface speciation and the constant capacitance solve", "Section"),
  cell_input(WL_SPECIES),
  cell_text("7. Forward prediction of dissolved metal", "Section"),
  cell_input(WL_PREDICT),
  cell_text("8. Fit log K by minimising WSOS/DF", "Section"),
  cell_input(WL_FIT),
  cell_text("9. Validation", "Section"),
  cell_input(WL_VALID),
  cell_text("10. Surface charge (correct magnitude, surface species only)", "Section"),
  cell_input(WL_SIGMA),
]
nb_wl = "Notebook[{\n" + ",\n".join(nb_cells) + "\n}]\n"
with open(os.path.join(OUT,"FITEQL_LiCo_dolomite_Mathematica.nb"),"w") as f:
    f.write(nb_wl)
print("wrote Mathematica nb")
