# %% [markdown]
# # Intrinsic log K for Li⁺ and Co²⁺ on dolomite — carbonate-site model, step by step
#
# Reformulated model: the metal cation binds the **carbonate site** (Pokrovsky 1999).
#
# ```
# >CO3H0 + Li+   ⇌  >CO3Li0 + H+        K_Li      (ΔZ = 0, no Boltzmann term)
# >CO3H0 + Co2+  ⇌  >CO3Co+ + H+        K_Co      (ΔZ = +1, carries exp(FΨ/RT))
# ```
#
# Two unknowns, one per metal. Everything else is fixed from Pokrovsky (1999).
#
# **Intrinsic vs apparent.**  K_app = ({>CO3Co+}{H+}) / ({>CO3H0}{Co2+}) straight from
# concentrations;  K_int = K_app · exp(ΔZ F Ψ / RT) strips out the electrostatic work.
# Intrinsic transfers to other brines; apparent does not.
#
# Method follows Belova et al. (2014, Ni on calcite/chalk, FITEQL with a dummy adsorbed
# component) and Pokrovsky et al. (1999, 2002). Every number and datum is embedded below;
# nothing is read from a file. Run the cells top to bottom.

# %%
import numpy as np, pandas as pd
import matplotlib.pyplot as plt
from scipy.optimize import brentq, least_squares
pd.set_option("display.width", 130); pd.set_option("display.precision", 4)

def banner(t): print("\n" + "="*88 + f"\n{t}\n" + "="*88)

# publication figure style: one palette, readable fonts, 300 dpi PNG + vector PDF in ./figures
import os
os.makedirs("figures", exist_ok=True)
plt.rcParams.update({"font.size": 10, "axes.labelsize": 10.5, "axes.titlesize": 10.5,
    "legend.fontsize": 8.5, "legend.frameon": False, "xtick.labelsize": 9.5, "ytick.labelsize": 9.5,
    "axes.spines.top": False, "axes.spines.right": False, "figure.dpi": 110,
    "savefig.dpi": 300, "savefig.bbox": "tight", "axes.grid": True, "grid.alpha": 0.22,
    "lines.linewidth": 1.6, "axes.titlelocation": "left", "axes.titleweight": "bold"})
# one palette for every figure: cobalt blue, lithium rust, model green, excluded grey
COL = {"Co": "#1f5f8b", "Li": "#c8542a", "model": "#3a7d5d", "grey": "#8a8a8a", "pH2": "#7a5c99"}
MK  = {"Co": "o", "Li": "s"}
def panel(ax, letter): ax.set_title(f"({letter})", loc="left", fontsize=10, pad=6)
def savefig(fig, name):
    """300 dpi PNG and vector PDF into ./figures; captions live in the paper, not in the figure."""
    fig.savefig(f"figures/{name}.png"); fig.savefig(f"figures/{name}.pdf"); plt.show()

# %% [markdown]
# ## Step 0 — the fixed parameter block (built once, never touched during fitting)
# Every value here traces to a published table or to the experimental set-up.

# %%
banner("STEP 0  Fixed parameter block")
F, R, T = 96485.0, 8.314, 298.15
FRT = F/(R*T)                                    # 38.92 V^-1

# --- experimental set-up (spreadsheet "Set up") ---
DOLOMITE_GL = 60.0          # 6 g per 100 mL
V_L, M_G   = 0.100, 6.0     # litres of solution, grams of dolomite per bottle
SSA        = 0.76           # m2/g (BET; 0.76 to 0.84 reported)
S_AREA     = DOLOMITE_GL*SSA                     # 45.6 m2 per litre
I_BATCH    = 0.684          # mol/L, 40 g/L NaCl
NACL_M     = 40.0/58.44
MM = dict(Li=6.941, Co=58.933, Ca=40.078, Mg=24.305)

# --- site density: Pokrovsky 1999, 1:1:2 Ca:Mg:CO3 (Belova used 8.22 umol/m2 for calcite) ---
SITE_DENS = dict(CO3=14.0e-6, Ca=7.0e-6, Mg=7.0e-6)      # mol/m2
S_T = {k: v*S_AREA for k, v in SITE_DENS.items()}         # mol/L of suspension

# --- capacitance: carbonate formula C = sqrt(I)/alpha, alpha = 0.006 (Pokrovsky 2002; Belova 2014) ---
ALPHA = 0.006
C_CAP = np.sqrt(I_BATCH)/ALPHA                            # F/m2

# --- surface constants: Pokrovsky, Schott & Thomas (1999) Table 3, dolomite columns ---
# Transcribed from the page. Reaction : calcite (Ca) | magnesite (Mg) | dolomite Ca | dolomite Mg
#  1. >CO3H0 = >CO3- + H+            : -5.1  | -4.65 +/- 0.15 | -4.8 +/- 0.2 | -4.8 +/- 0.2
#  2. >CO3H0 + Me2+ = >CO3Me+ + H+   : -1.7  | -2.2  +/- 0.15 | -1.8 +/- 0.2 | -2.0 +/- 0.2
#  3. >MeOH0 = >MeO- + H+            : -12   | -12   +/- 1    | -12  +/- 2   | -12  +/- 2
#  4. >MeOH0 + H+ = >MeOH2+          : 11.5  | 10.6  +/- 0.15 | 11.5 +/- 0.2 | 10.6 +/- 0.2
#  5. >MeOH0 + CO3-2 + 2H+ = >MeHCO3 : -3.5  | -2.4  +/- 0.5  | -4.0 +/- 0.5 | -3.5 +/- 0.5
#  6. >MeOH0 + CO3-2 + H+ = >MeCO3-  : 17.1  | 14.4  +/- 0.15 | 16.6 +/- 0.2 | 15.4 +/- 0.2
# Site density, same paper, section 3.3, verbatim: "The site density was assumed to be 7 umol/m2
# for Ca and Mg and 14 umol/m2 for carbonate. This value corresponds to eight metal or carbonate
# sites per square nanometer" (1:1:2 stoichiometry). The values are the paper's, to one decimal.
LOGK_UNC = {"CO3_deprot":0.2,"CO3Ca":0.2,"CO3Mg":0.2,"CaOH2":0.2,"CaO":2.0,"CaHCO3":0.5,
            "CaCO3":0.2,"MgOH2":0.2,"MgO":2.0,"MgHCO3":0.5,"MgCO3":0.2}
LOGK_SURF = {
    "CO3_deprot": -4.8,   # >CO3H0 = >CO3- + H+
    "CO3Ca":      -1.8,   # >CO3H0 + Ca2+ = >CO3Ca+ + H+      (competitor)
    "CO3Mg":      -2.0,   # >CO3H0 + Mg2+ = >CO3Mg+ + H+      (competitor)
    "CaOH2":      11.5,   # >CaOH0 + H+ = >CaOH2+             (carries charge, no metal)
    "CaO":       -12.0,   # >CaOH0 = >CaO- + H+
    "CaHCO3":     -4.0,   # >CaOH0 + CO3-2 + 2H+ = >CaHCO3(0) + H2O
    "CaCO3":      16.6,   # >CaOH0 + CO3-2 + H+  = >CaCO3- + H2O
    "MgOH2":      10.6,
    "MgO":       -12.0,
    "MgHCO3":     -3.5,
    "MgCO3":      15.4,
}
# --- open-system carbonate ---
PCO2, LOGKH = 10**-3.5, -1.469
A_CO2 = 10**LOGKH*PCO2
LK_H2CO3, LK_HCO3 = -6.345, 10.329            # CO2+H2O=HCO3-+H+ ; CO3-2+H+=HCO3-

# --- aqueous complex FORMATION log beta from free ions (25 C, I = 0; NIST/MINTEQ values) ---
# entries: (name, log beta, {component: stoich}, charge)   components: Me, Ca, Mg, Na, Cl, CO3, H
AQ = {
  "Ca": [("CaCl+",-0.70,{"Ca":1,"Cl":1},1),("CaCO3",3.22,{"Ca":1,"CO3":1},0),
         ("CaHCO3+",11.43,{"Ca":1,"CO3":1,"H":1},1),("CaOH+",-12.70,{"Ca":1,"H":-1},1)],
  "Mg": [("MgCl+",-0.14,{"Mg":1,"Cl":1},1),("MgCO3",2.98,{"Mg":1,"CO3":1},0),
         ("MgHCO3+",11.40,{"Mg":1,"CO3":1,"H":1},1),("MgOH+",-11.79,{"Mg":1,"H":-1},1)],
  "Na": [("NaCl",-0.78,{"Na":1,"Cl":1},0),("NaCO3-",1.27,{"Na":1,"CO3":1},-1),
         ("NaHCO3",10.48,{"Na":1,"CO3":1,"H":1},0),("NaOH",-14.18,{"Na":1,"H":-1},0)],
  "Co": [("CoCl+",0.30,{"Co":1,"Cl":1},1),("CoCl2",-0.20,{"Co":1,"Cl":2},0),
         ("CoCO3",4.23,{"Co":1,"CO3":1},0),("CoHCO3+",12.20,{"Co":1,"CO3":1,"H":1},1),
         ("CoOH+",-9.65,{"Co":1,"H":-1},1),("Co(OH)2",-18.80,{"Co":1,"H":-2},0)],
  "Li": [("LiCl",-0.50,{"Li":1,"Cl":1},0),("LiOH",-13.64,{"Li":1,"H":-1},0),
         ("LiCO3-",0.90,{"Li":1,"CO3":1},-1)],
}
Z_ION = dict(Li=1, Co=2, Ca=2, Mg=2, Na=1)
# --- solubility products (log Ksp, 25 C) for the saturation check ---
LOGKSP = {"CoCO3 (sphaerocobaltite)": -9.98, "Calcite": -8.48, "Dolomite (ordered)": -17.09,
          "Li2CO3": -2.5}   # Li2CO3 approximate; it is far undersaturated regardless

table2 = pd.DataFrame([
  ("Surface area per litre St", f"{S_AREA:.1f} m2/L", "60 g/L x 0.76 m2/g (BET)"),
  ("Carbonate site density", "14 umol/m2", "Pokrovsky 1999 (1:1:2 Ca:Mg:CO3)"),
  ("Total carbonate sites S_T", f"{S_T['CO3']:.3e} mol/L", "density x St"),
  ("Capacitance C = sqrt(I)/alpha", f"{C_CAP:.0f} F/m2", "alpha = 0.006 (Pokrovsky 2002; Belova 2014)"),
  ("Ionic strength", f"{I_BATCH:.3f} M", "40 g/L NaCl"),
  ("pCO2", "10^-3.5 atm", "open system"),
  ("Temperature", "25 C", ""),
] + [(k, f"{v:+.1f} +/- {LOGK_UNC[k]}", "Pokrovsky et al. 1999, Table 3, dolomite column") for k, v in LOGK_SURF.items()],
  columns=["parameter", "value", "source"])
print(table2.to_string(index=False))

# %% [markdown]
# ## Step 0b — embedded experimental data
# Single-ion batches, dolomite 6 g/100 mL in 40 g/L NaCl, sampled at days 0, 2, 4, 6.
# The equilibrium point is day 6. Note that the batches labelled "pH 6" and "pH 2" buffered
# up as the dolomite dissolved: the **measured** equilibrium pH is what enters the model.

# %%
banner("STEP 0b  Embedded data")
raw = [  # metal, batch, day, Ca_ppm, Mg_ppm, Me_ppm, pH_measured
 ("Co","pH6",0, 0.100, 0.100, 80.513, 6.320),("Co","pH6",2, 0.100, 0.100, 77.905, 7.400),
 ("Co","pH6",4, 0.100, 8.143, 76.411, 7.477),("Co","pH6",6, 0.126,12.941, 74.393, 8.090),
 ("Co","pH2",0, 0.100, 0.100, 80.513, 2.350),("Co","pH2",2,57.477,17.461, 79.661, 4.123),
 ("Co","pH2",4,24.178,45.592, 79.795, 5.660),("Co","pH2",6,21.658,51.379, 80.398, 7.110),
 ("Li","pH6",0, 0.663, 0.100,138.430, 6.420),("Li","pH6",2,37.633, 0.100,134.657, 7.647),
 ("Li","pH6",4,54.820, 0.100,136.722, 7.923),("Li","pH6",6,49.892, 0.870,135.262, 8.257),
 ("Li","pH2",0, 0.663, 0.100,138.430, 2.250),("Li","pH2",2,39.190, 0.100,136.263, 4.257),
 ("Li","pH2",4,79.452,21.316,135.956, 5.727),("Li","pH2",6,87.738,20.547,135.592, 7.047),
]
data = pd.DataFrame(raw, columns=["metal","batch","day","Ca_ppm","Mg_ppm","Me_ppm","pH"])
data["C0_ppm"] = data.groupby(["metal","batch"])["Me_ppm"].transform("first")
data["removal_%"] = 100*(data.C0_ppm - data.Me_ppm)/data.C0_ppm
print(data.to_string(index=False))
eq = data[data.day == 6].reset_index(drop=True)      # the four equilibrium points
print("\nEquilibrium points (day 6). Note measured pH vs the nominal label:")
print(eq[["metal","batch","pH","Ca_ppm","Mg_ppm","C0_ppm","Me_ppm","removal_%"]].to_string(index=False))

# %% [markdown]
# ## The capacity check — done before any fitting
# A surface complexation constant describes a **monolayer** reaction. Removal that exceeds
# the number of surface sites cannot be adsorption, whatever constant is fitted. The check is
# Γ = ΔC / (surface area per litre) against the crystallographic site density (14 µmol/m² for
# the dolomite carbonate site, Pokrovsky 1999; 8.22 µmol/m² for calcite {10.4}, Belova 2014).
#
# It is applied to (a) this single-ion dataset and (b) the highest-recovery runs reported in
# Elshebli et al. (fine particles, 0.84 m²/g, 6 g/100 mL: Co ≈ 78 → 8 mg/L; Li 8.5 % of
# 138 mg/L). Also note: the methods section quotes the surface concentration as
# 0.046–0.049 m²/L; with 0.84 m²/g × 60 g/L it is ≈ 50 m²/L, so the published figure is
# almost certainly m²/mL. Every SCM normalisation divides by this number.

# %%
banner("CAPACITY CHECK  Can surface complexation account for the removal?")
def capacity_row(label, dC_mg_L, metal, ssa, load_gL):
    area = ssa*load_gL                                   # m2/L
    G = dC_mg_L/MM[metal]/1e3/area                       # mol/m2
    return dict(case=label, removed_mg_L=dC_mg_L, area_m2_L=area, Gamma_umol_m2=G*1e6,
                x_carbonate_sites=G/SITE_DENS["CO3"], x_calcite_8p22=G/8.22e-6)
cap = pd.DataFrame([
    capacity_row("THIS DATA  Co, coarse, pH6 batch, day 6", 80.513-74.393, "Co", SSA, DOLOMITE_GL),
    capacity_row("THIS DATA  Li, coarse, pH6 batch, day 6", 138.43-135.262, "Li", SSA, DOLOMITE_GL),
    capacity_row("MANUSCRIPT Co, fine 0.84 m2/g, pH6 25C, ~70 mg/L removed", 70.0, "Co", 0.84, DOLOMITE_GL),
    capacity_row("MANUSCRIPT Li, fine 0.84 m2/g, 8.5 % of 138 mg/L", 0.085*138.43, "Li", 0.84, DOLOMITE_GL),
])
print(cap.to_string(index=False, float_format=lambda x: f"{x:.3g}"))
print("\nx_carbonate_sites = coverage divided by the carbonate site density: > 1 is physically")
print("unreachable by adsorption, and that is not a marginal exceedance a better constant could absorb.")
print("Ceiling on sorbable metal at 50.4 m2/L if EVERY carbonate site held metal:",
      f"{SITE_DENS['CO3']*0.84*DOLOMITE_GL*1e3:.2f} mmol/L  (Co removed in the best run: {70/MM['Co']:.2f} mmol/L).")
print("\nConclusion: the high-recovery runs in the manuscript are dominated by mineralisation")
print("(XRD: zabuyelite up to 8 %, sphaerocobaltite up to 4 %). Only this low-uptake single-ion")
print("dataset is inside the monolayer bound, so it is the only one where a sorption constant is")
print("even admissible, and Step 1 still has to clear the saturation test point by point.")

# %% [markdown]
# ## Step 1 — aqueous speciation at every data point
# Inputs: measured equilibrium pH, total Na, Cl, Ca, Mg, metal, fixed pCO₂. Davies activities
# at I = 0.68 M. Outputs: free-ion activities {Li⁺}, {Co²⁺}, {Ca²⁺}, {Mg²⁺}, {H⁺}, the
# **fraction of metal that is free**, and saturation indices for CoCO₃, calcite, dolomite.
# This is where part of the total cobalt disappears into CoCl⁺, CoCO₃⁰ and CoHCO₃⁺.

# %%
banner("STEP 1  Aqueous speciation")
def davies(I):
    f = -0.509*(np.sqrt(I)/(1+np.sqrt(I)) - 0.3*I)
    return {0:1.0, 1:10**f, 2:10**(4*f)}

def speciate(pH, tot, metal):
    """tot: dict of totals mol/L for Ca, Mg, Na, Cl and the metal (key = metal name)."""
    aH = 10**-pH
    aHCO3 = 10**LK_H2CO3*A_CO2/aH
    aCO3  = aHCO3/(10**LK_HCO3*aH)
    cations = ["Ca","Mg","Na",metal]
    I = 0.5*(tot["Na"]+tot["Cl"]+sum(Z_ION[c]**2*tot[c] for c in cations if c!="Na"))
    a = {c: tot[c]*0.5 for c in cations}; a["Cl"] = tot["Cl"]*0.7
    fixed = {"CO3":aCO3,"H":aH}
    for _ in range(300):
        g = davies(I)
        def conc(name, lb, st, z, act):
            v = 10**lb
            for comp, n in st.items():
                x = act[comp] if comp in act else fixed[comp]
                v *= x**n
            return v/g[abs(z)]
        # cations
        for c in cations:
            D = 1/g[Z_ION[c]]
            for name, lb, st, z in AQ[c]:
                D += conc(name, lb, {k:v for k,v in st.items() if k!=c}, z, {**a, **fixed})
            a[c] = tot[c]/D
        # chloride
        Dcl = 1/g[1]
        for c in cations:
            for name, lb, st, z in AQ[c]:
                if "Cl" in st:
                    Dcl += st["Cl"]*conc(name, lb, {k:v for k,v in st.items() if k!="Cl"}, z, {**a, **fixed})
        aCl_new = tot["Cl"]/Dcl
        # ionic strength
        species = []
        for c in cations:
            for name, lb, st, z in AQ[c]:
                species.append((z, conc(name, lb, st, z, {**a, "Cl":aCl_new, **fixed})))
        Inew = 0.5*(sum(Z_ION[c]**2*a[c]/g[Z_ION[c]] for c in cations) + aCl_new/g[1]
                    + aH/g[1] + 1e-14/aH/g[1] + aHCO3/g[1] + 4*aCO3/g[2]
                    + sum(z*z*cc for z, cc in species))
        done = abs(aCl_new - a["Cl"]) < 1e-12 and abs(Inew - I) < 1e-9
        a["Cl"] = aCl_new; I = 0.5*I + 0.5*Inew
        if done: break
    g = davies(I)
    # metal distribution
    free = a[metal]/g[Z_ION[metal]]
    dist = {metal+"2+" if Z_ION[metal]==2 else metal+"+": free}
    for name, lb, st, z in AQ[metal]:
        dist[name] = conc(name, lb, st, z, {**a, **fixed})
    frac_free = free/tot[metal]
    SI = {"CoCO3 (sphaerocobaltite)": np.log10(a["Co"]*aCO3) - LOGKSP["CoCO3 (sphaerocobaltite)"] if metal=="Co" else np.nan,
          "Calcite": np.log10(a["Ca"]*aCO3) - LOGKSP["Calcite"],
          "Dolomite (ordered)": np.log10(a["Ca"]*a["Mg"]*aCO3**2) - LOGKSP["Dolomite (ordered)"],
          "Li2CO3": np.log10(a["Li"]**2*aCO3) - LOGKSP["Li2CO3"] if metal=="Li" else np.nan}
    return dict(I=I, aH=aH, aCO3=aCO3, aHCO3=aHCO3, aCa=a["Ca"], aMg=a["Mg"], aMe=a[metal],
                aCl=a["Cl"], g=g, frac_free=frac_free, dist=dist, SI=SI)

def totals_for(row):
    m = row.metal
    return {"Ca":row.Ca_ppm/MM["Ca"]/1e3, "Mg":row.Mg_ppm/MM["Mg"]/1e3, "Na":NACL_M,
            "Cl":NACL_M, m:row.Me_ppm/MM[m]/1e3}

spec_rows = []
for r in data.itertuples():
    if r.day == 0: continue
    s = speciate(r.pH, totals_for(r), r.metal)
    spec_rows.append(dict(metal=r.metal, batch=r.batch, day=r.day, pH=r.pH, I=s["I"],
        a_Me=s["aMe"], frac_free=s["frac_free"], a_Ca=s["aCa"], a_Mg=s["aMg"], a_CO3=s["aCO3"],
        SI_CoCO3=s["SI"]["CoCO3 (sphaerocobaltite)"], SI_calcite=s["SI"]["Calcite"],
        SI_dolomite=s["SI"]["Dolomite (ordered)"], SI_Li2CO3=s["SI"]["Li2CO3"]))
spec = pd.DataFrame(spec_rows)
print(spec.to_string(index=False, float_format=lambda x: f"{x:.3g}"))
print("\nMetal distribution at the four equilibrium points (fraction of total dissolved):")
for r in eq.itertuples():
    s = speciate(r.pH, totals_for(r), r.metal)
    tot = sum(s["dist"].values())
    print(f"  {r.metal} {r.batch} (pH {r.pH:.2f}): " +
          ", ".join(f"{k} {100*v/tot:.1f}%" for k, v in s["dist"].items()))
print("\nReading: fraction free is the number to report. Any point with SI_CoCO3 > 0 is at or")
print("above CoCO3 saturation and is not usable for a pure adsorption constant.")

# %% [markdown]
# ## Step 2 — measured uptake → surface coverage Γ (mol/m²), with propagated error
# Γ = (C_total − C_eq) · V / (m · SSA).  ICP precision taken as ±3 % on both concentrations,
# so the error on the *difference* is large when removal is small. Those points are
# down-weighted, not dropped. No blank correction is available (record this in the methods).

# %%
banner("STEP 2  Surface coverage and its uncertainty")
ICP_REL = 0.03
def coverage(C0_ppm, Ceq_ppm, metal):
    dC = (C0_ppm - Ceq_ppm)/MM[metal]/1e3                      # mol/L removed
    G  = dC*V_L/(M_G*SSA)                                         # mol/m2  (= dC/S_AREA)
    sdC = np.hypot(ICP_REL*C0_ppm, ICP_REL*Ceq_ppm)/MM[metal]/1e3
    return G, sdC*V_L/(M_G*SSA)
cov_rows = []
for r in data[data.day>0].itertuples():
    G, sG = coverage(r.C0_ppm, r.Me_ppm, r.metal)
    cov_rows.append(dict(metal=r.metal, batch=r.batch, day=r.day, pH=r.pH,
                         removal_pct=100*(r.C0_ppm-r.Me_ppm)/r.C0_ppm,
                         Gamma_mol_m2=G, sigma_Gamma=sG, rel_err_pct=100*sG/G if G>0 else np.inf,
                         site_frac_pct=100*G/SITE_DENS["CO3"]))
covdf = pd.DataFrame(cov_rows)
print(covdf.to_string(index=False, float_format=lambda x: f"{x:.3g}"))
print("\nsite_frac_pct = coverage as a percentage of the carbonate site density (14 umol/m2).")
for metal in ["Li","Co"]:
    C0 = data[data.metal==metal].C0_ppm.iloc[0]/MM[metal]/1e3
    print(f"  Site ceiling for {metal}: S_T / C0 = {S_T['CO3']*1e3:.3f} mmol/L / {C0*1e3:.2f} mmol/L "
          f"= {100*S_T['CO3']/C0:.1f} % of the added metal can be on the surface at most.")
print("  Lithium at 20 mmol/L against 0.64 mmol/L of sites is SITE-LIMITED: measured 2.3 % removal is")
print("  71 % of the ceiling. In that regime the constant is insensitive; the 8.5 % Li recovery in")
print("  Elshebli et al. (2025) is 2.4x the ceiling and cannot be sorption alone.")

# %% [markdown]
# ## Step 3 — the carbonate-site mass balance
# S_T = {>CO3H0} + {>CO3⁻} + {>CO3Ca⁺} + {>CO3Mg⁺} + {>CO3Li0} + {>CO3Co⁺}
#
# Every term is proportional to {>CO3H0}, e.g. {>CO3Ca⁺} = K_Ca · {>CO3H0} · {Ca²⁺}/{H⁺} · exp(−FΨ/RT),
# so the balance collapses to one equation in one unknown. Ca and Mg are **competitors**, not
# spectators. The >CaOH and >MgOH sites are solved too because they carry surface charge.
#
# ## Step 4 — the electrostatic loop
# σ = (F/SSA) Σ zᵢ{speciesᵢ},  Ψ = σ/C.  Adsorbed Co carries charge, so Γ changes σ, which
# changes Ψ, which changes Γ: iterate until Ψ stops moving. A non-electrostatic (NEM, Ψ = 0)
# run is made alongside; if the two constants agree, the result is thermodynamic rather than
# model dependent.

# %%
banner("STEP 3 + 4  Site balance and the electrostatic loop")
def surface(s, metal, Gamma_known=None, logK_Me=None, electrostatic=True):
    """Solve the three site balances and the CCM loop at one data point.
    Either the metal coverage is KNOWN (point-by-point mode, Gamma_known in mol/m2) or a
    trial log K is given (forward mode). Returns all surface concentrations (mol/L), psi, sigma."""
    aH, aCO3, aCa, aMg, aMe = s["aH"], s["aCO3"], s["aCa"], s["aMg"], s["aMe"]
    z = Z_ION[metal]; dZ = z - 1                      # charge of >CO3Me
    def state(psi):
        b = np.exp(-FRT*psi)                          # a species of charge q carries b**q
        rm  = 10**LOGK_SURF["CO3_deprot"]/aH/b        # >CO3-   / >CO3H0
        rCa = 10**LOGK_SURF["CO3Ca"]*aCa/aH*b         # >CO3Ca+ / >CO3H0
        rMg = 10**LOGK_SURF["CO3Mg"]*aMg/aH*b
        if Gamma_known is not None:
            cMe = Gamma_known*S_AREA
            cH0 = (S_T["CO3"] - cMe)/(1 + rm + rCa + rMg)
        else:
            rMe = 10**logK_Me*aMe/aH*b**dZ
            cH0 = S_T["CO3"]/(1 + rm + rCa + rMg + rMe)
            cMe = rMe*cH0
        cm, cCa, cMg = rm*cH0, rCa*cH0, rMg*cH0
        # hydroxyl sites (no metal)
        dCa = 1 + 10**LOGK_SURF["CaOH2"]*aH*b + 10**LOGK_SURF["CaO"]/aH/b \
                + 10**LOGK_SURF["CaHCO3"]*aCO3*aH**2 + 10**LOGK_SURF["CaCO3"]*aCO3*aH/b
        dMg = 1 + 10**LOGK_SURF["MgOH2"]*aH*b + 10**LOGK_SURF["MgO"]/aH/b \
                + 10**LOGK_SURF["MgHCO3"]*aCO3*aH**2 + 10**LOGK_SURF["MgCO3"]*aCO3*aH/b
        CaOH, MgOH = S_T["Ca"]/dCa, S_T["Mg"]/dMg
        sp = dict(CO3H0=cH0, CO3m=cm, CO3Ca=cCa, CO3Mg=cMg, CO3Me=cMe,
                  CaOH2=10**LOGK_SURF["CaOH2"]*aH*b*CaOH, CaO=10**LOGK_SURF["CaO"]/aH/b*CaOH,
                  CaCO3=10**LOGK_SURF["CaCO3"]*aCO3*aH/b*CaOH,
                  MgOH2=10**LOGK_SURF["MgOH2"]*aH*b*MgOH, MgO=10**LOGK_SURF["MgO"]/aH/b*MgOH,
                  MgCO3=10**LOGK_SURF["MgCO3"]*aCO3*aH/b*MgOH)
        sigma = F/S_AREA*(sp["CaOH2"]+sp["MgOH2"]+sp["CO3Ca"]+sp["CO3Mg"]+dZ*sp["CO3Me"]
                          - sp["CO3m"]-sp["CaO"]-sp["MgO"]-sp["CaCO3"]-sp["MgCO3"])
        return sp, sigma
    if electrostatic:
        psi = brentq(lambda p: p - state(p)[1]/C_CAP, -0.5, 0.5, xtol=1e-7)   # 0.1 mV -> 1e-4; use tighter
    else:
        psi = 0.0
    sp, sigma = state(psi)
    return dict(sp=sp, psi=psi, sigma=sigma, dZ=dZ)

# occupancy of the carbonate site at the equilibrium points (metal coverage taken from data)
occ_rows = []
for r in eq.itertuples():
    s = speciate(r.pH, totals_for(r), r.metal)
    G, _ = coverage(r.C0_ppm, r.Me_ppm, r.metal)
    st = surface(s, r.metal, Gamma_known=G)
    sp = st["sp"]; tot = S_T["CO3"]
    occ_rows.append(dict(metal=r.metal, batch=r.batch, pH=r.pH,
        CO3H0=100*sp["CO3H0"]/tot, CO3minus=100*sp["CO3m"]/tot, CO3Ca=100*sp["CO3Ca"]/tot,
        CO3Mg=100*sp["CO3Mg"]/tot, CO3Me=100*sp["CO3Me"]/tot,
        psi_mV=1e3*st["psi"], sigma_C_m2=st["sigma"], boltz=np.exp(-FRT*st["psi"])))
occ = pd.DataFrame(occ_rows)
print("Carbonate-site occupancy (% of S_T) and the electrostatic state:")
print(occ.to_string(index=False, float_format=lambda x: f"{x:.3g}"))
print(f"\nC = {C_CAP:.0f} F/m2 keeps psi at a few mV, so the Boltzmann factor is close to 1:")
print("the electrostatic and non-electrostatic constants will nearly agree (Belova found 0.01 log units).")

# %% [markdown]
# ## Step 5a — point-by-point K (whiteboard algebra, no optimizer)
# With Γ measured, the mass-action expression is solved directly at each point:
#
# K_Me = ( {>CO3Me} · {H⁺} ) / ( {>CO3H0} · {Me} ) · exp(ΔZ F Ψ / RT)
#
# One K per data point. Day-6 points are the equilibrium values; day-2 and day-4 are shown
# as well because they widen the pH range for the diagnostic in Step 6, but they are
# pre-equilibrium and are marked as such.

# %%
banner("STEP 5a  Point-by-point log K")
def logK_point(row, electrostatic=True):
    s = speciate(row.pH, totals_for(row), row.metal)
    G, sG = coverage(row.C0_ppm, row.Me_ppm, row.metal)
    if G <= 0: return np.nan, np.nan, s, G, sG
    def lk_at(Gx):
        st = surface(s, row.metal, Gamma_known=Gx, electrostatic=electrostatic)
        sp = st["sp"]
        return np.log10(sp["CO3Me"]*s["aH"]/(sp["CO3H0"]*s["aMe"]) * np.exp(st["dZ"]*FRT*st["psi"])), st["psi"]
    lk, psi = lk_at(G)
    # propagate the coverage error through the mass-action inversion (asymmetric in log K)
    Gmax = 0.999*S_T["CO3"]/S_AREA
    lo = lk_at(max(G - sG, 1e-3*G))[0]; hi = lk_at(min(G + sG, Gmax))[0]
    return lk, psi, s, G, sG, (lo, hi)
pp = []
for r in data[data.day>0].itertuples():
    out = logK_point(r, True)
    if len(out) == 5:      # G <= 0
        lk_ccm, psi, s, G, sG = out; lo = hi = np.nan
    else:
        lk_ccm, psi, s, G, sG, (lo, hi) = out
    out2 = logK_point(r, False); lk_nem = out2[0]
    SI_mc = s["SI"]["CoCO3 (sphaerocobaltite)"] if r.metal=="Co" else s["SI"]["Li2CO3"]
    rel = 100*sG/G if G>0 else np.inf
    reason = ("SI(MeCO3) > 0: at or above saturation" if SI_mc > 0 else
              "removal within analytical noise (rel. error > 300 %)" if rel > 300 else "")
    pp.append(dict(metal=r.metal, batch=r.batch, day=r.day, equilibrium=(r.day==6), pH=r.pH,
                   Gamma=G, rel_err_pct=rel, SI_MeCO3=SI_mc,
                   logK_CCM=lk_ccm, logK_lo=lo, logK_hi=hi, logK_NEM=lk_nem,
                   psi_mV=1e3*psi if psi==psi else np.nan, usable=(reason==""), reason=reason))
ppdf = pd.DataFrame(pp)
print(ppdf.drop(columns=["reason"]).to_string(index=False, float_format=lambda x: f"{x:.3g}"))
print("\nlogK_lo / logK_hi: the constant recomputed at Gamma -/+ one propagated standard deviation.")
print("Exclusions, with the reason stated:")
for r in ppdf[~ppdf.usable].itertuples():
    print(f"  {r.metal} {r.batch} day {r.day}: {r.reason}")
print("\nFor scale: Pokrovsky's >CO3H0 + Ca2+ = >CO3Ca+ + H+ is log K = -1.8, Mg -2.0.")

# %% [markdown]
# ## Step 5b — global weighted least squares (the reported value and its uncertainty)
# The forward model predicts Γ from a trial log K; residuals are weighted by the propagated
# Γ error from Step 2. Fitted on the day-6 equilibrium points, CCM and NEM side by side.
# This is what Belova did in FITEQL with the `Niad` dummy component.

# %%
banner("STEP 5b  Global least-squares fit (equilibrium points)")
def fit_metal(metal, electrostatic, rows):
    pts = []
    for r in rows.itertuples():
        s = speciate(r.pH, totals_for(r), r.metal)
        G, sG = coverage(r.C0_ppm, r.Me_ppm, r.metal)
        pts.append((s, G, sG))
    def resid(theta):
        out = []
        for s, G, sG in pts:
            st = surface(s, metal, logK_Me=theta[0], electrostatic=electrostatic)
            out.append((st["sp"]["CO3Me"]/S_AREA - G)/sG)
        return np.array(out)
    x0 = np.nanmean(ppdf[(ppdf.metal==metal)&ppdf.equilibrium]["logK_CCM"])
    sol = least_squares(resid, [x0 if np.isfinite(x0) else -1.0], bounds=([-10],[10]), diff_step=1e-4)
    r = sol.fun; dof = max(len(r)-1, 1)
    s2 = np.sum(r**2)/dof
    se = float(np.sqrt(s2/np.sum(sol.jac**2))) if np.sum(sol.jac**2)>0 else np.nan
    return sol.x[0], se, float(np.sum(r**2)/dof), len(r)
usable_keys = set(zip(ppdf[ppdf.usable].metal, ppdf[ppdf.usable].batch, ppdf[ppdf.usable].day))
fit_rows = []
for metal in ["Li","Co"]:
    rows = data[(data.metal==metal)&(data.day>0)]
    rows = rows[[ (m,b,d) in usable_keys for m,b,d in zip(rows.metal, rows.batch, rows.day) ]]
    for ele, tag in [(True,"CCM"),(False,"NEM")]:
        lk, se, wsos, n = fit_metal(metal, ele, rows)
        col = "logK_CCM" if ele else "logK_NEM"
        u = ppdf[(ppdf.metal==metal)&ppdf.usable]
        spread = u[col].std(ddof=1) if len(u)>1 else np.nan
        # propagated per-point width, averaged: the honest floor on the uncertainty
        prop = np.nanmean(0.5*(u.logK_hi - u.logK_lo)) if ele else np.nan
        fit_rows.append(dict(metal=metal, model=tag, logK=lk, se_fit=se, spread_pbp=spread,
                             propagated=prop, WSOS_DF=wsos, n_points=n))
fits = pd.DataFrame(fit_rows)
print(fits.to_string(index=False, float_format=lambda x: f"{x:.3f}"))
print("\nn_points = usable points (all sampling days, weighted by the propagated coverage error).")
print("se_fit is the least-squares standard error; spread_pbp the standard deviation of the usable")
print("point-by-point values; propagated the mean half-width from the coverage error. The reported")
print("uncertainty is the LARGEST of the three. Lithium is site-limited (Step 2), so a small se_fit")
print("there is insensitivity, not precision: the propagated width is the number to quote.")

# %% [markdown]
# ## Step 6 — the diagnostic that tells you the model is right
# Point-by-point log K against pH and against surface loading Γ. A flat line means the pH
# dependence of uptake is fully explained by speciation, site protonation and the
# electrostatic term — a transferable constant. Drift with pH means a wrong surface species
# or a missing coulombic term; drift with loading means surface precipitation.

# %%
banner("STEP 6  Diagnostic: is log K flat?")
fig, axes = plt.subplots(1, 2, figsize=(9.2, 4.0), constrained_layout=True)
for ax, xcol, xlab, letter in [(axes[0],"pH","measured pH at sampling","a"),(axes[1],"Gamma","surface coverage Γ (mol m$^{-2}$)","b")]:
    for metal in ["Co","Li"]:
        d = ppdf[(ppdf.metal==metal)&np.isfinite(ppdf.logK_CCM)]
        u = d[d.usable]; x = d[~d.usable]
        ax.errorbar(u[xcol], u.logK_CCM, yerr=[u.logK_CCM-u.logK_lo, u.logK_hi-u.logK_CCM],
                    fmt=MK[metal], ms=6.5, color=COL[metal], ecolor=COL[metal], elinewidth=0.8, capsize=2.5,
                    alpha=0.95, ls="none", label=f"{metal}, usable point (bar: propagated interval)")
        ax.plot(x[xcol], x.logK_CCM, "x", color=COL["grey"], ms=6.5, mew=1.3, ls="none",
                label=f"{metal}, excluded (noise or SI > 0)")
        fitv = float(fits[(fits.metal==metal)&(fits.model=="CCM")].logK.iloc[0])
        ax.axhline(fitv, color=COL[metal], ls="--", lw=1.0, alpha=0.7, label=f"{metal}, weighted fit {fitv:+.2f}")
    ax.axhline(-1.8, color=COL["grey"], ls=":", lw=1.1, label="Pokrovsky >CO$_3$Ca$^+$ (−1.8)")
    ax.set_xlabel(xlab); ax.set_ylabel("point-by-point log K (CCM)"); ax.set_ylim(-7, 1.2); panel(ax, letter)
    if xcol == "Gamma": ax.set_xscale("log")
h, l = axes[0].get_legend_handles_labels()
fig.legend(h, l, loc="outside lower center", ncol=3, fontsize=8)
savefig(fig, "Fig6_logK_diagnostic")
for metal in ["Li","Co"]:
    d = ppdf[(ppdf.metal==metal)&ppdf.equilibrium&np.isfinite(ppdf.logK_CCM)]
    if len(d)>1:
        slope = np.polyfit(d.pH, d.logK_CCM, 1)[0]
        print(f"  {metal}: slope of log K vs pH over the equilibrium points = {slope:+.2f} per pH unit "
              f"({'flat within noise' if abs(slope)<0.5 else 'a trend - check the surface species / precipitation'})")

# %% [markdown]
# ## Figure 1 — kinetics (percent removal vs time, inset Γ)
# Justifies the equilibration time and shows whether a slow second stage (precipitation or
# recrystallisation) is present.

# %%
fig, (ax, bx) = plt.subplots(1, 2, figsize=(9.2, 3.8), constrained_layout=True)
for (metal, batch), d in data.groupby(["metal","batch"]):
    filled = batch == "pH6"
    kw = dict(marker=MK[metal], color=COL[metal], ms=6, mfc=COL[metal] if filled else "white", mew=1.4,
              ls="-" if filled else "--", label=f"{metal}, initial {batch.replace('pH','pH ')}")
    ax.plot(d.day, d["removal_%"], **kw)
    G = [coverage(c0, c, metal)[0]*1e6 for c0, c in zip(d.C0_ppm, d.Me_ppm)]
    bx.plot(d.day, G, **kw)
bx.axhline(SITE_DENS["CO3"]*1e6, color=COL["grey"], ls=":", lw=1.2, label="carbonate site density")
ax.set_xlabel("time (days)"); ax.set_ylabel("removal (% of initial)")
bx.set_xlabel("time (days)"); bx.set_ylabel("surface coverage Γ (µmol m$^{-2}$)")
ax.set_xticks([0,2,4,6]); bx.set_xticks([0,2,4,6]); bx.set_ylim(0, 15.5)
ax.legend(loc="upper left"); bx.legend(loc="upper left", frameon=True, framealpha=0.92, edgecolor="none")
panel(ax, "a"); panel(bx, "b")
savefig(fig, "Fig1_kinetics")

# %% [markdown]
# ## Figure 2 — adsorption edges (percent removal vs pH) with the model
# The model curve is the forward prediction with the fitted log K, using each batch's
# initial metal and its measured day-6 Ca and Mg. The position of the edge carries the
# proton stoichiometry.

# %%
def forward_removal(metal, pH, C0_ppm, Ca_ppm, Mg_ppm, logK, electrostatic=True):
    """Forward model with the metal mass balance closed: dissolved = total - adsorbed."""
    C0 = C0_ppm/MM[metal]/1e3; D = C0
    for _ in range(60):
        tot = {"Ca":Ca_ppm/MM["Ca"]/1e3, "Mg":Mg_ppm/MM["Mg"]/1e3, "Na":NACL_M, "Cl":NACL_M, metal:max(D,1e-12)}
        s = speciate(pH, tot, metal)
        st = surface(s, metal, logK_Me=logK, electrostatic=electrostatic)
        Dn = max(C0 - st["sp"]["CO3Me"], 1e-12)
        if abs(Dn-D) < 1e-12*C0: D = Dn; break
        D = 0.5*D + 0.5*Dn
    return 100*(C0-D)/C0
best = {m: float(fits[(fits.metal==m)&(fits.model=="CCM")].logK.iloc[0]) for m in ["Li","Co"]}
fig, axes = plt.subplots(1, 2, figsize=(9.2, 3.9), constrained_layout=True)
grid = np.linspace(3, 10, 141)
for ax, metal, letter in zip(axes, ["Co","Li"], "ab"):
    for batch, ls, lab in [("pH6","-","model, pH 6 batch chemistry"),("pH2","--","model, pH 2 batch chemistry")]:
        e = eq[(eq.metal==metal)&(eq.batch==batch)].iloc[0]
        curve = [forward_removal(metal, p, e.C0_ppm, e.Ca_ppm, e.Mg_ppm, best[metal]) for p in grid]
        ax.plot(grid, curve, ls, color=COL["model"], label=lab)
    d = data[(data.metal==metal)&(data.day>0)]
    ax.plot(d[d.day==6].pH, d[d.day==6]["removal_%"], MK[metal], color=COL[metal], ms=7.5, ls="none",
            label="measured, day 6", zorder=3)
    ax.plot(d[d.day<6].pH, d[d.day<6]["removal_%"], MK[metal], color=COL[metal], mfc="white", mew=1.3,
            ms=6, ls="none", label="measured, days 2 and 4 (pre-equilibrium)")
    if metal == "Li":
        C0 = data[data.metal=="Li"].C0_ppm.iloc[0]/MM["Li"]/1e3
        ax.axhline(100*S_T["CO3"]/C0, color=COL["grey"], ls=":", lw=1.2, label="site ceiling (all carbonate sites)")
    ax.set_xlabel("pH"); ax.set_ylabel(f"{metal} removal (% of initial)")
    ax.legend(loc="upper left" if metal=="Co" else "lower right",
              title=f"apparent log K = {best[metal]:+.2f}", title_fontsize=8.5); panel(ax, letter)
savefig(fig, "Fig2_adsorption_edges")

# %% [markdown]
# ## Figures 3 and 4 — isotherms (log Γ vs log C_eq) with the model
# The model line is the forward prediction at the batch's equilibrium pH and chemistry while
# the initial metal is varied. A slope of one is low-coverage adsorption; a steepening at high
# concentration is precipitation, not adsorption. At 100 mg/L this is where it would show.

# %%
fig, axes = plt.subplots(1, 2, figsize=(9.2, 3.9), constrained_layout=True)
C0_grid = np.logspace(-6, -1, 60)
for ax, metal, letter in zip(axes, ["Co","Li"], "ab"):
    for batch, ls, filled in [("pH6","-",True),("pH2","--",False)]:
        e = eq[(eq.metal==metal)&(eq.batch==batch)].iloc[0]
        Ceq, Gam = [], []
        for C0 in C0_grid:
            rem = forward_removal(metal, e.pH, C0*MM[metal]*1e3, e.Ca_ppm, e.Mg_ppm, best[metal])
            Ceq.append(C0*(1-rem/100)); Gam.append(C0*rem/100/S_AREA)
        ax.plot(np.log10(Ceq), np.log10(np.maximum(Gam,1e-30)), ls, color=COL["model"], label=f"model, pH {e.pH:.2f}")
        G, _ = coverage(e.C0_ppm, e.Me_ppm, metal)
        if G > 0:
            ax.plot(np.log10(e.Me_ppm/MM[metal]/1e3), np.log10(G), MK[metal], color=COL[metal], ms=7.5, ls="none",
                    mfc=COL[metal] if filled else "white", mew=1.4,
                    label=f"measured day 6, {batch.replace('pH','pH ')} batch", zorder=3)
    ax.axhline(np.log10(SITE_DENS["CO3"]), color=COL["grey"], ls=":", lw=1.2, label="carbonate site density")
    ax.set_xlabel("log C$_{eq}$ (mol L$^{-1}$)"); ax.set_ylabel("log Γ (mol m$^{-2}$)")
    ax.set_ylim(-9.4, -4.6); ax.legend(loc="lower right"); panel(ax, letter)
savefig(fig, "Fig3_4_isotherms")
print("The unit-slope region is adsorption; the flattening is site saturation (Langmuir behaviour).")

# %% [markdown]
# ## Figure 5 — distribution coefficient, log K_d vs pH
# K_d = (C_total − C_eq)/C_eq · V/m, in L/kg. The applied figure: it converts directly into a
# column capacity or bed-volume estimate.

# %%
fig, ax = plt.subplots(figsize=(6.4, 4.0), constrained_layout=True)
grid = np.linspace(3, 10, 141)
for metal in ["Co","Li"]:
    d = data[(data.metal==metal)&(data.day>0)].copy()
    d["Kd"] = (d.C0_ppm - d.Me_ppm)/d.Me_ppm*V_L/M_G*1e3
    noise = (d.metal=="Co")&(d.batch=="pH2")            # within analytical noise (Step 5a)
    e = d[(d.day==6)&~noise]; k = d[(d.day<6)&~noise&(d.Kd>0)]; x = d[noise&(d.Kd>0)]
    ax.plot(e.pH, np.log10(e.Kd), MK[metal], color=COL[metal], ms=7.5, ls="none", label=f"{metal}, day 6", zorder=3)
    ax.plot(k.pH, np.log10(k.Kd), MK[metal], color=COL[metal], mfc="white", mew=1.3, ms=6, ls="none",
            label=f"{metal}, days 2 and 4")
    if len(x): ax.plot(x.pH, np.log10(x.Kd), "x", color=COL["grey"], ms=7, mew=1.4, ls="none",
                       label="Co, pH 2 batch (within noise, excluded)")
    e0 = eq[(eq.metal==metal)&(eq.batch=="pH6")].iloc[0]
    rem = np.array([forward_removal(metal, p, e0.C0_ppm, e0.Ca_ppm, e0.Mg_ppm, best[metal]) for p in grid])
    Kd = rem/(100-rem)*V_L/M_G*1e3
    ax.plot(grid, np.log10(np.maximum(Kd,1e-3)), "-", color=COL[metal], lw=1.4, alpha=0.9, label=f"{metal}, model")
ax.set_xlabel("pH"); ax.set_ylabel("log $K_\\mathrm{d}$ (L kg$^{-1}$)"); ax.set_ylim(-2.1, 0.5)
ax.legend(loc="lower right", ncol=2)
savefig(fig, "Fig5_Kd_vs_pH")

# %% [markdown]
# ## Surface speciation plot — why uptake rises with pH
# Fractional occupancy of the carbonate site vs pH at the Co batch chemistry. It makes the
# Ca and Mg competition visible instead of asserting it.

# %%
fig, ax = plt.subplots(figsize=(6.8, 4.0), constrained_layout=True)
e0 = eq[(eq.metal=="Co")&(eq.batch=="pH6")].iloc[0]
grid = np.linspace(3, 10, 141)
fr = {k: [] for k in ["CO3H0","CO3m","CO3Ca","CO3Mg","CO3Me"]}
for p in grid:
    tot = totals_for(e0); s_ = speciate(p, tot, "Co")
    st = surface(s_, "Co", logK_Me=best["Co"])
    for k in fr: fr[k].append(100*st["sp"][k]/S_T["CO3"])
labels = {"CO3H0":">CO$_3$H$^0$","CO3m":">CO$_3^-$","CO3Ca":">CO$_3$Ca$^+$","CO3Mg":">CO$_3$Mg$^+$","CO3Me":">CO$_3$Co$^+$"}
pal = {"CO3H0":"#4c4c4c","CO3m":"#3a7d5d","CO3Ca":"#c8a44a","CO3Mg":"#7a5c99","CO3Me":COL["Co"]}
for k, v in fr.items(): ax.plot(grid, v, color=pal[k], lw=2.4 if k=="CO3Me" else 1.6, label=labels[k])
ax.set_xlabel("pH"); ax.set_ylabel("fraction of carbonate sites (%)"); ax.set_ylim(0, 100)
ax.legend(loc="center right", title=f"Co pH 6 batch: Ca {e0.Ca_ppm:.2f}, Mg {e0.Mg_ppm:.1f} mg L$^{{-1}}$", title_fontsize=8)
savefig(fig, "FigS1_surface_speciation")

# %% [markdown]
# ## Step 7 — report honestly
# Table 1: the constants, electrostatic and non-electrostatic, with uncertainty. Then the
# comparison against Pokrovsky's aqueous/surface analogy: on carbonates the surface constant
# tracks the aqueous MeCO₃⁰ constant (Van Cappellen et al. 1993; Pokrovsky et al. 1999), so
# from >CO₃Ca⁺ (−1.8) and log β(CaCO₃⁰) = 3.22 the expected cobalt value is
# −1.8 + [log β(CoCO₃⁰) − log β(CaCO₃⁰)].

# %%
banner("STEP 7  Table 1 and the literature comparison")
t1 = fits.pivot(index="metal", columns="model", values=["logK","se_fit","spread_pbp","WSOS_DF"])
print("Table 1  Intrinsic stability constants for the carbonate-site complexes, 25 C, I = 0.68 M NaCl")
for metal in ["Co","Li"]:
    c = fits[(fits.metal==metal)&(fits.model=="CCM")].iloc[0]; n = fits[(fits.metal==metal)&(fits.model=="NEM")].iloc[0]
    rxn = ">CO3H0 + Co2+ = >CO3Co+ + H+" if metal=="Co" else ">CO3H0 + Li+ = >CO3Li0 + H+"
    unc_c = np.nanmax([c.se_fit, c.spread_pbp, c.propagated]); unc_n = np.nanmax([n.se_fit, n.spread_pbp, c.propagated])
    print(f"  {rxn:32}  CCM log K = {c.logK:+.2f} +/- {unc_c:.2f}   NEM log K = {n.logK:+.2f} +/- {unc_n:.2f}")
print("\nUncertainty = largest of the fit standard error, the point-by-point spread and the propagated width.")
nco = int(fits[(fits.metal=="Co")&(fits.model=="CCM")].n_points.iloc[0])
print(f"Cobalt rests on {nco} usable point(s); that is an apparent value, not an intrinsic constant.")
lb_CaCO3, lb_CoCO3 = 3.22, 4.23
expected_Co = LOGK_SURF["CO3Ca"] + (lb_CoCO3 - lb_CaCO3)
got_Co = float(fits[(fits.metal=="Co")&(fits.model=="CCM")].logK.iloc[0])
print(f"\nPokrovsky / Van Cappellen aqueous-surface analogy for Co:")
print(f"  expected log K(>CO3Co+) = log K(>CO3Ca+) + [log b(CoCO3) - log b(CaCO3)] = {LOGK_SURF['CO3Ca']} + ({lb_CoCO3} - {lb_CaCO3}) = {expected_Co:+.2f}")
print(f"  fitted  log K(>CO3Co+) = {got_Co:+.2f}   ->  difference {got_Co-expected_Co:+.2f} log units")
print("  If the fitted value sits below the analogy line, the usual reason on a brine is Ca/Mg")
print("  competition for the sites (Step 3 occupancy) or a part of the removal being precipitation")
print("  rather than adsorption (Step 1 SI, Step 6 loading trend). Both are results, not caveats.")
print("\nJudgment calls to record in the methods: (1) total BET area vs a reduced reactive area")
print("(Belova's 70 % case raises log K); (2) >CO3Ca/>CO3Mg fractions follow the measured Ca and Mg")
print("point by point (done here) rather than being fixed; (3) no blank correction was available;")
print("(4) alkalinity is set by fixed pCO2 rather than measured DIC.")

# %% [markdown]
# ## Step 8 — turn the problem around: predict, then split the mechanism
# Because cobalt in the high-recovery runs is beyond the monolayer bound and this dataset's
# Co "pH 6" point is above sphaerocobaltite saturation, a *fitted* Co constant would be a lumped
# parameter absorbing the CoCO₃ solubility product, carbonate supply and nucleation. Instead:
#
# * **Cobalt**: take the constant from the aqueous–surface analogy (Van Cappellen et al. 1993;
#   Pokrovsky & Schott 2002): log K(>CO₃Co⁺) = log K(>CO₃Ca⁺) + [log β(CoCO₃⁰) − log β(CaCO₃⁰)]
#   = −1.8 + (4.23 − 3.22) = **−0.79**. Run the model forward with the fixed constants at each
#   time point, compute the sorbed amount, and read the difference from the measured removal
#   as the mineralised fraction, checked against the saturation index.
# * **Lithium**: the one metal whose sorption signal is not swamped, so it is the single
#   adjustable parameter, constrained by all the Li points (uptake plateaus by day 2), and
#   reported as conditional on the borrowed site parameters.

# %%
banner("STEP 8  Predict sorption from literature constants; split sorption from mineralisation")
LOGK_CO_PRED = LOGK_SURF["CO3Ca"] + (4.23 - 3.22)
print(f"Cobalt constant from the aqueous-surface analogy: log K(>CO3Co+) = {LOGK_CO_PRED:+.2f}")

def sorbed_forward(metal, row, logK, electrostatic=True):
    """predicted sorbed metal (mmol/L) and % of initial, forward model with mass balance closed"""
    rem = forward_removal(metal, row.pH, row.C0_ppm, row.Ca_ppm, row.Mg_ppm, logK, electrostatic)
    return rem/100*row.C0_ppm/MM[metal], rem
ceiling_mM = S_T["CO3"]*1e3
rows8 = []
for r in data[(data.metal=="Co")&(data.day>0)].itertuples():
    s = speciate(r.pH, totals_for(r), "Co")
    sorb_mM, sorb_pct = sorbed_forward("Co", r, LOGK_CO_PRED)
    meas_mM = (r.C0_ppm - r.Me_ppm)/MM["Co"]
    rows8.append(dict(batch=r.batch, day=r.day, pH=r.pH, SI_CoCO3=s["SI"]["CoCO3 (sphaerocobaltite)"],
                      measured_removed_mM=meas_mM, measured_pct=100*meas_mM*MM["Co"]/r.C0_ppm,
                      predicted_sorbed_mM=sorb_mM, predicted_sorbed_pct=sorb_pct,
                      residual_meas_minus_pred_mM=meas_mM - sorb_mM,
                      pred_over_meas=sorb_mM/meas_mM if meas_mM>0 else np.nan))
split = pd.DataFrame(rows8)
print("\nCobalt, this dataset (log K predicted, not fitted). Residual is SIGNED: negative means the")
print("model puts more cobalt on the surface than left solution.")
print(split.to_string(index=False, float_format=lambda x: f"{x:.3g}"))
ratio = split[split.batch=="pH6"].pred_over_meas
print(f"\nDIAGNOSIS. The analogy constant over-predicts sorption by {ratio.min():.0f}x to {ratio.max():.0f}x in the")
print("pH 6 batch (and by far more where measured removal is within noise). Something in the chain is")
print("wrong by more than an order of magnitude. The candidates, with the size of correction each needs:")
u = ppdf[(ppdf.metal=="Co")&ppdf.usable]
gap = LOGK_CO_PRED - u.logK_CCM.mean() if len(u) else np.nan
print(f"  (a) the borrowed constant: the usable points give an apparent log K of {u.logK_CCM.mean():+.2f}, "
      f"i.e. {gap:.1f} log units weaker than the analogy;")
f_area = (split[split.batch=="pH6"].measured_removed_mM/split[split.batch=="pH6"].predicted_sorbed_mM)
print(f"  (b) reactive area smaller than BET: with the analogy constant, a reactive fraction of "
      f"{f_area.min():.2f} to {f_area.max():.2f} of the BET area reproduces the pH 6 batch "
      f"(Belova's 70 % case is the same move, in the same direction);")
print("  (c) these batches are not at equilibrium: removal is still rising at day 6 (Figure 1),")
print("      so the day-6 coverage is a lower bound on the equilibrium coverage.")
print("None of these can be separated with two pH values and one concentration; the isotherm and")
print("pH-edge experiment described in the Discussion is what separates them.")
print(f"\nAbsolute ceiling if EVERY carbonate site held Co at 45.6 m2/L: {ceiling_mM:.3f} mmol/L "
      f"({100*ceiling_mM/(80.513/MM['Co']):.0f} % of the 80.5 mg/L added).")
print(f"At the manuscript's 0.84 m2/g the ceiling is {SITE_DENS['CO3']*0.84*DOLOMITE_GL*1e3:.2f} mmol/L "
      f"against {70/MM['Co']:.2f} mmol/L removed: at least {100*(1-SITE_DENS['CO3']*0.84*DOLOMITE_GL*1e3/(70/MM['Co'])):.0f} % of that")
print("removal is mineralisation even at the ceiling, and far more once Ca, Mg and H+ occupy sites.")

# lithium: the single adjustable parameter, from the usable points (Step 5b), with the honest width
fl = fits[(fits.metal=="Li")&(fits.model=="CCM")].iloc[0]
print(f"\nLithium (fitted, conditional on the borrowed site parameters):")
print(f"  usable points (n = {int(fl.n_points)}): log K(>CO3Li0) = {fl.logK:+.2f};  internal spread {fl.spread_pbp:.2f};")
print(f"  propagated interval at +/-3 % ICP precision: about {fl.logK-fl.propagated:+.1f} to {fl.logK+fl.propagated:+.1f}.")
# how much a better ICP precision would help: recompute the propagated half-width at 1 %
_saved = ICP_REL; ICP_REL = 0.01; w1 = []
for r in data[(data.metal=="Li")&(data.day>0)].itertuples():
    o = logK_point(r, True)
    if len(o) == 6: w1.append(0.5*(o[5][1]-o[5][0]))
ICP_REL = _saved
print(f"  at +/-1 % ICP precision the propagated half-width would fall to about {np.nanmean(w1):.1f} log units.")
print("  Lithium is site-limited (Step 2): the value is an order of magnitude, not a constant, until")
print("  an isotherm at lower Li concentration puts the surface below saturation.")

fig, (ax, bx) = plt.subplots(1, 2, figsize=(9.2, 3.8), constrained_layout=True)
bcol = {"pH6": COL["Co"], "pH2": COL["pH2"]}
for batch in ["pH6","pH2"]:
    d = split[split.batch==batch]; lab = batch.replace("pH","pH ")
    ax.plot(d.day, d.measured_pct, "o-", color=bcol[batch], ms=6.5, label=f"measured removal, {lab} batch")
    ax.plot(d.day, d.predicted_sorbed_pct, "o--", color=bcol[batch], mfc="white", mew=1.3, ms=6,
            label=f"predicted sorption (log K = {LOGK_CO_PRED:+.2f}), {lab} batch")
    bx.plot(d.day, d.SI_CoCO3, "o-", color=bcol[batch], ms=6.5, label=f"{lab} batch")
bx.axhline(0, color=COL["grey"], lw=1); bx.axhspan(0, 2, color=COL["grey"], alpha=0.12)
bx.text(2.05, 0.15, "supersaturated", fontsize=8.5, color=COL["grey"])
ax.set_xlabel("time (days)"); ax.set_ylabel("Co removed or sorbed (% of initial)")
bx.set_xlabel("time (days)"); bx.set_ylabel("saturation index, sphaerocobaltite CoCO$_3$")
ax.set_xticks([2,4,6]); bx.set_xticks([2,4,6]); ax.set_ylim(-1, 30); bx.set_ylim(-7.5, 1.5)
ax.legend(loc="center left", fontsize=8); bx.legend(loc="lower right"); panel(ax, "a"); panel(bx, "b")
savefig(fig, "Fig8_Co_sorption_vs_removal_SI")
print("Read the two curves together with the SI: predicted sorption sits ABOVE measured removal at")
print("every point, so in this coarse, low-uptake data set precipitation is not needed to explain the")
print("removal; the model is over-predicting sorption. The SI crossing zero at day 6 marks where")
print("precipitation would begin to add to it.")


# %% [markdown]
# ## Data check — the Ca/Mg release is opposite between the two batches
# The same dolomite in the same brine should release Ca and Mg in roughly the same proportion.
# In the Co pH 6 batch at day 6, Ca = 0.13 and Mg = 12.9 mg/L; in the Li pH 6 batch, Ca = 49.9 and
# Mg = 0.87 mg/L. These are opposite by two orders of magnitude. Because Ca and Mg are the
# competitors that set how many sites the metal can reach, this must be checked against the raw
# ICP run before it propagates. The sensitivity below shows how much it matters.

# %%
banner("DATA CHECK  Ca/Mg release pattern and its effect on the site balance")
chk = data[data.day>0][["metal","batch","day","Ca_ppm","Mg_ppm"]].copy()
chk["Ca_over_Mg_molar"] = (chk.Ca_ppm/MM["Ca"])/(chk.Mg_ppm/MM["Mg"])
print(chk.to_string(index=False, float_format=lambda x: f"{x:.3g}"))
print("\nDolomite dissolving congruently gives Ca/Mg near 1. Values of ~0.005 and ~90 in the two")
print("pH 6 batches cannot both be right for the same solid. Sensitivity: swap Ca and Mg in the")
print("Co pH 6 day-6 point and recompute the point-by-point constant.")
r0 = eq[(eq.metal=="Co")&(eq.batch=="pH6")].iloc[0]
sw = r0.copy(); sw["Ca_ppm"], sw["Mg_ppm"] = r0.Mg_ppm, r0.Ca_ppm
for lab, rr in [("as recorded", r0), ("Ca and Mg swapped", sw)]:
    out = logK_point(rr, True)
    print(f"  Co pH6 day 6, {lab:18}: log K = {out[0]:+.2f}   (Ca {rr.Ca_ppm:.2f}, Mg {rr.Mg_ppm:.2f} mg/L)")
print("The constant moves by the difference shown; if the raw run confirms the values, the difference")
print("is real incongruent dissolution and should be reported as such.")

# %% [markdown]
# ## Independent check — one point set up for Visual MINTEQ
# The speciation above is hand-coded with Davies at I = 0.68 M, near the edge of the equation's
# range. The direct answer to "is the code right" is to run one point through Visual MINTEQ
# (Gustafsson, 2014) and compare. The block below is the input for the Co pH 6 day-6 point and
# the numbers this notebook produces for it. An activity-model sensitivity (Davies vs the
# extended Debye–Hückel / B-dot form used by EQ3/6) is given alongside so the size of that
# uncertainty is known before the comparison.

# %%
banner("VISUAL MINTEQ CROSS-CHECK  input for one point and the values to compare")
r0 = eq[(eq.metal=="Co")&(eq.batch=="pH6")].iloc[0]; tot = totals_for(r0); s0 = speciate(r0.pH, tot, "Co")
print("Visual MINTEQ set-up (Main menu):")
print(f"  pH: fixed at {r0.pH:.2f}      Temperature 25 C      Activity model: Davies (then repeat with")
print("  Debye-Huckel)              Ionic strength: calculated")
print("  Gases: CO2(g) fixed, log partial pressure = -3.5")
print("  Components and totals (mol/L):")
for c in ["Na","Cl","Ca","Mg","Co"]:
    print(f"     {c:3}  {tot[c]:.4e}")
print("  Adsorption: none (this check is for the aqueous speciation only)")
print("\nValues this notebook gives for that point (compare with the MINTEQ output page):")
tot_d = sum(s0["dist"].values())
print(f"  ionic strength           {s0['I']:.4f} M")
print(f"  free Co2+ activity       {s0['aMe']:.4e}     fraction of total Co free {100*s0['frac_free']:.1f} %")
for k, v in s0["dist"].items(): print(f"  {k:9} {100*v/tot_d:6.2f} % of dissolved Co")
print(f"  CO3-2 activity           {s0['aCO3']:.4e}")
print(f"  SI sphaerocobaltite      {s0['SI']['CoCO3 (sphaerocobaltite)']:+.2f}   SI calcite {s0['SI']['Calcite']:+.2f}   SI dolomite {s0['SI']['Dolomite (ordered)']:+.2f}")
# activity-model sensitivity: B-dot (extended Debye-Huckel with ion-size and b-dot term)
def bdot_gammas(I):
    A, B, bdot = 0.5092, 0.3283, 0.041
    def g(z, a): return 10**(-A*z*z*np.sqrt(I)/(1+B*a*np.sqrt(I)) + bdot*I)
    return {0:1.0, 1:g(1, 4.0), 2:g(2, 6.0)}
gD, gB = davies(s0["I"]), bdot_gammas(s0["I"])
print(f"\nActivity-coefficient sensitivity at I = {s0['I']:.2f} M:  gamma(2+) Davies {gD[2]:.3f} vs B-dot {gB[2]:.3f};")
print(f"  gamma(1+) Davies {gD[1]:.3f} vs B-dot {gB[1]:.3f}. Free Co2+ activity would change by a factor of")
print(f"  {gB[2]/gD[2]:.2f} between the two models, i.e. {np.log10(gB[2]/gD[2]):+.2f} in log K. That is the floor")
print("  on how well any constant can be known at this ionic strength, whichever code computes it.")
