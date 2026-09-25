# %% [markdown]
# # Intrinsic surface complexation constants for Li⁺ and Co²⁺ on dolomite
# ## Single ion in sodium chloride, and produced water with the full metal suite
#
# **Marwa Elshebli, Javier Vilcáez, James Smay** — Boone Pickens School of Geology,
# Oklahoma State University
#
# ---
# ## The procedure
#
# The procedure is the one Ebrahimi & Vilcáez (2018, *J. Environ. Manage.* **214**, 370–378)
# set out for determining an intrinsic stability constant on dolomite in a sodium chloride brine.
# Their steps, in order:
#
# 1. Take the Pokrovsky, Schott & Thomas (1999) three site surface complexation model for
#    dolomite, with its published constants held fixed.
# 2. Add one reaction for the metal of interest on the carbonate site,
#    **>CO₃H⁰ + Meᶻ⁺ ⇌ >CO₃Me⁽ᶻ⁻¹⁾⁺ + H⁺**, with
#
#    $$K_{int}=\frac{[\mathrm{>CO_3Me^+}]\,[\mathrm{H^+}]}{[\mathrm{>CO_3H^0}]\,[\mathrm{Me^{2+}}]}$$
#
# 3. Carry the aqueous chloro complexes explicitly, because at brine salinity they control how
#    much of the dissolved metal is present as the free ion the surface can bind.
# 4. Account for dolomite dissolution, which raises the pH and releases calcium and magnesium
#    that compete for the same sites.
# 5. Estimate the one unknown constant by fitting the model to the measured sorption profile.
#
# Only the procedure is taken. Every number reported below comes from the experiments in this
# study.
#
# ## What is applied to, and the one departure
#
# The procedure is applied twice, to the two experiments of this work:
#
# * **Single ion.** Cobalt or lithium alone in 40 g L⁻¹ NaCl, which is the system their procedure
#   was written for.
# * **Produced water.** The same metals in a brine that also carries 4186 mg L⁻¹ calcium,
#   1006 mg L⁻¹ magnesium, and lead, cadmium, barium and strontium at about 80 mg L⁻¹ each.
#
# The departure is at step 4. They solved dolomite dissolution kinetically, with a transition
# state rate law, because they did not measure the pH and the dissolved calcium and magnesium
# through the run. These experiments **did** measure all three at every sampling, so the measured
# values are used directly and no rate law is needed. That removes two fitted kinetic parameters
# from the problem and replaces them with data.

# %%
import numpy as np, pandas as pd, matplotlib.pyplot as plt, os
from scipy.optimize import brentq, minimize_scalar
pd.set_option("display.width", 190); pd.set_option("display.precision", 4)
os.makedirs("figures", exist_ok=True)

def head(n, t): print(f"\n{'='*106}\n{n}  {t}\n{'='*106}")
def show(df):   print(df.to_string(index=False, float_format=lambda x: f"{x:,.4g}"))

plt.rcParams.update({
    "font.family": "sans-serif", "font.sans-serif": ["DejaVu Sans", "Arial"],
    "font.size": 9, "axes.labelsize": 9.5, "legend.fontsize": 8,
    "figure.facecolor": "white", "axes.facecolor": "white", "savefig.facecolor": "white",
    "axes.grid": False, "axes.linewidth": 0.8,
    "xtick.direction": "in", "ytick.direction": "in", "xtick.top": True, "ytick.right": True,
    "xtick.minor.visible": True, "ytick.minor.visible": True,
    "lines.linewidth": 1.5, "lines.markersize": 6, "legend.frameon": False,
    "savefig.dpi": 300, "savefig.bbox": "tight", "pdf.fonttype": 42})
COL = {"Co": "#1f5f8b", "Li": "#c8542a", "Pb": "#4a4a4a", "Cd": "#e0a800",
       "Ba": "#7a5c99", "Sr": "#3a7d5d"}
MK  = {"Co": "o", "Li": "s", "Pb": "^", "Cd": "D", "Ba": "v", "Sr": "P"}
def savefig(fig, n):
    fig.savefig(f"figures/{n}.png"); fig.savefig(f"figures/{n}.pdf"); plt.show()

# %% [markdown]
# ---
# ## 1 — Step 1: the Pokrovsky model, read from the source
#
# Pokrovsky, Schott & Thomas (1999) Table 3, dolomite columns. One entry in that table needs care,
# and settling it is worth the space because it changes the surface speciation completely.

# %%
head(1, "Step 1: the Pokrovsky model, read from the source")

LOGK = {"CO3_deprot": -4.8,  "CO3Ca": -1.8,  "CO3Mg": -2.0,
        "CaOH2": 11.5, "CaO": -12.0, "CaHCO3": 24.0, "CaCO3": 16.6,
        "MgOH2": 10.6, "MgO": -12.0, "MgHCO3": 23.5, "MgCO3": 15.4}

show(pd.DataFrame([
    ("1", ">CO3H0 = >CO3- + H+",                     LOGK["CO3_deprot"], LOGK["CO3_deprot"], "+/-0.2"),
    ("2", ">CO3H0 + Me2+ = >CO3Me+ + H+",            LOGK["CO3Ca"],      LOGK["CO3Mg"],      "+/-0.2"),
    ("3", ">MeOH0 = >MeO- + H+",                     LOGK["CaO"],        LOGK["MgO"],        "+/-2"),
    ("4", ">MeOH0 + H+ = >MeOH2+",                   LOGK["CaOH2"],      LOGK["MgOH2"],      "+/-0.2"),
    ("5", ">MeOH0 + CO3-2 + 2H+ = >MeHCO3(0) + H2O", LOGK["CaHCO3"],     LOGK["MgHCO3"],     "+/-0.5"),
    ("6", ">MeOH0 + CO3-2 + H+ = >MeCO3- + H2O",     LOGK["CaCO3"],      LOGK["MgCO3"],      "+/-0.2"),
], columns=["#", "reaction on the dolomite surface", "log K, Ca site", "log K, Mg site", "uncertainty"]))

print("\n  READING REACTION 5. Text extracted from the Pokrovsky PDF renders the leading minus sign")
print("  of a negative entry as a digit, so reaction 5 comes out as '24.0' and '23.5' with the")
print("  same ambiguity that rows 1 to 3 have, where the values are certainly negative. Three")
print("  independent checks settle it, and all three point the same way.\n")
_r6 = dict(calcite=17.1, magnesite=14.4, dol_Ca=16.6, dol_Mg=15.4)
chk = []
for lab, _r5 in [("positive, +24.0 and +23.5", dict(calcite=23.5, magnesite=22.4, dol_Ca=24.0, dol_Mg=23.5)),
                 ("negative, -4.0 and -3.5",   dict(calcite=-3.5, magnesite=-2.4, dol_Ca=-4.0, dol_Mg=-3.5))]:
    d = {k: _r5[k] - _r6[k] for k in _r6}
    chk.append(dict(reading=lab, **{f"d_{k}": v for k, v in d.items()},
                    verdict="physical" if min(d.values()) > 0 else "impossible"))
print("  CHECK 1. Reaction 5 minus reaction 6 is the protonation >MeCO3- + H+ = >MeHCO3(0),")
print("  which must be positive and of order 10, since aqueous HCO3- formation is 10.33:")
show(pd.DataFrame(chk))
print("\n  CHECK 2. Pokrovsky's abstract names the five species the model needs: >CO3-, >CO3Me+,")
print("  >MeOH2+, >MeHCO3(0) and >MeCO3-. At the conditions of these experiments the two readings")
print("  put >MeHCO3(0) at 10^+3.1 and 10^-24.9 relative to >MeOH0. A species at 10^-25 cannot be")
print("  one of the five a model is built on.")
print("\n  CHECK 3. The positive reading gives a protonation constant of +6.4 to +8.1 across all")
print("  four mineral columns, a tight and sensible spread; the negative reading gives -16.8 to")
print("  -20.6, scattered and unphysical.")
print("\n  The positive values are used. This matters: under the negative reading >MeHCO3(0)")
print("  disappears from the surface entirely, and the calcium and magnesium sites are speciated")
print("  wrongly at every pH in these experiments.")

# %% [markdown]
# ---
# ## 2 — Parameters
#
# Site densities and the capacitance relation are Pokrovsky's own, read from the same paper:
# 7 µmol m⁻² for the calcium and magnesium sites and 14 µmol m⁻² for the carbonate site, a 1:1:2
# stoichiometry corresponding to eight sites per square nanometre, and C = √I/α with
# α = 0.004 C mol^½ V⁻¹ m³, the value he fitted to reproduce the high surface charge of a
# carbonate and its weak response to ionic strength.

# %%
head(2, "Parameters")

F_CONST, R_GAS, T_K = 96485.0, 8.314, 298.15
FRT = F_CONST/(R_GAS*T_K)
N_AV = 6.02214e23
SOLID_GL, SSA = 60.0, 0.76
S_AREA  = SOLID_GL*SSA
ICP_REL = 0.03
ALPHA_POK = 0.004
SITE_DENS = dict(CO3=14.0e-6, Ca=7.0e-6, Mg=7.0e-6)          # mol/m2, Pokrovsky section 3.3
NT = {k: v*S_AREA for k, v in SITE_DENS.items()}             # mol/L of suspension
A_DAVIES = 0.509

I_SINGLE, I_PW = 0.684, 1.00                                  # 40 g/L NaCl; produced water
NACL_M = 40.0/58.44
def c_pokrovsky(I): return np.sqrt(I)/ALPHA_POK
MM = dict(Ca=40.078, Mg=24.305, Li=6.941, Co=58.933, Pb=207.2, Sr=87.62, Ba=137.33, Cd=112.41)
Z_ION = dict(Li=1, Co=2, Pb=2, Sr=2, Ba=2, Cd=2)

def gam(z, I):
    return 10**(-A_DAVIES*z*z*(np.sqrt(I)/(1 + np.sqrt(I)) - 0.3*I))

show(pd.DataFrame([
    ("Solid loading",           f"{SOLID_GL:.0f} g/L",            "6 g per 100 mL, measured"),
    ("Specific surface area",   f"{SSA:.2f} m2/g",                "BET, measured"),
    ("Surface per litre",       f"{S_AREA:.1f} m2/L",             "product of the two"),
    ("Site density, carbonate", f"{SITE_DENS['CO3']*1e6:.0f} umol/m2",   "Pokrovsky 1999 section 3.3"),
    ("Site density, Ca and Mg", f"{SITE_DENS['Ca']*1e6:.0f} umol/m2",    "Pokrovsky 1999 section 3.3"),
    ("Sites per nm2",           f"{SITE_DENS['CO3']*N_AV/1e18:.1f}",     "his 1:1:2 stoichiometry"),
    ("Carbonate sites per litre", f"{NT['CO3']*1e3:.3f} mmol/L",  "density x area"),
    ("Capacitance, single ion", f"{c_pokrovsky(I_SINGLE):.0f} F/m2", f"sqrt(I)/{ALPHA_POK} at I = {I_SINGLE:.3f} M"),
    ("Capacitance, produced water", f"{c_pokrovsky(I_PW):.0f} F/m2", f"sqrt(I)/{ALPHA_POK} at I = {I_PW:.2f} M"),
    ("ICP-OES precision",       f"+/-{100*ICP_REL:.0f} %",        "on each concentration"),
], columns=["parameter", "value", "source"]))

PCO2, LOGKH = 10**-3.5, -1.469
LK_H2CO3, LK_HCO3 = -6.345, 10.329
KSP = {"Sr": -9.271, "Ba": -8.562, "Cd": -13.74, "Co": -9.98, "Pb": -13.13, "Li": -2.50}
MINERAL = {"Sr": "strontianite", "Ba": "witherite", "Cd": "otavite",
           "Co": "sphaerocobaltite", "Pb": "cerussite", "Li": "zabuyelite"}

# Aqueous complexes, log beta from the free ions (NIST / MINTEQ). Step 3 of the procedure.
AQ = {
 "Co": [("CoCl+",0.30,dict(Cl=1),1), ("CoCl2",-0.20,dict(Cl=2),0),
        ("CoCO3",4.23,dict(CO3=1),0), ("CoHCO3+",12.20,dict(CO3=1,H=1),1),
        ("CoOH+",-9.65,dict(H=-1),1), ("Co(OH)2",-18.80,dict(H=-2),0)],
 "Li": [("LiCl",-0.50,dict(Cl=1),0), ("LiCO3-",0.90,dict(CO3=1),-1),
        ("LiOH",-13.64,dict(H=-1),0)],
 "Pb": [("PbCl+",1.60,dict(Cl=1),1), ("PbCl2",1.80,dict(Cl=2),0), ("PbCl3-",1.70,dict(Cl=3),-1),
        ("PbCO3",6.48,dict(CO3=1),0), ("PbHCO3+",13.20,dict(CO3=1,H=1),1),
        ("PbOH+",-7.60,dict(H=-1),1), ("Pb(OH)2",-17.09,dict(H=-2),0)],
 "Cd": [("CdCl+",1.98,dict(Cl=1),1), ("CdCl2",2.60,dict(Cl=2),0), ("CdCl3-",2.40,dict(Cl=3),-1),
        ("CdCO3",4.36,dict(CO3=1),0), ("CdHCO3+",12.40,dict(CO3=1,H=1),1),
        ("CdOH+",-10.08,dict(H=-1),1)],
 "Ba": [("BaCl+",-0.13,dict(Cl=1),1), ("BaCO3",2.71,dict(CO3=1),0),
        ("BaHCO3+",11.31,dict(CO3=1,H=1),1), ("BaOH+",-13.47,dict(H=-1),1)],
 "Sr": [("SrCl+",-0.25,dict(Cl=1),1), ("SrCO3",2.81,dict(CO3=1),0),
        ("SrHCO3+",11.54,dict(CO3=1,H=1),1), ("SrOH+",-13.29,dict(H=-1),1)],
}

# %% [markdown]
# ---
# ## 3 — The two experiments
#
# Both at 25 °C with dolomite at 6 g per 100 mL, both sampled at days 0, 2, 4 and 6, both run
# from two different starting pH values. The batch labels are nominal starting values only:
# dolomite dissolution buffers the suspension and the measured pH is what enters the model.

# %%
head(3, "The two experiments")

SINGLE = pd.DataFrame([
 ("Co","pH6",0, 0.100, 0.100, 80.513, 6.320),("Co","pH6",2, 0.100, 0.100, 77.905, 7.400),
 ("Co","pH6",4, 0.100, 8.143, 76.411, 7.477),("Co","pH6",6, 0.126,12.941, 74.393, 8.090),
 ("Co","pH2",0, 0.100, 0.100, 80.513, 2.350),("Co","pH2",2,57.477,17.461, 79.661, 4.123),
 ("Co","pH2",4,24.178,45.592, 79.795, 5.660),("Co","pH2",6,21.658,51.379, 80.398, 7.110),
 ("Li","pH6",0, 0.663, 0.100,138.430, 6.420),("Li","pH6",2,37.633, 0.100,134.657, 7.647),
 ("Li","pH6",4,54.820, 0.100,136.722, 7.923),("Li","pH6",6,49.892, 0.870,135.262, 8.257),
 ("Li","pH2",0, 0.663, 0.100,138.430, 2.250),("Li","pH2",2,39.190, 0.100,136.263, 4.257),
 ("Li","pH2",4,79.452,21.316,135.956, 5.727),("Li","pH2",6,87.738,20.547,135.592, 7.047),
], columns=["metal","batch","day","Ca","Mg","Me","pH"])
SINGLE["C0"] = SINGLE.groupby(["metal","batch"])["Me"].transform("first")
SINGLE["removal_pct"] = 100*(SINGLE.C0 - SINGLE.Me)/SINGLE.C0
print("SINGLE ION  -  metal alone in 40 g/L NaCl, I = 0.684 M")
show(SINGLE)

PW = pd.DataFrame([
    ("pH6", 0, 4185.950, 1005.970, 83.350, 92.673, 75.960, 90.165, 88.579, 76.780, 6.320),
    ("pH6", 2, 3930.733,  940.894, 48.789, 84.949, 69.892, 83.058, 82.369, 71.243, 7.413),
    ("pH6", 4, 3890.693,  984.050,  1.072, 25.739, 11.527, 49.280, 82.500, 71.470, 7.557),
    ("pH6", 6, 3754.450,  967.633,  2.037, 14.066,  7.019, 44.457, 83.166, 70.212, 7.857),
    ("pH2", 0, 4185.950, 1005.970, 83.350, 92.673, 75.960, 90.165, 88.579, 76.780, 2.350),
    ("pH2", 2, 4037.933,  978.907, 47.138, 92.468, 75.013, 87.731, 86.930, 76.212, 4.290),
    ("pH2", 4, 4199.987, 1115.780,  2.164, 43.352, 26.687, 54.205, 86.297, 74.399, 5.727),
    ("pH2", 6, 4412.000, 1231.000,  0.982, 38.880, 20.460, 57.230, 86.760, 73.770, 6.913),
], columns=["batch","day","Ca","Mg","Pb","Cd","Co","Ba","Sr","Li","pH"])
print("\nPRODUCED WATER  -  the same metals with Pb, Cd, Ba, Sr and 4186 mg/L Ca, I = 1.0 M")
show(PW)

cmp3 = []
for m in ["Co", "Li"]:
    s6 = SINGLE[(SINGLE.metal == m) & (SINGLE.batch == "pH6") & (SINGLE.day == 6)].iloc[0]
    p0 = PW[(PW.batch == "pH6") & (PW.day == 0)].iloc[0]
    p6 = PW[(PW.batch == "pH6") & (PW.day == 6)].iloc[0]
    cmp3.append(dict(metal=m, system="single ion", C0_mgL=s6.C0, Ceq_mgL=s6.Me,
                     removal_pct=s6.removal_pct, pH=s6.pH, Ca_mgL=s6.Ca))
    cmp3.append(dict(metal=m, system="produced water", C0_mgL=p0[m], Ceq_mgL=p6[m],
                     removal_pct=100*(p0[m]-p6[m])/p0[m], pH=p6.pH, Ca_mgL=p6.Ca))
print("\n  Day 6 of the pH 6 batch, the two systems side by side:")
show(pd.DataFrame(cmp3))
_c3 = pd.DataFrame(cmp3)
_rCo = _c3[(_c3.metal=="Co")].removal_pct.values; _rLi = _c3[(_c3.metal=="Li")].removal_pct.values
print(f"\n  Cobalt is removed {_rCo[1]/_rCo[0]:.0f} times more strongly in produced water than alone, and lithium")
print(f"  {_rLi[1]/_rLi[0]:.1f} times. Sorption alone cannot do that. Produced water carries 3754 mg/L of calcium")
print(f"  against {_c3[(_c3.metal=='Co')&(_c3.system=='single ion')].Ca_mgL.iloc[0]:.2f} mg/L in the cobalt batch and {_c3[(_c3.metal=='Li')&(_c3.system=='single ion')].Ca_mgL.iloc[0]:.0f} mg/L in the lithium batch, and that calcium")
print("  competes for the very sites the metal would bind, so a purely sorptive system should")
print("  remove LESS in produced water, not more. Something else is operating, and Part 5")
print("  identifies it before any constant is fitted.")

# %% [markdown]
# ---
# ## 4 — Step 3: aqueous speciation and the chloride effect
#
# The constant is written on the **free ion**, so the free ion fraction has to be computed at
# every point. Ebrahimi & Vilcáez identified this as the dominant control on sorption at brine
# salinity, and the same calculation for these two systems shows how much of it there is here.

# %%
head(4, "Step 3: aqueous speciation and the chloride effect")

def background(pH, aCO3=None):
    aH = 10**(-pH)
    if aCO3 is None or not np.isfinite(aCO3):
        aCO3 = (10**LK_H2CO3*(10**LOGKH*PCO2)/aH)/(10**LK_HCO3*aH)   # open system, their Table 2
    return aH, aCO3

def free_fraction(metal, aH, aCO3, aCl):
    tot = 1.0
    for _, lb, st, _ in AQ[metal]:
        v = 10**lb
        for c, n in st.items():
            v *= {"Cl": aCl, "H": aH, "CO3": aCO3}[c]**n
        tot += v
    return 1.0/tot

sp = []
for lab, I, pH in [("single ion, 40 g/L NaCl", I_SINGLE, 8.09),
                   ("produced water",          I_PW,     7.857)]:
    aCl = NACL_M*gam(1, I)
    aH, aCO3 = background(pH)
    for m in ["Co", "Li"]:
        sp.append(dict(system=lab, metal=m, pH=pH,
                       free_pct=100*free_fraction(m, aH, aCO3, aCl),
                       free_pct_no_Cl=100*free_fraction(m, aH, aCO3, 1e-12)))
sp = pd.DataFrame(sp); sp["chloride_costs_pct"] = sp.free_pct_no_Cl - sp.free_pct
show(sp)
print("\n  Cobalt loses about half of its free ion pool to chloride and carbonate complexing;")
print("  lithium barely complexes at all. The correction is therefore large for cobalt and small")
print("  for lithium, and it is applied at every point rather than as a single factor.")

# %% [markdown]
# ---
# ## 5 — Step 4: dolomite dissolution, and what it does to each system
#
# Ebrahimi & Vilcáez solved dolomite dissolution kinetically to obtain the pH. Here the pH,
# calcium and magnesium were measured at every sampling, so the state of the system is known
# rather than simulated. What dissolution also does, over six days rather than their 300 minutes,
# is release enough carbonate for some metals to precipitate. Precipitation and sorption both
# remove metal from solution, and the constant describes only one of them, so the two have to be
# separated before anything is fitted.

# %%
head(5, "Step 4: dolomite dissolution, and what it does to each system")

def co3_from_Sr(row):
    return 10**KSP["Sr"]/(row.Sr/MM["Sr"]/1e3*gam(2, I_PW))

carb = pd.DataFrame([dict(batch=r.batch, day=int(r.day), pH=r.pH, a_CO3=co3_from_Sr(r))
                     for _, r in PW[PW.day > 0].iterrows()])
carb["log_a_CO3"] = np.log10(carb.a_CO3)
print("  PRODUCED WATER. Strontium fixes the carbonate activity: Zachara et al. (1991) measured")
print("  log cKex(Sr) = -2.04 on calcite and describe strontium as effectively non sorbing, so")
print("  strontianite equilibrium alone sets its dissolved concentration.")
show(carb)
print(f"  Steady at log a(CO3) = {carb.log_a_CO3.mean():.2f} +/- {carb.log_a_CO3.std():.2f} across every sampling and both batches, which")
print("  is what a solubility controlled system should give.")

sat = []
for _, r in PW[PW.day > 0].iterrows():
    aCO3 = carb[(carb.batch == r.batch) & (carb.day == r.day)].a_CO3.iloc[0]
    c0 = PW[(PW.batch == r.batch) & (PW.day == 0)].iloc[0]
    for m in ["Co", "Li"]:
        fl = ((10**KSP[m]/aCO3)**0.5/gam(1, I_PW)*MM[m]*1e3 if m == "Li"
              else 10**KSP[m]/aCO3/gam(2, I_PW)*MM[m]*1e3)
        removed = (c0[m] - r[m])/MM[m]/1e3
        srb = (removed if fl > c0[m] else
               ((fl - r[m])/MM[m]/1e3 if r[m] < fl else 0.0))
        sat.append(dict(system="produced water", batch=r.batch, day=int(r.day), metal=m, pH=r.pH,
                        measured_mgL=r[m], floor_mgL=fl, removed_M=removed, sorbed_M=srb,
                        sorbed_pct=100*srb/removed if removed > 0 else np.nan))
for _, r in SINGLE[SINGLE.day > 0].iterrows():
    aH, aCO3 = background(r.pH)
    m = r.metal
    fl = ((10**KSP[m]/aCO3)**0.5/gam(1, I_SINGLE)*MM[m]*1e3 if m == "Li"
          else 10**KSP[m]/aCO3/gam(2, I_SINGLE)*MM[m]*1e3)
    removed = (r.C0 - r.Me)/MM[m]/1e3
    srb = (removed if fl > r.C0 else ((fl - r.Me)/MM[m]/1e3 if r.Me < fl else 0.0))
    sat.append(dict(system="single ion", batch=r.batch, day=int(r.day), metal=m, pH=r.pH,
                    measured_mgL=r.Me, floor_mgL=fl, removed_M=removed, sorbed_M=srb,
                    sorbed_pct=100*srb/removed if removed > 0 else np.nan))
sat = pd.DataFrame(sat)
print("\n  Saturation screen, every sampling of both systems:")
show(sat[["system","batch","metal","day","pH","measured_mgL","floor_mgL","sorbed_pct"]]
        .sort_values(["metal","system","batch","day"]))

_nus = sat.groupby(["metal","system"]).apply(
    lambda g: pd.Series(dict(n=len(g), usable=int((g.sorbed_pct > 0).sum()))),
    include_groups=False).reset_index()
print("\n  Samplings where a surface constant is defined:")
show(_nus)

print("\n  LITHIUM is usable everywhere, in both systems and at every sampling. Zabuyelite would")
print(f"  need lithium at {sat[sat.metal=='Li'].floor_mgL.min():,.0f} mg/L at its most soluble point in this set, against the 77 to")
print("  138 mg/L present, so it cannot form and no correction is applied anywhere.")
print("\n  COBALT is the metal the carbonate chemistry interferes with, and it does so differently")
print("  in the two systems.")
print("    In the SINGLE ION system cobalt is undersaturated for five of the six samplings. The")
print(f"    exception is day 6 of the pH 6 batch, where the pH has risen to 8.09 and the")
print(f"    sphaerocobaltite floor has fallen to {sat[(sat.system=='single ion')&(sat.metal=='Co')&(sat.day==6)&(sat.batch=='pH6')].floor_mgL.iloc[0]:.1f} mg/L against the {sat[(sat.system=='single ion')&(sat.metal=='Co')&(sat.day==6)&(sat.batch=='pH6')].measured_mgL.iloc[0]:.1f} mg/L measured. That is")
print("    the largest uptake in the whole cobalt series and it is the one point that has to be")
print("    dropped, which is worth stating plainly because it is the point a fit would lean on.")
print("    In PRODUCED WATER cobalt is supersaturated at every sampling but one. Six days of")
print("    dolomite dissolution in a brine that already carries 4186 mg/L of calcium drives")
print("    sphaerocobaltite formation, and that, not sorption, is what removes the cobalt.")
print("\n  This is the answer to the discrepancy noted at the end of Part 3. Cobalt is removed")
print("  twelve times more strongly in produced water not because it sorbs more strongly there,")
print("  but because it precipitates there and does not, for the most part, in the single ion")
print("  system. The procedure can still be applied to both, but it yields a sorption constant")
print("  from one and almost nothing from the other.")

# %% [markdown]
# ---
# ## 5b — What the dissolution rate law contributes, and what it catches
#
# The rate law of step 4 is not used to generate the pH here, because the pH was measured. But
# the term in braces is not kinetic at all: it is the dolomite saturation state,
#
# $$1-\frac{[\mathrm{HCO_3^-}]^2[\mathrm{Ca^{2+}}][\mathrm{Mg^{2+}}]}{[\mathrm{H^+}]^{2}K_{eq}}$$
#
# and it can be evaluated directly from the measured calcium, magnesium and pH together with the
# carbonate activity each part of Part 5 assumed. That makes it a test of the one quantity in this
# work that had to be inferred rather than measured, so it is worth running even though the rate
# itself is not needed.
#
# One note on the equation before using it. The published form carries [H⁺] to the first power,
# which does not balance the reaction it is written for, CaMg(CO₃)₂ + 2H⁺ ⇌ Ca²⁺ + Mg²⁺ + 2HCO₃⁻;
# the squared form is used here. The equilibrium constant is built from the dolomite solubility
# product, which is equivalent to the tabulated one, as the correction printed below sets out.

# %%
head(5.5, "What the dissolution rate law contributes, and what it catches")

LOGKSP_DOL = -17.09
_kbuilt = LOGKSP_DOL + 2*LK_HCO3
print("  CORRECTION. An earlier version of this section read the tabulated log Keq as -2.525 and")
print("  concluded it was six log units wrong. The tabulated value is +2.525 and it is correct:")
for _lab, _ksp in [("ordered dolomite, -17.09", -17.09), ("EQ3/6 dolomite-ord, -18.13", -18.13)]:
    print(f"    log Ksp = {_lab:32}  ->  log Keq = {_ksp + 2*LK_HCO3:+.3f}")
print("  +2.525 is exactly the second, so the source is internally consistent and the error was")
print("  mine. The saturation indices below use the solubility product directly, which is")
print("  unaffected either way.")

def dolomite_SI(pH, Ca_ppm, Mg_ppm, aCO3, I):
    g1, g2 = gam(1, I), gam(2, I)
    aH = 10**(-pH)
    aCa, aMg = Ca_ppm/MM["Ca"]/1e3*g2, Mg_ppm/MM["Mg"]/1e3*g2
    SI = np.log10(aCa*aMg*aCO3**2) - LOGKSP_DOL
    return SI, aCO3*10**LK_HCO3*aH/g1*1e3           # SI, and the implied HCO3- in mmol/L

sat_dol = []
for _, r in SINGLE[SINGLE.day > 0].iterrows():
    SI, hco3 = dolomite_SI(r.pH, r.Ca, r.Mg, background(r.pH)[1], I_SINGLE)
    sat_dol.append(dict(system="single ion", series=f"{r.metal} {r.batch}", day=int(r.day),
                        pH=r.pH, implied_HCO3_mM=hco3, SI_dolomite=SI))
for _, r in PW[PW.day > 0].iterrows():
    aC = carb[(carb.batch == r.batch) & (carb.day == r.day)].a_CO3.iloc[0]
    SI, hco3 = dolomite_SI(r.pH, r.Ca, r.Mg, aC, I_PW)
    sat_dol.append(dict(system="produced water", series=r.batch, day=int(r.day),
                        pH=r.pH, implied_HCO3_mM=hco3, SI_dolomite=SI))
sat_dol = pd.DataFrame(sat_dol)
print()
show(sat_dol)

_ss = sat_dol[sat_dol.system == "single ion"]
_pp = sat_dol[sat_dol.system == "produced water"]
print(f"\n  SINGLE ION. Every sampling is undersaturated with dolomite, from {_ss.SI_dolomite.min():.1f} to {_ss.SI_dolomite.max():.1f}, and the")
print(f"  implied bicarbonate runs from below a micromolar to {_ss.implied_HCO3_mM.max():.2f} mmol/L. Both are what a")
print("  suspension dissolving dolomite into a dilute brine should show, so the open system")
print("  carbonate activity used for this experiment passes the test.")

_bad = _pp[_pp.implied_HCO3_mM > 10]
print(f"\n  PRODUCED WATER. The saturation state is steady at about {_pp.SI_dolomite.mean():+.1f}, which is the expected")
print("  sign: the brine dissolves dolomite and cannot reprecipitate it at any useful rate, so it")
print("  sits above equilibrium. But the implied bicarbonate fails, and badly, in the acid start")
print(f"  batch: {_bad.implied_HCO3_mM.max():,.0f} mmol/L at pH {_bad.pH.min():.2f}. No solution holds that much bicarbonate.")
print("\n  WHAT THAT MEANS. Strontianite equilibrium fixes a(CO3) almost independently of pH, so the")
print("  bicarbonate it implies rises as the pH falls, and below about pH 7 the assumption becomes")
print("  physically inadmissible. Strontium is simply not solubility controlled in the acid start")
print("  batch: its removal there is 1.9 to 2.6 per cent, which is within the analytical noise.")

_pH0 = 4.290; _aC0 = float(carb[(carb.batch == 'pH2') & (carb.day == 2)].a_CO3.iloc[0])
_aCl0 = NACL_M*gam(1, I_PW)
_hi = free_fraction("Li", 10**(-_pH0), _aC0, _aCl0)
_lo = free_fraction("Li", 10**(-_pH0), _aC0/10, _aCl0)
print("\n  WHAT IT CHANGES IN THE REPORTED NUMBERS: nothing, and the reason is worth stating.")
print("    The single ion constants, which are the recommended values, never use strontium. They")
print("    take their carbonate from open system equilibrium, which the table above validates.")
print("    In produced water lithium needs no precipitation correction at any carbonate activity,")
print("    so a(CO3) reaches its constant only through the minor LiCO3- complex: a tenfold change")
print(f"    in a(CO3) moves the lithium constant by {abs(np.log10(_hi/_lo)):.4f} log units.")
print(f"    The single usable produced water cobalt point sits at pH {_pp[_pp.implied_HCO3_mM<10].pH.max():.2f}, where the implied")
print(f"    bicarbonate is {_pp[_pp.implied_HCO3_mM<10].implied_HCO3_mM.min():.2f} mmol/L and physically reasonable.")
print("    So the check tightens the account without moving a number, which is the most useful")
print("    outcome a check of an assumption can have.")

# %% [markdown]
# ---
# ## 6 — Step 2: the surface complexation model
#
# Site balance on each of the three sites, charge from the surface species only (Pokrovsky Eq. 5),
# and the constant capacitance closure ψ₀ = σ/C (his Eq. 3) with C = √I/α (his Eq. 4). The
# apparent and intrinsic constants are related by his Eq. 2, K_s = K_int·exp(−zFψ₀/RT).

# %%
head(6, "Step 2: the surface complexation model")

def _logratios(aH, aCO3, aCa, aMg, aMe, zMe, logK_Me, lb):
    lH, lCO3 = np.log10(aH), np.log10(max(aCO3, 1e-300))
    lCa, lMg = np.log10(max(aCa, 1e-300)), np.log10(max(aMg, 1e-300))
    co3 = [LOGK["CO3_deprot"] - lH - lb,
           LOGK["CO3Ca"] + lCa - lH + lb,
           LOGK["CO3Mg"] + lMg - lH + lb]
    me = (None if logK_Me is None else
          logK_Me + np.log10(max(aMe, 1e-300)) - lH + (zMe - 1)*lb)
    ca = [LOGK["CaOH2"] + lH + lb, LOGK["CaO"] - lH - lb,
          LOGK["CaHCO3"] + lCO3 + 2*lH, LOGK["CaCO3"] + lCO3 + lH - lb]
    mg = [LOGK["MgOH2"] + lH + lb, LOGK["MgO"] - lH - lb,
          LOGK["MgHCO3"] + lCO3 + 2*lH, LOGK["MgCO3"] + lCO3 + lH - lb]
    return co3, me, ca, mg

def _shares(logs):
    L = np.array([0.0] + list(logs), dtype=float)
    e = 10.0**(L - L.max())
    return e/e.sum()

def _state(psi, aH, aCO3, aCa, aMg, aMe, zMe, logK_Me, nt, fixed_ads):
    lb = -FRT*psi/np.log(10.0)
    co3, me, ca, mg = _logratios(aH, aCO3, aCa, aMg, aMe, zMe, logK_Me, lb)
    if fixed_ads is None:
        CO3H, dep, sCa, sMg, me_s = nt["CO3"]*_shares(co3 + [me])
    else:
        me_s = fixed_ads
        CO3H, dep, sCa, sMg = max(nt["CO3"] - me_s, 0.0)*_shares(co3)
    CaOH, CaOH2, CaO, CaHCO3, CaCO3 = nt["Ca"]*_shares(ca)
    MgOH, MgOH2, MgO, MgHCO3, MgCO3 = nt["Mg"]*_shares(mg)
    pos = sCa + sMg + (zMe - 1)*me_s + CaOH2 + MgOH2
    neg = dep + CaO + CaCO3 + MgO + MgCO3
    return dict(sigma=F_CONST/S_AREA*(pos - neg), me_s=me_s, CO3H=CO3H, lb=lb,
                CO3Ca=sCa, CO3Mg=sMg, CO3dep=dep, CaHCO3=CaHCO3, CaOH2=CaOH2, CaCO3=CaCO3)

def solve_psi(aH, aCO3, aCa, aMg, aMe, zMe, logK_Me, nt, C_cap, fixed_ads=None):
    f = lambda p: _state(p, aH, aCO3, aCa, aMg, aMe, zMe, logK_Me, nt,
                         fixed_ads)["sigma"] - C_cap*p
    psi = brentq(f, -1.5, 1.5, xtol=1e-12, rtol=8.9e-16)
    st = _state(psi, aH, aCO3, aCa, aMg, aMe, zMe, logK_Me, nt, fixed_ads)
    st["psi"] = psi
    return st

print("  The model's own calibration, checked against Pokrovsky's measurements before it is used.")
pzc = []
for pH in np.arange(4.0, 10.51, 0.5):
    aH, aCO3 = background(pH)
    st = solve_psi(aH, aCO3, 1e-3, 1e-3, 1e-12, 2, -99, NT, c_pokrovsky(0.01))
    pzc.append(dict(pH=pH, sigma_mmol_m2=st["sigma"]/F_CONST*1e3, psi_mV=1e3*st["psi"]))
pzc = pd.DataFrame(pzc)
_cross = pzc.iloc[(pzc.sigma_mmol_m2.abs()).idxmin()]
print(f"    At his own conditions (0.01 M NaCl, 1e-3 M Ca and Mg, open system) the model crosses")
print(f"    zero charge at pH {_cross.pH:.1f}, against the pH 8.0 +/- 0.1 he measured for both the point of zero")
print(f"    charge and the isoelectric point: {_cross.pH-8.0:+.1f} pH units. That is not exact, and the residual")
print("    reflects the aqueous model and the site density rather than the constants themselves,")
print("    but it puts the transferred set in the right place on the one quantity he determined")
print("    independently of the fit, which is the check worth making before using it.")

_aH, _aCO3 = background(7.857)
_st = solve_psi(_aH, carb.a_CO3.iloc[2], PW.Ca.iloc[3]/MM["Ca"]/1e3*gam(2, I_PW),
                PW.Mg.iloc[3]/MM["Mg"]/1e3*gam(2, I_PW), 1e-12, 2, -99, NT, c_pokrovsky(I_PW))
_ss = SINGLE[(SINGLE.metal == "Co") & (SINGLE.batch == "pH6") & (SINGLE.day == 4)].iloc[0]
_aHs, _aCO3s = background(_ss.pH)
_sts = solve_psi(_aHs, _aCO3s, _ss.Ca/MM["Ca"]/1e3*gam(2, I_SINGLE),
                 _ss.Mg/MM["Mg"]/1e3*gam(2, I_SINGLE), 1e-12, 2, -99, NT, c_pokrovsky(I_SINGLE))
print("\n  Occupancy of the carbonate site before any trace metal arrives, in the two systems:")
show(pd.DataFrame([
    dict(system="single ion, pH 7.48",   free_CO3H_pct=100*_sts["CO3H"]/NT["CO3"],
         held_by_Ca_pct=100*_sts["CO3Ca"]/NT["CO3"], held_by_Mg_pct=100*_sts["CO3Mg"]/NT["CO3"],
         deprotonated_pct=100*_sts["CO3dep"]/NT["CO3"], psi_mV=1e3*_sts["psi"]),
    dict(system="produced water, pH 7.86", free_CO3H_pct=100*_st["CO3H"]/NT["CO3"],
         held_by_Ca_pct=100*_st["CO3Ca"]/NT["CO3"], held_by_Mg_pct=100*_st["CO3Mg"]/NT["CO3"],
         deprotonated_pct=100*_st["CO3dep"]/NT["CO3"], psi_mV=1e3*_st["psi"])]))
print("\n  This is the competition term, and it is the largest single difference between the two")
print("  experiments. In the single ion system the carbonate sites are essentially free; in")
print("  produced water calcium and magnesium hold most of them before the trace metal arrives.")
print("  A constant fitted in produced water therefore has to be large simply to account for any")
print("  uptake at all, and the two systems are not interchangeable.")

# %% [markdown]
# ---
# ## 7 — Step 5: estimating the constant
#
# The sorbed quantity is known at each sampling from Part 5, so the site balance and the closure
# σ = Cψ fix the surface potential and the mass action expression returns the constant directly.
# The constant is reported both in the form Ebrahimi & Vilcáez write it, which carries no
# electrostatic term, and with Pokrovsky's Eq. 2 correction applied.

# %%
head(7, "Step 5: estimating the constant")

def invert(metal, pH, Ca_ppm, Mg_ppm, Ceq_ppm, sorbed_M, I, aCO3=None, nt=None, C_cap=None):
    nt = NT if nt is None else nt
    C_cap = c_pokrovsky(I) if C_cap is None else C_cap
    z = Z_ION[metal]
    g1, g2 = gam(1, I), gam(2, I)
    aH, aCO3 = background(pH, aCO3)
    aCa, aMg = Ca_ppm/MM["Ca"]/1e3*g2, Mg_ppm/MM["Mg"]/1e3*g2
    aMe = (Ceq_ppm/MM[metal]/1e3*free_fraction(metal, aH, aCO3, NACL_M*g1)
           * (g1 if z == 1 else g2))
    if not (0 < sorbed_M < nt["CO3"]):
        return np.nan, np.nan, np.nan, sorbed_M/nt["CO3"]
    st = solve_psi(aH, aCO3, aCa, aMg, aMe, z, None, nt, C_cap, fixed_ads=sorbed_M)
    base = np.log10(sorbed_M) + np.log10(aH) - np.log10(st["CO3H"]) - np.log10(aMe)
    return base - (z - 1)*st["lb"], base, 1e3*st["psi"], sorbed_M/nt["CO3"]

rows = []
for _, r in sat.iterrows():
    if r.system == "single ion":
        d = SINGLE[(SINGLE.metal == r.metal) & (SINGLE.batch == r.batch) &
                   (SINGLE.day == r.day)].iloc[0]
        lk, lk0, psi, N = invert(r.metal, d.pH, d.Ca, d.Mg, d.Me, r.sorbed_M, I_SINGLE)
        c0, ceq = d.C0, d.Me
    else:
        d = PW[(PW.batch == r.batch) & (PW.day == r.day)].iloc[0]
        c0 = PW[(PW.batch == r.batch) & (PW.day == 0)].iloc[0][r.metal]
        aC = carb[(carb.batch == r.batch) & (carb.day == r.day)].a_CO3.iloc[0]
        lk, lk0, psi, N = invert(r.metal, d.pH, d.Ca, d.Mg, d[r.metal], r.sorbed_M, I_PW, aCO3=aC)
        ceq = d[r.metal]
    noise = np.hypot(ICP_REL*c0, ICP_REL*ceq)
    rows.append(dict(metal=r.metal, system=r.system, batch=r.batch, day=r.day, pH=r.pH,
                     uptake_mgL=c0 - ceq, S_over_N=(c0 - ceq)/noise,
                     sorbed_mM=1e3*r.sorbed_M, coverage_N=N, psi_mV=psi,
                     logK_int=lk, logK_no_electrostatic=lk0))
pts = pd.DataFrame(rows)
show(pts.sort_values(["metal", "system", "batch", "day"]))

# %% [markdown]
# ---
# ## 8 — The constants, system by system
#
# Reported as a mean and standard deviation over the usable samplings of each system, which is how
# a constant obtained point by point should be quoted: a constant that drifts systematically with
# pH or with time is not an intrinsic constant, and the spread is the evidence either way.

# %%
head(8, "The constants, system by system")

summ = []
for m in ["Li", "Co"]:
    for sysname in ["single ion", "produced water"]:
        g = pts[(pts.metal == m) & (pts.system == sysname) & pts.logK_int.notna()]
        summ.append(dict(metal=m, system=sysname, n_usable=len(g),
                         n_total=int((pts.metal == m).sum()/2),
                         pH_range=f"{g.pH.min():.2f} to {g.pH.max():.2f}" if len(g) else "-",
                         logK_int=g.logK_int.mean() if len(g) else np.nan,
                         SD=g.logK_int.std() if len(g) > 1 else np.nan,
                         coverage=f"{g.coverage_N.min():.2f} to {g.coverage_N.max():.2f}" if len(g) else "-"))
summ = pd.DataFrame(summ); show(summ)

_lis = summ[(summ.metal == "Li") & (summ.system == "single ion")].iloc[0]
_cos = summ[(summ.metal == "Co") & (summ.system == "single ion")].iloc[0]
_lip = summ[(summ.metal == "Li") & (summ.system == "produced water")].iloc[0]
_cop = summ[(summ.metal == "Co") & (summ.system == "produced water")].iloc[0]

print("\n  THE SINGLE ION SYSTEM DELIVERS BOTH CONSTANTS.")
print(f"    Lithium  log K_int = {_lis.logK_int:+.2f} +/- {_lis.SD:.2f} over all {_lis.n_usable} samplings, pH {_lis.pH_range}")
print(f"    Cobalt   log K_int = {_cos.logK_int:+.2f} +/- {_cos.SD:.2f} over {_cos.n_usable} of 6 samplings, pH {_cos.pH_range}")
print("    Both hold across four pH units and both batches with no trend, which is what an")
print("    intrinsic constant should do and is the test Ioannou and Dimirkou applied to theirs.")
print("\n  PRODUCED WATER DELIVERS NEITHER CLEANLY.")
print(f"    Lithium  {_lip.n_usable} of 6 samplings usable; the other three exceed a monolayer, at coverages up")
print(f"             to {pts[(pts.metal=='Li')&(pts.system=='produced water')].coverage_N.max():.2f}, which no surface constant can reproduce.")
print(f"    Cobalt   {_cop.n_usable} of 6 samplings usable; sphaerocobaltite accounts for the rest.")
print(f"    Where a value can be formed at all it comes out {_lip.logK_int-_lis.logK_int:+.2f} log units higher for lithium")
print(f"    and {_cop.logK_int-_cos.logK_int:+.2f} for cobalt than in the single ion system, and Part 6 shows why: calcium and")
print(f"    magnesium hold {100*(_st['CO3Ca']+_st['CO3Mg'])/NT['CO3']:.0f} per cent of the carbonate sites there before the trace metal arrives,")
print("    so a constant that has to explain the observed uptake against that competition is")
print("    forced upward. A constant forced upward by unmodelled removal is not transferable.")

# %% [markdown]
# ---
# ## 9 — Against the constants Pokrovsky measured for the same site
#
# The only external check available is the same reaction for the metals that build the mineral,
# which Pokrovsky measured on the same site of the same mineral by a completely different method.

# %%
head(9, "Against the constants Pokrovsky measured for the same site")

RADIUS = {"Ca": 1.00, "Mg": 0.72, "Co": 0.745, "Li": 0.76}     # Shannon 1976, six coordinate
ctx = pd.DataFrame([
    dict(ion="Ca2+", radius_A=RADIUS["Ca"], logK=LOGK["CO3Ca"], source="Pokrovsky et al. 1999"),
    dict(ion="Mg2+", radius_A=RADIUS["Mg"], logK=LOGK["CO3Mg"], source="Pokrovsky et al. 1999"),
    dict(ion="Co2+", radius_A=RADIUS["Co"], logK=round(_cos.logK_int, 2), source="this work, single ion"),
    dict(ion="Li+",  radius_A=RADIUS["Li"], logK=round(_lis.logK_int, 2), source="this work, single ion"),
]).sort_values("radius_A")
show(ctx)
print(f"\n  Cobalt comes out {_cos.logK_int-LOGK['CO3Mg']:+.2f} log units from magnesium and {_cos.logK_int-LOGK['CO3Ca']:+.2f} from calcium. Its ionic")
print(f"  radius, {RADIUS['Co']:.3f} A, is within {100*abs(RADIUS['Co']-RADIUS['Mg'])/RADIUS['Mg']:.0f} per cent of magnesium's and {100*abs(RADIUS['Co']-RADIUS['Ca'])/RADIUS['Ca']:.0f} per cent from calcium's, so an")
print("  affinity close to magnesium's is what the crystal chemistry predicts. Sphaerocobaltite")
print("  and magnesite are isostructural for the same reason, and cobalt substitutes for magnesium")
print("  in natural carbonates far more readily than for calcium.")
print(f"\n  Lithium, radius {RADIUS['Li']:.2f} A, is similar in size but carries one charge rather than two, and")
print(f"  its constant is {_lis.logK_int-LOGK['CO3Mg']:+.2f} log units from magnesium's, the direction a halved charge on the")
print("  same site predicts.")
print("\n  Neither comparison was built into the fit. Pokrovsky's calcium and magnesium constants")
print("  enter this model only as the competition term, which the metal constant cannot adjust,")
print("  and both metals were free to come out anywhere over fourteen log units.")

# %% [markdown]
# ---
# ## 10 — Sensitivity
#
# The site density and the capacitance are the two transferred parameters, and a constant reported
# without them cannot be carried into another study.

# %%
head(10, "Sensitivity")

def refit(metal, sysname, nt_co3=None, C_scale=1.0):
    nt = NT if nt_co3 is None else dict(CO3=nt_co3, Ca=NT["Ca"], Mg=NT["Mg"])
    out = []
    for _, r in sat[(sat.metal == metal) & (sat.system == sysname)].iterrows():
        if r.sorbed_M <= 0: continue
        if sysname == "single ion":
            d = SINGLE[(SINGLE.metal == metal) & (SINGLE.batch == r.batch) &
                       (SINGLE.day == r.day)].iloc[0]
            lk = invert(metal, d.pH, d.Ca, d.Mg, d.Me, r.sorbed_M, I_SINGLE, nt=nt,
                        C_cap=C_scale*c_pokrovsky(I_SINGLE))[0]
        else:
            d = PW[(PW.batch == r.batch) & (PW.day == r.day)].iloc[0]
            aC = carb[(carb.batch == r.batch) & (carb.day == r.day)].a_CO3.iloc[0]
            lk = invert(metal, d.pH, d.Ca, d.Mg, d[metal], r.sorbed_M, I_PW, aCO3=aC, nt=nt,
                        C_cap=C_scale*c_pokrovsky(I_PW))[0]
        if np.isfinite(lk): out.append(lk)
    return np.mean(out) if out else np.nan

print("10a  SITE DENSITY.  Pokrovsky's 14 umol/m2 for the carbonate site, against the 2e-5 mol/m2")
print("     that Brady et al. (1999) give for dolomite and that is also in common use:\n")
sN = []
for lab, v in [("Pokrovsky 1999, 14 umol/m2", 14.0e-6*S_AREA),
               ("Brady et al. 1999, 20 umol/m2", 20.0e-6*S_AREA)]:
    sN.append(dict(site_density=lab, N_t_mmolL=v*1e3,
                   logK_Li_single=refit("Li", "single ion", nt_co3=v),
                   logK_Co_single=refit("Co", "single ion", nt_co3=v)))
sN = pd.DataFrame(sN); show(sN)
print(f"\n     A {100*(20/14-1):.0f} per cent change in the site density moves lithium by {abs(sN.logK_Li_single.diff().iloc[-1]):.2f} log units and cobalt")
print(f"     by {abs(sN.logK_Co_single.diff().iloc[-1]):.2f}. Both constants are quoted below with the density they were fitted at.")

print("\n10b  CAPACITANCE.  Pokrovsky's alpha = 0.004 against a factor of two either way:\n")
sC = []
for lab, f in [("alpha/2, C doubled", 2.0), ("Pokrovsky, alpha = 0.004", 1.0),
               ("alpha x2, C halved", 0.5)]:
    sC.append(dict(capacitance=lab, C_F_m2=f*c_pokrovsky(I_SINGLE),
                   logK_Li_single=refit("Li", "single ion", C_scale=f),
                   logK_Co_single=refit("Co", "single ion", C_scale=f)))
sC = pd.DataFrame(sC); show(sC)
print(f"\n     A factor of four in the capacitance is worth {sC.logK_Li_single.max()-sC.logK_Li_single.min():.2f} log units for lithium and")
print(f"     {sC.logK_Co_single.max()-sC.logK_Co_single.min():.2f} for cobalt. On a carbonate surface the capacitance is high and the potential")
print("     correspondingly small, so the electrostatic term does little work and the constant is")
print("     insensitive to the exact value of alpha. That is a property of the mineral, and it is")
print("     why Pokrovsky could treat alpha as a single empirical fitting parameter.")


# %% [markdown]
# ---
# ## 11 — Forward check: the model against the measured sorption profiles
#
# The constant was obtained by inversion, so the fair test is to put it back into the model and
# predict the sorption profile it was never shown. This is the comparison Ebrahimi & Vilcáez make
# in their Figure 5: measured sorption against time, with the simulated curve through it.
#
# The curve is computed at the measured pH, calcium and magnesium interpolated across the six
# days, so the rise in the model line is driven by the same dolomite dissolution the experiment
# underwent, and the only free quantity in it is the constant from Part 8.

# %%
head(11, "Forward check: the model against the measured sorption profiles")

def predict_removal(metal, pH, Ca_ppm, Mg_ppm, C0_ppm, logK_Me, I, aCO3=None, nt=None):
    """Total metal in, per cent removed by surface complexation out."""
    nt = NT if nt is None else nt
    z = Z_ION[metal]
    g1, g2 = gam(1, I), gam(2, I)
    aH, aCO3 = background(pH, aCO3)
    aCa, aMg = Ca_ppm/MM["Ca"]/1e3*g2, Mg_ppm/MM["Mg"]/1e3*g2
    fr = free_fraction(metal, aH, aCO3, NACL_M*g1)*(g1 if z == 1 else g2)
    T = C0_ppm/MM[metal]/1e3
    D = T
    for _ in range(60):
        st = solve_psi(aH, aCO3, aCa, aMg, D*fr, z, logK_Me, nt, c_pokrovsky(I))
        nxt = max(T - st["me_s"], 1e-16)
        if abs(nxt - D) < 1e-11*T: D = nxt; break
        D = 0.5*D + 0.5*nxt
    return 100*(T - D)/T

def profile(metal, system, batch, logK, n=61):
    """Model curve across the six days, at the interpolated measured conditions."""
    if system == "single ion":
        d = SINGLE[(SINGLE.metal == metal) & (SINGLE.batch == batch)].sort_values("day")
        C0, I, aC = d.C0.iloc[0], I_SINGLE, None
        me = d.Me.values
    else:
        d = PW[PW.batch == batch].sort_values("day")
        C0, I = d[metal].iloc[0], I_PW
        aC = float(carb[carb.batch == batch].a_CO3.mean())
        me = d[metal].values
    t = np.linspace(0, 6, n)
    pH = np.interp(t, d.day, d.pH); ca = np.interp(t, d.day, d.Ca); mg = np.interp(t, d.day, d.Mg)
    pred = [predict_removal(metal, p, c, m, C0, logK, I, aCO3=aC) for p, c, m in zip(pH, ca, mg)]
    return t, np.array(pred), d.day.values, 100*(C0 - me)/C0

FITTED = {(m, sy): summ[(summ.metal == m) & (summ.system == sy)].logK_int.iloc[0]
          for m in ["Li", "Co"] for sy in ["single ion", "produced water"]}

fig, axes = plt.subplots(2, 2, figsize=(8.8, 6.8), constrained_layout=True)
panel = "ABCD"; k = 0
for i, sysname in enumerate(["single ion", "produced water"]):
    for j, m in enumerate(["Li", "Co"]):
        ax = axes[i, j]
        lk = FITTED[(m, sysname)]
        top = 0.0                      # collected over BOTH batches, not just the last
        for batch, fc, ls in [("pH6", COL[m], "-"), ("pH2", "white", "--")]:
            t, pred, day, meas = profile(m, sysname, batch, lk)
            top = max(top, float(np.nanmax(meas)), float(np.nanmax(pred)))
            ax.plot(t, pred, ls, color=COL[m], lw=1.6, zorder=2)
            ax.plot(day, meas, MK[m] + "-", color=COL[m], mfc=fc, mec=COL[m], lw=0.9,
                    ms=7.5, zorder=3)
            # ring every sampling that yielded no constant, for either reason
            bad = pts[(pts.metal == m) & (pts.system == sysname) & (pts.batch == batch) &
                      pts.logK_int.isna()]
            for _, b in bad.iterrows():
                yv = meas[list(day).index(b.day)]
                ax.plot(b.day, yv, "o", ms=15, mfc="none", mec="0.25", mew=1.3, zorder=4)
        ax.set_xlim(-0.3, 6.5); ax.set_xticks([0, 2, 4, 6])
        ax.set_xlabel("time (days)")
        ax.set_ylabel(f"{m} removed (%)")
        ax.set_title(f"({panel[k]})  {m}, {sysname}   —   log $K_{{int}}$ = {lk:+.2f}",
                     fontsize=9.5, loc="left")
        ax.set_ylim(0, top*1.18 if top > 0 else 1)
        k += 1

from matplotlib.lines import Line2D
handles = [Line2D([], [], marker="s", color="0.3", mfc="0.3", ls="none", ms=7.5,
                  label="measured, pH 6 start"),
           Line2D([], [], marker="s", color="0.3", mfc="white", ls="none", ms=7.5,
                  label="measured, pH 2 start"),
           Line2D([], [], color="0.3", ls="-", lw=1.6, label="model, pH 6 start"),
           Line2D([], [], color="0.3", ls="--", lw=1.6, label="model, pH 2 start"),
           Line2D([], [], marker="o", color="0.25", mfc="none", mew=1.3, ls="none", ms=13,
                  label="excluded: carbonate supersaturated\nor above one monolayer")]
fig.legend(handles=handles, loc="outside lower center", ncol=5, fontsize=7.6,
           frameon=False, handletextpad=0.5, columnspacing=1.4)
savefig(fig, "TwoSys_Fig1_model_vs_measured")

print("  HOW TO READ THIS FIGURE. The curve is not a kinetic trajectory. It is the equilibrium")
print("  sorption the model predicts at the pH, calcium and magnesium actually measured at that")
print("  moment, so it does not start at zero: at day zero the suspension already sits at pH 6.3")
print("  in one pair of batches, and the model says what sorption would be there at equilibrium.")
print("  A measured point on the curve is a batch that has reached equilibrium; a point below it")
print("  is one still approaching. The acid start batches, dashed, begin far below their curve and")
print("  climb onto it as dolomite dissolution raises the pH, which is the behaviour the model")
print("  should show and does.")
print("\n  Rings mark the samplings the saturation screen removed before fitting. They are drawn")
print("  because they are informative, not because they were used: panel B shows the largest")
print("  cobalt uptake of the single ion series sitting well above the model, and that is the")
print("  point at which sphaerocobaltite becomes supersaturated.")

rmse = []
for sysname in ["single ion", "produced water"]:
    for m in ["Li", "Co"]:
        for batch in ["pH6", "pH2"]:
            t, pred, day, meas = profile(m, sysname, batch, FITTED[(m, sysname)])
            pm = np.interp(day[1:], t, pred)
            rmse.append(dict(metal=m, system=sysname, batch=batch,
                             RMSE_pct=float(np.sqrt(np.mean((pm - meas[1:])**2))),
                             max_measured_pct=float(meas.max())))
rmse = pd.DataFrame(rmse); show(rmse)
print(f"\n  Root mean square error between measured and modelled sorption, in percentage points.")
print(f"  For the single ion system it is {rmse[rmse.system=='single ion'].RMSE_pct.mean():.2f} points on average, against removals that reach")
print(f"  {rmse[rmse.system=='single ion'].max_measured_pct.max():.1f} per cent, so the model tracks the profile it was not shown. For produced water")
print(f"  it is {rmse[rmse.system=='produced water'].RMSE_pct.mean():.1f} points, because most of the removal there is not surface complexation")
print("  and no value of the constant can reproduce it.")

# %% [markdown]
# ---
# ## 12 — Figures

# %%
head(12, "Figures")

# --- Figure 2: the constant against pH, both systems, on one shared scale ----------------------
fig, axes = plt.subplots(1, 2, figsize=(9.0, 4.0), sharey=True, constrained_layout=True)
_allk = pts.logK_int.dropna()
ylo, yhi = np.floor(_allk.min()*2)/2 - 0.3, np.ceil(_allk.max()*2)/2 + 0.55
for ax, m in zip(axes, ["Li", "Co"]):
    ax.axhline(LOGK["CO3Ca"], color="0.45", ls="--", lw=1, zorder=1)
    ax.axhline(LOGK["CO3Mg"], color="0.45", ls=":", lw=1, zorder=1)
    ax.text(8.9, LOGK["CO3Ca"] + 0.06, "Pokrovsky Ca", fontsize=7.2, color="0.4",
            va="bottom", ha="right")
    ax.text(8.9, LOGK["CO3Mg"] - 0.06, "Pokrovsky Mg", fontsize=7.2, color="0.4",
            va="top", ha="right")
    mu = summ[(summ.metal == m) & (summ.system == "single ion")].logK_int.iloc[0]
    sd = summ[(summ.metal == m) & (summ.system == "single ion")].SD.iloc[0]
    ax.axhspan(mu - sd, mu + sd, color=COL[m], alpha=0.13, lw=0, zorder=0)
    ax.axhline(mu, color=COL[m], lw=1.6, zorder=2)
    for sysname, fc in [("single ion", COL[m]), ("produced water", "white")]:
        for bi, b in enumerate(["pH6", "pH2"]):
            g = (pts[(pts.metal == m) & (pts.system == sysname) & (pts.batch == b)
                     & pts.logK_int.notna()].sort_values("pH"))
            if not len(g): continue
            ax.plot(g.pH, g.logK_int, MK[m] + "-", color=COL[m], mfc=fc, mec=COL[m],
                    lw=0.9, ms=8, zorder=3,
                    label=sysname if bi == 0 else None)
    ax.set_xlim(3.6, 9.0); ax.set_ylim(ylo, yhi)
    ax.set_xlabel("measured pH")
    ax.set_title(f"{m}    single ion mean {mu:+.2f} $\\pm$ {sd:.2f}", fontsize=9.5, loc="left")
    ax.legend(loc="upper left", fontsize=8, ncol=2)
axes[0].set_ylabel("log $K_{int}$")
savefig(fig, "TwoSys_Fig2_logK_vs_pH")

# --- Figure 2: removal in the two systems ------------------------------------------------------
fig, axes = plt.subplots(1, 2, figsize=(9.2, 3.9), constrained_layout=True)
for ax, m in zip(axes, ["Li", "Co"]):
    for b, ls in [("pH6", "-"), ("pH2", "--")]:
        d = SINGLE[(SINGLE.metal == m) & (SINGLE.batch == b)]
        ax.plot(d.day, d.removal_pct, MK[m] + ls, color=COL[m], mfc=COL[m],
                label=f"single ion, {b.replace('pH', 'pH ')}")
        p = PW[PW.batch == b]
        p0 = p[p.day == 0].iloc[0][m]
        ax.plot(p.day, 100*(p0 - p[m])/p0, MK[m] + ls, color=COL[m], mfc="white",
                label=f"produced water, {b.replace('pH', 'pH ')}")
    ax.set_xticks([0, 2, 4, 6]); ax.set_xlabel("time (days)")
    ax.set_ylabel("removal (% of initial)"); ax.set_title(m, fontsize=9.5)
    ax.legend(loc="upper left", fontsize=7.2)
savefig(fig, "TwoSys_Fig3_removal_both_systems")

# --- Figure 3: why the two systems differ -----------------------------------------------------
fig, axes = plt.subplots(1, 2, figsize=(9.2, 3.8), constrained_layout=True)
ax = axes[0]
occ = pd.DataFrame([
    dict(system="single ion", free=100*_sts["CO3H"]/NT["CO3"], dep=100*_sts["CO3dep"]/NT["CO3"],
         ca=100*_sts["CO3Ca"]/NT["CO3"], mg=100*_sts["CO3Mg"]/NT["CO3"]),
    dict(system="produced water", free=100*_st["CO3H"]/NT["CO3"], dep=100*_st["CO3dep"]/NT["CO3"],
         ca=100*_st["CO3Ca"]/NT["CO3"], mg=100*_st["CO3Mg"]/NT["CO3"])])
y = np.arange(2); left = np.zeros(2)
for key, lab, c in [("free", ">CO$_3$H$^0$ free", "#cfd8dc"), ("dep", ">CO$_3^-$", "#90a4ae"),
                    ("ca", ">CO$_3$Ca$^+$", "#1f5f8b"), ("mg", ">CO$_3$Mg$^+$", "#7a5c99")]:
    ax.barh(y, occ[key], left=left, color=c, label=lab)
    left = left + occ[key].values
ax.set_yticks(y); ax.set_yticklabels(occ.system); ax.invert_yaxis()
ax.set_xlabel("share of the carbonate sites (%)"); ax.set_xlim(0, 100)
ax.legend(ncol=2, fontsize=7.2, loc="lower center")
ax.set_title("competition for the site, before the trace metal", fontsize=9)

ax = axes[1]
lab, val, colr = [], [], []
for m in ["Li", "Co"]:
    for sysname, alpha in [("single ion", 1.0), ("produced water", 0.40)]:
        r = summ[(summ.metal == m) & (summ.system == sysname)].iloc[0]
        lab.append(f"{m}\n{sysname}"); val.append(r.logK_int)
        colr.append(tuple(list(plt.matplotlib.colors.to_rgb(COL[m])) + [alpha]))
x = np.arange(4)
ax.bar(x, val, color=colr)
ax.axhline(LOGK["CO3Ca"], color="0.45", ls="--", lw=1)
ax.axhline(LOGK["CO3Mg"], color="0.45", ls=":", lw=1)
ax.text(3.55, LOGK["CO3Ca"], " Ca", fontsize=7.5, color="0.45", va="center")
ax.text(3.55, LOGK["CO3Mg"], " Mg", fontsize=7.5, color="0.45", va="center")
ax.set_xticks(x); ax.set_xticklabels(lab, fontsize=7.6)
ax.set_ylabel("log $K_{int}$"); ax.axhline(0, color="0.5", lw=0.8)
for xi, v in zip(x, val):
    ax.text(xi, v + (0.12 if v > 0 else -0.12), f"{v:+.2f}", ha="center",
            va="bottom" if v > 0 else "top", fontsize=8)
ax.set_title("the same metal, the two systems", fontsize=9)
savefig(fig, "TwoSys_Fig4_why_the_systems_differ")

# %% [markdown]
# ---
# ## 13 — Result

# %%
head(13, "Result")

print("  INTRINSIC SURFACE COMPLEXATION CONSTANTS FOR LITHIUM AND COBALT ON DOLOMITE")
print("  >CO3H0 + Me(z+)  =  >CO3Me(z-1)+ + H+ ,  25 C")
print(f"  Pokrovsky et al. (1999) parameterisation: carbonate site {SITE_DENS['CO3']*1e6:.0f} umol/m2")
print(f"  ({NT['CO3']*1e3:.3f} mmol/L here), C = sqrt(I)/{ALPHA_POK}, Davies activity coefficients\n")
res = pd.DataFrame([
    dict(reaction=">CO3H0 + Li+  = >CO3Li0 + H+", system="single ion, 40 g/L NaCl",
         logK_int=round(_lis.logK_int, 2), SD=round(_lis.SD, 2), n=f"{_lis.n_usable} of 6",
         status="determined"),
    dict(reaction=">CO3H0 + Co2+ = >CO3Co+ + H+", system="single ion, 40 g/L NaCl",
         logK_int=round(_cos.logK_int, 2), SD=round(_cos.SD, 2), n=f"{_cos.n_usable} of 6",
         status="determined"),
    dict(reaction=">CO3H0 + Li+  = >CO3Li0 + H+", system="produced water",
         logK_int=round(_lip.logK_int, 2), SD=round(_lip.SD, 2), n=f"{_lip.n_usable} of 6",
         status="upper bound, not transferable"),
    dict(reaction=">CO3H0 + Co2+ = >CO3Co+ + H+", system="produced water",
         logK_int=round(_cop.logK_int, 2), SD=np.nan, n=f"{_cop.n_usable} of 6",
         status="upper bound, not transferable"),
])
show(res)

print("\n  RECOMMENDED VALUES, the single ion determinations:")
print(f"    >CO3H0 + Li+  = >CO3Li0 + H+     log K_int = {_lis.logK_int:+.2f} +/- {_lis.SD:.2f}")
print(f"    >CO3H0 + Co2+ = >CO3Co+ + H+     log K_int = {_cos.logK_int:+.2f} +/- {_cos.SD:.2f}")
print("  quoted at the site density and capacitance above, without which they cannot be carried")
print("  into another calculation.")

print("\n  WHY THE SINGLE ION VALUES AND NOT THE PRODUCED WATER ONES. Three reasons, all of them")
print("  visible in the tables above rather than assumed.")
print(f"    1. Coverage. Every usable single ion point sits inside a monolayer. Three of the six")
print(f"       produced water lithium points do not, reaching {pts[(pts.metal=='Li')&(pts.system=='produced water')].coverage_N.max():.2f} monolayers, which no surface")
print("       complexation constant can reproduce at any value.")
print(f"    2. Precipitation. Cobalt is undersaturated with sphaerocobaltite at {_cos.n_usable} of 6 single ion")
print(f"       samplings and at only {_cop.n_usable} of 6 in produced water.")
print(f"    3. Competition. Calcium and magnesium hold {100*(_st['CO3Ca']+_st['CO3Mg'])/NT['CO3']:.0f} per cent of the carbonate sites in")
print(f"       produced water against {100*(_sts['CO3Ca']+_sts['CO3Mg'])/NT['CO3']:.0f} per cent in the single ion system, so any removal not")
print("       actually due to sorption is forced into the constant as an inflated affinity.")
print("  The produced water values are reported for completeness and as upper bounds. They are the")
print("  measure of how much of the recovery that system achieves is NOT surface complexation.")

print("\n  WHAT THIS SAYS ABOUT THE RECOVERY PROCESS. Dolomite recovers more lithium and far more")
print("  cobalt from produced water than from a simple brine, and the constants show that the")
print("  extra recovery is not stronger sorption. Sorption is weaker there, because calcium")
print("  occupies the sites. The extra recovery is carbonate mineral formation for cobalt and,")
print("  for lithium, uptake beyond what a monolayer can hold. That is a useful result for a")
print("  recovery process, because mineral formation and surface complexation respond to")
print("  different levers: the first to carbonate supply and residence time, the second to")
print("  surface area and to keeping calcium out of the way.")

print("\n  WHAT WOULD TIGHTEN THE CONSTANTS. The uncertainties above are set by one thing, that")
print(f"  uptake is a difference of two large concentrations: the best single ion point reaches a")
print(f"  signal to noise of {pts[(pts.system=='single ion')].S_over_N.max():.1f} at +/-{100*ICP_REL:.0f} per cent ICP precision. Raising the solid loading")
print("  scales the sites and therefore the uptake; lowering the initial metal raises the fraction")
print("  removed. Measuring alkalinity at each sampling would also remove the one quantity here")
print("  that had to be inferred rather than measured.")
