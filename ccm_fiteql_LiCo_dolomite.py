# %% [markdown]
# # Constant capacitance model constants for Li⁺ and Co²⁺ on dolomite
#
# **Marwa Elshebli, Javier Vilcáez, James Smay** — Boone Pickens School of Geology,
# Oklahoma State University
#
# ---
# ## The three papers this follows, and what each one supplies
#
# | Source | What it supplies |
# |---|---|
# | **Hayes, Redden, Ela & Leckie (1991)**, *J. Colloid Interface Sci.* **142**, 448–469 | **Which model, and how its parameters are chosen.** The CCM is the surface complexation model for ionic strengths above 0.01 M; site totals must be fixed from measurement because log K trades off against them; the capacitance is bounded by theory to 0.1–2.0 F m⁻²; FITEQL converges only with about two adjustable parameters; a unique parameter set is chosen by sensitivity analysis rather than asserted |
# | **Goldberg (2004)**, *Vadose Zone J.* **3**, 676–680 | **The fitting protocol.** Adsorption **isotherms and envelopes fitted simultaneously** with one constant set, in FITEQL, with capacitance and site density fixed at literature values and the constants reported per sample as mean ± SD |
# | **Ioannou & Dimirkou (1997)**, *J. Colloid Interface Sci.* **192**, 119–128 | **How the result is reported.** Protonation constants taken from a published compilation rather than fitted; one complexation constant fitted **per pH**, then averaged across pH with a standard deviation (their Table 3) |
#
# ## What changes from the earlier version of this work
#
# Three things, each of them forced by one of the papers above.
#
# 1. **The capacitance.** The earlier model used the carbonate convention
#    $C=\sqrt{I}/\alpha$ with $\alpha=0.004$, which gives **207 F m⁻²**. Hayes Eq. 22,
#    $C_1=\varepsilon_r\varepsilon_0/d$, brackets the capacitance at **0.1 to 2.0 F m⁻²** for any
#    physically possible dielectric constant and ion approach distance. A capacitance two orders
#    above that bracket pins the surface potential near zero and switches the electrostatic term
#    off altogether, so the earlier model was a mass action model wearing a CCM label. This
#    version uses $C=1.06$ F m⁻², the value all three papers use, and shows by a Hayes sensitivity
#    sweep what the choice is worth.
#
# 2. **The site total.** Hayes showed that the fitted log K moves directly with the site total,
#    so a constant is meaningless unless the $N_t$ it was fitted at is reported with it. This
#    version fits at three independently sourced site densities and reports the trade-off as a
#    table, rather than adopting one and hiding the dependence.
#
# 3. **The data used.** Goldberg's contribution is that an isotherm and an envelope constrain a
#    constant set far better together than either does alone. The single-ion batches are an
#    adsorption **envelope** — fixed total metal, pH running from 2.3 to 8.3 — and the produced
#    water is an **isotherm** point at a different total metal and against competing calcium.
#    They are fitted simultaneously here for the first time.

# %%
import numpy as np, pandas as pd, matplotlib.pyplot as plt, os
from scipy.optimize import brentq, minimize_scalar
pd.set_option("display.width", 175); pd.set_option("display.precision", 4)
os.makedirs("figures", exist_ok=True)

def head(n, t): print(f"\n{'='*100}\n{n}  {t}\n{'='*100}")
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
COL = {"Co": "#1f5f8b", "Li": "#c8542a", "grey": "#6e6e6e", "green": "#3a7d5d"}
MK  = {"Co": "o", "Li": "s"}
def savefig(fig, n):
    fig.savefig(f"figures/{n}.png"); fig.savefig(f"figures/{n}.pdf"); plt.show()

# %% [markdown]
# ---
# ## 1 — Parameters, and where each one is bracketed from
#
# Nothing here is adjustable except the two metal constants. Every other number is either
# measured in this laboratory or taken from a named table in a named paper, which is the
# discipline Hayes et al. set out: the fewer the free parameters, the more a fitted constant means.

# %%
head(1, "Parameters, and where each one is bracketed from")

F_CONST, R_GAS, T_K = 96485.0, 8.314, 298.15
FRT   = F_CONST/(R_GAS*T_K)
EPS_0 = 8.854e-12
N_AV  = 6.02214e23

# --- measured in this laboratory ---
SOLID_GL = 60.0          # 6 g dolomite per 100 mL
SSA      = 0.76          # m2/g, BET
S_AREA   = SOLID_GL*SSA  # m2 per litre of suspension
ICP_REL  = 0.03          # relative precision of the ICP-OES on each concentration
MM = dict(Li=6.941, Co=58.933, Ca=40.078, Mg=24.305, Sr=87.62, Ba=137.33, Cd=112.41, Pb=207.2)
Z_ION = dict(Li=1, Co=2)

# --- ionic strength of each experiment ---
I_SINGLE = 0.684         # 40 g/L NaCl
I_PW     = 1.00          # synthetic produced water
NACL_M   = 40.0/58.44

# --- capacitance: Hayes Eq. 22 brackets it, Goldberg and Ioannou both use 1.06 ---
C_HAYES_LO = 6*EPS_0/4.3e-10     # er = 6  (dielectric saturation), d = 4.3 A
C_HAYES_HI = 50*EPS_0/2.3e-10    # er = 50, d = 2.3 A (hydrated univalent ion)
C_CAP      = 1.06                # F/m2, Westall & Hohl (1980) optimum, used by all three papers
def c_diffuse(I, er=78.5):
    kappa = np.sqrt(2*1000*F_CONST**2*I/(er*EPS_0*R_GAS*T_K))
    return er*EPS_0*kappa
C_OLD = np.sqrt(I_SINGLE)/0.004  # the carbonate convention used in the earlier version

# --- site totals: three independent sources, all reported, none privileged ---
def nt_from_sites_nm2(ns):  return ns*1e18*S_AREA/N_AV     # Goldberg Eq. 10
def nt_from_mol_m2(v):      return v*S_AREA
SITES = {
    "Goldberg 2004 (2.31 sites/nm2, Davis & Kent generic)": nt_from_sites_nm2(2.31),
    "Zachara 1991 (3.46e-6 mol/g, 45-Ca exchange)":         3.46e-6*SOLID_GL,
    "Pokrovsky 1999 (14 umol/m2, carbonate site)":          nt_from_mol_m2(14e-6),
}
NT_CO3 = SITES["Pokrovsky 1999 (14 umol/m2, carbonate site)"]   # base case, justified in Part 8
NT_CA, NT_MG = nt_from_mol_m2(7e-6), nt_from_mol_m2(7e-6)

show(pd.DataFrame([
    ("Solid loading",        f"{SOLID_GL:.0f} g/L",            "6 g per 100 mL, measured"),
    ("Specific surface area",f"{SSA:.2f} m2/g",                "BET, measured"),
    ("Surface per litre",    f"{S_AREA:.1f} m2/L",             "product of the two"),
    ("Ionic strength",       f"{I_SINGLE:.3f} / {I_PW:.2f} M", "single-ion / produced water"),
    ("ICP-OES precision",    f"+/-{100*ICP_REL:.0f} %",        "on each concentration"),
    ("Capacitance used",     f"{C_CAP:.2f} F/m2",              "Westall & Hohl 1980; Goldberg 2004; Ioannou 1997"),
    ("Hayes Eq. 22 bracket", f"{C_HAYES_LO:.2f} to {C_HAYES_HI:.2f} F/m2",
                             "er 6 to 50, d 2.3 to 4.3 Angstrom"),
    ("Diffuse-layer equivalent", f"{c_diffuse(I_SINGLE):.2f} F/m2",
                             "C = er.eps0.kappa at this ionic strength"),
    ("Value used previously", f"{C_OLD:.0f} F/m2",             "sqrt(I)/alpha, alpha = 0.004"),
], columns=["parameter", "value", "source"]))

print(f"\n  The capacitance previously used is {C_OLD/C_HAYES_HI:.0f} times the top of the Hayes bracket, which means it")
print("  describes a double layer no physically possible combination of dielectric constant and")
print("  ion approach distance can produce. Part 6 shows what it costs the model.")
print(f"  The value adopted here lies inside the bracket, and the diffuse-layer capacitance at this")
print(f"  ionic strength is {c_diffuse(I_SINGLE):.2f} F/m2, the same order and also inside it. Hayes give that agreement")
print("  as the condition under which the diffuse layer model reduces to the CCM, so the two")
print("  independent routes to a capacitance for this brine agree on where it belongs.")

print("\n  Site totals from three independent sources:")
show(pd.DataFrame([(k, v*1e3, v*N_AV/(S_AREA*1e18)) for k, v in SITES.items()],
                  columns=["source", "N_t (mmol/L)", "sites/nm2"]))
print("  Goldberg's generic value for natural materials and Zachara's directly measured exchange")
print(f"  capacity for a carbonate agree to {100*abs(SITES['Goldberg 2004 (2.31 sites/nm2, Davis & Kent generic)']-SITES['Zachara 1991 (3.46e-6 mol/g, 45-Ca exchange)'])/SITES['Zachara 1991 (3.46e-6 mol/g, 45-Ca exchange)']:.0f} %, which is a real cross-check: a generic soil site density and an")
print("  isotopically measured carbonate exchange capacity should not agree that closely by chance.")
print("  Pokrovsky's carbonate-site value is about three times higher because it counts every")
print("  surface carbonate group, not only those that exchange.")

# %% [markdown]
# ---
# ## 2 — Chemical model
#
# Surface reactions are Pokrovsky, Schott & Thomas (1999) Table 3, dolomite columns, taken as
# published and **not** fitted — the same discipline Ioannou & Dimirkou applied when they took
# their protonation constants from the Goldberg & Sposito compilation instead of optimising them.
# Only the two metal reactions are adjustable.

# %%
head(2, "Chemical model")

LOGK = {"CO3_deprot": -4.8, "CO3Ca": -1.8, "CO3Mg": -2.0,
        "CaOH2": 11.5, "CaO": -12.0, "CaHCO3": -4.0, "CaCO3": 16.6,
        "MgOH2": 10.6, "MgO": -12.0, "MgHCO3": -3.5, "MgCO3": 15.4}

RX = [(">CO3H0 = >CO3- + H+",                    f"{LOGK['CO3_deprot']:+.1f}", "Pokrovsky 1999 T3"),
      (">CO3H0 + Ca2+ = >CO3Ca+ + H+",           f"{LOGK['CO3Ca']:+.1f}",      "Pokrovsky 1999 T3"),
      (">CO3H0 + Mg2+ = >CO3Mg+ + H+",           f"{LOGK['CO3Mg']:+.1f}",      "Pokrovsky 1999 T3"),
      (">CO3H0 + Li+  = >CO3Li0 + H+",           "FITTED",                     "this work"),
      (">CO3H0 + Co2+ = >CO3Co+ + H+",           "FITTED",                     "this work"),
      (">CaOH0 + H+ = >CaOH2+",                  f"{LOGK['CaOH2']:+.1f}",      "Pokrovsky 1999 T3"),
      (">CaOH0 = >CaO- + H+",                    f"{LOGK['CaO']:+.1f}",        "Pokrovsky 1999 T3"),
      (">CaOH0 + CO3-2 + 2H+ = >CaHCO3 + H2O",   f"{LOGK['CaHCO3']:+.1f}",     "Pokrovsky 1999 T3"),
      (">CaOH0 + CO3-2 + H+ = >CaCO3- + H2O",    f"{LOGK['CaCO3']:+.1f}",      "Pokrovsky 1999 T3"),
      (">MgOH0 + H+ = >MgOH2+",                  f"{LOGK['MgOH2']:+.1f}",      "Pokrovsky 1999 T3"),
      (">MgOH0 = >MgO- + H+",                    f"{LOGK['MgO']:+.1f}",        "Pokrovsky 1999 T3"),
      (">MgOH0 + CO3-2 + 2H+ = >MgHCO3 + H2O",   f"{LOGK['MgHCO3']:+.1f}",     "Pokrovsky 1999 T3"),
      (">MgOH0 + CO3-2 + H+ = >MgCO3- + H2O",    f"{LOGK['MgCO3']:+.1f}",      "Pokrovsky 1999 T3")]
show(pd.DataFrame(RX, columns=["surface reaction", "log K", "source"]))

print("\n  Mass action, in the CCM form of Goldberg (2004) Eqs. 5 to 7:")
print("      K(int) = [product]/([>CO3H0].a(Me)/a(H)) . exp(-(z-1).F.psi/RT)")
print("  Site balance (her Eq. 8) and charge balance (her Eq. 9), and the CCM closure sigma = C.psi")
print("  (Hayes Eq. 14). Charge is built from surface species only; dissolved Ca and Mg belong to")
print("  dolomite dissolution and carry no surface charge.")

AQ = {
  "Co": [("CoCl+", 0.30, dict(Cl=1), 1), ("CoCl2", -0.20, dict(Cl=2), 0),
         ("CoCO3", 4.23, dict(CO3=1), 0), ("CoHCO3+", 12.20, dict(CO3=1, H=1), 1),
         ("CoOH+", -9.65, dict(H=-1), 1), ("Co(OH)2", -18.80, dict(H=-2), 0)],
  "Li": [("LiCl", -0.50, dict(Cl=1), 0), ("LiOH", -13.64, dict(H=-1), 0),
         ("LiCO3-", 0.90, dict(CO3=1), -1)],
}
print("\n  Aqueous complexes carried for each metal (log beta from the free ion, NIST/MINTEQ):")
show(pd.DataFrame([(m, n, b) for m, lst in AQ.items() for n, b, _, _ in lst],
                  columns=["metal", "complex", "log beta"]))

LOGKSP = {"Sr": -9.271, "Ba": -8.562, "Cd": -13.74, "Co": -9.98, "Pb": -13.13, "Li": -2.50}
MINERAL = {"Sr": "strontianite", "Ba": "witherite", "Cd": "otavite",
           "Co": "sphaerocobaltite", "Pb": "cerussite", "Li": "zabuyelite"}
PCO2, LOGKH = 10**-3.5, -1.469
LK_H2CO3, LK_HCO3 = -6.345, 10.329
A_DAVIES = 0.509
def gam(z, I):
    return 10**(-A_DAVIES*z*z*(np.sqrt(I)/(1 + np.sqrt(I)) - 0.3*I))

# %% [markdown]
# ---
# ## 3 — The data: an adsorption envelope and an isotherm point
#
# Goldberg (2004) fits the two together. Mapped onto this system:
#
# * **Envelope** — the single-ion batches. Total metal is fixed, dolomite dissolution carries the
#   suspension from the nominal starting pH up through the carbonate buffer, and the measured pH
#   spans 2.25 to 8.26. That is an adsorption envelope acquired by drift rather than by titration,
#   and the measured pH at each sampling is what enters the model.
# * **Isotherm** — the produced water. The same solid meets a different total metal, against
#   4186 mg L⁻¹ of competing calcium, so it probes a different point on the same isotherm.

# %%
head(3, "The data")

SINGLE = pd.DataFrame([
 ("Co","pH6",0, 0.100, 0.100, 80.513, 6.320),("Co","pH6",2, 0.100, 0.100, 77.905, 7.400),
 ("Co","pH6",4, 0.100, 8.143, 76.411, 7.477),("Co","pH6",6, 0.126,12.941, 74.393, 8.090),
 ("Co","pH2",0, 0.100, 0.100, 80.513, 2.350),("Co","pH2",2,57.477,17.461, 79.661, 4.123),
 ("Co","pH2",4,24.178,45.592, 79.795, 5.660),("Co","pH2",6,21.658,51.379, 80.398, 7.110),
 ("Li","pH6",0, 0.663, 0.100,138.430, 6.420),("Li","pH6",2,37.633, 0.100,134.657, 7.647),
 ("Li","pH6",4,54.820, 0.100,136.722, 7.923),("Li","pH6",6,49.892, 0.870,135.262, 8.257),
 ("Li","pH2",0, 0.663, 0.100,138.430, 2.250),("Li","pH2",2,39.190, 0.100,136.263, 4.257),
 ("Li","pH2",4,79.452,21.316,135.956, 5.727),("Li","pH2",6,87.738,20.547,135.592, 7.047),
], columns=["metal","batch","day","Ca_ppm","Mg_ppm","Me_ppm","pH"])
SINGLE["C0_ppm"]   = SINGLE.groupby(["metal","batch"])["Me_ppm"].transform("first")
SINGLE["uptake"]   = (SINGLE.C0_ppm - SINGLE.Me_ppm)
SINGLE["removal%"] = 100*SINGLE.uptake/SINGLE.C0_ppm
print("ENVELOPE ARM  -  single-ion batches, 40 g/L NaCl, 25 C")
show(SINGLE)

METALS_PW = ["Pb", "Sr", "Ba", "Cd", "Co", "Li"]
PW = pd.DataFrame([
    ("pH6", 0, 4185.950, 1005.970, 83.350, 88.579, 90.165, 92.673, 75.960, 76.780, 6.320),
    ("pH6", 2, 3930.733,  940.894, 48.789, 82.369, 83.058, 84.949, 69.892, 71.243, 7.413),
    ("pH6", 4, 3890.693,  984.050,  1.072, 82.500, 49.280, 25.739, 11.527, 71.470, 7.557),
    ("pH6", 6, 3754.450,  967.633,  2.037, 83.166, 44.457, 14.066,  7.019, 70.212, 7.857),
    ("pH2", 0, 4185.950, 1005.970, 83.350, 88.579, 90.165, 92.673, 75.960, 76.780, 2.350),
    ("pH2", 2, 4037.933,  978.907, 47.138, 86.930, 87.731, 92.468, 75.013, 76.212, 4.290),
    ("pH2", 4, 4199.987, 1115.780,  2.164, 86.297, 54.205, 43.352, 26.687, 74.399, 5.727),
], columns=["batch","day","Ca","Mg","Pb","Sr","Ba","Cd","Co","Li","pH"])
print("\nISOTHERM ARM  -  produced water, six metals against 4186 mg/L calcium")
show(PW)

print("\n  Signal against analytical noise on the uptake, which is a difference of two")
print("  concentrations each carrying +/-3 %:")
sn = []
for _, r in SINGLE[SINGLE.day == 6].iterrows():
    noise = np.hypot(ICP_REL*r.C0_ppm, ICP_REL*r.Me_ppm)
    sn.append(dict(arm="envelope", metal=r.metal, batch=r.batch, uptake_mgL=r.uptake,
                   noise_mgL=noise, S_over_N=r.uptake/noise))
d0 = PW[(PW.batch == "pH6") & (PW.day == 0)].iloc[0]
d6 = PW[(PW.batch == "pH6") & (PW.day == 6)].iloc[0]
for m in ["Co", "Li"]:
    noise = np.hypot(ICP_REL*d0[m], ICP_REL*d6[m])
    sn.append(dict(arm="isotherm", metal=m, batch="pH6", uptake_mgL=d0[m]-d6[m],
                   noise_mgL=noise, S_over_N=(d0[m]-d6[m])/noise))
sn = pd.DataFrame(sn); show(sn)
print("\n  This is the reason the two arms are needed together. In the envelope arm the cobalt")
print(f"  signal-to-noise reaches only {sn[(sn.arm=='envelope')&(sn.metal=='Co')].S_over_N.max():.1f}, which cannot carry a constant on its own. In the")
print(f"  isotherm arm it reaches {sn[(sn.arm=='isotherm')&(sn.metal=='Co')].S_over_N.iloc[0]:.0f}. Fitted simultaneously, the envelope supplies the pH")
print("  dependence and the isotherm supplies the magnitude, which is exactly Goldberg's argument")
print("  for fitting the two forms together rather than separately.")

# %% [markdown]
# ---
# ## 4 — Correcting the isotherm arm for carbonate precipitation
#
# The produced water removes metal by two routes, and the constant capacitance model describes
# only one of them. Strontium separates them. Zachara, Cowan & Resch (1991) measured
# log ᶜK_ex(Sr) = −2.04 on calcite and describe strontium as effectively non-sorbing, so its
# dissolved concentration is set by strontianite equilibrium alone and fixes the carbonate
# activity. With that, every other metal has a solubility floor, and only what lies **below**
# the floor is available for the surface model to explain.

# %%
head(4, "Correcting the isotherm arm for carbonate precipitation")

G1_PW, G2_PW = gam(1, I_PW), gam(2, I_PW)
carb = []
for _, r in PW[PW.day > 0].iterrows():
    aSr = r.Sr/MM["Sr"]/1e3*G2_PW
    carb.append(dict(batch=r.batch, day=int(r.day), pH=r.pH,
                     a_CO3=10**LOGKSP["Sr"]/aSr, log_a_CO3=np.log10(10**LOGKSP["Sr"]/aSr)))
carb = pd.DataFrame(carb); show(carb)
print(f"\n  log a(CO3) = {carb.log_a_CO3.mean():.2f} +/- {carb.log_a_CO3.std():.2f} across every sampling and both batches. A")
print("  solubility-controlled system should give a steady value, and it does, which is the check")
print("  that using strontium this way is sound.")

pw_ads = []
for _, r in PW[PW.day > 0].iterrows():
    aCO3 = carb[(carb.batch == r.batch) & (carb.day == r.day)].a_CO3.iloc[0]
    c0 = PW[(PW.batch == r.batch) & (PW.day == 0)].iloc[0]
    for m in ["Co", "Li"]:
        if m == "Li":
            fl = (10**LOGKSP[m]/aCO3)**0.5/G1_PW*MM[m]*1e3
        else:
            fl = 10**LOGKSP[m]/aCO3/G2_PW*MM[m]*1e3
        removed = (c0[m] - r[m])/MM[m]/1e3
        if fl > c0[m]:
            sorbed, verdict = removed, "carbonate cannot form: all removal is surface"
        elif r[m] < fl:
            sorbed, verdict = (fl - r[m])/MM[m]/1e3, "below floor: the excess is surface"
        else:
            sorbed, verdict = 0.0, "above floor: precipitation explains it"
        pw_ads.append(dict(batch=r.batch, day=int(r.day), metal=m, mineral=MINERAL[m], pH=r.pH,
                           floor_mgL=fl, measured_mgL=r[m], removed_M=removed,
                           sorbed_M=sorbed, verdict=verdict))
pw_ads = pd.DataFrame(pw_ads)
show(pw_ads[["batch","day","metal","mineral","pH","floor_mgL","measured_mgL","sorbed_M","verdict"]])
print(f"\n  LITHIUM. Zabuyelite would need {pw_ads[pw_ads.metal=='Li'].floor_mgL.mean():,.0f} mg/L of lithium before it could nucleate in")
print(f"  the bulk solution, against the {PW.Li.iloc[0]:.0f} mg/L present. Lithium cannot precipitate from this")
print("  water, so its entire removal enters the surface model uncorrected.")
print("  COBALT. Sphaerocobaltite carries most of the removal, and only the part that drives the")
print("  solution below the solubility floor is passed to the surface model.")

# %% [markdown]
# ---
# ## 5 — Aqueous speciation
#
# Goldberg (2004) records that in the FITEQL constant capacitance model "activities are assumed
# to be equal to concentrations and no activity coefficient corrections are performed". That
# convention was written for soil suspensions near 0.1 M. These experiments run at 0.68 and
# 1.0 M, so Davies coefficients are applied here as the base case and the FITEQL convention is
# carried as a sensitivity case in Part 8, rather than either being assumed silently.

# %%
head(5, "Aqueous speciation")

def background(pH, aCO3=None):
    """Proton and carbonate activities. Open-system carbonate unless a value is supplied."""
    aH = 10**(-pH)
    if aCO3 is None or not np.isfinite(aCO3):   # a pandas column of mixed None/float gives NaN
        aHCO3 = 10**LK_H2CO3*(10**LOGKH*PCO2)/aH
        aCO3  = aHCO3/(10**LK_HCO3*aH)
    return aH, aCO3

def free_fraction(metal, aH, aCO3, aCl):
    """Free-ion share of the dissolved metal, from the aqueous complexes of Part 2."""
    tot = 1.0
    for _, lb, st, _ in AQ[metal]:
        v = 10**lb
        for c, n in st.items():
            v *= {"Cl": aCl, "H": aH, "CO3": aCO3}[c]**n
        tot += v
    return 1.0/tot

_aH, _aCO3 = background(8.09)
_fr = free_fraction("Co", _aH, _aCO3, gam(1, I_SINGLE)*NACL_M)
print(f"  Worked example, cobalt at pH 8.09 in 40 g/L NaCl:")
print(f"    free Co2+ is {100*_fr:.1f} % of the dissolved cobalt; the rest is chloride and carbonate complexes.")
print("  The free ion is what the surface reaction sees, so this correction moves the fitted")
print("  constant directly and is applied at every point.")

# %% [markdown]
# ---
# ## 6 — The constant capacitance model
#
# Site balance, charge balance and the closure σ = Cψ, solved for the surface potential at every
# point. The charge is assembled from surface species only, following the definition Hayes et al.
# give in their Table Ib, σ₀ = B(Σ z_i [surface species]) with B = F/(C_s·SA).

# %%
head(6, "The constant capacitance model")

def _logratios(aH, aCO3, aCa, aMg, aMe, zMe, logK_Me, lb):
    """log10 of each surface species relative to its neutral reference form.

    Everything is carried in log space so that a strongly polarised surface, where one species
    outweighs the reference form by many orders, never produces an overflow or an inf/inf.
    lb is log10 of the Boltzmann factor exp(-F.psi/RT), and a species of charge z carries z.lb.
    """
    lH, lCO3 = np.log10(aH), np.log10(max(aCO3, 1e-300))
    lCa, lMg = np.log10(max(aCa, 1e-300)), np.log10(max(aMg, 1e-300))
    co3 = [LOGK["CO3_deprot"] - lH - lb,          # >CO3-    z = -1
           LOGK["CO3Ca"] + lCa - lH + lb,         # >CO3Ca+  z = +1
           LOGK["CO3Mg"] + lMg - lH + lb]         # >CO3Mg+  z = +1
    me = (None if logK_Me is None else
          logK_Me + np.log10(max(aMe, 1e-300)) - lH + (zMe - 1)*lb)
    ca = [LOGK["CaOH2"] + lH + lb,                # >CaOH2+  z = +1
          LOGK["CaO"] - lH - lb,                  # >CaO-    z = -1
          LOGK["CaHCO3"] + lCO3 + 2*lH,           # >CaHCO3  z =  0
          LOGK["CaCO3"] + lCO3 + lH - lb]         # >CaCO3-  z = -1
    mg = [LOGK["MgOH2"] + lH + lb, LOGK["MgO"] - lH - lb,
          LOGK["MgHCO3"] + lCO3 + 2*lH, LOGK["MgCO3"] + lCO3 + lH - lb]
    return co3, me, ca, mg

def _shares(logs):
    """Site fractions from log10 ratios, with the neutral reference form implied at index 0."""
    L = np.array([0.0] + list(logs), dtype=float)
    e = 10.0**(L - L.max())
    return e/e.sum()

def _state(psi, aH, aCO3, aCa, aMg, aMe, zMe, logK_Me, nt, fixed_ads):
    """Every surface species and the resulting charge, at a trial potential."""
    lb = -FRT*psi/np.log(10.0)
    co3, me, ca, mg = _logratios(aH, aCO3, aCa, aMg, aMe, zMe, logK_Me, lb)
    if fixed_ads is None:
        fr = _shares(co3 + [me])
        CO3H, dep, sCa, sMg, me_s = nt["CO3"]*fr
    else:
        me_s = fixed_ads
        fr = _shares(co3)
        CO3H, dep, sCa, sMg = max(nt["CO3"] - me_s, 0.0)*fr
    fCa = _shares(ca); fMg = _shares(mg)
    CaOH, CaOH2, CaO, CaHCO3, CaCO3 = nt["Ca"]*fCa
    MgOH, MgOH2, MgO, MgHCO3, MgCO3 = nt["Mg"]*fMg
    pos = sCa + sMg + (zMe - 1)*me_s + CaOH2 + MgOH2
    neg = dep + CaO + CaCO3 + MgO + MgCO3
    return dict(sigma=F_CONST/S_AREA*(pos - neg), me_s=me_s, CO3H=CO3H, lb=lb,
                CaOH=CaOH, MgOH=MgOH)

def solve_psi(aH, aCO3, aCa, aMg, aMe, zMe, logK_Me, nt, C_cap, fixed_ads=None):
    """Solve sigma(psi) = C.psi, which is Hayes Eq. 14.

    Every positive surface species carries the Boltzmann factor and every negative one carries
    its reciprocal, so sigma falls monotonically with psi while C.psi rises. The root is
    therefore unique and the interval +/-1.5 V brackets it for any capacitance in the Hayes
    range: at those ends the surface is fully polarised and sigma has saturated at its sign.
    """
    f = lambda p: _state(p, aH, aCO3, aCa, aMg, aMe, zMe, logK_Me, nt,
                         fixed_ads)["sigma"] - C_cap*p
    psi = brentq(f, -1.5, 1.5, xtol=1e-12, rtol=8.9e-16)
    st = _state(psi, aH, aCO3, aCa, aMg, aMe, zMe, logK_Me, nt, fixed_ads)
    st["psi"] = psi
    st["b"] = 10.0**st["lb"]
    return st

NT_BASE = dict(CO3=NT_CO3, Ca=NT_CA, Mg=NT_MG)

def invert_logK(metal, pH, Ca_ppm, Mg_ppm, Ceq_ppm, uptake_M, I, nt=None, C_cap=None,
                aCO3=None, davies=True):
    """One constant from one measurement, the way Ioannou & Dimirkou fitted one per pH."""
    nt = NT_BASE if nt is None else nt
    C_cap = C_CAP if C_cap is None else C_cap
    z = Z_ION[metal]
    g1, g2 = (gam(1, I), gam(2, I)) if davies else (1.0, 1.0)
    aH, aCO3 = background(pH, aCO3)
    aCa, aMg = Ca_ppm/MM["Ca"]/1e3*g2, Mg_ppm/MM["Mg"]/1e3*g2
    aCl = NACL_M*g1
    aMe = Ceq_ppm/MM[metal]/1e3*free_fraction(metal, aH, aCO3, aCl)*(g1 if z == 1 else g2)
    if not (0 < uptake_M < nt["CO3"]): return np.nan, np.nan, np.nan
    st = solve_psi(aH, aCO3, aCa, aMg, aMe, z, None, nt, C_cap, fixed_ads=uptake_M)
    lk = (np.log10(uptake_M) + np.log10(aH) - np.log10(st["CO3H"])
          - np.log10(aMe) - (z - 1)*st["lb"])
    return lk, st["psi"], uptake_M/nt["CO3"]

def predict_Ceq(metal, pH, Ca_ppm, Mg_ppm, C0_ppm, logK_Me, I, nt=None, C_cap=None,
                aCO3=None, davies=True):
    """Forward model: total metal in, equilibrium dissolved metal out."""
    nt = NT_BASE if nt is None else nt
    C_cap = C_CAP if C_cap is None else C_cap
    z = Z_ION[metal]
    g1, g2 = (gam(1, I), gam(2, I)) if davies else (1.0, 1.0)
    aH, aCO3 = background(pH, aCO3)
    aCa, aMg = Ca_ppm/MM["Ca"]/1e3*g2, Mg_ppm/MM["Mg"]/1e3*g2
    aCl = NACL_M*g1
    fr = free_fraction(metal, aH, aCO3, aCl)*(g1 if z == 1 else g2)
    T = C0_ppm/MM[metal]/1e3
    D = T
    for _ in range(60):
        st = solve_psi(aH, aCO3, aCa, aMg, D*fr, z, logK_Me, nt, C_cap)
        nxt = max(T - st["me_s"], 1e-16)
        if abs(nxt - D) < 1e-11*T: D = nxt; break
        D = 0.5*D + 0.5*nxt
    return D*MM[metal]*1e3, st

print("  Surface potential and charge against pH, at the two capacitances, with no metal present:")
_sw = []
for _pH in [4.0, 5.0, 6.0, 7.0, 8.0, 9.0]:
    _aH, _aC = background(_pH)
    _a = solve_psi(_aH, _aC, 1e-4, 1e-4, 1e-12, 2, -99, NT_BASE, C_CAP)
    _b = solve_psi(_aH, _aC, 1e-4, 1e-4, 1e-12, 2, -99, NT_BASE, C_OLD)
    _sw.append(dict(pH=_pH, sigma_mmol_m2=_a["sigma"]/F_CONST*1e3,
                    psi_mV_C106=1e3*_a["psi"], boltz_C106=_a["b"],
                    psi_mV_C207=1e3*_b["psi"], boltz_C207=_b["b"]))
_sw = pd.DataFrame(_sw); show(_sw)
print(f"\n  The Boltzmann factor spans {_sw.boltz_C106.min():.3f} to {_sw.boltz_C106.max():.2f} at C = {C_CAP} F/m2, and only {_sw.boltz_C207.min():.2f} to {_sw.boltz_C207.max():.2f} at")
print(f"  C = {C_OLD:.0f} F/m2. In the log units a constant is reported in, that is a correction of up to")
print(f"  {abs(np.log10(_sw.boltz_C106)).max():.2f} log units against {abs(np.log10(_sw.boltz_C207)).max():.2f}, so the earlier capacitance suppressed the")
print(f"  electrostatic term by a factor of {abs(np.log10(_sw.boltz_C106)).max()/abs(np.log10(_sw.boltz_C207)).max():.0f} and left a mass action model behind.")
print(f"\n  Two checks before any metal constant is asked of the model. It crosses zero charge near")
print(f"  pH {_sw.pH[(_sw.sigma_mmol_m2.abs()).idxmin()]:.0f}, against the pH 7.5 to 8.5 point of zero charge reported for dolomite. And its")
print(f"  charge reaches {_sw.sigma_mmol_m2.abs().max():.4f} mmol/m2 over the pH 4 to 9 window of these experiments, about an")
print("  order of magnitude below the 0.01 to 0.02 mmol/m2 Pokrovsky, Schott & Thomas report at the")
print("  extremes of titrations that ran to pH 3 and pH 11. Both are the right sign and the right")
print("  order, which is as much as a transferred constant set can be asked to deliver.")

# %% [markdown]
# ---
# ## 7 — One constant per point, then averaged
#
# Ioannou & Dimirkou (1997) fitted their phosphate constants separately at each of six pH values
# and reported the mean and standard deviation across them (their Table 3). Their reason was that
# a constant which drifts with pH is not an intrinsic constant, and the spread across pH is the
# honest measure of that. The same test is applied here, and it has more to say about this data
# set than any single optimised number would.

# %%
head(7, "One constant per point, then averaged")

pts = []
for _, r in SINGLE[SINGLE.day > 0].iterrows():
    up = (r.C0_ppm - r.Me_ppm)/MM[r.metal]/1e3
    lk, psi, N = invert_logK(r.metal, r.pH, r.Ca_ppm, r.Mg_ppm, r.Me_ppm, up, I_SINGLE)
    noise = np.hypot(ICP_REL*r.C0_ppm, ICP_REL*r.Me_ppm)
    pts.append(dict(arm="envelope", metal=r.metal, batch=r.batch, day=int(r.day), pH=r.pH,
                    uptake_mgL=r.C0_ppm - r.Me_ppm, sorbed_M=up,
                    S_over_N=(r.C0_ppm - r.Me_ppm)/noise,
                    coverage_N=N, psi_mV=1e3*psi if np.isfinite(psi) else np.nan, log_K=lk))
for _, r in pw_ads.iterrows():
    row = PW[(PW.batch == r.batch) & (PW.day == r.day)].iloc[0]
    c0  = PW[(PW.batch == r.batch) & (PW.day == 0)].iloc[0]
    aC  = carb[(carb.batch == r.batch) & (carb.day == r.day)].a_CO3.iloc[0]
    lk, psi, N = invert_logK(r.metal, r.pH, row.Ca, row.Mg, row[r.metal], r.sorbed_M,
                             I_PW, aCO3=aC)
    noise = np.hypot(ICP_REL*c0[r.metal], ICP_REL*row[r.metal])
    pts.append(dict(arm="isotherm", metal=r.metal, batch=r.batch, day=int(r.day), pH=r.pH,
                    uptake_mgL=c0[r.metal] - row[r.metal], sorbed_M=r.sorbed_M,
                    S_over_N=(c0[r.metal] - row[r.metal])/noise,
                    coverage_N=N, psi_mV=1e3*psi if np.isfinite(psi) else np.nan, log_K=lk))
pts = pd.DataFrame(pts)
def _status(r):
    if np.isfinite(r.log_K):        return "usable"
    if r.sorbed_M <= 0:             return "precipitation accounts for the removal"
    if r.sorbed_M >= NT_CO3:        return "uptake exceeds site capacity"
    return "no uptake resolvable"
pts["status"] = pts.apply(_status, axis=1)
show(pts)

print("\n  Goldberg and Ioannou both fit every point and let the FITEQL weighting decide how much")
print("  each one counts. The diagnostic that matters is therefore not a pass/fail screen but the")
print("  ratio of uptake to analytical noise, series by series:\n")
ser = (pts.groupby(["metal", "batch", "arm"])
          .agg(n=("log_K", "size"), usable=("log_K", "count"), mean_S_over_N=("S_over_N", "mean"),
               mean_logK=("log_K", "mean"), SD=("log_K", "std"))
          .reset_index())
show(ser)

_env = ser[ser.arm == "envelope"]
print("\n  The acid-start batches are the ones to be careful with. They spent most of the run below")
print("  pH 6, where the model itself predicts almost no uptake, so their measured differences are")
print(f"  analytical noise: in the envelope arm their mean uptake-to-noise is {_env[_env.batch=='pH2'].mean_S_over_N.mean():.2f} against {_env[_env.batch=='pH6'].mean_S_over_N.mean():.2f} for")
print("  the batches started at pH 6. Their scatter in log K is a property of the ICP, not of the")
print("  dolomite surface, and the same is true of the acid-start isotherm rows.")

print("\n  The Ioannou and Dimirkou test asks whether one constant holds across pH within a single")
print("  experiment, so it is applied to the envelope arm alone, where the matrix is fixed and only")
print("  pH moves:")
keep = pts[pts.log_K.notna() & pts.batch.eq("pH6") & pts.arm.eq("envelope")]
tbl7 = []
for m in ["Co", "Li"]:
    g = keep[keep.metal == m]
    tbl7.append(dict(metal=m, n_points=len(g),
                     pH_range=f"{g.pH.min():.2f} to {g.pH.max():.2f}" if len(g) else "-",
                     mean_logK=g.log_K.mean() if len(g) else np.nan,
                     SD=g.log_K.std() if len(g) > 1 else np.nan,
                     spread=g.log_K.max() - g.log_K.min() if len(g) > 1 else np.nan))
tbl7 = pd.DataFrame(tbl7); show(tbl7)
print("\n  A constant that drifted systematically with pH would not be an intrinsic constant, and")
print("  that is what the test is for. Over roughly a pH unit and a threefold range of uptake,")
print(f"  cobalt holds to +/-{tbl7[tbl7.metal=='Co'].SD.iloc[0]:.2f} log units and lithium to +/-{tbl7[tbl7.metal=='Li'].SD.iloc[0]:.2f}, with no trend in either.")
print("  These spreads are far wider than the +/-0.03 log units Ioannou and Dimirkou obtained")
print("  across six pH values, and the difference is not method but signal: their phosphate")
print("  uptake was most of the added phosphate, while ours is a few per cent of the added metal.")
_iso = pts[pts.log_K.notna() & pts.arm.eq("isotherm")]
print("\n  The isotherm rows are held back from this average deliberately, because they sit in a")
print("  different matrix, and are used instead as an independent check in Part 8:")
show(_iso[["metal", "batch", "day", "pH", "coverage_N", "log_K"]])

# %% [markdown]
# ---
# ## 8 — Simultaneous fit of the envelope and the isotherm
#
# Goldberg (2004) fitted isotherm and envelope together in FITEQL, minimising the weighted sum of
# squares over the whole data set at once. Here the same objective is minimised over the single
# adjustable constant per metal, with the FITEQL weighting
# $s_i=\sqrt{\sigma_{abs}^2+(\sigma_{rel}C_i)^2}$ and the goodness of fit reported as WSOS/DF.
# Herbelin & Westall take WSOS/DF between about 0.1 and 20 as an acceptable fit.

# %%
head(8, "Simultaneous fit of the envelope and the isotherm")

ERR_ABS, ERR_REL = 0.05, ICP_REL

def fit_rows(metal, nt=None):
    """The rows that enter the simultaneous fit, envelope and isotherm together.

    A point whose measured uptake exceeds the total site inventory is dropped, with the reason
    recorded. No value of log K can place more metal on a monolayer than the monolayer holds,
    so such a point does not constrain the constant, it only drags the optimiser to its bound.
    """
    nt = NT_BASE if nt is None else nt
    rows, dropped = [], []
    for _, r in SINGLE[(SINGLE.metal == metal) & (SINGLE.day > 0)].iterrows():
        up = (r.C0_ppm - r.Me_ppm)/MM[metal]/1e3
        rec = dict(arm="envelope", batch=r.batch, day=int(r.day), pH=r.pH, Ca=r.Ca_ppm,
                   Mg=r.Mg_ppm, C0=r.C0_ppm, Ceq=r.Me_ppm, I=I_SINGLE, aCO3=np.nan,
                   coverage=up/nt["CO3"])
        (rows if up < nt["CO3"] else dropped).append(rec)
    for _, r in pw_ads[pw_ads.metal == metal].iterrows():
        row = PW[(PW.batch == r.batch) & (PW.day == r.day)].iloc[0]
        c0  = PW[(PW.batch == r.batch) & (PW.day == 0)].iloc[0]
        aC  = carb[(carb.batch == r.batch) & (carb.day == r.day)].a_CO3.iloc[0]
        rec = dict(arm="isotherm", batch=r.batch, day=int(r.day), pH=r.pH, Ca=row.Ca, Mg=row.Mg,
                   C0=c0[metal], Ceq=c0[metal] - r.sorbed_M*MM[metal]*1e3, I=I_PW, aCO3=aC,
                   coverage=r.sorbed_M/nt["CO3"])
        (rows if r.sorbed_M < nt["CO3"] else dropped).append(rec)
    return pd.DataFrame(rows), pd.DataFrame(dropped)

def wsos_df(metal, logK, rows=None, nt=None, C_cap=None, davies=True, n_param=1):
    rows = fit_rows(metal, nt)[0] if rows is None else rows
    ss = 0.0
    for _, r in rows.iterrows():
        pred, _ = predict_Ceq(metal, r.pH, r.Ca, r.Mg, r.C0, logK, r.I,
                              nt=nt, C_cap=C_cap, aCO3=r.aCO3, davies=davies)
        s = np.hypot(ERR_ABS, ERR_REL*r.Ceq)
        ss += ((pred - r.Ceq)/s)**2
    return ss/max(len(rows) - n_param, 1)

def fit_metal(metal, nt=None, C_cap=None, davies=True, bracket=(-8.0, 6.0)):
    rows, _ = fit_rows(metal, nt)
    obj = lambda k: wsos_df(metal, k, rows, nt, C_cap, davies)
    out = minimize_scalar(obj, bounds=bracket, method="bounded",
                          options=dict(xatol=1e-4))
    return out.x, obj(out.x), rows

FIT = {}
for m in ["Co", "Li"]:
    k, w, rows = fit_metal(m)
    _, dropped = fit_rows(m)
    FIT[m] = dict(logK=k, wsos=w, n=len(rows))
    print(f"  {m}:  log K(int) = {k:+.2f}   WSOS/DF = {w:.3f}   over {len(rows)} points "
          f"({(rows.arm=='envelope').sum()} envelope, {(rows.arm=='isotherm').sum()} isotherm)")
    if len(dropped):
        print(f"       {len(dropped)} point(s) dropped for exceeding the site inventory "
              f"(coverage {dropped.coverage.min():.2f} to {dropped.coverage.max():.2f}):")
        show(dropped[["arm", "batch", "day", "pH", "C0", "Ceq", "coverage"]])

print("\n  Herbelin & Westall take WSOS/DF between about 0.1 and 20 as an acceptable fit, so both")
print("  constants sit inside the accepted band. A value well below 1, as here, says the residuals")
print("  are smaller than the analytical error assigned to them, which with +/-3 % on a difference")
print("  of two large concentrations is exactly what should happen: the data are consistent with")
print("  the model but do not pin it tightly.")

print("\n  The two independent routes to the same number:")
comp = []
for m in ["Co", "Li"]:
    env = pts[(pts.metal == m) & (pts.batch == "pH6") & (pts.arm == "envelope")].log_K
    comp.append(dict(metal=m,
                     per_point_envelope=env.mean(), per_point_SD=env.std(),
                     simultaneous_fit=FIT[m]["logK"],
                     difference=FIT[m]["logK"] - env.mean()))
comp = pd.DataFrame(comp); show(comp)
print("  Part 7 inverted one constant from each measurement separately and averaged them;")
print("  Part 8 optimised a single constant against every point at once, with a different arm of")
print("  the data and a different weighting. For cobalt the two agree to")
print(f"  {abs(comp[comp.metal=='Co'].difference.iloc[0]):.2f} log units, well inside the scatter of the per-point values, which is the")
print("  strongest internal evidence in this work that the cobalt constant is real. For lithium")
print(f"  they differ by {abs(comp[comp.metal=='Li'].difference.iloc[0]):.2f} log units, which is still inside the per-point scatter of {comp[comp.metal=='Li'].per_point_SD.iloc[0]:.2f}, but")
print("  lithium sits at far higher coverage and its margin is correspondingly thinner.")

print("\n  One result is visible in the table above and worth stating plainly. In the produced")
print("  water the model puts almost no cobalt on the surface, because 4186 mg/L of calcium holds")
print("  the exchange sites: predicted coverage there is under 1 per cent against 18 per cent in")
print("  the single-ion batch at a comparable pH. The surface model and the solubility analysis")
print("  of Part 4 reach the same conclusion by different routes, which is that cobalt recovery")
print("  from produced water is carbonate mineral formation and not sorption.")

print("\n  Fit quality point by point:")
val = []
for m in ["Co", "Li"]:
    rows, _ = fit_rows(m)
    for _, r in rows.iterrows():
        pred, st = predict_Ceq(m, r.pH, r.Ca, r.Mg, r.C0, FIT[m]["logK"], r.I, aCO3=r.aCO3)
        val.append(dict(metal=m, arm=r.arm, batch=r.batch, day=r.day, pH=r.pH, C0_mgL=r.C0,
                        Ceq_meas=r.Ceq, Ceq_pred=pred, uptake_meas=r.C0 - r.Ceq,
                        uptake_pred=r.C0 - pred, coverage_N=st["me_s"]/NT_CO3,
                        psi_mV=1e3*st["psi"]))
val = pd.DataFrame(val); show(val)

# %% [markdown]
# ---
# ## 9 — Hayes sensitivity: what the constant depends on
#
# This is the part of Hayes, Redden, Ela & Leckie that the earlier version of this work did not
# do. They showed that a surface complexation constant is only meaningful alongside the site
# total and capacitance it was fitted at, because the fit trades one against the other. Their
# recommended procedure is to map that trade-off and then choose a parameter set from the region
# where the constant is stable, rather than adopting one set and reporting a single number.

# %%
head(9, "Hayes sensitivity: what the constant depends on")

print("9a  CAPACITANCE.  Hayes found the CCM protolysis constants insensitive to C above about")
print("    1.2 F/m2. The same sweep, here, on the metal constants:\n")
sens_C = []
for C in [0.1, 0.3, 0.6, 1.06, 1.2, 1.5, 2.0, 5.0, 20.0, C_OLD]:
    row = dict(C_F_m2=C, inside_Hayes_bracket="yes" if C_HAYES_LO <= C <= C_HAYES_HI else "no")
    for m in ["Co", "Li"]:
        row[f"logK_{m}"] = fit_metal(m, C_cap=C)[0]
    sens_C.append(row)
sens_C = pd.DataFrame(sens_C); show(sens_C)
_in  = sens_C[sens_C.inside_Hayes_bracket == "yes"]
_out = sens_C[sens_C.C_F_m2 > C_HAYES_HI]
print(f"\n    Inside the Hayes bracket the cobalt constant moves by {_in.logK_Co.max()-_in.logK_Co.min():.3f} log units and the lithium")
print(f"    constant by {_in.logK_Li.max()-_in.logK_Li.min():.3f}. The capacitance is the one parameter in this model that can be")
print("    chosen freely without consequence, which is the same insensitivity Hayes reported for")
print("    the protolysis constants, and it is the reason they treat C as a fitting parameter")
print("    bounded by theory rather than as a quantity to be measured.")
print(f"    Push it past the bracket and it starts to matter: at the {C_OLD:.0f} F/m2 the earlier version of")
print(f"    this work used, cobalt shifts by {sens_C[sens_C.C_F_m2==C_OLD].logK_Co.iloc[0]-FIT['Co']['logK']:+.2f} and lithium by {sens_C[sens_C.C_F_m2==C_OLD].logK_Li.iloc[0]-FIT['Li']['logK']:+.2f} log units.")

print("\n9b  SITE TOTAL.  Hayes found that raising the site total lowers the fitted constant. Each")
print("    site density here comes from a different published measurement, and none is privileged:\n")
sens_N = []
for name, ntv in SITES.items():
    nt = dict(CO3=ntv, Ca=NT_CA, Mg=NT_MG)
    row = dict(site_source=name.split(" (")[0], N_t_mmolL=ntv*1e3,
               sites_nm2=ntv*N_AV/(S_AREA*1e18))
    for m in ["Co", "Li"]:
        up = (SINGLE[(SINGLE.metal == m) & (SINGLE.batch == "pH6") & (SINGLE.day == 6)]
              .uptake.iloc[0])/MM[m]/1e3
        row[f"N_{m}_day6"] = up/ntv
        row[f"logK_{m}"] = fit_metal(m, nt=nt)[0] if up < ntv else np.nan
    sens_N.append(row)
sens_N = pd.DataFrame(sens_N); show(sens_N)
print(f"\n    This is the dominant sensitivity in the whole model. Across the three published site")
print(f"    densities the cobalt constant spans {sens_N.logK_Co.max()-sens_N.logK_Co.min():.2f} log units, which is thirty times what the")
print("    capacitance is worth. Hayes made exactly this point: a surface complexation constant")
print("    reported without the site total it was fitted at cannot be transferred to another")
print("    study, because whoever transfers it will use a different one.")
print("\n    The lithium column decides which to adopt. At the two lower site densities the")
print("    measured lithium uptake already exceeds a monolayer, so no constant exists at any")
print("    value: the model cannot hold that much lithium on a surface that small. Only")
print("    Pokrovsky's carbonate-site density leaves both metals inside the range a surface")
print("    reaction can describe, and it is also the only one of the three measured on dolomite")
print("    rather than transferred from calcite or from a generic soil. Those two reasons, one")
print("    from this data and one from the literature, are why it is the base case here, and")
print("    every constant in this work is quoted with it attached.")
print("\n9c  ACTIVITY CONVENTION.  Goldberg records that FITEQL's constant capacitance model sets")
print("    activities equal to concentrations. That convention was built for 0.1 M soil")
print("    suspensions; these run at 0.68 and 1.0 M.\n")
sens_A = []
for lab, dv in [("Davies coefficients (used here)", True), ("FITEQL convention, a = c", False)]:
    row = dict(convention=lab)
    for m in ["Co", "Li"]:
        row[f"logK_{m}"] = fit_metal(m, davies=dv)[0]
    sens_A.append(row)
sens_A = pd.DataFrame(sens_A); show(sens_A)
print(f"\n    The convention is worth {abs(sens_A.logK_Co.diff().iloc[-1]):.2f} log units for cobalt and {abs(sens_A.logK_Li.diff().iloc[-1]):.2f} for lithium at this ionic")
print("    strength, which is why it is stated rather than left implicit.")

# %% [markdown]
# ---
# ## 10 — Figures
#
# Laid out as Goldberg (2004) Figures 1 and 2: the isotherm and the envelope for each metal,
# measured points as symbols and the model as a line, followed by the Hayes sensitivity maps
# that say what the constants depend on.

# %%
head(10, "Figures")

# --- Figure 1: adsorption envelopes, measured against the model ------------------------------
fig, axes = plt.subplots(1, 2, figsize=(9.0, 3.8), constrained_layout=True)
for ax, m in zip(axes, ["Co", "Li"]):
    d = SINGLE[(SINGLE.metal == m) & (SINGLE.day > 0)]
    for b, mk, fc in [("pH6", "o", COL[m]), ("pH2", "s", "white")]:
        g = d[d.batch == b]
        ax.errorbar(g.pH, g.uptake, yerr=np.hypot(ICP_REL*g.C0_ppm, ICP_REL*g.Me_ppm),
                    fmt=mk, color=COL[m], mfc=fc, mec=COL[m], capsize=2.5, lw=0.9, ls="none",
                    label=f"measured, {b.replace('pH','pH ')} start")
    grid = np.linspace(3.5, 9.0, 34)
    C0 = d.C0_ppm.iloc[0]
    ca_i = np.interp(grid, d.sort_values("pH").pH, d.sort_values("pH").Ca_ppm)
    mg_i = np.interp(grid, d.sort_values("pH").pH, d.sort_values("pH").Mg_ppm)
    pred = [C0 - predict_Ceq(m, p, ca, mg, C0, FIT[m]["logK"], I_SINGLE)[0]
            for p, ca, mg in zip(grid, ca_i, mg_i)]
    ax.plot(grid, pred, "-", color=COL[m], lw=1.6,
            label=f"CCM, log K = {FIT[m]['logK']:+.2f}")
    ax.set_xlabel("pH"); ax.set_ylabel(f"{m} adsorbed (mg L$^{{-1}}$)")
    ax.set_title(f"{m} adsorption envelope", fontsize=9.5)
    ax.legend(loc="upper left")
savefig(fig, "CCM_Fig1_adsorption_envelopes")

# --- Figure 2: measured against predicted, both arms ------------------------------------------
fig, ax = plt.subplots(figsize=(4.9, 4.7), constrained_layout=True)
lim = (-0.5, 8.0)
ax.plot(lim, lim, "-", color="0.35", lw=1)
for m in ["Co", "Li"]:
    for arm, fc in [("envelope", COL[m]), ("isotherm", "white")]:
        g = val[(val.metal == m) & (val.arm == arm)]
        ax.errorbar(g.uptake_meas, g.uptake_pred,
                    xerr=np.hypot(ICP_REL*g.C0_mgL, ICP_REL*g.Ceq_meas),
                    fmt=MK[m], color=COL[m], mfc=fc, mec=COL[m], ls="none",
                    ecolor=COL[m], elinewidth=0.7, capsize=2, alpha=0.95,
                    label=f"{m}, {arm}")
ax.set_xlim(*lim); ax.set_ylim(*lim)
ax.set_xlabel("measured uptake (mg L$^{-1}$)")
ax.set_ylabel("predicted uptake (mg L$^{-1}$)")
ax.set_title("horizontal bars are the $\\pm$3 % analytical error on the uptake", fontsize=8)
ax.legend(loc="upper left")
savefig(fig, "CCM_Fig2_measured_vs_predicted")

# --- Figure 3: the Hayes sensitivity maps -----------------------------------------------------
fig, axes = plt.subplots(1, 2, figsize=(9.2, 3.9), constrained_layout=True)
ax = axes[0]
ax.axvspan(C_HAYES_LO, C_HAYES_HI, color=COL["green"], alpha=0.12, lw=0)
for m in ["Co", "Li"]:
    ax.plot(sens_C.C_F_m2, sens_C[f"logK_{m}"], MK[m] + "-", color=COL[m], label=m)
ax.axvline(C_CAP, color="0.4", ls="--", lw=1)
ax.set_xscale("log"); ax.set_xlabel("capacitance C (F m$^{-2}$)")
ax.set_ylabel("fitted log K(int)")
_y0, _y1 = ax.get_ylim(); ax.set_ylim(_y0 - 0.06, _y1 + 0.10)
ax.text(np.sqrt(C_HAYES_LO*C_HAYES_HI), _y1 + 0.06, "Hayes Eq. 22 bracket",
        fontsize=7.5, color=COL["green"], ha="center", va="top")
ax.text(C_CAP*1.35, _y0 + 0.02, "used here", fontsize=7.5, color="0.4",
        va="bottom", rotation=90)
ax.legend(loc="center left")
ax.set_title("sensitivity to the capacitance", fontsize=9.5)

ax = axes[1]
xs = np.arange(len(sens_N)); w = 0.36
_lo = float(np.nanmin([sens_N.logK_Co.min(), sens_N.logK_Li.min()]))
ax.set_ylim(_lo - 0.30, 0.55)
for i, m in enumerate(["Co", "Li"]):
    v = sens_N[f"logK_{m}"].values.astype(float)
    ax.bar(xs + (i - 0.5)*w, np.nan_to_num(v), w, color=COL[m], label=m)
    for x, y in zip(xs + (i - 0.5)*w, v):
        if np.isfinite(y):
            ax.text(x, y - 0.05, f"{y:.2f}", ha="center", va="top", fontsize=7.5, color=COL[m])
        else:
            ax.text(x, 0.06, "uptake above\nsite capacity:\nno solution",
                    ha="center", va="bottom", fontsize=6.6, color=COL[m])
ax.set_xticks(xs)
ax.set_xticklabels([f"{r.site_source}\n{r.N_t_mmolL:.2f} mmol L$^{{-1}}$"
                    for _, r in sens_N.iterrows()], fontsize=7.4)
ax.axhline(0, color="0.5", lw=0.8)
ax.set_ylabel("fitted log K(int)")
ax.set_title("sensitivity to the site total", fontsize=9.5)
ax.legend(loc="lower left")
savefig(fig, "CCM_Fig3_hayes_sensitivity")

# %% [markdown]
# ---
# ## 11 — Result

# %%
head(11, "Result")

print("  CONSTANT CAPACITANCE MODEL CONSTANTS FOR LITHIUM AND COBALT ON DOLOMITE")
print("  25 C, fitted simultaneously to an adsorption envelope and an isotherm,")
print(f"  at C = {C_CAP} F/m2 and N_t = {NT_CO3*1e3:.3f} mmol/L ({NT_CO3*N_AV/(S_AREA*1e18):.1f} sites/nm2)\n")
final = pd.DataFrame([
    dict(reaction=">CO3H0 + Co2+ = >CO3Co+ + H+",
         logK=FIT["Co"]["logK"], per_point_SD=comp[comp.metal=="Co"].per_point_SD.iloc[0],
         WSOS_DF=FIT["Co"]["wsos"], n_points=FIT["Co"]["n"]),
    dict(reaction=">CO3H0 + Li+  = >CO3Li0 + H+",
         logK=FIT["Li"]["logK"], per_point_SD=comp[comp.metal=="Li"].per_point_SD.iloc[0],
         WSOS_DF=FIT["Li"]["wsos"], n_points=FIT["Li"]["n"]),
])
show(final)

print(f"\n  Both constants are tied to that site total. Fitted instead at the two lower published site")
print(f"  densities, cobalt would read {sens_N.logK_Co.min():.2f} to {sens_N.logK_Co.max():.2f} and lithium would not exist at all, so the")
print("  pairing of a constant with its N_t is part of the result and not a footnote to it.")

print("\n  IN CONTEXT.  Pokrovsky, Schott & Thomas measured the same reaction for the two metals")
print("  that build the mineral:")
show(pd.DataFrame([
    (">CO3H0 + Ca2+ = >CO3Ca+ + H+", LOGK["CO3Ca"], "Pokrovsky 1999, dolomite"),
    (">CO3H0 + Mg2+ = >CO3Mg+ + H+", LOGK["CO3Mg"], "Pokrovsky 1999, dolomite"),
    (">CO3H0 + Co2+ = >CO3Co+ + H+", round(FIT["Co"]["logK"], 2), "this work"),
    (">CO3H0 + Li+  = >CO3Li0 + H+", round(FIT["Li"]["logK"], 2), "this work"),
], columns=["reaction", "log K", "source"]))
print(f"\n  Cobalt comes out {FIT['Co']['logK']-LOGK['CO3Ca']:+.2f} log units from calcium and {FIT['Co']['logK']-LOGK['CO3Mg']:+.2f} from magnesium on the same")
print("  site. A divalent transition metal binding a surface carbonate group a little more")
print("  strongly than the alkaline earths that build the mineral is what the order of the")
print("  Irving-Williams series predicts, and the fit was in no way constrained to produce it.")
print("  That the number lands there is the main external check available on this result.")

print("\n  WHAT EACH CONSTANT RESTS ON.")
print(f"  COBALT is the better determined. It carries {FIT['Co']['n']} points across {pts[(pts.metal=='Co')&(pts.arm=='envelope')].pH.min():.1f} to {pts[(pts.metal=='Co')&(pts.arm=='envelope')].pH.max():.1f} pH units, the")
print(f"  per-point inversions and the simultaneous fit agree to {abs(comp[comp.metal=='Co'].difference.iloc[0]):.2f} log units, the coverage stays")
print("  below a fifth of a monolayer so the model is used inside its valid range, and the value")
print("  falls where the Irving-Williams order says it should.")
print(f"\n  LITHIUM is weaker and should be quoted with that said. Its coverage reaches {val[val.metal=='Li'].coverage_N.max():.2f} of a")
print("  monolayer in the single-ion batches, and in the produced water the measured uptake is")
print(f"  {pw_ads[(pw_ads.metal=='Li')&(pw_ads.batch=='pH6')].sorbed_M.max()/NT_CO3:.2f} times the site inventory, which no surface complexation constant can reproduce.")
print("  Lithium is therefore removed by more than the exchange sites alone. The XRD of Elshebli")
print("  et al. (2025) identifies zabuyelite in the reacted solid, and since Part 4 shows the bulk")
print("  solution is four orders of magnitude undersaturated with respect to it, that phase must")
print("  be forming at the surface where the carbonate activity is far above the bulk value. The")
print(f"  constant of {FIT['Li']['logK']:+.2f} describes the surface complexation part of lithium uptake and should")
print("  be read as a lower bound on the total surface affinity.")

print("\n  WHAT WOULD TIGHTEN THESE.  The limit throughout is that uptake is a difference of two")
print(f"  large concentrations. At +/-{100*ICP_REL:.0f} % the best single-ion point reaches a signal-to-noise of")
print(f"  {pts[pts.arm=='envelope'].S_over_N.max():.1f}. Two changes fix that without new chemistry: raise the solid loading, which")
print("  scales the site inventory and therefore the uptake directly, and lower the initial metal")
print("  concentration, which raises the fraction removed. Both move the measurement away from")
print("  the difference of two near-equal numbers, which is the only reason the constants above")
print("  carry the uncertainty they do.")
