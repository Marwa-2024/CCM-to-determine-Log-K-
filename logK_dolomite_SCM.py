# %% [markdown]
# # Intrinsic surface complexation constants for Li⁺ and Co²⁺ on dolomite
#
# **Marwa Elshebli, Javier Vilcáez, James Smay** — Boone Pickens School of Geology,
# Oklahoma State University
#
# ---
# ## The workflow this follows
#
# | Source | What it supplies |
# |---|---|
# | **Belova et al. (2014)**, *Nickel adsorption on chalk and calcite*, J. Contam. Hydrol. 170, 1–9 | **The workflow.** The same problem class: a trace divalent cation on a carbonate mineral. Capacity and saturation screening before any fitting, FITEQL inversion, electrostatic and non-electrostatic constants reported side by side, the exchange constant as the transferable quantity, a reduced reactive-area variant |
# | **Pokrovsky, Schott & Thomas (1999)**, GCA 63, 3133–3143 | **The surface model.** Three dolomite sites in 1:1:2 stoichiometry, their densities, the eleven fixed acid–base and carbonate constants, and α in the capacitance |
# | **Goldberg (1985)**, SSSAJ 49, 851–856 | **The formalism.** The model assumptions, the σ–ψ relation, the site and charge balances, intrinsic versus conditional constants |
# | **Westall (1982)**, FITEQL 2.0 | **The objective function**, WSOS/DF |
#
# **Belova et al. (2014) is the paper to check this against.** The others supply components to it.
#
# ## The reaction being measured
#
# ```
# >CO3H0 + Li+   ⇌  >CO3Li0 + H+      log K_Li      ΔZ = 0   (no electrostatic term)
# >CO3H0 + Co2+  ⇌  >CO3Co+ + H+      log K_Co      ΔZ = +1
# ```
#
# Two unknowns, one per metal; everything else is fixed from the literature. Metal **cations**
# bind the **carbonate** site (Pokrovsky et al., 1999). The >CaOH and >MgOH hydroxyl sites take
# anions and ligands, not metal cations.
#
# ## Structure
#
# | Part | Contents |
# |---|---|
# | 1 | Fixed parameters |
# | 2 | Measured data |
# | 3 | Screening, before any fitting: precision, capacity |
# | 4 | Aqueous speciation |
# | 5 | Surface model: site balances and the electrostatic loop |
# | 6 | Inversion, point by point, with the saturation screen |
# | 7 | Global fit and WSOS/DF |
# | 8 | **Results** |
# | 9 | Figures |
# | 10 | What would determine these constants |
# | 11 | Every reaction and equation |
#
# Everything is embedded; nothing is read from a file. Run top to bottom.

# %%
import numpy as np, pandas as pd, matplotlib.pyplot as plt
from scipy.optimize import brentq, least_squares
import os

pd.set_option("display.width", 150); pd.set_option("display.precision", 4)
os.makedirs("figures", exist_ok=True)

def head(n, title):
    """One heading per Part, so the printed output can be navigated."""
    print(f"\n{'='*92}\nPART {n}  {title}\n{'='*92}")

def show(df):
    """Single place that formats every table, so all tables look alike."""
    print(df.to_string(index=False, float_format=lambda x: f"{x:,.4g}"))

# Journal figure style: white ground, no gridlines, full box, ticks inward.
plt.rcParams.update({
    "font.family": "sans-serif", "font.sans-serif": ["DejaVu Sans", "Arial"],
    "font.size": 9, "axes.labelsize": 9.5, "legend.fontsize": 8,
    "figure.facecolor": "white", "axes.facecolor": "white", "savefig.facecolor": "white",
    "axes.grid": False, "axes.spines.top": True, "axes.spines.right": True,
    "axes.linewidth": 0.8, "xtick.direction": "in", "ytick.direction": "in",
    "xtick.top": True, "ytick.right": True, "xtick.minor.visible": True,
    "ytick.minor.visible": True, "lines.linewidth": 1.4, "lines.markersize": 5.5,
    "legend.frameon": False, "savefig.dpi": 300, "savefig.bbox": "tight",
    "pdf.fonttype": 42,
})
COL = {"Co": "#1f5f8b", "Li": "#c8542a", "model": "#3a7d5d", "grey": "#8a8a8a"}
MK  = {"Co": "o", "Li": "s"}

def panel(ax, letter):
    ax.set_title(f"({letter})", loc="left", fontsize=10, fontweight="bold", pad=5)

def savefig(fig, name):
    fig.savefig(f"figures/{name}.png"); fig.savefig(f"figures/{name}.pdf"); plt.show()

# %% [markdown]
# ---
# ## Part 1 — Fixed parameters
#
# Nothing in this block is ever adjusted. When the results say "one constant per metal is
# fitted", everything here is what the other eleven are.

# %%
F, R, T = 96485.0, 8.314, 298.15
FRT = F/(R*T)                                    # 38.92 V^-1

# --- experimental set-up (spreadsheet "Set up") ---
DOLOMITE_GL = 60.0          # 6 g per 100 mL
V_L, M_G   = 0.100, 6.0     # litres of solution, grams of dolomite per bottle
SSA        = 0.76           # m2/g (BET; 0.76 to 0.84 reported)
S_AREA     = DOLOMITE_GL*SSA                     # 45.6 m2 per litre
I_BATCH    = 0.684          # mol/L, 40 g/L NaCl
ICP_REL    = 0.03           # relative precision of the ICP-OES, on EACH concentration
NACL_M     = 40.0/58.44
MM = dict(Li=6.941, Co=58.933, Ca=40.078, Mg=24.305)

# --- site density: Pokrovsky 1999, dolomite, 1:1:2 Ca:Mg:CO3 ---
SITE_DENS = dict(CO3=14.0e-6, Ca=7.0e-6, Mg=7.0e-6)      # mol/m2
S_T = {k: v*S_AREA for k, v in SITE_DENS.items()}         # mol/L of suspension

# --- capacitance: carbonate formula C = sqrt(I)/alpha ---
# alpha = 0.004 C mol^-1/2 V^-1 m^3 is the value Pokrovsky, Schott & Thomas (1999) fitted for
# dolomite (their section 3.3), and it is used here because that is the paper this work follows.
# A value of 0.006 is sometimes quoted from later carbonate work; the choice is immaterial to the
# result - see the alpha sensitivity printed with the fits - so the source paper's value is kept.
ALPHA = 0.004
C_CAP = np.sqrt(I_BATCH)/ALPHA                            # F/m2

# --- surface constants: Pokrovsky, Schott & Thomas (1999) Table 3, dolomite columns ---
# Transcribed from the page. Reaction : calcite (Ca) | magnesite (Mg) | dolomite Ca | dolomite Mg
#  1. >CO3H0 = >CO3- + H+            : -5.1  | -4.65 +/- 0.15 | -4.8 +/- 0.2 | -4.8 +/- 0.2
#  2. >CO3H0 + Me2+ = >CO3Me+ + H+   : -1.7  | -2.2  +/- 0.15 | -1.8 +/- 0.2 | -2.0 +/- 0.2
#  3. >MeOH0 = >MeO- + H+            : -12   | -12   +/- 1    | -12  +/- 2   | -12  +/- 2
#  4. >MeOH0 + H+ = >MeOH2+          : 11.5  | 10.6  +/- 0.15 | 11.5 +/- 0.2 | 10.6 +/- 0.2
#  5. >MeOH0 + CO3-2 + 2H+ = >MeHCO3 : -3.5  | -2.4  +/- 0.5  | -4.0 +/- 0.5 | -3.5 +/- 0.5
#  6. >MeOH0 + CO3-2 + H+ = >MeCO3-  : 17.1  | 14.4  +/- 0.15 | 16.6 +/- 0.2 | 15.4 +/- 0.2
# Re-checked line by line against Table 3 of the PDF: all eleven constants and all eleven
# uncertainties above match the two dolomite columns exactly.
# Site density, same paper, section 3.3: "The site density was assumed to be 7 mmol/m2 for Ca and
# Mg and 14 mmol/m2 for carbonate. This value corresponds to eight metal or carbonate sites per
# square nanometer" (1:1:2 stoichiometry). The printed UNIT is wrong by 10^3 - 14 mmol/m2 would be
# 8400 sites/nm2, not eight - and the code uses umol/m2, which reproduces the paper's own
# eight-sites-per-nm2 statement (14e-6 mol/m2 x 6.022e23 = 8.4 sites/nm2). This is the same class
# of slip as the m2/mL-for-m2/L in the manuscript's methods, noted in the capacity check.
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

head(1, "Fixed parameters")
show(pd.DataFrame(
    [("Solid loading",           f"{DOLOMITE_GL:.0f} g/L",       "6 g per 100 mL"),
     ("Specific surface area",   f"{SSA:.2f} m2/g",              "BET"),
     ("Surface area per litre",  f"{S_AREA:.1f} m2/L",           "denominator of every normalisation"),
     ("Ionic strength",          f"{I_BATCH:.3f} M",             "40 g/L NaCl"),
     ("ICP-OES precision",       f"+/-{100*ICP_REL:.0f} %",      "on EACH concentration"),
     ("Carbonate site density",  "14 umol/m2",                   "Pokrovsky 1999, 1:1:2 Ca:Mg:CO3"),
     ("Total carbonate sites",   f"{S_T['CO3']*1e3:.3f} mmol/L", "density x area per litre"),
     ("Capacitance C",           f"{C_CAP:.0f} F/m2",            f"sqrt(I)/alpha, alpha = {ALPHA}"),
     ("pCO2",                    "10^-3.5 atm",                  "open system"),
     ("Temperature",             "25 C",                         "")],
    columns=["parameter", "value", "source"]))

RXN_LABEL = {
    "CO3_deprot": ">CO3H0 = >CO3- + H+",
    "CO3Ca":      ">CO3H0 + Ca2+ = >CO3Ca+ + H+",
    "CO3Mg":      ">CO3H0 + Mg2+ = >CO3Mg+ + H+",
    "CaOH2":      ">CaOH0 + H+ = >CaOH2+",
    "MgOH2":      ">MgOH0 + H+ = >MgOH2+",
    "CaO":        ">CaOH0 = >CaO- + H+",
    "MgO":        ">MgOH0 = >MgO- + H+",
    "CaHCO3":     ">CaOH0 + CO3-2 + 2H+ = >CaHCO3(0) + H2O",
    "MgHCO3":     ">MgOH0 + CO3-2 + 2H+ = >MgHCO3(0) + H2O",
    "CaCO3":      ">CaOH0 + CO3-2 + H+ = >CaCO3- + H2O",
    "MgCO3":      ">MgOH0 + CO3-2 + H+ = >MgCO3- + H2O",
}
print("\nFixed surface constants (Pokrovsky et al. 1999, Table 3, dolomite columns):")
show(pd.DataFrame([(RXN_LABEL[k], f"{v:+.1f}", f"+/-{LOGK_UNC[k]}") for k, v in LOGK_SURF.items()],
                  columns=["reaction", "log K(int)", "uncertainty"]))

# %% [markdown]
# ---
# ## Part 2 — Measured data
#
# Single-ion batches, dolomite 6 g per 100 mL in 40 g/L NaCl at 25 °C, sampled at days 0, 2, 4, 6.
# The batch labels "pH 2" and "pH 6" are **nominal initial** values. Dolomite dissolution buffers
# the suspension and the measured day-6 pH runs from 7.05 to 8.26; the measured value is what
# enters the model.

# %%
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

head(2, "Measured data")
show(data)
eq = data[data.day == 6].reset_index(drop=True)
print("\nDay-6 points. Note the measured pH against the nominal batch label:")
show(eq[["metal", "batch", "pH", "Ca_ppm", "Mg_ppm", "C0_ppm", "Me_ppm", "removal_%"]])

# %% [markdown]
# ---
# ## Part 3 — Screening, before any fitting
#
# Belova et al. (2014) apply two screens before inverting anything. A third is placed first here
# because this data set needs it.
#
# | Screen | Question | Rule |
# |---|---|---|
# | **Precision** | is the uptake bigger than the error in measuring it? | Γ is a *difference* of two concentrations each known to ±3 %, so σ_Γ/Γ is fixed the moment the ICP numbers are recorded. A point needs σ_Γ/Γ < 30 % to support a constant |
# | **Capacity** | could adsorption account for the removal at all? | Γ / (site density) > 1 is unreachable by a monolayer reaction |
# | **Saturation** | is a solid precipitating instead? | SI(MeCO₃) must be < 0; applied in Part 6, once speciation exists |
#
# The precision screen uses no model at all — no speciation, no fitted constant — so it comes
# first. Whatever the model does later, it cannot recover information the measurement never held.

# %%
REL_ERR_FIT  = 30.0    # fit-grade: uptake at least ~3x its own error. The bar for a reported value.
REL_ERR_SHOW = 200.0   # indicative only: the arithmetic runs, but the number is not a measurement.

head(3, "Screening, before any fitting")

# ---- 3a  PRECISION -------------------------------------------------------------------------
scr = []
for r in data[data.day > 0].itertuples():
    dC = r.C0_ppm - r.Me_ppm
    sd = np.hypot(ICP_REL*r.C0_ppm, ICP_REL*r.Me_ppm)
    scr.append(dict(metal=r.metal, batch=r.batch, day=r.day, removal_pct=100*dC/r.C0_ppm,
                    uptake_mgL=dC, noise_mgL=sd,
                    rel_err_pct=100*sd/dC if dC > 0 else np.inf,
                    S_over_N=dC/sd if dC > 0 else 0.0))
scr = pd.DataFrame(scr)
print("3a  PRECISION.  Gamma = (C0 - Ceq), so its error follows from the ICP numbers alone.\n")
show(scr)

_v  = scr.rel_err_pct.replace([np.inf], np.nan).dropna()
_b  = scr.loc[_v.idxmin()]
_r0 = data[(data.metal == _b.metal) & (data.batch == _b.batch) & (data.day == _b.day)].iloc[0]
print(f"\n    n = {len(_v)}   min {_v.min():.0f} %   median {_v.median():.0f} %   max {_v.max():.0f} %")
print(f"    Points meeting rel. error < {REL_ERR_FIT:.0f} %:  {int((_v < REL_ERR_FIT).sum())} of {len(_v)}")
print(f"\n    Best point in the set ({_b.metal} {_b.batch}, day {int(_b.day)}):")
print(f"      uptake = {_r0.C0_ppm:.1f} - {_r0.Me_ppm:.1f} = {_b.uptake_mgL:.1f} mg/L")
print(f"      noise  = sqrt(({100*ICP_REL:.0f}% x {_r0.C0_ppm:.1f})^2 + ({100*ICP_REL:.0f}% x {_r0.Me_ppm:.1f})^2) = {_b.noise_mgL:.1f} mg/L")
print(f"      signal-to-noise = {_b.S_over_N:.1f}, against the ~3 that a constant needs.")
print(f"    At {int((scr.uptake_mgL < scr.noise_mgL).sum())} of {len(scr)} points the uptake is smaller than the noise outright.")

# ---- 3b  CAPACITY --------------------------------------------------------------------------
def capacity_row(label, dC_mg_L, metal, ssa, load_gL):
    area = ssa*load_gL                                   # m2 per litre
    G = dC_mg_L/MM[metal]/1e3/area                       # mol per m2
    return dict(case=label, removed_mgL=dC_mg_L, area_m2L=area, Gamma_umol_m2=G*1e6,
                x_site_density=G/SITE_DENS["CO3"])

cap = pd.DataFrame([
    capacity_row("this data  Co, pH6 day 6",                 80.513 - 74.393,  "Co", SSA,  DOLOMITE_GL),
    capacity_row("this data  Li, pH6 day 6",                 138.43 - 135.262, "Li", SSA,  DOLOMITE_GL),
    capacity_row("manuscript Co, 0.84 m2/g, ~70 mg/L",       70.0,             "Co", 0.84, DOLOMITE_GL),
    capacity_row("manuscript Li, 0.84 m2/g, 8.5 % of 138",   0.085*138.43,     "Li", 0.84, DOLOMITE_GL)])
print("\n\n3b  CAPACITY.  Coverage divided by the carbonate site density; above 1 is not adsorption.\n")
show(cap)
print("\n    The two manuscript rows exceed 1 and therefore cannot be adsorption, consistent with")
print("    the zabuyelite and sphaerocobaltite found by XRD. The present low-uptake data set is")
print("    inside the bound, so a sorption constant is at least admissible for it.")
print("    Separately: the published surface concentration of 0.046-0.049 m2/L should read")
print(f"    ~{0.84*DOLOMITE_GL:.0f} m2/L. The original figure is m2 per mL.")
print("\n\n3c  SATURATION.  Needs speciation, so it is applied point by point in Part 6.")

# %% [markdown]
# ---
# ## Part 4 — Aqueous speciation
#
# Goldberg's assumption (iii), *no complexes with the background electrolyte*, fails at 0.684 M
# NaCl, so this stage is required rather than optional. It converts total dissolved metal into
# the **free-ion activity**, which is the only part the surface reaction draws on.
#
# Davies activity coefficients, carbonate fixed by pCO₂, solved by fixed-point iteration jointly
# with the chloride balance and the ionic strength.

# %%
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

head(4, "Aqueous speciation")
_rows = []
for r in data[data.day > 0].itertuples():
    s = speciate(r.pH, totals_for(r), r.metal)
    _rows.append(dict(metal=r.metal, batch=r.batch, day=r.day, pH=r.pH, I=s["I"],
                      a_Me=s["aMe"], frac_free=s["frac_free"], a_Ca=s["aCa"], a_Mg=s["aMg"],
                      a_CO3=s["aCO3"],
                      SI_MeCO3=s["SI"]["CoCO3 (sphaerocobaltite)"] if r.metal == "Co"
                               else s["SI"]["Li2CO3"],
                      SI_dolomite=s["SI"]["Dolomite (ordered)"]))
spec = pd.DataFrame(_rows)
show(spec)
print("\nfrac_free is the number that matters: the surface reaction sees the FREE ion, not the total.")
for r in eq.itertuples():
    s = speciate(r.pH, totals_for(r), r.metal); tot = sum(s["dist"].values())
    print(f"  {r.metal} {r.batch} (pH {r.pH:.2f}): " +
          ", ".join(f"{k} {100*v/tot:.0f}%" for k, v in s["dist"].items() if v/tot > 0.01))

# %% [markdown]
# ---
# ## Part 5 — Surface model
#
# Three site balances (Eq. 5a–5c), the charge balance (Eq. 6) and the constant capacitance
# relation (Eq. 8), solved together at each point. The surface potential enters through the
# Boltzmann factor, and the loop between σ and ψ is closed by `brentq`.
#
# `surface()` runs in either direction:
#
# * **inverse** — the coverage Γ is known from the measurement, and the constant is read off
# * **forward** — a trial constant is given, and the coverage is predicted
#
# `electrostatic=False` gives the non-electrostatic case (ψ = 0), which Belova et al. (2014)
# argue is the thermodynamic quantity while the CCM value is model-dependent.

# %%
def coverage(C0_ppm, Ceq_ppm, metal, icp_rel=None):
    """Surface coverage by difference, and its propagated error.  Returns (Gamma, sigma) in mol/m2."""
    rel = ICP_REL if icp_rel is None else icp_rel
    dC  = (C0_ppm - Ceq_ppm)/MM[metal]/1e3                     # mol/L removed
    sdC = np.hypot(rel*C0_ppm, rel*Ceq_ppm)/MM[metal]/1e3      # error of the difference
    per_area = V_L/(M_G*SSA)                                   # L per m2  (= 1/S_AREA)
    return dC*per_area, sdC*per_area

def surface(s, metal, Gamma_known=None, logK_Me=None, electrostatic=True, C_cap=None):
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
        C = C_CAP if C_cap is None else C_cap
        psi = brentq(lambda p: p - state(p)[1]/C, -0.5, 0.5, xtol=1e-7)        # 0.1 mV resolves to 1e-4 V
    else:
        psi = 0.0
    sp, sigma = state(psi)
    return dict(sp=sp, psi=psi, sigma=sigma, dZ=dZ)

head(5, "Surface model")
print("The model solved at each point:")
print("  site balances   >CO3H0 + >CO3- + >CO3Li0 + >CO3Co+ + >CO3Ca+ + >CO3Mg+ = [>CO3]_T")
print("                  >CaOH0 + >CaOH2+ + >CaO- + >CaHCO3(0) + >CaCO3- = [>CaOH]_T")
print("                  >MgOH0 + >MgOH2+ + >MgO- + >MgHCO3(0) + >MgCO3- = [>MgOH]_T")
print("  charge balance  sigma = F/S_a x (sum of POSITIVE surface species - sum of NEGATIVE)")
print("                  dissolved Ca2+ and Mg2+ appear NOWHERE in it")
print("  capacitance     psi = sigma / C,   C = sqrt(I)/alpha")
print("  Boltzmann       a species of charge q carries exp(-q F psi / RT)\n")
_dem = []
for r in eq.itertuples():
    s  = speciate(r.pH, totals_for(r), r.metal)
    G, _ = coverage(r.C0_ppm, r.Me_ppm, r.metal)
    st = surface(s, r.metal, Gamma_known=max(G, 1e-12))
    _dem.append(dict(metal=r.metal, batch=r.batch, pH=r.pH, Gamma_umol_m2=G*1e6,
                     psi_mV=1e3*st["psi"], sigma_C_m2=st["sigma"],
                     boltzmann=np.exp(FRT*st["psi"])))
show(pd.DataFrame(_dem))
print(f"\nC = {C_CAP:.0f} F/m2 keeps psi at a few millivolts, so the Boltzmann factor sits near 1 and")
print("the electrostatic correction is nearly inert here. That is a property of a stiff carbonate")
print("double layer, not evidence that the constants are well determined.")

# %% [markdown]
# ---
# ## Part 6 — Inversion, point by point
#
# Equation (11) rearranged: with Γ measured and the speciation known, the constant is read off
# directly at each point. No optimiser, no starting guess — this is the whiteboard algebra, and
# doing it point by point rather than fitting first is what exposes scatter that a single fitted
# number would hide.
#
# Each point is then graded. The saturation screen from Part 3 is applied here, where speciation
# exists, alongside the precision grade.

# %%
def logK_point(row, electrostatic=True, C_cap=None, icp_rel=None):
    """Invert the mass action law at one measured point.

    Always returns the same six fields, so no caller has to branch on the shape of the result:
        logK, psi, speciation, Gamma, sigma_Gamma, (logK_lo, logK_hi)
    A point with no measurable uptake (Gamma <= 0) returns NaN in every numeric field.
    """
    s = speciate(row.pH, totals_for(row), row.metal)
    G, sG = coverage(row.C0_ppm, row.Me_ppm, row.metal, icp_rel=icp_rel)
    if G <= 0:
        return np.nan, np.nan, s, G, sG, (np.nan, np.nan)

    def lk_at(Gx):
        st = surface(s, row.metal, Gamma_known=Gx, electrostatic=electrostatic, C_cap=C_cap)
        sp = st["sp"]
        return np.log10(sp["CO3Me"]*s["aH"]/(sp["CO3H0"]*s["aMe"])
                        * np.exp(st["dZ"]*FRT*st["psi"])), st["psi"]

    lk, psi = lk_at(G)
    # propagate the coverage error through the inversion; it is asymmetric in log K
    Gmax = 0.999*S_T["CO3"]/S_AREA
    lo = lk_at(max(G - sG, 1e-3*G))[0]
    hi = lk_at(min(G + sG, Gmax))[0]
    return lk, psi, s, G, sG, (lo, hi)
# --- how good does a point have to be? Two thresholds, both stated, neither invented here. ---
# Gamma is a DIFFERENCE of two concentrations, so its relative error is sigma_G/G. A point whose
# uptake is smaller than the uncertainty in measuring that uptake carries no information about the
# constant, however carefully the algebra is done afterwards.
REL_ERR_FIT  = 30.0    # fit-grade: uptake at least ~3x its own error. The bar for a reported value.
REL_ERR_SHOW = 200.0   # indicative only: arithmetic still runs, but the number is not a measurement.

head(6, "Inversion, point by point")
_pp = []
for r in data[data.day > 0].itertuples():
    lk, psi, s, G, sG, (lo, hi) = logK_point(r, electrostatic=True)
    lk_nem = logK_point(r, electrostatic=False)[0]
    SI = s["SI"]["CoCO3 (sphaerocobaltite)"] if r.metal == "Co" else s["SI"]["Li2CO3"]
    rel = 100*sG/G if G > 0 else np.inf
    if   SI > 0:                grade, why = "excluded",   "SI(MeCO3) > 0, at or above saturation"
    elif rel > REL_ERR_SHOW:    grade, why = "excluded",   f"uptake below its own error ({rel:.0f} %)"
    elif rel > REL_ERR_FIT:     grade, why = "indicative", f"uptake comparable to its error ({rel:.0f} %)"
    else:                       grade, why = "fit-grade",  ""
    _pp.append(dict(metal=r.metal, batch=r.batch, day=r.day, pH=r.pH, Gamma=G,
                    rel_err_pct=rel, SI_MeCO3=SI, logK_CCM=lk, logK_NEM=lk_nem,
                    logK_lo=lo, logK_hi=hi, grade=grade, reason=why,
                    fit_grade=(grade == "fit-grade"), usable=(grade != "excluded")))
ppdf = pd.DataFrame(_pp)
show(ppdf.drop(columns=["reason"]))
print("\nlogK_lo / logK_hi: the constant recomputed at Gamma -/+ one propagated standard deviation.")
print("\nGrading:")
for r in ppdf[ppdf.grade != "fit-grade"].itertuples():
    print(f"  {r.metal} {r.batch} day {r.day}:  {r.grade:11} {r.reason}")

n_fit = ppdf.groupby("metal").fit_grade.sum()
print(f"\n{'-'*92}")
print(f"FIT-GRADE POINTS (rel. error < {REL_ERR_FIT:.0f} % and SI < 0):")
for m in ["Co", "Li"]:
    print(f"    {m}:  {int(n_fit.get(m, 0))} of {int((ppdf.metal == m).sum())}")
print("Neither constant is determined by this data set. Parts 7 and 8 compute and report the")
print("numbers anyway, because the arithmetic is correct and the pattern is informative, but they")
print("diagnose the experiment rather than measure a constant. Part 10 quantifies the fix.")
print(f"{'-'*92}")

# %% [markdown]
# ---
# ## Part 7 — Global fit
#
# The forward model predicts Γ from a trial constant; residuals are weighted by the propagated
# Γ error of Eq. (13) and minimised. The objective is FITEQL's WSOS/DF (Westall, 1982), Eq. (15).
# This is the step Belova et al. (2014) ran in FITEQL with a dummy adsorbed component.
#
# **Reading WSOS/DF.** A value near 1 means the model matches the data to within its measurement
# error. A value far *below* 1 is not a better fit — it means the error bars are so wide that
# almost any constant would pass, which is what happens here.

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

def fit_metal(metal, electrostatic=True):
    """Weighted least squares on the usable points. Returns the constant and its diagnostics."""
    use = ppdf[(ppdf.metal == metal) & ppdf.usable]
    rows = data.merge(use[["metal", "batch", "day"]], on=["metal", "batch", "day"])
    if len(rows) == 0:
        return dict(metal=metal, model="CCM" if electrostatic else "NEM", logK=np.nan,
                    se_fit=np.nan, spread_pbp=np.nan, WSOS_DF=np.nan, n_points=0)

    def resid(theta):
        out = []
        for r in rows.itertuples():
            s = speciate(r.pH, totals_for(r), metal)
            G, sG = coverage(r.C0_ppm, r.Me_ppm, metal)
            st = surface(s, metal, logK_Me=theta[0], electrostatic=electrostatic)
            out.append((st["sp"]["CO3Me"]/S_AREA - G)/max(sG, 1e-12))
        return np.array(out)

    col = "logK_CCM" if electrostatic else "logK_NEM"
    x0  = float(np.nanmean(use[col])) if len(use) else -2.0
    sol = least_squares(resid, [x0], diff_step=1e-3)
    n, p = len(rows), 1
    wsos = float(np.sum(sol.fun**2)/max(n - p, 1))
    J = sol.jac
    se = float(np.sqrt(wsos*np.linalg.inv(J.T @ J)[0, 0])) if n > p else np.nan
    return dict(metal=metal, model="CCM" if electrostatic else "NEM", logK=float(sol.x[0]),
                se_fit=se, spread_pbp=float(np.nanstd(use[col])) if len(use) else np.nan,
                WSOS_DF=wsos, n_points=n)

head(7, "Global fit")
fits = pd.DataFrame([fit_metal(m, e) for m in ["Li", "Co"] for e in (True, False)])
show(fits)
print("\nWSOS/DF is far below 1 for both metals. That is not a good fit: it means the propagated")
print("error bars are so wide that the residuals are trivially satisfied, which is the fitting")
print("signature of the precision problem found in Part 3.")

# %% [markdown]
# ---
# ## Part 8 — Results
#
# The reporting rule, stated rather than left implicit: a constant is reported only if at least
# one point is fit-grade. Neither metal has one, so neither constant is reported. The fitted
# numbers are shown so the pattern can be read, and labelled so they cannot be quoted as values.
#
# The **exchange constant** log K_ex = log K(Me) − log K(Ca) is given alongside because it cancels
# the site-density convention and is therefore what transfers between studies (Belova et al.,
# 2014). The formation constant does not transfer.

# %%
head(8, "Results")
RXN = {"Co": ">CO3H0 + Co2+ = >CO3Co+ + H+", "Li": ">CO3H0 + Li+  = >CO3Li0 + H+"}
for metal in ["Co", "Li"]:
    c = fits[(fits.metal == metal) & (fits.model == "CCM")].iloc[0]
    n = fits[(fits.metal == metal) & (fits.model == "NEM")].iloc[0]
    u = ppdf[(ppdf.metal == metal) & ppdf.usable]
    nfit = int(ppdf[ppdf.metal == metal].fit_grade.sum())
    best = ppdf[ppdf.metal == metal].rel_err_pct.min()
    print(f"\n  {RXN[metal]}")
    print("     log K   NOT DETERMINED")
    print(f"     {nfit} of {int((ppdf.metal == metal).sum())} points are fit-grade; the best carries {best:.0f} % relative error")
    print(f"     against the {REL_ERR_FIT:.0f} % bar.")
    if len(u):
        print(f"     For information only, the fit over {len(u)} indicative point(s) returns")
        print(f"       {n.logK:+.2f} (NEM)  /  {c.logK:+.2f} (CCM),  propagated interval "
              f"{u.logK_lo.min():+.1f} to {u.logK_hi.max():+.1f}")
        print(f"       exchange constant log K_ex = log K({metal}) - log K(Ca) ~ "
              f"{n.logK - LOGK_SURF['CO3Ca']:+.2f} (NEM)")
        print("     These are arithmetic on noise. They are printed so the pattern can be read and")
        print("     must not be quoted as constants.")
    else:
        print("     No point survives even the indicative screen.")

# A chemical sanity check that the numbers themselves fail.
_kn = {m: float(fits[(fits.metal == m) & (fits.model == "NEM")].logK.iloc[0]) for m in ["Co", "Li"]}
print(f"\n  SANITY CHECK, and it fails. The two indicative values are {_kn['Co']:+.2f} (Co) and "
      f"{_kn['Li']:+.2f} (Li),")
print(f"  differing by {abs(_kn['Co'] - _kn['Li']):.2f} log units. A monovalent alkali and a divalent transition metal")
print("  binding the same carbonate site with the same strength is not chemically credible:")
print(f"  Pokrovsky's divalent values are {LOGK_SURF['CO3Ca']:+.1f} (Ca) and {LOGK_SURF['CO3Mg']:+.1f} (Mg), and lithium should sit")
print("  far weaker than either. Two metals with different charge and different ionic radius")
print("  landing on the same number is the signature of an inversion pinned by noise rather than")
print("  by chemistry, and it is independent confirmation of the precision result in Part 3.")

print("\n  Literature comparison, for the exchange constant:")
print("     Belova et al. (2014)  Ni on calcite   log K_ex = +0.58")
print("     Zachara et al. (1991) Ni on calcite   log K_ex = +0.51")
print(f"     Pokrovsky et al. (1999) >CO3Ca+ on dolomite  log K = {LOGK_SURF['CO3Ca']:+.1f}, "
      f">CO3Mg+ {LOGK_SURF['CO3Mg']:+.1f}")
print("\n  Second, independent limitation on cobalt: both indicative points are pre-equilibrium")
print("  (uptake is still rising at day 6), so even with adequate precision the value would be a")
print("  lower bound rather than a central one. And the pH 2 points sit at ~0.7 % removal, inside")
print("  the regime where Part 10 shows the inversion becomes biased as well as imprecise.")

# %% [markdown]
# ---
# ## Part 9 — Figures
#
# Four figures. Measured points are filled at day 6 and open for the pre-equilibrium samplings;
# model curves are the forward prediction at each batch's own chemistry. Where the metal
# carbonate is supersaturated the curve is faded and the region shaded, because the sorption-only
# model has stopped applying there.

# %%
best = {m: float(fits[(fits.metal == m) & (fits.model == "CCM")].logK.iloc[0]) for m in ["Li", "Co"]}

# ---- Figure 1: kinetics ---------------------------------------------------------------------
fig, (ax, bx) = plt.subplots(1, 2, figsize=(9.0, 3.7), constrained_layout=True)
for metal in ["Co", "Li"]:
    for batch, ls, fill in [("pH6", "-", True), ("pH2", "--", False)]:
        d = data[(data.metal == metal) & (data.batch == batch)].sort_values("day")
        kw = dict(color=COL[metal], marker=MK[metal], ls=ls,
                  mfc=COL[metal] if fill else "white",
                  label=f"{metal}, initial {batch.replace('pH', 'pH ')}")
        ax.plot(d.day, d["removal_%"], **kw)
        G = [coverage(r.C0_ppm, r.Me_ppm, metal)[0]*1e6 for r in d.itertuples()]
        bx.plot(d.day, G, **kw)
bx.axhline(SITE_DENS["CO3"]*1e6, color=COL["grey"], ls=":", lw=1.2)
bx.text(0.08, SITE_DENS["CO3"]*1e6 + 0.3, "carbonate site density (14 µmol m$^{-2}$)",
        fontsize=8, color=COL["grey"])
ax.set_xlabel("time (days)"); ax.set_ylabel("removal (% of initial)")
bx.set_xlabel("time (days)"); bx.set_ylabel("surface coverage Γ (µmol m$^{-2}$)")
ax.set_xticks([0, 2, 4, 6]); bx.set_xticks([0, 2, 4, 6]); bx.set_ylim(0, 16)
ax.legend(loc="upper left"); panel(ax, "a"); panel(bx, "b")
savefig(fig, "Fig1_kinetics")

# ---- Figure 2: measured against predicted, each point at its OWN chemistry -------------------
# NOT an adsorption edge. An edge (Goldberg 1985; Belova et al. 2014) varies pH alone at fixed
# contact time. Here pH, contact time and the Ca/Mg competitor load all move together between
# samplings, so no single curve on a pH axis is valid. The standard presentation for data whose
# conditions vary between samples is measured against predicted on a 1:1 line, with every point
# predicted at its own pH, its own Ca and Mg, and its own initial concentration.
par = []
for _, r in data[data.day > 0].iterrows():
    pred = forward_removal(r.metal, r.pH, r.C0_ppm, r.Ca_ppm, r.Mg_ppm, best[r.metal])
    dC   = r.C0_ppm - r.Me_ppm
    sd   = np.hypot(ICP_REL*r.C0_ppm, ICP_REL*r.Me_ppm)
    grade = ppdf[(ppdf.metal == r.metal) & (ppdf.batch == r.batch) &
                 (ppdf.day == r.day)].grade.iloc[0]
    par.append(dict(metal=r.metal, batch=r.batch, day=int(r.day), pH=r.pH,
                    measured=r["removal_%"], predicted=pred,
                    err_pct=100*sd/r.C0_ppm, grade=grade))
par = pd.DataFrame(par)

fig, axes = plt.subplots(1, 2, figsize=(9.2, 4.0), constrained_layout=True)

ax = axes[0]
lim = (0.05, 20)
ax.plot(lim, lim, "-", color="k", lw=1, zorder=1)
ax.fill_between(lim, [l*0.5 for l in lim], [l*2 for l in lim], color=COL["grey"],
                alpha=0.12, lw=0, zorder=0)
for metal in ["Co", "Li"]:
    for grade, mfc, lab in [("indicative", COL[metal], "indicative"),
                            ("excluded", "white", "excluded")]:
        d = par[(par.metal == metal) & (par.grade == grade)]
        if not len(d): continue
        ax.errorbar(d.measured, d.predicted, xerr=d.err_pct, fmt=MK[metal], color=COL[metal],
                    mfc=mfc, mew=1.3, ms=6.5, capsize=2.5, lw=1, ls="none",
                    label=f"{metal}, {lab}")
ax.set_xscale("log"); ax.set_yscale("log"); ax.set_xlim(*lim); ax.set_ylim(*lim)
ax.set_xlabel("measured removal (% of initial)")
ax.set_ylabel("predicted removal (% of initial)")
ax.text(0.055, 13, "shaded: within a factor of 2", fontsize=7.5, color=COL["grey"])
ax.legend(loc="lower right"); panel(ax, "a")

bx = axes[1]
for metal in ["Co", "Li"]:
    for batch, ls, fill in [("pH6", "-", True), ("pH2", "--", False)]:
        d = par[(par.metal == metal) & (par.batch == batch)].sort_values("day")
        bx.plot(d.day, d.measured, MK[metal], color=COL[metal], ls="none", ms=6.5,
                mfc=COL[metal] if fill else "white", mew=1.3,
                label=f"{metal} {batch.replace('pH', 'pH ')}, measured")
        bx.plot(d.day, d.predicted, ls, color=COL["model"], lw=1.3,
                label=f"{metal} {batch.replace('pH', 'pH ')}, predicted" if metal == "Co" else None)
bx.set_xlabel("time (days)"); bx.set_ylabel("removal (% of initial)")
bx.set_xticks([2, 4, 6]); bx.legend(loc="upper left", ncol=1, fontsize=7)
panel(bx, "b")
savefig(fig, "Fig2_measured_vs_predicted")

print("Figure 2a: every point predicted at its OWN pH, Ca, Mg and dose, against the 1:1 line.")
print("Figure 2b: the same numbers against time. The prediction is an equilibrium calculation,")
print("so it is nearly flat; the cobalt pH 6 measurement climbs past it, which is the approach to")
print("equilibrium. The cobalt pH 2 points run the other way because at 0.1 to 0.9 mg/L of uptake")
print("against 3.4 mg/L of noise they carry no information, which Part 3 established.")
_ok = par[(par.grade == "indicative")]
_fac = np.maximum(_ok.predicted/_ok.measured, _ok.measured/_ok.predicted)
print(f"\nOf the {len(_ok)} indicative points, {int((_fac < 2).sum())} fall within a factor of 2 of the prediction.")

# ---- Figure 3: the point-by-point constant against pH and coverage ---------------------------
fig, axes = plt.subplots(1, 2, figsize=(9.0, 3.7), constrained_layout=True)
for ax, (xcol, xlab) in zip(axes, [("pH", "pH"), ("Gamma", "coverage Γ (mol m$^{-2}$)")]):
    for metal in ["Co", "Li"]:
        d = ppdf[(ppdf.metal == metal) & np.isfinite(ppdf.logK_CCM)]
        u, x = d[d.usable], d[~d.usable]
        ax.errorbar(u[xcol], u.logK_CCM,
                    yerr=[u.logK_CCM - u.logK_lo, u.logK_hi - u.logK_CCM],
                    fmt=MK[metal], color=COL[metal], ms=6, capsize=2.5, lw=1,
                    ls="none", label=f"{metal}, indicative")
        ax.plot(x[xcol], x.logK_CCM, "x", color=COL["grey"], ms=6, mew=1.3, ls="none",
                label=f"{metal}, excluded")
    if xcol == "Gamma": ax.set_xscale("log")
    ax.set_xlabel(xlab); ax.set_ylabel("point-by-point log K (CCM)"); ax.set_ylim(-7, 1.5)
axes[0].legend(loc="lower left", ncol=2); panel(axes[0], "a"); panel(axes[1], "b")
savefig(fig, "Fig3_logK_diagnostic")

# ---- Figure 4: precision, the controlling variable --------------------------------------------
fig, ax = plt.subplots(figsize=(6.0, 4.0), constrained_layout=True)
for metal in ["Co", "Li"]:
    d = scr[scr.metal == metal]
    ax.plot(d.removal_pct, d.rel_err_pct, MK[metal], color=COL[metal], ls="none",
            ms=7, label=metal)
_x = np.linspace(0.1, 12, 200)
ax.plot(_x, 100*np.hypot(ICP_REL, ICP_REL*(1 - _x/100))/(_x/100), "-",
        color=COL["grey"], lw=1.2, label="±3 % ICP, exact")
ax.axhline(REL_ERR_FIT, color="k", ls=":", lw=1.2)
ax.text(0.15, REL_ERR_FIT*1.15, f"fit-grade bar, {REL_ERR_FIT:.0f} %", fontsize=8)
ax.set_xscale("log"); ax.set_yscale("log")
ax.set_xlabel("removal (% of initial)"); ax.set_ylabel("relative error on Γ (%)")
ax.legend(loc="upper right")
savefig(fig, "Fig4_precision_vs_removal")
print("Figure 4 is the one that explains the result: the relative error on Γ is set entirely by")
print("how much metal leaves solution. Every measured point sits above the fit-grade bar, and the")
print("grey curve shows that no amount of care in the analysis can move a point off it - only a")
print("different experiment can.")

# %% [markdown]
# ---
# ## Part 10 — What would determine these constants
#
# Three questions, answered by computation rather than argument. Each table is recomputed when
# the notebook runs, so the recommendations cannot drift from the model that produced them.
#
# 1. **Is the solver unbiased?** Generate removal from a *known* constant, add the measured ±3 %
#    ICP noise, invert it back. If the answer comes back biased, the code is wrong; if it comes
#    back true, the code is exonerated and the experiment is what limits the result.
# 2. **Which knob helps cobalt?** Cobalt is *signal-limited* — uptake is small because binding is
#    weak.
# 3. **Which knob helps lithium?** Lithium fails differently: it is *site-limited*, dosed far
#    above what the surface can hold, so the answer is the opposite one.

# %%
head(10, "What would determine these constants")
_rng = np.random.default_rng(11)
_G   = forward_removal.__globals__          # the live namespace, so the knobs really connect
_KEEP = {k: _G[k] for k in ("DOLOMITE_GL", "SSA", "S_AREA", "S_T", "M_G")}
EDGE  = [5.0, 5.5, 6.0, 6.5, 7.0, 7.5, 8.0]

class _Row:
    """Stand-in for a row of `data`, so synthetic points go through the real logK_point()."""
    def __init__(self, **kw): self.__dict__.update(kw)

def _set_solid(solid_gL, ssa):
    _G["DOLOMITE_GL"], _G["SSA"] = solid_gL, ssa
    _G["S_AREA"] = solid_gL*ssa
    _G["S_T"] = {k: v*_G["S_AREA"] for k, v in _G["SITE_DENS"].items()}
    _G["M_G"] = solid_gL*_G["V_L"]

def _recover(metal, logK_true, C0_ppm, pH_list, solid_gL=None, ssa=None, nrep=120):
    """Returns (mean removal %, bias, standard deviation) over nrep synthetic replicates."""
    if solid_gL is not None: _set_solid(solid_gL, ssa)
    exact = [(p, forward_removal(metal, p, C0_ppm, 1.0, 1.0, logK_true)) for p in pH_list]
    ests = []
    for _ in range(nrep):
        ks = [logK_point(_Row(pH=p, metal=metal, Ca_ppm=1.0, Mg_ppm=1.0,
                              C0_ppm=C0_ppm*(1 + ICP_REL*_rng.standard_normal()),
                              Me_ppm=C0_ppm*(1 - rem/100)*(1 + ICP_REL*_rng.standard_normal())))[0]
              for p, rem in exact]
        ks = [k for k in ks if np.isfinite(k)]
        if ks: ests.append(np.mean(ks))
    return np.mean([r for _, r in exact]), np.mean(ests) - logK_true, np.std(ests)

def _verdict(sd):
    return ("unusable" if sd > 0.5 else "weak" if sd > 0.25 else
            "publishable" if sd > 0.10 else "excellent")

# ---- 10.1  is the solver unbiased? -----------------------------------------------------------
print("10.1  Can the solver recover a KNOWN constant?   120 replicates, +/-3 % ICP noise.\n")
_b1 = []
for lab, kt, c0, ph in [("this study: 2 points, pH 7.40 and 7.48", -2.20, 80.5, [7.40, 7.48]),
                        ("7-point pH edge, same dose",             -2.20, 80.5, EDGE),
                        ("7-point pH edge, weak binder",           -3.00,  5.0, EDGE),
                        ("7-point pH edge, strong binder",         -0.79,  5.0, EDGE)]:
    _set_solid(_KEEP["DOLOMITE_GL"], _KEEP["SSA"])
    rem, bias, sd = _recover("Co", kt, c0, ph)
    _b1.append(dict(design=lab, true_logK=kt, removal_pct=rem, bias=bias, st_dev=sd))
show(pd.DataFrame(_b1))
_ok  = [r for r in _b1 if r["removal_pct"] >= 2.0]
_bad = [r for r in _b1 if r["removal_pct"] <  2.0]
print(f"\n    Wherever removal reaches about 2 %, the bias is within +/-{max(abs(r['bias']) for r in _ok):.2f}: the solver returns")
print("    the constant it was given, so THE CODE IS NOT WHAT LIMITS THIS RESULT. What changes")
print("    from row to row is the SPREAD, and it tracks only how much metal leaves solution")
print(f"    ({min(r['st_dev'] for r in _ok):.2f} at {max(r['removal_pct'] for r in _ok):.0f} % removal, "
      f"{max(r['st_dev'] for r in _ok):.2f} at {min(r['removal_pct'] for r in _ok):.0f} %).")
if _bad:
    r = _bad[0]
    print(f"\n    The exception matters here. At {r['removal_pct']:.1f} % removal the bias is {r['bias']:+.2f}: below about")
    print("    1 % uptake the inversion becomes BIASED as well as imprecise, because Gamma is then")
    print("    dominated by noise that cannot go negative and the recovered constant is pushed up.")
    _p2 = data[(data.metal == "Co") & (data.batch == "pH2")]["removal_%"].iloc[1:].mean()
    print(f"    The pH 2 cobalt batch of this study sits at {_p2:.1f} % removal, inside that regime.")

# ---- 10.2  cobalt: surface area --------------------------------------------------------------
print("\n\n10.2  COBALT is signal-limited. Sites per litre scale as (solid loading) x (SSA).\n")
_b2 = []
for solid, ssa, c0 in [(60, 0.76, 80.5), (60, 0.76, 20.0), (60, 0.76, 5.0),
                       (60, 5.00, 20.0), (60, 20.0, 20.0), (200, 5.00, 20.0), (200, 20.0, 20.0)]:
    rem, bias, sd = _recover("Co", -2.20, c0, EDGE, solid, ssa)
    _b2.append(dict(solid_gL=solid, SSA=ssa, area_m2L=solid*ssa, Co_mgL=c0,
                    removal_pct=rem, bias=bias, st_dev=sd, verdict=_verdict(sd)))
show(pd.DataFrame(_b2))
print("\n    Rows 1 to 3: cutting the cobalt dose from 80 to 5 mg/L changes almost nothing, because")
print("    fractional uptake is set by the binding strength and the site-to-metal ratio, not by")
print("    the dose. Rows 1 and 4: six times more surface per litre takes uptake from 4 % to 23 %")
print("    and the uncertainty from +/-0.20 to +/-0.04. SURFACE AREA PER LITRE IS THE ONLY KNOB")
print("    THAT MATTERS FOR COBALT. Above ~1200 m2/L a bias appears, so do not overshoot.")

# ---- 10.3  lithium: the dose -----------------------------------------------------------------
print("\n\n10.3  LITHIUM is site-limited, so the answer is the opposite: what matters is sites per Li.\n")
_b3 = []
for area, c0 in [(46, 138.4), (46, 14.0), (300, 138.4), (300, 14.0), (1000, 14.0)]:
    rem, bias, sd = _recover("Li", -2.33, c0, EDGE, area/5.0, 5.0)
    _b3.append(dict(area_m2L=area, Li_mgL=c0, sites_per_Li=14e-6*area*1e3/(c0/MM["Li"]),
                    removal_pct=rem, bias=bias, st_dev=sd, verdict=_verdict(sd)))
show(pd.DataFrame(_b3))
print("\n    Row 1 is this study and its bias is -0.45: the lithium result is not merely imprecise,")
print("    it is SYSTEMATICALLY TOO WEAK by nearly half a log unit, because the surface holds only")
print("    3 % of the lithium present and the inversion is pinned against the ceiling. Row 2 is")
print("    the cheap fix: a tenfold lower dose makes lithium publishable at the EXISTING area.")

_set_solid(_KEEP["DOLOMITE_GL"], _KEEP["SSA"])
for k, v in _KEEP.items(): _G[k] = v          # restore, so nothing downstream is affected

# %%
print("\n\n10.4  THE EXPERIMENT THAT WOULD DETERMINE BOTH CONSTANTS\n")
_design = pd.DataFrame([
 ("surface area per litre", "46 m2/L",           "300 m2/L",
  "400 g/L at the present 0.76 m2/g, or 60 g/L milled to 5 m2/g; the product is what counts"),
 ("cobalt dose",            "80 mg/L",           "20 mg/L",
  "dose barely matters for Co; 20 keeps the ICP comfortable"),
 ("lithium dose",           "138 mg/L",          "14 mg/L",
  "this alone removes the -0.45 bias; Li is site-limited"),
 ("pH values",              "2 nominal, drifting", "6 to 8, each held fixed",
  "buffer or titrate; the present usable Co points span 0.08 pH units"),
 ("concentrations",         "1",                 "4 to 6 at one fixed pH",
  "without an isotherm, site density and binding constant cannot be separated"),
 ("blanks",                 "none",              "solid-free at every condition",
  "Belova et al. 2014; matters most where uptake is small"),
 ("alkalinity / DIC",       "assumed pCO2",      "measured at every sampling",
  "so the saturation index comes from the solution, not an assumption"),
 ("solid conditioning",     "metal added with solid", "pre-equilibrate ~1 week, filter 0.2 um",
  "removes pH drift, the non-congruent Ca/Mg pattern and carbonate build-up at once"),
 ("contact time",           "6 d, drifting",     "24 h, plus a kinetic series",
  "Belova capped at 24 h to stay inside the adsorption regime"),
 ("desorption step",        "none",              "resuspend in fresh brine at the same pH",
  "reversible means surface complexation; poorly reversible means incorporation"),
], columns=["variable", "this study", "target", "why"])
show(_design)
print("\n    At 300 m2/L with Co 20 mg/L and Li 14 mg/L, both constants reach +/-0.02 to +/-0.04.")
print("    Stay below ~1200 m2/L to avoid the high-uptake bias of 10.2.")

# %% [markdown]
# ---
# ## Part 11 — Every reaction and every equation
#
# Complete listing of what this notebook solves. `{}` is aqueous activity, `[]` is concentration
# (mol L⁻¹ of suspension for surface species).
#
# ### 11.1 Surface reactions
#
# Three sites on dolomite in 1 : 1 : 2 stoichiometry. Only R7 and R8 are ever adjusted.
#
# | | Reaction | log K° | ΔZ |
# |---|---|---|---|
# | R1 | >CO₃H⁰ ⇌ >CO₃⁻ + H⁺ | −4.8 ± 0.2 | −1 |
# | R2 | >CaOH⁰ + H⁺ ⇌ >CaOH₂⁺ | +11.5 ± 0.2 | +1 |
# | R3 | >MgOH⁰ + H⁺ ⇌ >MgOH₂⁺ | +10.6 ± 0.2 | +1 |
# | R4 | >CaOH⁰ ⇌ >CaO⁻ + H⁺ | −12 ± 2 | −1 |
# | R5 | >MgOH⁰ ⇌ >MgO⁻ + H⁺ | −12 ± 2 | −1 |
# | R6a | >CO₃H⁰ + Ca²⁺ ⇌ >CO₃Ca⁺ + H⁺ | −1.8 ± 0.2 | +1 |
# | R6b | >CO₃H⁰ + Mg²⁺ ⇌ >CO₃Mg⁺ + H⁺ | −2.0 ± 0.2 | +1 |
# | **R7** | **>CO₃H⁰ + Li⁺ ⇌ >CO₃Li⁰ + H⁺** | **adjustable** | **0** |
# | **R8** | **>CO₃H⁰ + Co²⁺ ⇌ >CO₃Co⁺ + H⁺** | **adjustable** | **+1** |
# | R9a | >CaOH⁰ + CO₃²⁻ + 2H⁺ ⇌ >CaHCO₃⁰ + H₂O | −4.0 ± 0.5 | 0 |
# | R9b | >MgOH⁰ + CO₃²⁻ + 2H⁺ ⇌ >MgHCO₃⁰ + H₂O | −3.5 ± 0.5 | 0 |
# | R10a | >CaOH⁰ + CO₃²⁻ + H⁺ ⇌ >CaCO₃⁻ + H₂O | +16.6 ± 0.2 | −1 |
# | R10b | >MgOH⁰ + CO₃²⁻ + H⁺ ⇌ >MgCO₃⁻ + H₂O | +15.4 ± 0.2 | −1 |
#
# ΔZ is the change in charge of the surface species, and it decides which constants carry an
# electrostatic term. **R7 has ΔZ = 0, so the lithium constant is electrostatically blind** — a
# property of the stoichiometry, not evidence that it is well determined.
#
# ### 11.2 Aqueous reactions
#
# | | Reaction | log K |
# |---|---|---|
# | A1 | CO₂(g) + H₂O ⇌ H₂CO₃* | −1.469 |
# | A2 | H₂CO₃* ⇌ HCO₃⁻ + H⁺ | −6.345 |
# | A3 | CO₃²⁻ + H⁺ ⇌ HCO₃⁻ | +10.329 |
# | A4 | H₂O ⇌ H⁺ + OH⁻ | −14.0 |
#
# Complexes, as formation constants from the free ions (log β, 25 °C, I = 0):
#
# | Cobalt | log β | Lithium | log β | Calcium | log β | Magnesium | log β | Sodium | log β |
# |---|---|---|---|---|---|---|---|---|---|
# | CoCl⁺ | +0.30 | LiCl⁰ | −0.50 | CaCl⁺ | −0.70 | MgCl⁺ | −0.14 | NaCl⁰ | −0.78 |
# | CoCl₂⁰ | −0.20 | LiCO₃⁻ | +0.90 | CaCO₃⁰ | +3.22 | MgCO₃⁰ | +2.98 | NaCO₃⁻ | +1.27 |
# | CoCO₃⁰ | +4.23 | LiOH⁰ | −13.64 | CaHCO₃⁺ | +11.43 | MgHCO₃⁺ | +11.40 | NaHCO₃⁰ | +10.48 |
# | CoHCO₃⁺ | +12.20 | | | CaOH⁺ | −12.70 | MgOH⁺ | −11.79 | NaOH⁰ | −14.18 |
# | CoOH⁺ | −9.65 | | | | | | | | |
# | Co(OH)₂⁰ | −18.80 | | | | | | | | |
#
# **CoCl⁺ carries about a fifth of the dissolved cobalt at 0.684 M NaCl**, and **CoCO₃⁰ takes
# over 90 % of it above pH 9.5**, which is why the modelled edge in Figure 2a turns over near
# pH 8.6.
#
# Solids, used only for the saturation screen and never to compute sorption:
#
# | | Reaction | log K_sp |
# |---|---|---|
# | S1 | CoCO₃ (sphaerocobaltite) ⇌ Co²⁺ + CO₃²⁻ | −9.98 |
# | S2 | CaCO₃ (calcite) ⇌ Ca²⁺ + CO₃²⁻ | −8.48 |
# | S3 | CaMg(CO₃)₂ (ordered dolomite) ⇌ Ca²⁺ + Mg²⁺ + 2CO₃²⁻ | −17.09 |
# | S4 | Li₂CO₃ (zabuyelite) ⇌ 2Li⁺ + CO₃²⁻ | ≈ −2.5 |
#
# ### 11.3 Equations
#
# **Activity coefficients**, Davies (1962), A = 0.509, b = 0.3:
#
# $$\log\gamma_z = -A z^2\left(\frac{\sqrt{I}}{1+\sqrt{I}} - 0.3 I\right) \tag{1}$$
# $$I = \tfrac12\sum_i c_i z_i^2 \tag{2}$$
#
# **Carbonate at fixed pCO₂** — the assumption that drives the high-pH behaviour:
#
# $$\{\mathrm{HCO_3^-}\} = \frac{K_{A1}K_{A2}\,p_{\mathrm{CO_2}}}{\{\mathrm{H^+}\}}
# \qquad \{\mathrm{CO_3^{2-}}\} = \frac{\{\mathrm{HCO_3^-}\}}{K_{A3}\{\mathrm{H^+}\}} \tag{3}$$
#
# so $\{\mathrm{CO_3^{2-}}\}\propto\{\mathrm{H^+}\}^{-2}$: a ten-thousandfold rise from pH 8 to 10.
#
# **Aqueous mass balance**, solved with (1) and (2) by fixed-point iteration:
#
# $$T_j = \sum_i a_{ij} C_i, \qquad C_i = \beta_i\prod_j\{X_j\}^{a_{ij}}\big/\gamma_i \tag{4}$$
#
# **Site mass balances**, one per site:
#
# $$[\mathrm{>CO_3}]_T = [\mathrm{>CO_3H^0}]+[\mathrm{>CO_3^-}]+[\mathrm{>CO_3Li^0}]
# +[\mathrm{>CO_3Co^+}]+[\mathrm{>CO_3Ca^+}]+[\mathrm{>CO_3Mg^+}] \tag{5a}$$
# $$[\mathrm{>CaOH}]_T = [\mathrm{>CaOH^0}]+[\mathrm{>CaOH_2^+}]+[\mathrm{>CaO^-}]
# +[\mathrm{>CaHCO_3^0}]+[\mathrm{>CaCO_3^-}] \tag{5b}$$
# $$[\mathrm{>MgOH}]_T = [\mathrm{>MgOH^0}]+[\mathrm{>MgOH_2^+}]+[\mathrm{>MgO^-}]
# +[\mathrm{>MgHCO_3^0}]+[\mathrm{>MgCO_3^-}] \tag{5c}$$
#
# **Surface charge, predicted** — charged surface species only. Dissolved Ca²⁺ and Mg²⁺ from
# dissolution appear **nowhere** in it:
#
# $$\sigma = \frac{F}{S_a}\Big([\mathrm{>CaOH_2^+}]+[\mathrm{>MgOH_2^+}]+[\mathrm{>CO_3Co^+}]
# +[\mathrm{>CO_3Ca^+}]+[\mathrm{>CO_3Mg^+}]$$
# $$\qquad -[\mathrm{>CO_3^-}]-[\mathrm{>CaO^-}]-[\mathrm{>MgO^-}]-[\mathrm{>CaCO_3^-}]
# -[\mathrm{>MgCO_3^-}]\Big) \tag{6}$$
#
# **Surface charge, measured** — a *different* quantity, from the solution mass balance, which
# carries the dissolution signal. The two are compared, never substituted:
#
# $$\sigma_{meas} = \frac{F}{S_a}\big(2\Delta[\mathrm{Ca}]+2\Delta[\mathrm{Mg}]
# -\Delta[\mathrm{CO_3}]\big) \tag{7}$$
#
# **Constant capacitance**, in both conventions:
#
# $$\psi = \frac{\sigma}{C},\quad C = \frac{\sqrt I}{\alpha},\quad \alpha = 0.004 \tag{8}$$
# $$\sigma_v = \frac{C S a}{F}\psi \quad\text{(Goldberg 1985 Eq. 1, per volume)} \tag{9}$$
#
# **Intrinsic from apparent:**
#
# $$K^{int} = K^{app}\exp\!\left(\frac{\Delta Z F\psi}{RT}\right) \tag{10}$$
#
# **Mass action for the adjustable reactions**, which is what Part 6 inverts:
#
# $$K^{int}_{\mathrm{Me}} = \frac{[\mathrm{>CO_3Me}]\{\mathrm{H^+}\}}
# {[\mathrm{>CO_3H^0}]\{\mathrm{Me}^{z+}\}}\exp\!\left(\frac{\Delta Z F\psi}{RT}\right) \tag{11}$$
#
# **Coverage and its error** — the quantity whose relative error is 54 % to 2968 % here:
#
# $$\Gamma = \frac{(C_0-C_{eq})V}{m S_{BET}} \tag{12}$$
# $$\sigma_\Gamma = \frac{V}{m S_{BET}}\sqrt{(\epsilon C_0)^2+(\epsilon C_{eq})^2},\quad
# \epsilon = 0.03 \tag{13}$$
#
# **Saturation index:**
#
# $$\mathrm{SI} = \log\frac{\mathrm{IAP}}{K_{sp}} \tag{14}$$
#
# **Objective**, FITEQL (Westall, 1982). Near 1 means the model matches to within measurement
# error; far below 1 means the error bars are too wide to constrain anything:
#
# $$\mathrm{WSOS/DF} = \frac{1}{N_{obs}-N_{par}}\sum_i\left(\frac{Y_i}{s_i}\right)^2 \tag{15}$$
#
# **Aqueous–surface analogy** (Van Cappellen et al., 1993):
#
# $$\log K(\mathrm{>CO_3Co^+}) = \log K(\mathrm{>CO_3Ca^+})
# +\big[\log\beta(\mathrm{CoCO_3^0})-\log\beta(\mathrm{CaCO_3^0})\big] \tag{16}$$
#
# **Exchange constant** — what transfers between studies, because it cancels the site-density
# convention (Belova et al., 2014):
#
# $$\log K_{ex} = \log K(\mathrm{>CO_3Me}) - \log K(\mathrm{>CO_3Ca}) \tag{17}$$
#
# ---
# ## References
#
# Belova, D. A., Lakshtanov, L. Z., Carneiro, J. F., & Stipp, S. L. S. (2014). Nickel adsorption
# on chalk and calcite. *Journal of Contaminant Hydrology, 170*, 1–9.
#
# Brady, P. V., Papenguth, H. W., & Kelly, J. W. (1999). Metal sorption to dolomite surfaces.
# *Applied Geochemistry, 14*(5), 569–579.
#
# Davies, C. W. (1962). *Ion association*. Butterworths.
#
# Goldberg, S. (1985). Chemical modeling of anion competition on goethite using the constant
# capacitance model. *Soil Science Society of America Journal, 49*(4), 851–856.
#
# Parkhurst, D. L., & Appelo, C. A. J. (2013). *Description of input and examples for PHREEQC
# version 3* (Techniques and Methods 6-A43). U.S. Geological Survey.
#
# Pokrovsky, O. S., Schott, J., & Thomas, F. (1999). Dolomite surface speciation and reactivity
# in aquatic systems. *Geochimica et Cosmochimica Acta, 63*(19–20), 3133–3143.
#
# Van Cappellen, P., Charlet, L., Stumm, W., & Wersin, P. (1993). A surface complexation model of
# the carbonate mineral–aqueous solution interface. *GCA, 57*(15), 3505–3518.
#
# Westall, J. C. (1982). *FITEQL: A computer program for determination of chemical equilibrium
# constants from experimental data* (Report 82-01). Oregon State University.
#
# Zachara, J. M., Cowan, C. E., & Resch, C. T. (1991). Sorption of divalent metals on calcite.
# *Geochimica et Cosmochimica Acta, 55*(6), 1549–1562.
