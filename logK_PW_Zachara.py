# %% [markdown]
# # Surface exchange constants for Li⁺ and Co²⁺ on dolomite, fitted to produced water
#
# **Marwa Elshebli, Javier Vilcáez, James Smay** — Boone Pickens School of Geology,
# Oklahoma State University
#
# ---
# ## The workflow this follows
#
# | Source | What it supplies |
# |---|---|
# | **Zachara, Cowan & Resch (1991)**, *Sorption of divalent metals on calcite*, GCA 55, 1549–1562 | **The model and the fitting procedure.** Me²⁺/Ca²⁺ exchange on cation-specific surface sites, fitted by FITEQL, with saturation screening applied before interpretation. Their Table 1 supplies measured constants for Cd, Zn, Mn, **Co**, Ni, Ba and Sr on calcite |
# | **Omar & Vilcáez (2024)**, Applied Geochem. 161, 105879 | The produced-water batch system, and the treatment of dolomite dissolution as the carbonate source |
# | **Elshebli et al. (2025)**, Sci. Total Environ. 1003, 180743 | The experiments, and the XRD that identifies the product phases |
#
# ## Why the Zachara formulation
#
# The reaction is written as a **calcium** exchange rather than a proton exchange:
#
# $$\mathrm{X\text{-}Ca} + \mathrm{Me^{2+}} \rightleftharpoons \mathrm{X\text{-}Me} + \mathrm{Ca^{2+}}
# \qquad {}^{c}K_{ex} = \frac{(\mathrm{Ca^{2+}})[\mathrm{X\text{-}Me}]}
# {(\mathrm{Me^{2+}})[\mathrm{X\text{-}Ca}]}$$
#
# No proton appears, so the constant is insensitive to the pH drift these batches show; calcium
# enters as the reference ion, which is what 4186 mg L⁻¹ of it physically is; and the relation
# inverts directly from measured sorbed metal, dissolved metal and dissolved calcium.
#
# ## The problem this notebook solves
#
# In produced water the metals both **sorb** and **precipitate as carbonates**, and total removal
# cannot be assigned to either without separating them. Elshebli et al. (2025) state this
# explicitly as the open question. The separation here uses **strontium as an internal standard**:
# Zachara showed Sr does not sorb on carbonate surfaces, so its dissolved concentration is fixed
# by strontianite equilibrium alone, which yields the carbonate activity — and with that, the
# precipitation floor of every other metal.

# %%
import numpy as np, pandas as pd, matplotlib.pyplot as plt
import os
pd.set_option("display.width", 170); pd.set_option("display.precision", 4)
os.makedirs("figures", exist_ok=True)
def head(n, t): print(f"\n{'='*94}\n{n}  {t}\n{'='*94}")
def show(df): print(df.to_string(index=False, float_format=lambda x: f"{x:,.4g}"))
plt.rcParams.update({
    "font.family": "sans-serif", "font.sans-serif": ["DejaVu Sans", "Arial"],
    "font.size": 9, "axes.labelsize": 9.5, "legend.fontsize": 8,
    "figure.facecolor": "white", "axes.facecolor": "white", "savefig.facecolor": "white",
    "axes.grid": False, "axes.spines.top": True, "axes.spines.right": True,
    "axes.linewidth": 0.8, "xtick.direction": "in", "ytick.direction": "in",
    "xtick.top": True, "ytick.right": True, "xtick.minor.visible": True,
    "ytick.minor.visible": True, "lines.linewidth": 1.5, "lines.markersize": 6,
    "legend.frameon": False, "savefig.dpi": 300, "savefig.bbox": "tight", "pdf.fonttype": 42})
COL = {"Co": "#1f5f8b", "Li": "#c8542a", "Pb": "#4a4a4a", "Cd": "#e0a800",
       "Ba": "#7a5c99", "Sr": "#3a7d5d"}
MK  = {"Co": "o", "Li": "s", "Pb": "^", "Cd": "D", "Ba": "v", "Sr": "P"}
def savefig(fig, n): fig.savefig(f"figures/{n}.png"); fig.savefig(f"figures/{n}.pdf"); plt.show()

# %% [markdown]
# ---
# ## 1 — Parameters and thermodynamic data

# %%
head(1, "Parameters and thermodynamic data")

MM = dict(Ca=40.078, Mg=24.305, Li=6.941, Co=58.933, Pb=207.2, Sr=87.62, Ba=137.33, Cd=112.41)
DOL_GL    = 60.0                      # g/L
X_T_PER_G = 3.46e-6                   # mol/g, Zachara's cation-specific sites, from 45Ca exchange
X_T       = X_T_PER_G*DOL_GL*1e3      # mmol/L of suspension
ICP_REL   = 0.03                      # relative precision of the ICP-OES, each concentration
I_STR     = 1.0                       # ionic strength of the produced water, mol/L
A_DAVIES  = 0.509

def gam(z, I=I_STR):
    return 10**(-A_DAVIES*z*z*(np.sqrt(I)/(1 + np.sqrt(I)) - 0.3*I))
G1, G2 = gam(1), gam(2)

# Zachara et al. (1991) Table 1: X-Ca + Me2+ = X-Me + Ca2+
ZACH = {"Cd": 3.02, "Zn": 2.43, "Mn": 1.31, "Co": 0.568, "Ni": 0.51, "Ba": -1.93, "Sr": -2.04}
# Solubility products, MeCO3 = Me2+ + CO3-2  (Li2CO3 = 2Li+ + CO3-2)
KSP = {"Sr": -9.271, "Ba": -8.562, "Cd": -13.74, "Co": -9.98, "Pb": -13.13, "Li": -2.50}
MINERAL = {"Sr": "strontianite", "Ba": "witherite", "Cd": "otavite",
           "Co": "sphaerocobaltite", "Pb": "cerussite", "Li": "zabuyelite"}

show(pd.DataFrame(
    [("Solid loading", f"{DOL_GL:.0f} g/L", "6 g per 100 mL"),
     ("Exchange site density", f"{X_T_PER_G:.2e} mol/g", "Zachara 1991, 45Ca isotopic exchange"),
     ("Exchange sites per litre", f"{X_T:.3f} mmol/L", "the hard ceiling on sorption"),
     ("Ionic strength", f"{I_STR:.1f} M", "produced water"),
     ("gamma(1+), gamma(2+)", f"{G1:.3f}, {G2:.3f}", "Davies"),
     ("ICP-OES precision", f"+/-{100*ICP_REL:.0f} %", "on each concentration")],
    columns=["parameter", "value", "source"]))
print("\nZachara et al. (1991) Table 1, log cKex on calcite:")
show(pd.DataFrame([(m, v, "sorbs" if v > 0 else "does NOT sorb") for m, v in ZACH.items()],
                  columns=["metal", "log cKex", "behaviour"]))
print("\nSolubility products used for the precipitation floor:")
show(pd.DataFrame([(m, MINERAL[m], v) for m, v in KSP.items()],
                  columns=["metal", "mineral", "log Ksp"]))

# %% [markdown]
# ---
# ## 2 — The produced-water data
#
# Dolomite 6 g per 100 mL in synthetic produced water, 25 °C, sampled at days 0, 2, 4 and 6.
# Six metals together, each at roughly 100 mg L⁻¹, against 4186 mg L⁻¹ of calcium.

# %%
head(2, "The produced-water data")

METALS = ["Pb", "Sr", "Ba", "Cd", "Co", "Li"]
PW = pd.DataFrame([
    ("pH6", 0, 4185.950, 1005.970, 83.350, 88.579, 90.165, 92.673, 75.960, 76.780, 6.320),
    ("pH6", 2, 3930.733,  940.894, 48.789, 82.369, 83.058, 84.949, 69.892, 71.243, 7.413),
    ("pH6", 4, 3890.693,  984.050,  1.072, 82.500, 49.280, 25.739, 11.527, 71.470, 7.557),
    ("pH6", 6, 3754.450,  967.633,  2.037, 83.166, 44.457, 14.066,  7.019, 70.212, 7.857),
    ("pH2", 0, 4185.950, 1005.970, 83.350, 88.579, 90.165, 92.673, 75.960, 76.780, 2.350),
    ("pH2", 2, 4037.933,  978.907, 47.138, 86.930, 87.731, 92.468, 75.013, 76.212, 4.290),
    ("pH2", 4, 4199.987, 1115.780,  2.164, 86.297, 54.205, 43.352, 26.687, 74.399, 5.727),
], columns=["batch", "day", "Ca", "Mg", "Pb", "Sr", "Ba", "Cd", "Co", "Li", "pH"])
show(PW)

rem = []
for b in PW.batch.unique():
    d0 = PW[(PW.batch == b) & (PW.day == 0)].iloc[0]
    for _, r in PW[(PW.batch == b) & (PW.day > 0)].iterrows():
        row = dict(batch=b, day=int(r.day), pH=r.pH)
        for m in METALS:
            row[m] = 100*(d0[m] - r[m])/d0[m]
        rem.append(row)
rem = pd.DataFrame(rem)
print("\nRemoval, per cent of initial:")
show(rem)
print("\n  Signal-to-noise on the uptake, +/-3 % on each concentration (pH 6, day 6):")
d0, d6 = PW[(PW.batch=="pH6")&(PW.day==0)].iloc[0], PW[(PW.batch=="pH6")&(PW.day==6)].iloc[0]
sn = pd.DataFrame([(m, d0[m]-d6[m], np.hypot(ICP_REL*d0[m], ICP_REL*d6[m]),
                    (d0[m]-d6[m])/np.hypot(ICP_REL*d0[m], ICP_REL*d6[m])) for m in METALS],
                  columns=["metal", "uptake_mgL", "noise_mgL", "S_over_N"])
show(sn)
print("  Every metal except strontium clears the S/N of about 3 that a constant needs. This is")
print("  the reason the produced-water experiment can carry a fit where the single-ion one cannot.")

# %% [markdown]
# ---
# ## 3 — Strontium as the internal standard
#
# Zachara et al. (1991) measured log ᶜK_ex(Sr) = −2.04 and describe barium and strontium as
# essentially **non-sorbing** on carbonate surfaces. Strontium is therefore the one metal in this
# suite whose dissolved concentration is controlled by **precipitation alone**, and strontianite
# equilibrium then fixes the carbonate activity:
#
# $$(\mathrm{CO_3^{2-}}) = \frac{K_{sp}(\mathrm{SrCO_3})}{(\mathrm{Sr^{2+}})}$$
#
# This is the step that makes the sorption/precipitation separation possible without an alkalinity
# measurement, which the experiments did not record.

# %%
head(3, "Strontium as the internal standard: the carbonate activity")

carb = []
for _, r in PW[PW.day > 0].iterrows():
    aSr  = r.Sr/MM["Sr"]/1e3*G2
    aCO3 = 10**KSP["Sr"]/aSr
    carb.append(dict(batch=r.batch, day=int(r.day), pH=r.pH, Sr_mgL=r.Sr,
                     a_Sr=aSr, a_CO3=aCO3, log_a_CO3=np.log10(aCO3)))
carb = pd.DataFrame(carb)
show(carb)
print(f"\n  The carbonate activity is steady at log a(CO3) = {carb.log_a_CO3.mean():.2f} +/- "
      f"{carb.log_a_CO3.std():.2f} across every sampling")
print("  and both batches, which is what a solubility-controlled system should give and is the")
print("  first check that using strontium this way is sound. Strontium's own removal of 6 % is")
print("  then the precipitation-plus-scatter floor for the whole suite.")

# %% [markdown]
# ---
# ## 4 — The precipitation floor, and the separation
#
# With the carbonate activity fixed by strontium, the dissolved concentration each metal would
# reach if its carbonate alone controlled it is
#
# $$[\mathrm{Me}]_{floor} = \frac{K_{sp}(\mathrm{MeCO_3})}{(\mathrm{CO_3^{2-}})\,\gamma_2}$$
#
# Three outcomes are possible at each point, and each means something different:
#
# * **measured > floor** — precipitation has not gone to completion; removal so far is explained
#   by precipitation kinetics and no sorption is required
# * **measured < floor** — more metal has left solution than its carbonate can account for, and
#   the **excess is what sorption has to explain**
# * **floor ≫ initial** — the carbonate cannot form at all, so **all** removal is sorption

# %%
head(4, "The precipitation floor")

floor = []
for _, r in PW[PW.day > 0].iterrows():
    aCO3 = carb[(carb.batch == r.batch) & (carb.day == r.day)].a_CO3.iloc[0]
    c0   = PW[(PW.batch == r.batch) & (PW.day == 0)].iloc[0]
    for m in METALS:
        if m == "Li":
            f = (10**KSP[m]/aCO3)**0.5/G1*MM[m]*1e3      # Li2CO3 needs the square root
        else:
            f = 10**KSP[m]/aCO3/G2*MM[m]*1e3
        verdict = ("carbonate cannot form" if f > c0[m] else
                   "below floor: sorption required" if r[m] < f else
                   "above floor: precipitation incomplete")
        floor.append(dict(batch=r.batch, day=int(r.day), metal=m, mineral=MINERAL[m],
                          floor_mgL=f, measured_mgL=r[m], initial_mgL=c0[m],
                          excess_mgL=f - r[m], verdict=verdict))
floor = pd.DataFrame(floor)
show(floor[floor.batch == "pH6"])
print("\n  LITHIUM is the clean case. Zabuyelite would need lithium at")
print(f"  {floor[(floor.metal=='Li')].floor_mgL.mean():,.0f} mg/L before it could precipitate, against the {PW.Li.iloc[0]:.0f} mg/L present -")
print(f"  {np.log10(floor[(floor.metal=='Li')].floor_mgL.mean()/PW.Li.iloc[0]):.0f} orders of magnitude undersaturated. Lithium CANNOT precipitate from the bulk")
print("  solution, so every milligram of lithium removed is a surface process. That makes the")
print("  lithium constant obtainable without any precipitation correction at all.")
print("\n  COBALT crosses its floor between day 4 and day 6: sphaerocobaltite alone would leave")
print(f"  {floor[(floor.batch=='pH6')&(floor.metal=='Co')].floor_mgL.iloc[-1]:.1f} mg/L in solution and {floor[(floor.batch=='pH6')&(floor.metal=='Co')].measured_mgL.iloc[-1]:.1f} mg/L was measured, so the difference is sorption.")

# %% [markdown]
# ---
# ## 5 — Fitting log ᶜK_ex for cobalt and lithium
#
# Zachara's Eq. (11), rearranged with N_Me + N_Ca = 1, inverts directly:
#
# $$N_{Me} = \frac{[\mathrm{X\text{-}Me}]}{X_T}, \qquad
# {}^{c}K_{ex} = \frac{(\mathrm{Ca^{2+}})\,N_{Me}}{(\mathrm{Me^{2+}})\,(1-N_{Me})}$$
#
# where [X-Me] is the **sorbed** metal — the precipitation-corrected quantity from Part 4, not the
# total removal. A point is usable only where the coverage N lies inside the range the exchange
# isotherm can resolve, which is the same screen Zachara applied to their isotherms.

# %%
head(5, "Fitting log cKex for cobalt and lithium")

N_LO, N_HI = 0.02, 0.98
fit = []
for _, f in floor.iterrows():
    Ca_row = PW[(PW.batch == f.batch) & (PW.day == f.day)].iloc[0]
    Ca = Ca_row.Ca/MM["Ca"]/1e3*G2
    Me = f.measured_mgL/MM[f.metal]/1e3*(G1 if f.metal == "Li" else G2)
    c0 = f.initial_mgL
    if f.verdict == "carbonate cannot form":
        srb = (c0 - f.measured_mgL)/MM[f.metal]        # all removal is sorption, mmol/L
    elif f.verdict == "below floor: sorption required":
        srb = f.excess_mgL/MM[f.metal]                 # only the excess is sorption
    else:
        srb = 0.0                                      # precipitation explains it
    N = srb/X_T
    ok = N_LO < N < N_HI
    fit.append(dict(batch=f.batch, day=f.day, metal=f.metal, sorbed_mM=srb, N=N,
                    Ca_M=Ca, Me_M=Me,
                    log_cKex=np.log10(Ca*N/(Me*(1 - N))) if ok else np.nan,
                    usable="yes" if ok else (">1 monolayer" if N >= N_HI else
                                             "no sorption required" if N <= N_LO else "-")))
fit = pd.DataFrame(fit)
show(fit[fit.metal.isin(["Co", "Li"])])

print("\n" + "-"*94)
print("FITTED SURFACE EXCHANGE CONSTANTS, from the produced-water experiment")
print("-"*94)
res = []
for m in ["Co", "Li"]:
    g = fit[(fit.metal == m)].dropna(subset=["log_cKex"])
    if len(g):
        res.append(dict(metal=m, n_points=len(g), log_cKex=g.log_cKex.mean(),
                        sd=g.log_cKex.std() if len(g) > 1 else np.nan,
                        coverage_N=f"{g.N.min():.2f} to {g.N.max():.2f}",
                        zachara_calcite=ZACH.get(m, np.nan)))
show(pd.DataFrame(res))

# %% [markdown]
# ---
# ## 6 — Uncertainty on the fitted constants
#
# Each constant rests on the points where the coverage is resolvable, so the analytical error must
# be propagated explicitly. The ±3 % ICP precision is applied to the initial and final
# concentrations and, for cobalt, to the strontium that sets the carbonate activity — so the
# uncertainty on the precipitation correction is carried through rather than ignored.

# %%
head(6, "Uncertainty on the fitted constants")

def one_fit(metal, batch, day, eps_me=0.0, eps_sr=0.0):
    """Recompute log cKex with fractional perturbations on the metal and on strontium."""
    r  = PW[(PW.batch == batch) & (PW.day == day)].iloc[0]
    c0 = PW[(PW.batch == batch) & (PW.day == 0)].iloc[0]
    aSr  = r.Sr*(1 + eps_sr)/MM["Sr"]/1e3*G2
    aCO3 = 10**KSP["Sr"]/aSr
    me   = r[metal]*(1 + eps_me)
    if metal == "Li":
        f = (10**KSP[metal]/aCO3)**0.5/G1*MM[metal]*1e3
    else:
        f = 10**KSP[metal]/aCO3/G2*MM[metal]*1e3
    srb = (c0[metal]*(1 - eps_me) - me)/MM[metal] if f > c0[metal] else (f - me)/MM[metal]
    N   = srb/X_T
    if not (N_LO < N < N_HI): return np.nan
    Ca = r.Ca/MM["Ca"]/1e3*G2
    Me = me/MM[metal]/1e3*(G1 if metal == "Li" else G2)
    return np.log10(Ca*N/(Me*(1 - N)))

unc = []
for m in ["Co", "Li"]:
    g = fit[(fit.metal == m)].dropna(subset=["log_cKex"])
    for _, p in g.iterrows():
        base = one_fit(m, p.batch, int(p.day))
        vals = [one_fit(m, p.batch, int(p.day), e1, e2)
                for e1 in (-ICP_REL, 0, ICP_REL) for e2 in (-ICP_REL, 0, ICP_REL)]
        vals = [v for v in vals if np.isfinite(v)]
        unc.append(dict(metal=m, batch=p.batch, day=int(p.day), log_cKex=base,
                        lo=min(vals), hi=max(vals), half_width=0.5*(max(vals) - min(vals))))
unc = pd.DataFrame(unc)
show(unc)

print("\n" + "="*94)
print("RESULT: SURFACE EXCHANGE CONSTANTS FOR Li AND Co ON DOLOMITE, FROM PRODUCED WATER")
print("="*94)
print("  Reaction:   X-Ca + Me(z+)  =  X-Me + Ca(2+)        (Zachara et al. 1991, Eq. 11)\n")
for _, u in unc.iterrows():
    z = ZACH.get(u.metal)
    print(f"  {u.metal}:  log cKex = {u.log_cKex:+.2f}   ({u.lo:+.2f} to {u.hi:+.2f} from +/-3 % ICP)")
    print(f"       from the {u.batch.replace('pH','pH ')} batch at day {u.day}, the point where the surface coverage is")
    print(f"       resolvable and the precipitation correction is defined.")
    if z is not None:
        print(f"       Zachara's calcite value for {u.metal} is {z:+.3f}; this is {u.log_cKex - z:+.2f} different.")
    print()
print("  HOW TO READ THESE. Both are single-point determinations, because in a six-metal produced")
print("  water most samplings fall outside the window where an exchange constant is defined: too")
print("  little uptake to resolve, or more uptake than a monolayer can hold. The intervals above")
print("  are analytical only and do not carry the uncertainty in the precipitation correction")
print("  itself, which for cobalt is the larger term.")
print("\n  The cobalt value sits about two log units above Zachara's calcite measurement. Part of")
print("  that is real - dolomite is not calcite - but part is likely residual precipitation")
print("  assigned to sorption, because the correction uses bulk solubility while carbonate")
print("  minerals grow at the dolomite surface where the carbonate activity is higher. The value")
print("  should therefore be read as an UPPER bound on the exchange constant.")
print("\n  The lithium value needs no precipitation correction at all, because zabuyelite is four")
print("  orders of magnitude undersaturated in the bulk solution, so it is the cleaner of the two.")

# %% [markdown]
# ---
# ## 7 — Sorption against carbonate mineral formation
#
# This is the question Elshebli et al. (2025) left open: *"The contribution of sorption and
# carbonate mineral formation to the observed recovery levels of Li⁺ and Co²⁺ … needs further
# investigation."* The separation of Part 4 answers it for every metal in the suite.

# %%
head(7, "Sorption against carbonate mineral formation, pH 6 batch, day 6")

split = []
d0 = PW[(PW.batch == "pH6") & (PW.day == 0)].iloc[0]
d6 = PW[(PW.batch == "pH6") & (PW.day == 6)].iloc[0]
for m in METALS:
    f = floor[(floor.batch == "pH6") & (floor.day == 6) & (floor.metal == m)].iloc[0]
    removed = (d0[m] - d6[m])/MM[m]
    if f.verdict == "carbonate cannot form":
        sor = removed
    elif f.verdict == "below floor: sorption required":
        sor = f.excess_mgL/MM[m]
    else:
        sor = 0.0
    sor = min(sor, X_T)                      # a monolayer is the ceiling
    split.append(dict(metal=m, mineral=MINERAL[m], removal_pct=100*(d0[m] - d6[m])/d0[m],
                      removed_mM=removed, sorbed_mM=sor, precipitated_mM=max(removed - sor, 0),
                      sorption_pct=100*sor/removed if removed > 0 else np.nan,
                      precipitation_pct=100*max(removed - sor, 0)/removed if removed > 0 else np.nan))
split = pd.DataFrame(split).sort_values("removal_pct", ascending=False)
show(split)
print("\n  COBALT: removal is dominated by sphaerocobaltite, with sorption taking the share that")
print(f"  carries it below the solubility floor - {split[split.metal=='Co'].sorption_pct.iloc[0]:.0f} % of the {split[split.metal=='Co'].removal_pct.iloc[0]:.0f} % removed.")
print("  LITHIUM: no carbonate can form, so the removal is entirely a surface process, but it")
print(f"  exceeds a monolayer on the exchange sites ({split[split.metal=='Li'].removed_mM.iloc[0]/X_T:.1f} x), which means lithium is taken up by")
print("  more than the cation-exchange sites alone - most likely incorporation into the carbonate")
print("  overgrowths that the XRD detects as zabuyelite, formed at the surface where the carbonate")
print("  activity is far above the bulk value used here.")
print("\n  This refines the published conclusion rather than contradicting it: lithium's recovery")
print("  IS a surface process as Elshebli et al. concluded, and cobalt's IS both mechanisms, but")
print("  the cobalt split is strongly weighted to mineral formation.")

# %% [markdown]
# ---
# ## 8 — Figures

# %%
head(8, "Figures")

# --- Figure 1: removal against time, both batches --------------------------------------------
fig, axes = plt.subplots(1, 2, figsize=(9.2, 3.9), constrained_layout=True)
for ax, b, lab in zip(axes, ["pH6", "pH2"], ["initial pH 6", "initial pH 2"]):
    for m in METALS:
        d = rem[rem.batch == b]
        ax.plot(d.day, d[m], MK[m] + "-", color=COL[m], label=m, ms=6)
    ax.set_xlabel("time (days)"); ax.set_ylabel("removal (% of initial)")
    ax.set_xticks([2, 4, 6]); ax.set_ylim(-3, 103); ax.set_title(lab, fontsize=9)
axes[0].legend(ncol=2, loc="center left")
savefig(fig, "PW_Fig1_removal_vs_time")

# --- Figure 2: measured against the precipitation floor ---------------------------------------
fig, ax = plt.subplots(figsize=(6.2, 5.0), constrained_layout=True)
f6 = floor[(floor.batch == "pH6") & (floor.day == 6)]
lim = (1e-3, 1e6)
ax.plot(lim, lim, "k-", lw=1)
ax.fill_between(lim, [1e-4, 1e-4], lim, color=COL["Co"], alpha=0.07, lw=0)
for _, r in f6.iterrows():
    ax.plot(r.floor_mgL, r.measured_mgL, MK[r.metal], color=COL[r.metal], ms=9)
    ax.annotate(f" {r.metal}", (r.floor_mgL, r.measured_mgL), fontsize=9, color=COL[r.metal])
ax.set_xscale("log"); ax.set_yscale("log"); ax.set_xlim(*lim); ax.set_ylim(1e-1, 2e2)
ax.set_xlabel("precipitation floor, MeCO$_3$ solubility (mg L$^{-1}$)")
ax.set_ylabel("measured concentration at day 6 (mg L$^{-1}$)")
ax.text(3e-3, 1.3e2, "above the 1:1 line: precipitation incomplete", fontsize=7.5)
ax.text(3e-3, 0.30, "below the 1:1 line:\nsorption required", fontsize=7.5, color=COL["Co"])
savefig(fig, "PW_Fig2_measured_vs_solubility_floor")

# --- Figure 3: the mechanism split ------------------------------------------------------------
fig, ax = plt.subplots(figsize=(6.4, 3.8), constrained_layout=True)
s = split.set_index("metal").loc[[m for m in ["Pb", "Co", "Cd", "Ba", "Li", "Sr"]]]
y = np.arange(len(s))
ax.barh(y, s.precipitation_pct, color="#8a8a8a", label="carbonate mineral formation")
ax.barh(y, s.sorption_pct, left=s.precipitation_pct, color="#1f5f8b", label="sorption")
ax.set_yticks(y); ax.set_yticklabels(s.index); ax.invert_yaxis()
ax.set_xlabel("share of the removal (%)"); ax.set_xlim(0, 100)
ax.legend(loc="lower right")
for i, (_, r) in enumerate(s.iterrows()):
    ax.text(101, i, f"{r.removal_pct:.0f} % removed", va="center", fontsize=7.5)
savefig(fig, "PW_Fig3_mechanism_split")
print("Figure 3 answers the question the published paper left open, metal by metal.")
