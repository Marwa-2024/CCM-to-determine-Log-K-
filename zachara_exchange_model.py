# %% [markdown]
# # Sorption of Li⁺ and Co²⁺ on dolomite: the Zachara ion-exchange model
#
# **Following Zachara, Cowan & Resch (1991)**, *Sorption of divalent metals on calcite*,
# Geochimica et Cosmochimica Acta 55, 1549–1562 — the paper that fitted surface exchange
# constants for seven divalent metals on a carbonate, **cobalt among them**, using FITEQL.
#
# ## Why this model and not the proton-exchange one
#
# Every framework tried so far writes the surface reaction as a **proton** exchange,
# `>CO3H0 + Me2+ = >CO3Me+ + H+`, and treats calcium as a nuisance competitor. Zachara writes it
# as a **calcium** exchange:
#
# $$\mathrm{X\text{-}Ca} + \mathrm{Me^{2+}} \rightleftharpoons \mathrm{X\text{-}Me} + \mathrm{Ca^{2+}}
# \qquad {}^{c}K_{ex} = \frac{(\mathrm{Ca^{2+}})[\mathrm{X\text{-}Me}]}
# {(\mathrm{Me^{2+}})[\mathrm{X\text{-}Ca}]} \qquad (11)$$
#
# Three consequences, all of which matter here:
#
# 1. **No proton appears**, so the constant is not sensitive to the pH drift that makes these
#    batch experiments hard to interpret.
# 2. **Calcium is the reference ion, not an interference.** The produced water carries
#    4186 mg L⁻¹ of calcium; in this formulation that is the denominator of the exchange, which is
#    physically what it is.
# 3. **It inverts directly from the measurements.** With sorbed metal, dissolved metal and
#    dissolved calcium all measured, K_ex follows without an optimiser.
#
# Zachara fitted the constants in FITEQL by fixing the calcium half-reaction at log K = 15 so the
# sites begin calcium-saturated, fitting the metal half-reaction to the sorption edge, and taking
# log ᶜK_ex = log K_Me − log K_Ca. They worked deliberately **below saturation with the metal
# solid phases**, and marked the isotherm points where precipitation began (10⁻⁵ to 10⁻⁴ mol L⁻¹
# initial metal) so that sorption and precipitation were never conflated.

# %%
import numpy as np, pandas as pd, matplotlib.pyplot as plt
import os
pd.set_option("display.width", 165); pd.set_option("display.precision", 4)
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

# %% [markdown]
# ## 1 — Zachara's fitted constants, Table 1
#
# These are measured values on calcite at I = 0.1, fitted by FITEQL to sorption edges below
# metal-solid saturation. Three of the six metals in the produced water appear here directly.

# %%
head(1, "Zachara, Cowan & Resch (1991), Table 1")

ZACHARA = pd.DataFrame([
    ("Cd", 0.97, -0.02, 3.02,  "fitted"),
    ("Zn", 0.74, -0.25, 2.43,  "fitted"),
    ("Mn", 0.80, -0.19, 1.31,  "fitted"),
    ("Co", 0.72, -0.27, 0.568, "fitted; needs Co(CO3)2-2 above pH 8.5 (log K 7.46)"),
    ("Ni", 0.69, -0.30, 0.51,  "fitted; needs a second complex X-NiOH"),
    ("Ba", 1.34,  0.35, -1.93, "fitted; essentially does not sorb"),
    ("Sr", 1.12,  0.13, -2.04, "fitted; essentially does not sorb"),
    ("Ca", 0.99,  0.00, 0.0,   "reference ion, by definition"),
], columns=["metal", "ionic_radius_A", "delta_r_vs_Ca", "log_cKex", "note"])
show(ZACHARA)

X_T_PER_G = 3.46e-6          # mol/g, cation-specific surface sites, from 45Ca isotopic exchange
print(f"\n  Site density X_T = {X_T_PER_G:.2e} mol/g, measured by 45Ca isotopic exchange.")
print("  Zachara caution that only about 10 % of these may be accessible to trace sorbates")
print("  (their p.1555, citing Zachara et al. 1988 and Cowan et al. 1990), so the sorption")
print("  capacity computed below is an UPPER bound.")
print("\n  Selectivity: Cd > Zn > Mn > Co > Ni >> Ba = Sr, which is the order of ionic radius")
print("  relative to Ca: metals that fit the calcite lattice sorb, metals that do not, do not.")

# %% [markdown]
# ## 2 — The data, and the sorption capacity in Zachara's units
#
# The exchange model has a hard ceiling: sorption cannot exceed the cation-specific site
# inventory X_T, because each site holds one cation. That ceiling is evaluated first, for both
# experiments, exactly as Zachara evaluated saturation before interpreting an isotherm.

# %%
head(2, "The data and the sorption ceiling")

MM = dict(Ca=40.078, Mg=24.305, Li=6.941, Co=58.933, Pb=207.2, Sr=87.62, Ba=137.33, Cd=112.41)
DOL_GL = 60.0
X_T = X_T_PER_G*DOL_GL                       # mol/L of suspension
print(f"  Dolomite {DOL_GL:.0f} g/L  x  {X_T_PER_G:.2e} mol/g  =  X_T = {X_T*1e3:.3f} mmol/L of exchange sites")

SINGLE = pd.DataFrame([
    ("Co","pH6",0,   0.100,  0.100, 80.513, 6.320), ("Co","pH6",2,   0.100,  0.100, 77.905, 7.400),
    ("Co","pH6",4,   0.100,  8.143, 76.411, 7.477), ("Co","pH6",6,   0.126, 12.941, 74.393, 8.090),
    ("Co","pH2",0,   0.100,  0.100, 80.513, 2.350), ("Co","pH2",2,  57.477, 17.461, 79.661, 4.123),
    ("Co","pH2",4,  24.178, 45.592, 79.795, 5.660), ("Co","pH2",6,  21.658, 51.379, 80.398, 7.110),
    ("Li","pH6",0,   0.663,  0.100,138.430, 6.420), ("Li","pH6",2,  37.633,  0.100,134.657, 7.647),
    ("Li","pH6",4,  54.820,  0.100,136.722, 7.923), ("Li","pH6",6,  49.892,  0.870,135.262, 8.257),
    ("Li","pH2",0,   0.663,  0.100,138.430, 2.250), ("Li","pH2",2,  39.190,  0.100,136.263, 4.257),
    ("Li","pH2",4,  79.452, 21.316,135.956, 5.727), ("Li","pH2",6,  87.738, 20.547,135.592, 7.047),
], columns=["metal","batch","day","Ca_ppm","Mg_ppm","Me_ppm","pH"])
SINGLE["C0_ppm"] = SINGLE.groupby(["metal","batch"]).Me_ppm.transform("first")

PW_METALS = ["Pb","Sr","Ba","Cd","Co","Li"]
PW = pd.DataFrame([
    ("pH6",0, 4185.950, 1005.970, 83.350, 88.579, 90.165, 92.673, 75.960, 76.780, 6.320),
    ("pH6",2, 3930.733,  940.894, 48.789, 82.369, 83.058, 84.949, 69.892, 71.243, 7.413),
    ("pH6",4, 3890.693,  984.050,  1.072, 82.500, 49.280, 25.739, 11.527, 71.470, 7.557),
    ("pH6",6, 3754.450,  967.633,  2.037, 83.166, 44.457, 14.066,  7.019, 70.212, 7.857),
    ("pH2",0, 4185.950, 1005.970, 83.350, 88.579, 90.165, 92.673, 75.960, 76.780, 2.350),
    ("pH2",2, 4037.933,  978.907, 47.138, 86.930, 87.731, 92.468, 75.013, 76.212, 4.290),
    ("pH2",4, 4199.987, 1115.780,  2.164, 86.297, 54.205, 43.352, 26.687, 74.399, 5.727),
], columns=["batch","day","Ca_ppm","Mg_ppm","Pb","Sr","Ba","Cd","Co","Li","pH"])

rows = []
for m in PW_METALS:
    c0 = PW[(PW.batch=="pH6")&(PW.day==0)][m].iloc[0]
    c6 = PW[(PW.batch=="pH6")&(PW.day==6)][m].iloc[0]
    rem = (c0-c6)/MM[m]/1e3
    rows.append(dict(experiment="produced water", metal=m, removed_mM=rem*1e3,
                     x_site_ceiling=rem/X_T, removal_pct=100*(1-c6/c0)))
for m in ["Co","Li"]:
    d = SINGLE[(SINGLE.metal==m)&(SINGLE.batch=="pH6")]
    rem = (d.C0_ppm.iloc[0]-d.Me_ppm.iloc[-1])/MM[m]/1e3
    rows.append(dict(experiment="single ion", metal=m, removed_mM=rem*1e3,
                     x_site_ceiling=rem/X_T,
                     removal_pct=100*(1-d.Me_ppm.iloc[-1]/d.C0_ppm.iloc[0])))
cap = pd.DataFrame(rows)
show(cap)
print("\n  x_site_ceiling > 1 means more metal left solution than there are sites to hold it, so")
print("  sorption cannot be the mechanism however large the exchange constant. Zachara applied")
print("  exactly this test to their isotherms and excluded the points that failed it.")

# %% [markdown]
# ## 3 — How much sorption does Zachara's own cobalt constant predict here?
#
# Zachara measured log ᶜK_ex(Co) = 0.568 on calcite. Applying it to the produced water needs no
# fitting at all: rearranging Eq. (11) with N_Me + N_Ca = 1,
#
# $$\frac{N_{Me}}{1-N_{Me}} = {}^{c}K_{ex}\,\frac{(\mathrm{Me^{2+}})}{(\mathrm{Ca^{2+}})}
# \qquad\Longrightarrow\qquad
# [\mathrm{X\text{-}Me}] = X_T\,\frac{{}^{c}K_{ex}(\mathrm{Me^{2+}})/(\mathrm{Ca^{2+}})}
# {1 + {}^{c}K_{ex}(\mathrm{Me^{2+}})/(\mathrm{Ca^{2+}})}$$
#
# The measured calcium concentration goes straight into the denominator, which is the whole point
# of writing the reaction this way.

# %%
head(3, "Sorption predicted by Zachara's cobalt constant")

def sorbed_by_exchange(logKex, Me_M, Ca_M, X_T=X_T):
    """Zachara Eq. 11 rearranged. Returns mol/L of metal held on the exchange sites."""
    r = 10**logKex*Me_M/Ca_M
    return X_T*r/(1.0 + r)

pred = []
for _, r in PW[PW.batch == "pH6"].iterrows():
    Ca = r.Ca_ppm/MM["Ca"]/1e3
    for m in ["Co", "Cd", "Ba", "Sr"]:
        lk = float(ZACHARA[ZACHARA.metal == m].log_cKex.iloc[0])
        Me = r[m]/MM[m]/1e3
        s  = sorbed_by_exchange(lk, Me, Ca)
        c0 = PW[(PW.batch=="pH6")&(PW.day==0)][m].iloc[0]
        removed = (c0 - r[m])/MM[m]/1e3
        pred.append(dict(day=int(r.day), metal=m, log_cKex=lk, Ca_mM=Ca*1e3,
                         Me_aq_mM=Me*1e3, sorbed_mM=s*1e3, removed_mM=removed*1e3,
                         pct_of_removal_by_sorption=(100*s/removed if removed > 0 else np.nan)))
pred = pd.DataFrame(pred)
show(pred[pred.day > 0])

co6 = pred[(pred.metal == "Co") & (pred.day == 6)].iloc[0]
print(f"\n  COBALT AT DAY 6, produced water:")
print(f"    calcium in solution        {co6.Ca_mM:8.1f} mmol/L")
print(f"    cobalt in solution         {co6.Me_aq_mM:8.4f} mmol/L")
print(f"    cobalt held by exchange    {co6.sorbed_mM:8.5f} mmol/L   (Zachara's constant, no fitting)")
print(f"    cobalt actually removed    {co6.removed_mM:8.4f} mmol/L")
print(f"    sorption accounts for      {co6.pct_of_removal_by_sorption:8.2f} % of the removal")
print(f"\n  Calcium at {co6.Ca_mM:.0f} mmol/L outcompetes cobalt at {co6.Me_aq_mM:.3f} mmol/L by a factor of")
print(f"  {co6.Ca_mM/co6.Me_aq_mM:,.0f} on the exchange sites, and Zachara's selectivity of "
      f"{10**co6.log_cKex:.1f} cannot overcome that.")
print("  The exchange sites stay calcium-saturated throughout, which is the condition Zachara")
print("  assumed when they fixed log K_Ca at 15. Essentially all of the cobalt removal in")
print("  produced water is carbonate mineral formation, not sorption.")

# %% [markdown]
# ## 4 — The single-ion experiment, and why the two give opposite answers
#
# The produced water keeps the exchange sites calcium-saturated, which is Zachara's assumed
# condition and the reason sorption contributes almost nothing there. The single-ion experiment
# is the mirror image: calcium is 0.13 mg L⁻¹, essentially absent, so the sites are *not*
# calcium-saturated and cobalt has nothing to compete with.
#
# This is where a constant can be obtained, by inverting Eq. (11) directly on the measurement.

# %%
head(4, "The single-ion experiment: inverting Eq. 11 for the constant")

inv = []
for metal in ["Co", "Li"]:
    for batch in ["pH6", "pH2"]:
        d = SINGLE[(SINGLE.metal == metal) & (SINGLE.batch == batch)]
        for _, r in d[d.day > 0].iterrows():
            Ca  = r.Ca_ppm/MM["Ca"]/1e3
            Me  = r.Me_ppm/MM[metal]/1e3
            srb = (r.C0_ppm - r.Me_ppm)/MM[metal]/1e3
            N   = srb/X_T
            ok  = 0.02 < N < 0.98                      # away from both ends of the isotherm
            k   = np.log10(Ca*N/(Me*(1 - N))) if ok else np.nan
            inv.append(dict(metal=metal, batch=batch, day=int(r.day), pH=r.pH,
                            Ca_mM=Ca*1e3, Me_aq_mM=Me*1e3, sorbed_mM=srb*1e3,
                            N_Me=N, log_cKex=k,
                            usable="yes" if ok else ("N > 0.98, sites full" if N >= 0.98
                                                     else "N < 0.02, below detection")))
inv = pd.DataFrame(inv)
show(inv)

good = inv.dropna(subset=["log_cKex"])
print("\n  Inverting Eq. (11) point by point, where the surface coverage N is inside the range the")
print("  exchange model can resolve:")
for metal in ["Co", "Li"]:
    g = good[good.metal == metal]
    if len(g):
        print(f"    {metal}:  log cKex = {g.log_cKex.mean():+.2f} +/- {g.log_cKex.std():.2f}   "
              f"from {len(g)} point(s), coverage N = {g.N_Me.min():.2f} to {g.N_Me.max():.2f}")
    else:
        print(f"    {metal}:  no point falls inside the resolvable coverage range")

zco = float(ZACHARA[ZACHARA.metal == "Co"].log_cKex.iloc[0])
print(f"\n  Zachara's calcite value for cobalt is log cKex = {zco:+.3f}.")
print("\n  THE CAVEAT THAT HAS TO BE STATED. Zachara's formulation references the exchange to")
print("  calcium, and they enforced calcium-saturated sites by fixing log K_Ca at 15. The")
print(f"  single-ion batches carry {SINGLE[(SINGLE.metal=='Co')&(SINGLE.batch=='pH6')].Ca_ppm.iloc[-1]:.2f} mg/L of calcium, which is not a calcium-saturated")
print("  surface, so the denominator of Eq. (11) is being evaluated far outside the conditions")
print("  the constant was calibrated under. The number above is what the equation returns; it is")
print("  not on the same footing as Zachara's value, and the produced-water experiment - which IS")
print("  calcium-saturated - is the one that matches their conditions but is dominated by")
print("  precipitation. That is the bind this data set is in, stated plainly.")

# %% [markdown]
# ## 5 — The answer to the question the published paper left open
#
# Elshebli et al. (2025) concluded that *"the contribution of sorption and carbonate mineral
# formation to the observed recovery levels of Li⁺ and Co²⁺ … needs further investigation."*
# Zachara's model answers it for the produced water without any fitting, because the exchange
# constants for four of the six metals are already measured.

# %%
head(5, "Sorption versus carbonate mineral formation, produced water, day 6")

split = []
for m in ["Cd", "Co", "Ba", "Sr"]:
    p = pred[(pred.metal == m) & (pred.day == 6)].iloc[0]
    split.append(dict(metal=m, log_cKex=p.log_cKex, removal_pct=100*p.removed_mM /
                      (PW[(PW.batch=="pH6")&(PW.day==0)][m].iloc[0]/MM[m]),
                      removed_mM=p.removed_mM, sorbed_mM=p.sorbed_mM,
                      sorption_pct_of_removal=p.pct_of_removal_by_sorption,
                      precipitation_pct=100 - p.pct_of_removal_by_sorption))
show(pd.DataFrame(split))
print("\n  Cadmium is the exception that proves the model works: it has by far the strongest")
print("  exchange constant (3.02) and is the only metal for which sorption takes a real share.")
print("  Cobalt, with a constant 2.5 log units weaker, gets 0.08 %. Barium and strontium, which")
print("  Zachara showed do not sorb at all, get 0.003 % and 0.03 % - and strontium is the control,")
print("  because its measured removal of 6 % must then be precipitation and analytical scatter.")
print("\n  For COBALT in produced water the answer is therefore: carbonate mineral formation")
print("  accounts for essentially all of the removal, and sorption for well under one per cent.")
print("  The reason is not that cobalt binds weakly - it is that 94 mmol/L of calcium holds every")
print("  exchange site, and there are only 0.21 mmol/L of them against 1.17 mmol/L of cobalt to")
print("  be removed.")
