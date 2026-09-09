#!/usr/bin/env python3
"""
Quantitative synthesis of the lithium-recovery evidence table.

Because the source reviews report point estimates without dispersion, an
inverse-variance random-effects model is not identifiable.  The synthesis is
therefore a pooled descriptive meta-analysis: robust location and spread
statistics, bootstrap confidence intervals on the median, non-parametric
between-family tests, and rank correlations for temporal trends.
"""
import os, json, csv
import numpy as np
import pandas as pd
from scipy import stats
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import ScalarFormatter

HERE = os.path.dirname(os.path.abspath(__file__))
FIG = os.path.join(HERE, "figures")
os.makedirs(FIG, exist_ok=True)
rng = np.random.default_rng(20260909)

df = pd.read_csv(os.path.join(HERE, "li_recovery_dataset.csv"))
for c in ["capacity_mg_g", "energy_Wh_mol", "recovery_pct", "purity_pct",
          "sel_Li_Mg", "cycles", "retention_pct", "year"]:
    df[c] = pd.to_numeric(df[c], errors="coerce")

OUT = {}

# --------------------------------------------------------------------------
def boot_ci_median(x, n=10000, alpha=0.05):
    x = np.asarray(x, float)
    x = x[~np.isnan(x)]
    if len(x) < 3:
        return (np.nan, np.nan)
    b = rng.choice(x, size=(n, len(x)), replace=True)
    m = np.median(b, axis=1)
    return (float(np.percentile(m, 100 * alpha / 2)),
            float(np.percentile(m, 100 * (1 - alpha / 2))))

def describe(x, label):
    x = np.asarray(pd.to_numeric(pd.Series(x), errors="coerce").dropna(), float)
    if len(x) == 0:
        return None
    lo, hi = boot_ci_median(x)
    d = dict(
        label=label, n=int(len(x)),
        mean=float(np.mean(x)), sd=float(np.std(x, ddof=1)) if len(x) > 1 else 0.0,
        median=float(np.median(x)),
        q1=float(np.percentile(x, 25)), q3=float(np.percentile(x, 75)),
        min=float(np.min(x)), max=float(np.max(x)),
        ci_lo=lo, ci_hi=hi,
        cv=float(np.std(x, ddof=1) / np.mean(x)) if len(x) > 1 and np.mean(x) else np.nan,
    )
    return d

def table(rows, title):
    print("\n" + "=" * 96)
    print(title)
    print("=" * 96)
    hdr = f"{'group':<42}{'n':>4}{'median':>10}{'IQR':>20}{'mean':>10}{'SD':>10}{'range':>20}"
    print(hdr)
    print("-" * 96)
    for r in rows:
        if r is None:
            continue
        print(f"{r['label']:<42}{r['n']:>4}{r['median']:>10.2f}"
              f"{('%.2f-%.2f' % (r['q1'], r['q3'])):>20}{r['mean']:>10.2f}"
              f"{r['sd']:>10.2f}{('%.3g-%.3g' % (r['min'], r['max'])):>20}")
    return rows

# ==========================================================================
# 1. ADSORPTION CAPACITY
# ==========================================================================
ads = df[df.family == "Adsorption / ion sieve"].dropna(subset=["capacity_mg_g"])
hyb = df[df.family == "Hybrid membrane-adsorbent"].dropna(subset=["capacity_mg_g"])
els = df[(df.family == "Electrosorption / ion pumping")].dropna(subset=["capacity_mg_g"])

rows = [describe(ads.capacity_mg_g, "All adsorbents / ion sieves")]
for sub, g in ads.groupby("subclass"):
    if len(g) >= 2:
        rows.append(describe(g.capacity_mg_g, "   " + sub))
rows.append(describe(hyb.capacity_mg_g, "Hybrid adsorptive membranes"))
rows.append(describe(els.capacity_mg_g, "Electrosorption electrodes"))
OUT["capacity"] = table([r for r in rows if r], "TABLE A. Lithium uptake capacity (mg Li per g of material)")

# Ti- vs Mn-based sieves
ti = ads[ads.subclass.str.contains("Ti", na=False)].capacity_mg_g.dropna()
mn = ads[ads.subclass.str.contains("Mn", na=False)].capacity_mg_g.dropna()
if len(ti) >= 3 and len(mn) >= 3:
    u, p = stats.mannwhitneyu(ti, mn, alternative="two-sided")
    OUT["ti_vs_mn"] = dict(n_ti=int(len(ti)), n_mn=int(len(mn)),
                           med_ti=float(np.median(ti)), med_mn=float(np.median(mn)),
                           U=float(u), p=float(p))
    print(f"\nTi-based (n={len(ti)}, median {np.median(ti):.1f}) vs Mn-based "
          f"(n={len(mn)}, median {np.median(mn):.1f}): Mann-Whitney U={u:.1f}, p={p:.3f}")

# adsorbent vs electrosorption capacity
u, p = stats.mannwhitneyu(ads.capacity_mg_g, els.capacity_mg_g, alternative="two-sided")
OUT["ads_vs_els"] = dict(U=float(u), p=float(p),
                         med_ads=float(ads.capacity_mg_g.median()),
                         med_els=float(els.capacity_mg_g.median()))
print(f"Adsorbents vs electrosorption electrodes: U={u:.1f}, p={p:.4f}")

# ==========================================================================
# 2. ENERGY CONSUMPTION
# ==========================================================================
en_el = df[(df.family == "Electrosorption / ion pumping")].dropna(subset=["energy_Wh_mol"])
en_ed = df[(df.family == "Electrodialysis / membrane electro-process")].dropna(subset=["energy_Wh_mol"])
en_cdi = df[(df.family == "Capacitive deionisation")].dropna(subset=["energy_Wh_mol"])

rows = [describe(en_el.energy_Wh_mol, "Electrosorption / ion pumping"),
        describe(en_ed.energy_Wh_mol, "Electrodialysis family"),
        describe(en_cdi.energy_Wh_mol, "Capacitive deionisation")]
OUT["energy"] = table([r for r in rows if r], "TABLE B. Specific energy consumption (Wh per mol Li recovered)")

u, p = stats.mannwhitneyu(en_el.energy_Wh_mol, en_ed.energy_Wh_mol, alternative="two-sided")
OUT["energy_test"] = dict(U=float(u), p=float(p),
                          med_el=float(en_el.energy_Wh_mol.median()),
                          med_ed=float(en_ed.energy_Wh_mol.median()),
                          ratio=float(en_ed.energy_Wh_mol.median() / en_el.energy_Wh_mol.median()))
print(f"\nElectrosorption vs electrodialysis energy: U={u:.1f}, p={p:.4f}, "
      f"median ratio {OUT['energy_test']['ratio']:.1f}x")

# ==========================================================================
# 3. RECOVERY EFFICIENCY BY FAMILY
# ==========================================================================
rec_groups, rec_rows = {}, []
for fam in ["Precipitation", "Solvent extraction",
            "Electrodialysis / membrane electro-process", "Capacitive deionisation"]:
    v = df[df.family == fam].recovery_pct.dropna()
    if len(v) >= 3:
        rec_groups[fam] = v.values
        rec_rows.append(describe(v, fam))
OUT["recovery"] = table(rec_rows, "TABLE C. Reported lithium recovery / extraction efficiency (%)")
if len(rec_groups) >= 3:
    H, p = stats.kruskal(*rec_groups.values())
    OUT["recovery_kw"] = dict(H=float(H), p=float(p), k=len(rec_groups))
    print(f"\nKruskal-Wallis across families: H={H:.2f}, df={len(rec_groups)-1}, p={p:.4f}")
    # pairwise
    keys = list(rec_groups)
    pw = []
    for i in range(len(keys)):
        for j in range(i + 1, len(keys)):
            u, pp = stats.mannwhitneyu(rec_groups[keys[i]], rec_groups[keys[j]],
                                       alternative="two-sided")
            pw.append(dict(a=keys[i], b=keys[j], U=float(u), p=float(pp),
                           p_bonf=float(min(1.0, pp * (len(keys) * (len(keys) - 1) / 2)))))
    OUT["recovery_pairwise"] = pw
    print("\nPairwise (Bonferroni-adjusted):")
    for r in pw:
        print(f"  {r['a'][:28]:<30} vs {r['b'][:28]:<30} p={r['p']:.4f}  p_adj={r['p_bonf']:.4f}")

# ==========================================================================
# 4. PRODUCT PURITY
# ==========================================================================
pur_rows = []
for fam in ["Electrosorption / ion pumping", "Electrodialysis / membrane electro-process",
            "Precipitation"]:
    v = df[df.family == fam].purity_pct.dropna()
    if len(v) >= 3:
        pur_rows.append(describe(v, fam))
OUT["purity"] = table(pur_rows, "TABLE D. Purity of the recovered lithium stream (%)")

# ==========================================================================
# 5. Li/Mg SELECTIVITY
# ==========================================================================
sel_rows = []
sel_groups = {}
for fam in ["Electrosorption / ion pumping", "Nanofiltration / selective membrane",
            "Electrodialysis / membrane electro-process"]:
    v = df[df.family == fam].sel_Li_Mg.dropna()
    if len(v) >= 3:
        sel_rows.append(describe(v, fam))
        sel_groups[fam] = v.values
OUT["selectivity"] = table(sel_rows, "TABLE E. Li/Mg separation factor (dimensionless)")
if len(sel_groups) >= 2:
    keys = list(sel_groups)
    u, p = stats.mannwhitneyu(sel_groups[keys[0]], sel_groups[keys[1]], alternative="two-sided")
    OUT["sel_test"] = dict(a=keys[0], b=keys[1], U=float(u), p=float(p))
    print(f"\n{keys[0][:35]} vs {keys[1][:35]}: U={u:.1f}, p={p:.4f}")

# ==========================================================================
# 6. TEMPORAL TREND
# ==========================================================================
tt = df.dropna(subset=["year", "capacity_mg_g"])
tt = tt[tt.family.isin(["Adsorption / ion sieve", "Electrosorption / ion pumping"])]
rho, p = stats.spearmanr(tt.year, tt.capacity_mg_g)
OUT["trend_capacity"] = dict(n=int(len(tt)), rho=float(rho), p=float(p))
print(f"\nCapacity vs publication year (all materials): Spearman rho={rho:.3f}, "
      f"p={p:.4f}, n={len(tt)}")

te = df.dropna(subset=["year", "energy_Wh_mol"])
te = te[te.family == "Electrosorption / ion pumping"]
rho2, p2 = stats.spearmanr(te.year, te.energy_Wh_mol)
OUT["trend_energy"] = dict(n=int(len(te)), rho=float(rho2), p=float(p2))
print(f"Electrosorption energy vs year: Spearman rho={rho2:.3f}, p={p2:.4f}, n={len(te)}")

# log-linear fit for energy trend
if len(te) >= 5:
    x = te.year.values.astype(float)
    y = np.log10(te.energy_Wh_mol.values.astype(float))
    sl, ic, r, pv, se = stats.linregress(x, y)
    half_life = np.log10(2) / abs(sl) if sl else np.nan
    OUT["energy_regression"] = dict(slope_log10_per_year=float(sl), r2=float(r ** 2),
                                    p=float(pv), years_per_halving=float(half_life))
    print(f"log10(energy) vs year: slope={sl:.4f}/yr, R2={r**2:.3f}, p={pv:.4f}, "
          f"halving time {half_life:.1f} yr")

# ==========================================================================
# 7. CYCLING STABILITY
# ==========================================================================
cy = df.dropna(subset=["cycles", "retention_pct"])
if len(cy) >= 5:
    rho3, p3 = stats.spearmanr(cy.cycles, cy.retention_pct)
    OUT["cycling"] = dict(n=int(len(cy)), rho=float(rho3), p=float(p3),
                          median_retention=float(cy.retention_pct.median()),
                          median_cycles=float(cy.cycles.median()))
    # fade per cycle
    fade = (100 - cy.retention_pct) / cy.cycles
    OUT["fade"] = describe(fade, "capacity fade (% per cycle)")
    print(f"\nRetention vs cycle number: rho={rho3:.3f}, p={p3:.4f}, n={len(cy)}")
    print(f"Median retention {cy.retention_pct.median():.0f}% at median "
          f"{cy.cycles.median():.0f} cycles; median fade "
          f"{np.median(fade):.3f}% per cycle")

# ==========================================================================
# 8. REPORTING COMPLETENESS  (a reporting-bias analogue)
# ==========================================================================
metrics = ["capacity_mg_g", "energy_Wh_mol", "recovery_pct", "purity_pct",
           "sel_Li_Mg", "cycles", "retention_pct"]
comp = {}
for m in metrics:
    comp[m] = float(100.0 * df[m].notna().sum() / len(df))
OUT["completeness"] = comp
print("\n" + "=" * 96)
print("TABLE F. Reporting completeness across the 131 extracted records (%)")
print("=" * 96)
for m, v in sorted(comp.items(), key=lambda kv: -kv[1]):
    print(f"  {m:<20}{v:>8.1f}%")

# ==========================================================================
# FIGURES
# ==========================================================================
plt.rcParams.update({"font.size": 9, "font.family": "DejaVu Sans",
                     "axes.spines.top": False, "axes.spines.right": False,
                     "figure.dpi": 200, "savefig.bbox": "tight"})
C = ["#2F6C8F", "#C46A3B", "#4E8A6B", "#8B5E9E", "#B0473F", "#7A7A52"]

# Fig 1 - capacity distribution by material family
fig, ax = plt.subplots(figsize=(6.6, 3.6))
groups, labels = [], []
for name, sel in [("Mn oxide\nsieves", ads[ads.subclass.str.contains("Mn", na=False)]),
                  ("Ti oxide\nsieves", ads[ads.subclass.str.contains("Ti", na=False)]),
                  ("LDH\nsorbents", ads[ads.subclass.str.contains("hydroxide", na=False)]),
                  ("MOF\nsorbents", ads[ads.subclass.str.contains("organic", na=False)]),
                  ("Hybrid\nmembranes", hyb),
                  ("Electro-\nsorption", els)]:
    v = sel.capacity_mg_g.dropna().values
    if len(v):
        groups.append(v); labels.append(f"{name}\n(n={len(v)})")
bp = ax.boxplot(groups, tick_labels=labels, patch_artist=True, widths=0.55,
                medianprops=dict(color="black", lw=1.4), showfliers=False)
for patch, c in zip(bp["boxes"], C):
    patch.set_facecolor(c); patch.set_alpha(0.45); patch.set_edgecolor(c); patch.set_lw(1.2)
for i, v in enumerate(groups):
    ax.scatter(np.full(len(v), i + 1) + rng.normal(0, 0.055, len(v)), v,
               s=16, color=C[i % len(C)], edgecolor="white", lw=0.5, zorder=3)
ax.set_ylabel("Li uptake capacity (mg g$^{-1}$)")
ax.set_title("Figure 1. Lithium uptake capacity by material family", loc="left", fontsize=10)
ax.grid(axis="y", ls=":", alpha=0.5)
fig.savefig(os.path.join(FIG, "fig1_capacity.png")); plt.close(fig)

# Fig 2 - energy consumption, log scale
fig, ax = plt.subplots(figsize=(6.0, 3.4))
groups, labels = [], []
for name, v in [("Electrosorption /\nion pumping", en_el.energy_Wh_mol.dropna().values),
                ("Electrodialysis\nfamily", en_ed.energy_Wh_mol.dropna().values),
                ("Capacitive\ndeionisation", en_cdi.energy_Wh_mol.dropna().values)]:
    if len(v):
        groups.append(v); labels.append(f"{name}\n(n={len(v)})")
bp = ax.boxplot(groups, tick_labels=labels, patch_artist=True, widths=0.5,
                medianprops=dict(color="black", lw=1.4), showfliers=False)
for patch, c in zip(bp["boxes"], C):
    patch.set_facecolor(c); patch.set_alpha(0.45); patch.set_edgecolor(c); patch.set_lw(1.2)
for i, v in enumerate(groups):
    ax.scatter(np.full(len(v), i + 1) + rng.normal(0, 0.05, len(v)), v, s=18,
               color=C[i], edgecolor="white", lw=0.5, zorder=3)
ax.set_yscale("log"); ax.set_ylabel("Specific energy (Wh mol$^{-1}$ Li)")
ax.set_title("Figure 2. Specific energy consumption by process family",
             loc="left", fontsize=10)
ax.grid(axis="y", ls=":", alpha=0.5, which="both")
fig.savefig(os.path.join(FIG, "fig2_energy.png")); plt.close(fig)

# Fig 3 - recovery efficiency by family
fig, ax = plt.subplots(figsize=(6.4, 3.4))
groups, labels = [], []
short = {"Precipitation": "Precipitation", "Solvent extraction": "Solvent\nextraction",
         "Electrodialysis / membrane electro-process": "Electrodialysis\nfamily",
         "Capacitive deionisation": "Capacitive\ndeionisation"}
for fam, v in rec_groups.items():
    groups.append(v); labels.append(f"{short[fam]}\n(n={len(v)})")
bp = ax.boxplot(groups, tick_labels=labels, patch_artist=True, widths=0.5,
                medianprops=dict(color="black", lw=1.4), showfliers=False)
for patch, c in zip(bp["boxes"], C):
    patch.set_facecolor(c); patch.set_alpha(0.45); patch.set_edgecolor(c); patch.set_lw(1.2)
for i, v in enumerate(groups):
    ax.scatter(np.full(len(v), i + 1) + rng.normal(0, 0.05, len(v)), v, s=18,
               color=C[i], edgecolor="white", lw=0.5, zorder=3)
ax.set_ylabel("Reported Li recovery / extraction efficiency (%)")
ax.set_ylim(0, 105)
ax.set_title("Figure 3. Recovery efficiency by technology family", loc="left", fontsize=10)
ax.grid(axis="y", ls=":", alpha=0.5)
fig.savefig(os.path.join(FIG, "fig3_recovery.png")); plt.close(fig)

# Fig 4 - energy vs year for electrosorption
fig, ax = plt.subplots(figsize=(6.0, 3.4))
ax.scatter(te.year, te.energy_Wh_mol, s=40, color=C[0], edgecolor="white", lw=0.7, zorder=3)
if len(te) >= 5:
    xs = np.linspace(te.year.min(), te.year.max(), 50)
    ax.plot(xs, 10 ** (ic + sl * xs), color=C[1], lw=1.8,
            label=f"log-linear fit ($R^2$={r**2:.2f}, p={pv:.3f})")
    ax.legend(frameon=False, fontsize=8)
ax.set_yscale("log"); ax.set_xlabel("Publication year")
ax.set_ylabel("Specific energy (Wh mol$^{-1}$ Li)")
ax.set_title("Figure 4. Energy intensity of electrosorption systems over time",
             loc="left", fontsize=10)
ax.grid(ls=":", alpha=0.5, which="both")
fig.savefig(os.path.join(FIG, "fig4_energy_trend.png")); plt.close(fig)

# Fig 5 - selectivity, log scale
fig, ax = plt.subplots(figsize=(6.0, 3.4))
groups, labels = [], []
sshort = {"Electrosorption / ion pumping": "Electrosorption",
          "Nanofiltration / selective membrane": "Nanofiltration /\nselective membrane",
          "Electrodialysis / membrane electro-process": "Electrodialysis\nfamily"}
for fam, v in sel_groups.items():
    groups.append(v); labels.append(f"{sshort[fam]}\n(n={len(v)})")
bp = ax.boxplot(groups, tick_labels=labels, patch_artist=True, widths=0.45,
                medianprops=dict(color="black", lw=1.4), showfliers=False)
for patch, c in zip(bp["boxes"], C):
    patch.set_facecolor(c); patch.set_alpha(0.45); patch.set_edgecolor(c); patch.set_lw(1.2)
for i, v in enumerate(groups):
    ax.scatter(np.full(len(v), i + 1) + rng.normal(0, 0.05, len(v)), v, s=18,
               color=C[i], edgecolor="white", lw=0.5, zorder=3)
ax.set_yscale("log"); ax.set_ylabel("Li/Mg separation factor")
ax.set_title("Figure 5. Li/Mg selectivity by separation principle", loc="left", fontsize=10)
ax.grid(axis="y", ls=":", alpha=0.5, which="both")
fig.savefig(os.path.join(FIG, "fig5_selectivity.png")); plt.close(fig)

# Fig 6 - reporting completeness
fig, ax = plt.subplots(figsize=(6.0, 3.2))
names = {"capacity_mg_g": "Uptake capacity", "energy_Wh_mol": "Specific energy",
         "recovery_pct": "Recovery efficiency", "purity_pct": "Product purity",
         "sel_Li_Mg": "Li/Mg selectivity", "cycles": "Cycle number tested",
         "retention_pct": "Capacity retention"}
ks = sorted(comp, key=lambda k: comp[k])
ax.barh([names[k] for k in ks], [comp[k] for k in ks], color=C[0], alpha=0.75, height=0.62)
for i, k in enumerate(ks):
    ax.text(comp[k] + 1, i, f"{comp[k]:.0f}%", va="center", fontsize=8)
ax.set_xlim(0, 62); ax.set_xlabel("Records reporting the metric (% of 131)")
ax.set_title("Figure 6. Reporting completeness across the evidence base",
             loc="left", fontsize=10)
ax.grid(axis="x", ls=":", alpha=0.5)
fig.savefig(os.path.join(FIG, "fig6_completeness.png")); plt.close(fig)

# Fig 7 - capacity retention vs cycles
fig, ax = plt.subplots(figsize=(6.0, 3.4))
ax.scatter(cy.cycles, cy.retention_pct, s=44, color=C[2], edgecolor="white", lw=0.7, zorder=3)
for _, rr in cy.iterrows():
    ax.annotate(f"{rr['study'].split()[0]} {int(rr['year'])}", (rr.cycles, rr.retention_pct),
                fontsize=6, xytext=(3, 3), textcoords="offset points", alpha=0.75)
ax.set_xscale("log")
ax.set_xlabel("Number of cycles tested"); ax.set_ylabel("Capacity retention (%)")
ax.set_title("Figure 7. Cycling stability of lithium-capture electrodes",
             loc="left", fontsize=10)
ax.grid(ls=":", alpha=0.5, which="both")
fig.savefig(os.path.join(FIG, "fig7_cycling.png")); plt.close(fig)

with open(os.path.join(HERE, "results.json"), "w") as fh:
    json.dump(OUT, fh, indent=2, default=float)
print("\nFigures ->", FIG)
print("Results  ->", os.path.join(HERE, "results.json"))
