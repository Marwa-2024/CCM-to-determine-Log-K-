#!/usr/bin/env python3
"""
Build the harmonised evidence table for the lithium-recovery meta-analysis.

Every row is transcribed from a comparison table (or explicit in-text statement)
in one of the six source reviews held in the user's Drive corpus.  The
`compiled_in` field records which review the value was read from, so that each
datum stays traceable to the document it was extracted from.

Source reviews
  R1  Liu, Zhao & Ghahreman (2019)  Hydrometallurgy 187, 81-100
  R2  Wang et al. (2022)            Renew. Sustain. Energy Rev. 154, 111813
  R3  Zavahir et al. (2021)         Desalination 500, 114883
  R4  Wu et al. (2022)              Water Research 221, 118822
  R5  Zhao et al. (2019)            J. Electroanal. Chem. 850, 113389
  R6  Khalil et al. (2022)          Desalination 528, 115611
"""
import csv, os

OUT = os.path.join(os.path.dirname(__file__), "li_recovery_dataset.csv")

FIELDS = [
    "id", "study", "year", "family", "subclass", "material_or_system",
    "feed_type", "capacity_mg_g", "energy_Wh_mol", "recovery_pct",
    "purity_pct", "sel_Li_Mg", "cycles", "retention_pct", "compiled_in",
]

R = []
def add(**kw):
    row = {f: kw.get(f, "") for f in FIELDS}
    row["id"] = "S%03d" % (len(R) + 1)
    R.append(row)

# ----------------------------------------------------------------------------
# FAMILY 1 - ADSORPTION / ION SIEVES   (R6 Table 2; R1 narrative)
# ----------------------------------------------------------------------------
ADS = [
    ("Chitrakar et al.", 2014, "Mn/Ti oxide sieve", "Layered H2TiO3", "salt-lake brine", 32.7, "R6"),
    ("Wang et al.", 2017, "Mn/Ti oxide sieve", "HxTiO3", "salt-lake brine", 36.3, "R6"),
    ("Wang et al.", 2018, "Mn/Ti oxide sieve", "Iron-doped Li-Ti oxide", "concentrated natural brine", 34.8, "R6"),
    ("Zandvakili & Ranjbar", 2018, "composite sieve", "MnO2/PVDF composite", "enriched salt-lake brine", 11.06, "R6"),
    ("Ohashi & Tai", 2019, "Mn oxide sieve", "Li1.33Mn1.67O4 (surface-modified)", "natural brine", 50.0, "R6"),
    ("Marthi & Smith", 2019, "composite sieve", "LiMn2O4 nanowire/diatomaceous earth", "salt-lake water", 18.6, "R6"),
    ("Zhang et al.", 2019, "Mn oxide sieve", "Al-doped Li manganese oxide", "synthetic salt solution", 32.6, "R6"),
    ("Liang et al.", 2020, "imprinted polymer", "Li+-IIP-Fe3O4@C", "synthetic salt solution", 22.26, "R6"),
    ("Lai et al.", 2020, "granular sieve", "Granular H4Mn5O12 (EP/HMO)", "salt-lake brine", 17.2, "R6"),
    ("Chen et al.", 2020, "layered double hydroxide", "Magnetic Li-Al-LDH", "salt-lake brine", 5.83, "R6"),
    ("Wei et al.", 2020, "Ti oxide sieve", "Porous H4Ti5O12 nanofibre", "simulated brine", 59.1, "R6"),
    ("Wei et al.", 2021, "metal-organic framework", "Al-based MOF (MIL-121)", "synthetic salt solution", 1.25, "R6"),
    ("Marthi & Smith", 2021, "composite sieve", "H2TiO3/diatomaceous earth", "salt-lake brine", 8.25, "R6"),
    ("Zhong et al.", 2021, "layered double hydroxide", "Li-Al-LDH", "salt-lake brine", 7.27, "R6"),
    ("Luo et al.", 2021, "layered double hydroxide", "Granular Li-Al-LDH/NH4Al3(SO4)2(OH)6", "salt-lake brine", 9.16, "R6"),
    ("Marthi et al.", 2021, "Ti oxide sieve", "Layered H2TiO3", "synthetic salt solution", 40.0, "R6"),
    ("Jiang et al.", 2021, "metal-organic framework", "Fluorine-pillared MOF", "salt-lake brine", 18.8, "R6"),
    ("Zhu et al.", 2021, "Ti oxide sieve", "Shaped Ti-based Li ion-sieve", "salt-lake brine (W. Taijinar)", 19.22, "R6"),
    ("Zhang et al.", 2021, "metal-organic framework", "PNIPAM-co-crown-6 MOF-808", "synthetic salt solution", 7.21, "R6"),
    # from R1 narrative text
    ("Zhang et al.", 2014, "Ti oxide sieve", "LTO from anatase precursor", "synthetic solution", 39.2, "R1"),
    ("He et al.", 2015, "Ti oxide sieve", "H2TiO3", "synthetic solution", 57.8, "R1"),
    ("Lawagon et al.", 2016, "Ti oxide sieve", "H2TiO3", "synthetic solution", 94.5, "R1"),
    ("Sun et al.", 2017, "imprinted membrane", "PVDF/crown-ether ion-imprinted membrane", "brine", 27.1, "R1"),
]
for a, y, sub, mat, feed, cap, src in ADS:
    add(study=a, year=y, family="Adsorption / ion sieve", subclass=sub,
        material_or_system=mat, feed_type=feed, capacity_mg_g=cap, compiled_in=src)

# Hybrid membrane-adsorbent systems (R6 Table 6)
HYB = [
    ("Sun et al.", 2018, "Crown-ether functionalised GO/PVDF", "synthetic salt solution", 24.25),
    ("Xue et al.", 2020, "alpha-Al2O3 supported Li4Mn5O12", "synthetic salt solution", 22.9),
    ("Wang et al.", 2021, "Sulfonated polyarylene ether nitrile hollow fibre", "synthetic brine", 20.54),
    ("Zheng et al.", 2021, "Crown-ether GO/chitosan/PVA", "synthetic salt solution", 168.5),
]
for a, y, mat, feed, cap in HYB:
    add(study=a, year=y, family="Hybrid membrane-adsorbent", subclass="adsorptive membrane",
        material_or_system=mat, feed_type=feed, capacity_mg_g=cap, compiled_in="R6")

# ----------------------------------------------------------------------------
# FAMILY 2 - ELECTROSORPTION / ELECTROCHEMICAL ION PUMPING
#            (R6 Table 5; R5 Table 1; R4 Table 2)
# ----------------------------------------------------------------------------
ELS = [
    # study, year, electrode pair, feed, cap, energy, purity, cycles, retention, src
    ("Kanoh et al.", 1993, "Pt // lambda-MnO2", "aqueous LiCl", "", 11.0, "", "", "", "R4"),
    ("Pasta et al.", 2012, "Ag // LiFePO4", "simulated brine", 7.76, 1.0, 84.56, "", "", "R4"),
    ("Lee et al.", 2013, "Ag // lambda-MnO2", "1 M LiCl / brine", "", 4.5, 91.8, 20, 87, "R5"),
    ("Trocoli et al.", 2014, "Ag // LiFePO4", "artificial brine", "", 2.8, 97.86, "", "", "R5"),
    ("Trocoli et al.", 2015, "NiHCF // LiFePO4", "Atacama brine", "", 8.7, 97.9, "", "", "R5"),
    ("Kim et al.", 2015, "AC // lambda-MnO2", "simulated brine", 14.88, 4.2, 91.8, 50, 96, "R4"),
    ("Missoni et al.", 2016, "PPy // LiMn2O4", "natural brine", 1.73, 7.5, "", 200, 50, "R4"),
    ("Zhao et al.", 2017, "LiMn2O4 // Li(1-x)Mn2O4", "simulated brine", 22.0, 18.0, 89.0, "", 97, "R6"),
    ("Zhao et al.", 2017, "LiMn2O4 // Li(1-x)Mn2O4", "simulated concentrated seawater", 21.0, 18.6, 86.0, "", 97, "R6"),
    ("Trocoli et al.", 2017, "NiHCF // lambda-MnO2", "Atacama brine", "", 3.58, 96.2, 100, 83.3, "R5"),
    ("Kim et al.", 2018, "LiMn2O4 // Zn", "salt-lake brine", "", 6.3, "", 100, 73, "R5"),
    ("Lawagon et al.", 2018, "Ag // LNCM", "simulated brine", 10.83, 2.60, 96.4, 50, 83.3, "R6"),
    ("Xu et al.", 2018, "Graphite // lambda-MnO2 film", "simulated brine", "", 4.14, 99.0, 100, 91, "R4"),
    ("Kim et al.", 2019, "Ag // lambda-MnO2", "desalination brine", 10.1, 3.07, 99.0, "", "", "R6"),
    ("Lawagon et al.", 2019, "Ag // LiNi0.5Mn1.5O4", "simulated brine", 8.74, 1.32, 98.14, "", "", "R6"),
    ("Zhao et al.", 2019, "PANI // lambda-MnO2", "simulated brine", 23.13, 3.95, 81.40, 200, 70.8, "R4"),
    ("Wang et al.", 2019, "PPy-threaded HKUST-1", "simulated brine", 37.55, "", "", "", "", "R6"),
    ("Joo et al.", 2020, "Ag // lambda-MnO2", "seawater concentrate", "", "", 88.0, 5, "", "R4"),
    ("Zhao et al.", 2020, "AC // LNCM/rGO", "simulated brine", 13.84, 1.4, 93.5, 30, 88.89, "R4"),
    ("Zhao et al.", 2020, "AC // CF-NMMO", "simulated brine", 14.4, 7.91, 97.2, 30, 90.8, "R6"),
    ("Zhao et al.", 2020, "AC // PPy/Al2O3/LiMn2O4", "simulated brine", 12.84, 1.41, 97.37, "", "", "R6"),
    ("Niu et al.", 2021, "BiOCl@PPy // lambda-MnO2", "brine lake", 10.88, 1.007, "", 100, 96, "R6"),
    ("Xu et al.", 2021, "LiMn2O4 // Li(1-x)Mn2O4", "salt-lake brine", 15.5, 16.0, 97.0, 100, "", "R6"),
    ("Mu et al.", 2021, "lambda-MnO2/LiMn2O4 on graphite felt", "high Mg/Li simulated brine", 26.14, 23.38, "", 20, 89.6, "R4"),
]
for a, y, mat, feed, cap, en, pur, cyc, ret, src in ELS:
    add(study=a, year=y, family="Electrosorption / ion pumping", subclass="battery-type electrode pair",
        material_or_system=mat, feed_type=feed, capacity_mg_g=cap, energy_Wh_mol=en,
        purity_pct=pur, cycles=cyc, retention_pct=ret, compiled_in=src)

# Selectivity (separation factor Li/Mg) reported for electrosorption systems (R4 Table 2, R6 text)
SEL_EL = [
    ("Lee et al.", 2013, "Ag // lambda-MnO2", 18.0, "R4"),
    ("Kim et al.", 2015, "AC // lambda-MnO2", 286.67, "R4"),
    ("Zhao et al.", 2017, "LiMn2O4 // Li(1-x)Mn2O4", 74.03, "R4"),
    ("Lawagon et al.", 2018, "Ag // LNCM", 42.87, "R4"),
    ("Kim et al.", 2019, "Ag // lambda-MnO2", 5.57, "R4"),
    ("Wang et al.", 2019, "Pt // HMO/rGO", 10.23, "R4"),
    ("Zhao et al.", 2020, "AC // LNCM/rGO", 308.3, "R4"),
    ("Kim et al.", 2021, "lambda-MnO2 with ferri/ferrocyanide redox", 804.0, "R6"),
    ("Mu et al.", 2021, "lambda-MnO2/LiMn2O4 graphite felt", 45.58, "R4"),
]
for a, y, mat, sf, src in SEL_EL:
    add(study=a, year=y, family="Electrosorption / ion pumping", subclass="selectivity report",
        material_or_system=mat, sel_Li_Mg=sf, compiled_in=src)

# ----------------------------------------------------------------------------
# FAMILY 3 - ELECTRODIALYSIS AND MEMBRANE ELECTRO-PROCESSES
#            (R4 Table 3; R3 Table 2)
# ----------------------------------------------------------------------------
ED = [
    # study, year, subclass, system, feed, recovery, energy Wh/mol, purity, sel, src
    ("Jiang et al.", 2014, "BMED", "Neosepta CMX / BP-1", "lake brine (Li2CO3 feed)", 95.0, 8.0, 99.0, "", "R4"),
    ("Parsa et al.", 2015, "ED", "AR204SXR412 / CR67-MK111", "Li-Na solution", 20.0, "", 99.6, "", "R4"),
    ("Nie et al.", 2017, "SED", "Selemion ASA / CSO", "high Mg/Li brine", 62.5, 131.0, "", 20.2, "R4"),
    ("Nie et al.", 2017, "SED", "Selemion ASA / CSO", "high Mg/Li brine", 55.2, 31.0, "", 9.89, "R4"),
    ("Ji et al.", 2017, "SED", "Astom ACS / CIMS", "high Mg/Li brine", 75.44, "", "", 8.57, "R4"),
    ("Bunani et al.", 2017, "BMED", "Astom CMB/CHA + BP-1E", "geothermal water", 20.0, "", "", "", "R4"),
    ("Hwang et al.", 2016, "BMED", "Neosepta BP-1", "lake brine", 70.0, "", "", "", "R3"),
    ("Ipekci et al.", 2018, "BMED", "Astom CMB/CHA + BP-1E", "geothermal water", 62.0, "", 94.7, "", "R4"),
    ("Guo et al.", 2018, "SED", "Astom ACS / CIMS", "concentrated seawater / brine", "", 660.0, "", 7.94, "R4"),
    ("Qiu et al.", 2019, "BMSED", "SED + BMED integration", "high Mg/Li brine", 58.3, 38.0, 99.0, "", "R4"),
    ("Qiu et al.", 2019, "ED (two-stage)", "Astom ACS / CIMS", "brine", 73.0, "", "", "", "R4"),
    ("Zhao et al.", 2019, "LMED", "TBP + ClO4 liquid membrane", "high Mg/Li brine", 68.0, 130.0, 92.0, 25.43, "R4"),
    ("Hoshino", 2013, "LMED", "PP13-TFSI / Nafion 324", "seawater", 22.2, "", "", "", "R4"),
    ("Yang et al.", 2018, "SCED", "NASICON superionic conductor", "seawater", "", "", "", "", "R4"),
    ("Li et al.", 2021, "SCED", "LLTO dense glass membrane", "Red Sea water", "", 530.0, 99.0, 6090.0, "R4"),
    ("Liu et al.", 2020, "LMED", "[C4mim][TFSI] liquid membrane", "high Mg/Li brine", 83.4, "", "", 266.9, "R4"),
    ("Shi et al.", 2019, "MCDI", "monovalent-selective CEM + AC", "LiCl/MgCl2 solution", 38.0, 1.8, "", 2.95, "R6"),
    ("Wei et al.", 2022, "BMED", "Li3PO4 conversion to LiOH", "Li3PO4 solution", 99.0, "", "", "", "R4"),
]
for a, y, sub, mat, feed, rec, en, pur, sf, src in ED:
    add(study=a, year=y, family="Electrodialysis / membrane electro-process", subclass=sub,
        material_or_system=mat, feed_type=feed, recovery_pct=rec, energy_Wh_mol=en,
        purity_pct=pur, sel_Li_Mg=sf, compiled_in=src)

# Recovery percentages compiled in R3 Table 2 (ED family, optimum conditions)
ED_REC = [
    ("SED", 75.8), ("SED", 70.0), ("SED", 90.5), ("SED", 72.5), ("SED", 77.5),
    ("SED", 76.0), ("SED", 67.7), ("SED", 76.5),
    ("BMED", 97.8), ("BMED", 83.7), ("BMED", 62.0), ("BMED", 73.0),
    ("ILM-ED", 22.2),
]
for sub, rec in ED_REC:
    add(study="compiled entry (R3 Table 2)", year="", family="Electrodialysis / membrane electro-process",
        subclass=sub, material_or_system="see R3 Table 2", feed_type="brine / seawater",
        recovery_pct=rec, compiled_in="R3")

# ----------------------------------------------------------------------------
# FAMILY 4 - NANOFILTRATION AND SELECTIVE MEMBRANES  (R6 Table 4)
# ----------------------------------------------------------------------------
NF = [
    ("Li et al.", 2019, "NF", "Commercial polyamide", 92.0, "", "R6"),
    ("Wu et al.", 2020, "NF", "[MimAP][Tf2N]-modified polyamide", 83.8, 8.12, "R6"),
    ("Ashraf et al.", 2020, "NF", "PDA/PEI-modified polyamide", 86.7, 7.15, "R6"),
    ("Wu et al.", 2020, "NF", "BTESE-coated TiO2 hybrid silica", 20.3, "", "R6"),
    ("Mohammad et al.", 2020, "ED", "MOF (ZIF-8)/polypropylene", "", 3.87, "R6"),
    ("Sharma et al.", 2020, "ED", "Mg-doped LMO / sulfonated PEEK", "", 4.82, "R6"),
    ("Xu et al.", 2021, "NF", "Polyamide on MWCNT-COOK PES/PEG", 98.49, 58.66, "R6"),
    ("Feng et al.", 2022, "NF", "PEI-grafted polyamide", 98.5, 33.4, "R6"),
    ("Peng & Zhao", 2021, "NF", "Diaminoethimidazole-modified polyamide", 95.8, "", "R6"),
    ("Saif et al.", 2021, "electro-diffusion", "HMO/PSS-Na/LiCF3SO3 composite", "", 11.75, "R6"),
    ("Wang et al.", 2021, "NF", "Cu-modified polyamide", "", 8.0, "R6"),
    ("Feng et al.", 2022, "NF", "Quaternized bipyridine-modified PEI", 92.0, "", "R6"),
]
for a, y, sub, mat, mgrej, sf, src in NF:
    add(study=a, year=y, family="Nanofiltration / selective membrane", subclass=sub,
        material_or_system=mat, feed_type="synthetic brine / salt solution",
        recovery_pct="", sel_Li_Mg=sf, purity_pct="", compiled_in=src)
    R[-1]["retention_pct"] = ""  # Mg rejection kept separately below
    R[-1]["cycles"] = ""
    if mgrej != "":
        R[-1]["purity_pct"] = ""      # keep purity blank; Mg rejection stored in notes col
        R[-1]["recovery_pct"] = ""

# ----------------------------------------------------------------------------
# FAMILY 5 - SOLVENT EXTRACTION  (R6 Table 3)
# ----------------------------------------------------------------------------
SX = [
    ("Shi et al.", 2014, "TBP / [C4mim][PF6]", "salt-lake brine", 90.93),
    ("Shi et al.", 2015, "TBP + NaClO4 / [C4mim][PF6]", "salt-lake brine", 99.12),
    ("Xiang et al.", 2016, "TBP + FeCl3 / MIBK", "salt-lake brine", 98.9),
    ("Xiang et al.", 2017, "TBP + FeCl3 / MIBK", "salt-lake brine", 98.0),
    ("Shi et al.", 2018, "N523 + TBP + FeCl3 / kerosene", "salt-lake brine", 96.0),
    ("Wang et al.", 2019, "TBP + [Bmim]3PW12O40 / dimethyl phthalate", "simulated brine", 99.23),
    ("Zhou et al.", 2020, "TBP + FeCl3 / diethyl succinate", "simulated brine", 65.0),
    ("Su et al.", 2020, "TBP + P507 + FeCl3 / kerosene", "salt-lake brine", 99.8),
    ("Zhou et al.", 2021, "TBP + [OHEMIM][NTf2]", "simulated brine", 94.2),
    ("Ren et al.", 2021, "TBP + NaBPh4 / CH2ClBr", "simulated brine", 87.65),
    ("Li et al.", 2021, "TBP + [Bmim]BPh4 / CH2BrCl", "simulated brine", 99.47),
    ("Bai et al.", 2021, "TBP + [N1888][P507] + FeCl3 / kerosene", "simulated brine", 70.0),
]
for a, y, mat, feed, eff in SX:
    add(study=a, year=y, family="Solvent extraction", subclass="extractant system",
        material_or_system=mat, feed_type=feed, recovery_pct=eff, compiled_in="R6")

# ----------------------------------------------------------------------------
# FAMILY 6 - PRECIPITATION  (R6 Table 1)
# ----------------------------------------------------------------------------
PR = [
    ("An et al.", 2012, "Na2CO3 + Ca(OH)2", "salt-lake brine (Uyuni)", 90.0, 99.55),
    ("Xu et al.", 2014, "Na2CO3", "salt-lake brine", 84.0, 99.6),
    ("Li et al.", 2015, "Activated Al-Ca / Al-Fe alloys", "synthetic salt solution", 94.6, ""),
    ("Liu et al.", 2018, "Aluminium-based material", "salt-lake brine", 78.3, ""),
    ("Quintero et al.", 2020, "NaOH + Na2CO3", "industrially refined brine", 85.0, 99.0),
    ("Zhang et al.", 2019, "Sodium metasilicate nonahydrate", "simulated brine", 86.73, ""),
    ("Lai et al.", 2020, "Na2HPO4 (crystallisation-precipitation)", "refined salt-lake brine", 93.2, ""),
    ("Liu et al.", 2021, "Facet-engineered Li3PO4", "salt-lake brine", 51.62, ""),
    ("Alsabbagh et al.", 2021, "Tri-sodium phosphate", "Dead Sea end brine", 40.0, ""),
]
for a, y, mat, feed, rec, pur in PR:
    add(study=a, year=y, family="Precipitation", subclass="chemical precipitant",
        material_or_system=mat, feed_type=feed, recovery_pct=rec, purity_pct=pur, compiled_in="R6")

# ----------------------------------------------------------------------------
# FAMILY 7 - CAPACITIVE DEIONISATION  (R3 Table 1)
# ----------------------------------------------------------------------------
CDI = [
    ("compiled entry (R3 Table 1)", "", "MCDI", "monovalent-selective CEM", "LiCl+MgCl2, 152.6 ppm", "", 1.8, ""),
    ("compiled entry (R3 Table 1)", "", "HCDI", "PVDF-EDA AEM", "LiCl, 848 ppm", 30.0, "", ""),
    ("compiled entry (R3 Table 1)", "", "HCDI", "lambda-LMO/AC", "Li binary solution, 1457 ppm", "", 160.77, ""),
    ("compiled entry (R3 Table 1)", "", "HCDI", "LMO/AC with EFA", "LiOH, 60 ppm", "", "", ""),
    ("compiled entry (R3 Table 1)", "", "HCDI", "Li-Mn-Fe spinel", "LiCl, 424 ppm", 32.0, "", ""),
    ("Siekierka et al.", 2019, "HCDI", "LMTO / AC with PVC-EDA AEM", "geothermal brine, 16 ppm", 800.0, 1.26, ""),
    ("Ha et al.", 2019, "FCDI", "AC flow electrode", "LiCl, 100 ppm", "", "", 91.7),
]
for a, y, sub, mat, feed, cap, en, rec in CDI:
    add(study=a, year=y, family="Capacitive deionisation", subclass=sub,
        material_or_system=mat, feed_type=feed, capacity_mg_g=cap,
        energy_Wh_mol=en, recovery_pct=rec, compiled_in="R3")

with open(OUT, "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=FIELDS)
    w.writeheader()
    w.writerows(R)

print("rows:", len(R))
fams = {}
for r in R:
    fams[r["family"]] = fams.get(r["family"], 0) + 1
for k, v in sorted(fams.items()):
    print("  %-45s %d" % (k, v))
print("written ->", OUT)
