// Table data for the lithium recovery review and meta-analysis.
// Each table records the source review from which the values were transcribed.

module.exports = {

  // ---------------------------------------------------------------- Table 1
  brines: {
    title: "Table 1. Major cation composition of representative lithium-bearing brines (g L-1).",
    note: "Transcribed from Wu et al. (2022), Table 1, with Chinese salt-lake ranges from Liu et al. (2019), Table 1. Blank cells indicate values not reported in the source.",
    head: ["Resource", "Li+", "Mg2+", "Na+", "Ca2+", "K+", "SO4 2-"],
    widths: [2900, 1150, 1150, 1250, 1150, 1150, 1290],
    rows: [
      ["Atacama Salar, Chile", "3.02", "17.60", "61.90", "0.41", "28.20", "37.90"],
      ["Uyuni Salar, Bolivia", "0.84", "16.70", "105.40", "3.33", "15.70", "21.30"],
      ["East Taijinar, China", "0.14", "5.64", "117.03", "-", "0.43", "3.79"],
      ["West Taijinar, China", "0.26", "15.36", "102.40", "-", "0.19", "8.44"],
      ["Yiliping, China", "0.32", "23.20", "80.10", "0.15", "11.89", "15.27"],
      ["Chott Djerid, Tunisia", "0.06", "3.40", "80.00", "1.60", "5.60", "6.70"],
      ["North Arm, Great Salt Lake, USA", "0.04", "9.38", "100.80", "0.35", "5.50", "19.70"],
      ["Wairakei geothermal brine, NZ", "0.01", "-", "0.11", "-", "0.16", "-"],
      ["Qarhan, China (range)", "0.21-0.35", "Mg/Li 184-547", "-", "-", "-", "-"],
      ["Da Qaidam, China (range)", "0.10-1.30", "Mg/Li 6.9-1170", "-", "-", "-", "-"],
      ["Zabuye, China (range)", "0.42-1.61", "Mg/Li 0.024", "-", "-", "-", "-"],
      ["Seawater (typical)", "0.00017-0.00021", "1.29", "10.8", "0.41", "0.39", "2.71"],
    ],
  },

  // ---------------------------------------------------------------- Table 2
  precipitation: {
    title: "Table 2. Lithium recovery by precipitation reported in recent primary studies.",
    note: "Transcribed from Khalil et al. (2022), Table 1.",
    head: ["Year", "Feed", "Reagent / precipitant", "T (C)", "Recovery (%)", "Purity (%)"],
    widths: [800, 2400, 2900, 900, 1300, 1700],
    rows: [
      ["2012", "Salt-lake brine (Uyuni)", "Na2CO3 + Ca(OH)2", "90", ">90", "99.55"],
      ["2014", "Salt-lake brine", "Na2CO3", "80", "84", "99.6"],
      ["2015", "Synthetic salt solution", "Activated Al-Ca and Al-Fe alloys", "70", "94.6", "-"],
      ["2018", "Salt-lake brine", "Aluminium-based material", "90", "78.3", "-"],
      ["2019", "Industrially refined brine", "NaOH + Na2CO3", "60", ">85", "99"],
      ["2019", "Simulated brine", "Sodium metasilicate nonahydrate", "25", "86.73", "-"],
      ["2020", "Refined salt-lake brine", "Na2HPO4", "40", "93.2", "-"],
      ["2021", "Salt-lake brine", "Facet-engineered Li3PO4", "30", "51.62", "-"],
      ["2021", "Seawater brine", "Tri-sodium phosphate", "40", "40", "-"],
    ],
  },

  // ---------------------------------------------------------------- Table 3
  adsorption: {
    title: "Table 3. Lithium uptake capacity of adsorbents and ion sieves.",
    note: "Rows 1-19 transcribed from Khalil et al. (2022), Table 2; rows 20-23 from statements in the text of Liu et al. (2019).",
    head: ["Year", "Feed", "Adsorbent", "Capacity (mg g-1)"],
    widths: [900, 3100, 4200, 1800],
    rows: [
      ["2014", "Salt-lake brine", "Layered H2TiO3", "32.7"],
      ["2017", "Salt-lake brine", "HxTiO3", "36.3"],
      ["2018", "Concentrated natural brine", "Iron-doped Li-Ti oxide", "34.8"],
      ["2018", "Enriched salt-lake brine", "MnO2/PVDF composite", "11.06"],
      ["2019", "Natural brine", "Li1.33Mn1.67O4 (surface-modified)", ">50"],
      ["2019", "Salt-lake water", "LiMn2O4 nanowire / diatomaceous earth", "18.6"],
      ["2019", "Synthetic salt solution", "Aluminium-doped Li manganese oxide", "32.6"],
      ["2020", "Synthetic salt solution", "Li+-IIP-Fe3O4@C", "22.26"],
      ["2020", "Salt-lake brine", "Granular H4Mn5O12 (EP/HMO)", "17.2"],
      ["2020", "Salt-lake brine", "Magnetic Li-Al layered double hydroxide", "5.83"],
      ["2020", "Simulated brine", "Porous H4Ti5O12 nanofibre", "59.1"],
      ["2021", "Synthetic salt solution", "Al-based MOF (MIL-121)", "1.25"],
      ["2021", "Salt-lake brine", "H2TiO3 / diatomaceous earth", "8.25"],
      ["2021", "Salt-lake brine", "Li-Al layered double hydroxide", "7.27"],
      ["2021", "Salt-lake brine", "Granular Li-Al-LDH / NH4Al3(SO4)2(OH)6", "9.16"],
      ["2021", "Synthetic salt solution", "Layered H2TiO3", "40"],
      ["2021", "Salt-lake brine", "Fluorine-pillared MOF", "18.8"],
      ["2021", "Salt-lake brine (W. Taijinar)", "Shaped Ti-based Li ion sieve", "19.22"],
      ["2021", "Synthetic salt solution", "PNIPAM-co-crown-6 MOF-808", "7.21"],
      ["2014", "Synthetic solution", "LTO from anatase precursor", "39.2"],
      ["2015", "Synthetic solution", "H2TiO3", "57.8"],
      ["2016", "Synthetic solution", "H2TiO3", "94.5"],
      ["2017", "Brine", "PVDF / crown-ether imprinted membrane", "27.1"],
    ],
  },

  // ---------------------------------------------------------------- Table 4
  solvent: {
    title: "Table 4. Lithium extraction efficiency of solvent extraction systems.",
    head: ["Year", "Feed", "Extractant", "Co-extractant", "Diluent", "Efficiency (%)"],
    note: "Transcribed from Khalil et al. (2022), Table 3.",
    widths: [800, 2000, 1300, 2100, 2200, 1600],
    rows: [
      ["2014", "Salt-lake brine", "TBP", "-", "[C4mim][PF6]", "90.93"],
      ["2015", "Salt-lake brine", "TBP", "NaClO4", "[C4mim][PF6]", "99.12"],
      ["2016", "Salt-lake brine", "TBP", "FeCl3", "MIBK", "98.9"],
      ["2017", "Salt-lake brine", "TBP", "FeCl3", "MIBK", ">98"],
      ["2018", "Salt-lake brine", "N523, TBP", "FeCl3", "Kerosene", "96"],
      ["2019", "Simulated brine", "TBP", "[Bmim]3PW12O40", "Dimethyl phthalate", "99.23"],
      ["2020", "Simulated brine", "TBP", "FeCl3", "Diethyl succinate", "65"],
      ["2020", "Salt-lake brine", "TBP", "P507, FeCl3", "Kerosene", "99.8"],
      ["2021", "Simulated brine", "TBP", "[OHEMIM][NTf2]", "-", "94.2"],
      ["2021", "Simulated brine", "TBP", "NaBPh4", "CH2ClBr", "87.65"],
      ["2021", "Simulated brine", "TBP", "[Bmim]BPh4", "CH2BrCl", "99.47"],
      ["2021", "Simulated brine", "TBP", "[N1888][P507], FeCl3", "Kerosene", ">70"],
    ],
  },

  // ---------------------------------------------------------------- Table 5
  membranes: {
    title: "Table 5. Performance of nanofiltration and selective membranes for Li/Mg separation.",
    note: "Transcribed from Khalil et al. (2022), Table 4. LMH = L m-2 h-1.",
    head: ["Year", "Mode", "Membrane material", "Mg2+ rejection (%)", "Li/Mg selectivity"],
    widths: [800, 1500, 4200, 1900, 1600],
    rows: [
      ["2019", "NF", "Commercial polyamide", "92", "-"],
      ["2020", "NF", "[MimAP][Tf2N]-modified polyamide", "83.8", "8.12"],
      ["2020", "NF", "PDA/PEI-modified polyamide", "86.7", "7.15"],
      ["2020", "NF", "BTESE-coated TiO2 hybrid silica", "20.3", "-"],
      ["2020", "ED", "MOF (ZIF-8) / polypropylene", "-", "3.87"],
      ["2020", "ED", "Mg-doped LMO / sulfonated PEEK", "-", "4.82"],
      ["2021", "NF", "Polyamide on MWCNT-COOK PES/PEG", "98.49", "58.66"],
      ["2021", "NF", "PEI-grafted polyamide", "98.5", "33.4"],
      ["2021", "NF", "Diaminoethimidazole-modified polyamide", "95.8", "-"],
      ["2021", "Electro-diffusion", "HMO/PSS-Na/LiCF3SO3 composite", "-", "11.75"],
      ["2021", "NF", "Cu-modified polyamide", "-", "8.0"],
      ["2022", "NF", "Quaternized bipyridine-modified PEI", "92", "-"],
    ],
  },

  // ---------------------------------------------------------------- Table 6
  electro: {
    title: "Table 6. Performance of electrosorption and electrochemical ion-pumping systems.",
    note: "Compiled from Khalil et al. (2022), Table 5; Zhao et al. (2019), Table 1; and Wu et al. (2022), Table 2. Energy is expressed per mole of lithium recovered. Retention is given as percentage of initial capacity after the stated cycle count.",
    head: ["Year", "Electrode pair (counter // capture)", "Capacity (mg g-1)", "Energy (Wh mol-1)", "Purity (%)", "Retention / cycles"],
    widths: [700, 3500, 1500, 1500, 1200, 1800],
    rows: [
      ["1993", "Pt // lambda-MnO2", "-", "11.0", "-", "-"],
      ["2012", "Ag // LiFePO4", "7.76", "1.0", "84.6", "-"],
      ["2013", "Ag // lambda-MnO2", "-", "4.5", "91.8", "87 / 20"],
      ["2014", "Ag // LiFePO4", "-", "2.8", "97.9", "-"],
      ["2015", "NiHCF // LiFePO4", "-", "8.7", "97.9", "-"],
      ["2015", "AC // lambda-MnO2", "14.88", "4.2", "91.8", "96 / 50"],
      ["2016", "PPy // LiMn2O4", "1.73", "5-10", "-", "50 / 200"],
      ["2017", "LiMn2O4 // Li(1-x)Mn2O4 (brine)", "22.0", "18.0", "89", "97"],
      ["2017", "LiMn2O4 // Li(1-x)Mn2O4 (seawater)", "21.0", "18.6", "86", "97"],
      ["2017", "NiHCF // lambda-MnO2", "-", "3.58", "96.2", "83.3 / 100"],
      ["2018", "LiMn2O4 // Zn", "-", "6.3", "-", "73 / 100"],
      ["2018", "Ag // LNCM", "10.83", "2.60", "96.4", "83.3 / 50"],
      ["2018", "Graphite // lambda-MnO2 film", "-", "4.14", "99.0", "91 / 100"],
      ["2019", "Ag // lambda-MnO2", "10.1", "3.07", "99.0", "-"],
      ["2019", "Ag // LiNi0.5Mn1.5O4", "8.74", "1.32", "98.14", "-"],
      ["2019", "PANI // lambda-MnO2", "23.13", "3.95", "81.4", "70.8 / 200"],
      ["2019", "PPy-threaded HKUST-1", "37.55", "-", "-", "-"],
      ["2020", "Ag // lambda-MnO2 (seawater concentrate)", "-", "-", "88.0", "-"],
      ["2020", "AC // LNCM/rGO", "13.84", "1.4", "93.5", "88.9 / 30"],
      ["2020", "AC // CF-NMMO", "14.4", "7.91", "97.2", "90.8 / 30"],
      ["2020", "AC // PPy/Al2O3/LiMn2O4", "12.84", "1.41", "97.37", "-"],
      ["2021", "BiOCl@PPy // lambda-MnO2", "10.88", "1.007", "-", "96 / 100"],
      ["2021", "LiMn2O4 // Li(1-x)Mn2O4", "15-16", "16.0", "97.0", "- / 100"],
      ["2021", "lambda-MnO2/LiMn2O4 on graphite felt", "26.14", "23.38", "-", "89.6 / 20"],
    ],
  },

  // ---------------------------------------------------------------- Table 7
  ed: {
    title: "Table 7. Performance of electrodialysis and related membrane electro-processes.",
    note: "Compiled from Wu et al. (2022), Table 3 and Zavahir et al. (2021), Table 2. Energy values reported in kWh mol-1 in the source have been converted to Wh mol-1.",
    head: ["Year", "Mode", "Membrane / liquid phase", "Recovery (%)", "Energy (Wh mol-1)", "Li/Mg SF"],
    widths: [700, 1300, 4000, 1400, 1700, 1100],
    rows: [
      ["2013", "LMED", "PP13-TFSI / Nafion 324", "22.2", "-", "-"],
      ["2014", "BMED", "Neosepta CMX / BP-1", "95", "8", "-"],
      ["2015", "ED", "AR204SXR412 / CR67-MK111", "20", "-", "-"],
      ["2016", "BMED", "Neosepta BP-1", "70", "-", "-"],
      ["2017", "SED", "Selemion ASA / CSO", "62.5", "131", "20.2"],
      ["2017", "SED", "Selemion ASA / CSO", "55.2", "31", "9.89"],
      ["2017", "SED", "Astom ACS / CIMS", "75.4", "-", "8.57"],
      ["2017", "BMED", "Astom CMB/CHA + Neosepta BP-1E", "20", "-", "-"],
      ["2018", "SED", "Astom ACS / CIMS", "-", "660", "7.94"],
      ["2018", "BMED", "Astom CMB/CHA + Neosepta BP-1E", "62", "-", "-"],
      ["2018", "SCED", "NASICON superionic conductor", "-", "-", "-"],
      ["2019", "BMSED", "Integrated SED + BMED", "58.3", "38", "-"],
      ["2019", "ED (two-stage)", "Astom ACS / CIMS", "73", "-", "-"],
      ["2019", "LMED", "TBP + ClO4- liquid membrane", "68", "130", "25.4"],
      ["2019", "MCDI", "Monovalent-selective CEM + AC", "38", "1.8", "2.95"],
      ["2020", "LMED", "[C4mim][TFSI] liquid membrane", "83.4", "-", "266.9"],
      ["2021", "SCED", "LLTO dense glass membrane", "-", "530", "6090"],
      ["2022", "BMED", "Li3PO4 to LiOH conversion", "99", "-", "-"],
    ],
  },

  // ---------------------------------------------------------------- Table 8
  cdi: {
    title: "Table 8. Capacitive deionisation systems investigated for lithium capture.",
    note: "Transcribed from Zavahir et al. (2021), Table 1. Recovery ratio is expressed per gram of electrode or adsorbent material.",
    head: ["Mode", "Electrode / membrane", "Feed (ppm Li)", "Voltage (V)", "Recovery ratio (mg g-1)", "SEC (Wh mol-1)"],
    widths: [1100, 2900, 1800, 1200, 1900, 1300],
    rows: [
      ["MCDI", "Monovalent-selective CEM", "152.6 (LiCl+MgCl2)", "0.6-1.4", "-", "1.8"],
      ["HCDI", "PVDF-EDA anion exchange membrane", "848 (LiCl)", "1.23", "30", "-"],
      ["HCDI", "lambda-LMO / activated carbon", "1457", "1.0", "-", "160.77"],
      ["HCDI", "LMO/AC, electrostatic field assist", "60 (LiOH)", "3.5", "-", "-"],
      ["HCDI", "LMO/AC", "50 (LiOH)", "0.1-2.0", "-", "-"],
      ["HCDI", "Li-Mn-Fe spinel", "424 (LiCl)", "1.23", "32", "-"],
      ["HCDI", "LMTO / AC with PVC-EDA AEM", "16 (geothermal)", "2.0", "800", "1.26"],
      ["FCDI", "Activated carbon flow electrode", "100 (LiCl)", "1.2", "-", "-"],
    ],
  },

  // ---------------------------------------------------------------- Table 9
  capacity: {
    title: "Table 9. Pooled lithium uptake capacity by material family (mg g-1).",
    note: "Bootstrap 95% confidence intervals on the median use 10 000 resamples.",
    head: ["Material family", "n", "Median", "IQR", "Mean", "SD", "Range", "95% CI (median)"],
    widths: [2800, 600, 1000, 1600, 1000, 900, 1700, 1600],
    rows: [
      ["All adsorbents / ion sieves", "23", "22.26", "10.11-37.75", "28.27", "22.04", "1.25-94.5", "17.2-36.3"],
      ["  Titanium oxide sieves", "6", "48.90", "39.40-58.77", "51.64", "25.58", "19.2-94.5", "-"],
      ["  Mn/Ti mixed oxide sieves", "3", "34.80", "33.75-35.55", "34.60", "1.81", "32.7-36.3", "-"],
      ["  Manganese oxide sieves", "2", "41.30", "36.95-45.65", "41.30", "12.30", "32.6-50.0", "-"],
      ["  Composite sieves", "3", "11.06", "9.66-14.83", "12.64", "5.35", "8.25-18.6", "-"],
      ["  Layered double hydroxides", "3", "7.27", "6.55-8.21", "7.42", "1.67", "5.83-9.16", "-"],
      ["  Metal-organic frameworks", "3", "7.21", "4.23-13.01", "9.09", "8.92", "1.25-18.8", "-"],
      ["Hybrid adsorptive membranes", "4", "23.57", "22.31-60.31", "59.05", "72.98", "20.5-168.5", "-"],
      ["Electrosorption electrodes", "16", "14.12", "10.65-21.25", "15.71", "8.60", "1.73-37.6", "10.8-21.0"],
    ],
  },

  // ---------------------------------------------------------------- Table 10
  energy: {
    title: "Table 10. Pooled specific energy consumption by process family (Wh mol-1 Li).",
    note: "Electrosorption and electrodialysis differ significantly (Mann-Whitney U = 23.0, p = 0.0014).",
    head: ["Process family", "n", "Median", "IQR", "Mean", "SD", "Range"],
    widths: [3400, 700, 1200, 1900, 1200, 1300, 1500],
    rows: [
      ["Electrosorption / ion pumping", "22", "4.17", "2.65-8.50", "6.93", "6.50", "1.0-23.4"],
      ["Electrodialysis family", "8", "84.00", "25.25-230.75", "191.22", "256.49", "1.8-660"],
      ["Capacitive deionisation", "3", "1.80", "1.53-81.29", "54.61", "91.94", "1.26-160.8"],
    ],
  },

  // ---------------------------------------------------------------- Table 11
  recovery: {
    title: "Table 11. Pooled lithium recovery or extraction efficiency by technology family (%).",
    note: "Kruskal-Wallis across the three families: H = 15.22, df = 2, p = 0.00049.",
    head: ["Technology family", "n", "Median", "IQR", "Mean", "SD", "Range"],
    widths: [3400, 700, 1200, 1900, 1200, 1300, 1500],
    rows: [
      ["Precipitation", "9", "85.00", "78.30-90.00", "78.16", "19.21", "40-94.6"],
      ["Solvent extraction", "12", "97.00", "90.11-99.15", "91.52", "11.89", "65-99.8"],
      ["Electrodialysis family", "28", "71.25", "61.08-76.75", "65.97", "22.72", "20-99"],
    ],
  },

  // ---------------------------------------------------------------- Table 12
  purity: {
    title: "Table 12. Purity of the recovered lithium stream by process family (%).",
    note: "Purity is reported as the percentage of lithium among cations in the recovery solution or product, as defined in each source study.",
    head: ["Process family", "n", "Median", "IQR", "Mean", "SD", "Range"],
    widths: [3400, 700, 1200, 1900, 1200, 1300, 1500],
    rows: [
      ["Electrosorption / ion pumping", "18", "96.30", "89.70-97.74", "93.45", "5.51", "81.4-99.0"],
      ["Electrodialysis family", "6", "99.00", "95.78-99.00", "97.22", "3.12", "92.0-99.6"],
      ["Precipitation", "3", "99.55", "99.28-99.57", "99.38", "0.33", "99.0-99.6"],
    ],
  },

  // ---------------------------------------------------------------- Table 13
  selectivity: {
    title: "Table 13. Li/Mg separation factor by separation principle (dimensionless).",
    note: "Electrosorption and nanofiltration differ significantly (Mann-Whitney U = 59.0, p = 0.027). The electrodialysis mean is dominated by a single superionic-conductor result of 6090.",
    head: ["Separation principle", "n", "Median", "IQR", "Mean", "SD", "Range"],
    widths: [3400, 700, 1200, 1900, 1200, 1300, 1500],
    rows: [
      ["Electrosorption / ion pumping", "9", "45.58", "18.00-286.67", "177.25", "262.38", "5.57-804"],
      ["Nanofiltration / selective membrane", "8", "8.06", "6.57-17.16", "16.97", "19.31", "3.87-58.7"],
      ["Electrodialysis family", "8", "15.04", "8.41-85.80", "804.0", "2137.7", "2.95-6090"],
    ],
  },

  // ---------------------------------------------------------------- Table 14
  completeness: {
    title: "Table 14. Reporting completeness across the 131 extracted performance records.",
    note: "Percentages give the share of records for which the metric was reported in the source comparison table.",
    head: ["Performance metric", "Records reporting (n)", "Share of corpus (%)"],
    widths: [4600, 2600, 2600],
    rows: [
      ["Recovery / extraction efficiency", "50", "38.2"],
      ["Uptake or extraction capacity", "46", "35.1"],
      ["Specific energy consumption", "33", "25.2"],
      ["Product purity", "27", "20.6"],
      ["Li/Mg selectivity", "25", "19.1"],
      ["Cycle number tested", "14", "10.7"],
      ["Capacity retention", "14", "10.7"],
    ],
  },

  // ---------------------------------------------------------------- Table 15
  comparison: {
    title: "Table 15. Cross-cutting comparison of the seven technology families.",
    note: "Maturity classes follow the assessments given in Zavahir et al. (2021), Khalil et al. (2022) and Wu et al. (2022). TRL bands are indicative rather than formally assessed.",
    head: ["Family", "Typical capacity or recovery", "Energy", "Li/Mg selectivity", "Maturity", "Principal limitation"],
    widths: [1900, 2000, 1200, 1400, 1400, 2000],
    rows: [
      ["Evaporation + precipitation", "40-95% recovery", "Solar, low electrical", "Poor above Mg/Li 6", "Commercial", "12-24 month cycle; large water and land use"],
      ["Adsorption / ion sieve", "Median 22 mg g-1", "Low to moderate", "High intrinsic", "Pilot to early commercial", "Sorbent dissolution; slow kinetics; elution cost"],
      ["Solvent extraction", "Median 97% efficiency", "Moderate", "Moderate to high", "Pilot", "Organic inventory; acid stripping; solvent loss"],
      ["Nanofiltration", "80-98% Mg rejection", "Pressure-driven", "Median 8", "Pilot to commercial", "Single stage insufficient; fouling and scaling"],
      ["Electrosorption / ion pump", "Median 14 mg g-1", "Median 4.2 Wh mol-1", "Median 46", "Laboratory to pilot", "Electrode fade; counter-electrode cost; scale-up"],
      ["Electrodialysis family", "Median 71% recovery", "Median 84 Wh mol-1", "Median 15", "Bench to commercial", "Membrane cost and scaling; co-ion leakage"],
      ["Capacitive deionisation", "Variable", "1.3-161 Wh mol-1", "Low without sieve", "Laboratory", "Immature; inconsistent performance metrics"],
    ],
  },
};
