# Lithium Recovery: Systematic Review and Quantitative Meta-Analysis

A 63 page review manuscript synthesising the quantitative evidence on lithium recovery from
aqueous resources, together with the full analysis pipeline that produced it.

## Deliverable

`Lithium_Recovery_Review_and_Meta_Analysis.docx` — 63 pages, 15 tables, 7 figures,
283 references, plus a 131 row appendix reproducing the complete evidence table.

## Evidence base

Six full text review articles were read in their entirety and every comparison table they
contain was transcribed into a single harmonised evidence table:

| ID | Source review |
|----|---------------|
| R1 | Liu, Zhao & Ghahreman (2019), *Hydrometallurgy* 187, 81-100 |
| R2 | Wang et al. (2022), *Renewable and Sustainable Energy Reviews* 154, 111813 |
| R3 | Zavahir et al. (2021), *Desalination* 500, 114883 |
| R4 | Wu et al. (2022), *Water Research* 221, 118822 |
| R5 | Zhao et al. (2019), *Journal of Electroanalytical Chemistry* 850, 113389 |
| R6 | Khalil et al. (2022), *Desalination* 528, 115611 |

The resulting table holds 131 performance records spanning 8 technology families and roughly
120 distinct primary studies. Every row carries the identifier of the review it was
transcribed from, so any value traces back to its source document.

## Headline findings

- Electrosorption and electrochemical ion pumping consume a median of **4.17 Wh per mol Li**
  against **84.0 Wh per mol** for the electrodialysis family — a factor of 20
  (Mann-Whitney U = 23.0, p = 0.0014).
- Titanium oxide ion sieves reach the highest median capacity (48.9 mg/g), but their
  presumed capacity advantage over manganese sieves is **not** statistically supported
  (U = 30.5, p = 0.32). Their real advantage is stability.
- Reporting is the binding problem. No metric is reported by even 40% of records;
  durability by only 10.7%. Median capacity fade where reported is 0.26% per cycle.
- No evidence of performance improvement over time survives confounding
  (capacity vs year: rho = -0.32, p = 0.047).

## Reproducing the analysis

```bash
cd analysis
python3 build_dataset.py     # writes li_recovery_dataset.csv (131 records)
python3 analyse.py           # writes results.json and figures/
cd ../doc
npm install docx
node build.js                # writes the .docx
```

`analyse.py` uses a pooled descriptive meta-analysis rather than an inverse variance
random effects model: the source reviews report point estimates without dispersion, so
study level variances are not identifiable. It reports medians with bootstrap confidence
intervals, non parametric between family tests, and rank correlations for temporal trends.

## Layout

```
analysis/   build_dataset.py, analyse.py, li_recovery_dataset.csv, results.json
doc/        refs.js, tables.js, body1.js, body2.js, build.js
figures/    fig1..fig7 (PNG, 200 dpi)
```

## Provenance note

All pooled statistics derive exclusively from values transcribed from the six source
reviews. Work published after those reviews closed was identified by literature search and
is discussed qualitatively in Section 8; it is excluded from the quantitative synthesis
because full text could not be retrieved for verification. No performance value in the
manuscript was estimated or inferred — rows whose column assignment was ambiguous in the
source layout were excluded rather than guessed.
