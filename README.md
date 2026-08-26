# Positive-Case Count and Synthetic Data Fidelity: A CART-Based Analysis of MIMIC-IV
[![Full Report](https://img.shields.io/badge/📄_Read_Full_Report-PDF-blue?style=for-the-badge)](report/mimic_fidelity_report.pdf)

CART-based synthetic data generation (`synthpop`) on MIMIC-IV (n=70,954 ICU stays), testing whether fidelity inside a rare stratum is governed by the stratum's positive count or by the size of the table it sits in. Two sweeps hold the generator fixed and vary only the input: one varies the positive count `k` at fixed sample size, the other varies sample size `n` at fixed `k`.

Four independent measurements follow the count and none follows sample size. Divergence reaches exactly zero at `k = 80` in all six series; holding `k = 20` and raising `n` sixteen-fold never brings it near zero. The same contrast shows up in estimate spread, in whether CART spends a split on the variable at all, and in a plain logistic fit. But the mechanism does not rescue the analysis: at MIMIC-IV's effect sizes a 500-row subsample carries 6–8% power, and no count reachable inside one makes the estimate usable.

[See key results →](#key-results)

## Motivation

Synthetic data is usually evaluated with global fidelity metrics: does the synthetic table look like the real one overall? Those metrics are dominated by the common strata, which is exactly where nothing interesting happens. Rare subgroups carry little weight in the aggregate, yet they are often what clinical research is about, and they are where a generator is most likely to lose structure.

The harder problem is attribution. If a rare-stratum effect goes missing after synthesis, is that the generator's failure or was the effect never estimable from that many rows in the first place? Most evaluations cannot tell, because they never establish what the source data could support. This project builds the control in: every synthetic estimate is paired with a real one from the same rows, and the whole grid is read against a full-pool reference and a power calculation.

The answer is that a prior finding transfers — positive count, not sample size, sets the threshold — and that it transfers as a mechanism, not as a remedy.

## Design Decisions

**Why two sweeps?**

At fixed `n`, raising the positive count also raises prevalence, so a single sweep cannot say which one is responsible. Fixing `k` and varying `n` moves prevalence on its own. `k` was held at 20 because that is the last count where fits still fail often enough to matter (6–17% of replicates, against 0–4% at `k = 40` and none from `k = 80` up) — sample size gets its best chance to show an effect there.

**Why report failure separately from spread?**

A fit fails two ways: **inestimable** (the table holds no positives, the coefficient comes back `NA`) and **divergent** (the coefficient runs past |2|, a likelihood maximized at the boundary rather than a real effect). Averaging failures into the quantiles inflates apparent instability at low counts, since divergent fits push the quantiles outward, and it hides the sharpest result in the study — the common threshold at `k = 80`. Both are reported as shares of all `M` replicates, so converged + inestimable + divergent = 1.

**Why a permutation null for the tree diagnostic?**

"Does CART use the comorbidity as a split?" looks like a clean fidelity metric. It is not readable without a control: the same measurement on a shuffled column rises from 0.00 to 0.32 as class balance improves, because a balanced binary variable is structurally easy to split on whether or not it carries signal. Every configuration is re-run with the comorbidity column shuffled, which preserves the positive count and destroys the association.

**Why a GLM arm alongside CART?**

To separate "CART cannot see this" from "nobody can see this." A likelihood-ratio test on the identical subsamples answers the same yes/no question the tree is asked. It detects the comorbidity in 2–16% of replicates throughout, no consistent gain over shuffled labels — so the null result is a property of the design, not of the generator.

## Key Results

**Full-pool reference (n=70,954)**

| Comorbidity | Creatinine | Interaction | SE | p |
|---|---|---|---|---|
| CHF | 0.208 | −0.0497 | 0.0128 | 1.1 × 10⁻⁴ |
| Diabetes | 0.234 | −0.0709 | 0.0126 | 2.0 × 10⁻⁸ |
| COPD | 0.188 | +0.0668 | 0.0197 | 6.9 × 10⁻⁴ |

All three interactions are real and all three are small — roughly three to four times smaller than the creatinine main effect they modify.

**The two levers, matched at roughly eightfold**

| Measurement | count: k 20 → 150 | size: n 500 → 4000 |
|---|---|---|
| Divergence | 6–17% → **exactly zero** (all six series) | −29% to +25%, none near zero |
| Interquartile width | **−42% to −66%** (all six series) | −19% to +10%, no consistent sign |
| CART uses the comorbidity | 26–40% at k=250 | 4–14% at every n |
| GLM detection | — | 2–16%, no trend in n |

Two of the four never touch the synthesizer.

**The sharpest single comparison**

At `k = 5` and 1% prevalence the synthetic table loses the stratum 4–6% of the time. At `k = 20` and the same 1% prevalence (n = 2000) it never does. Prevalence held equal, count differs, outcome differs.

**Sign recovery (COPD, the only positive interaction, +0.0668)**

| | synthetic median |
|---|---|
| k = 20 (Experiment 1) | −0.177 |
| k = 20, n = 250…4000 (Experiment 2) | −0.077 to −0.159 |
| k = 40 (Experiment 1) | **+0.039** |

Raising the count fixes the sign. Raising `n` sixteen-fold does not.

**Tree diagnostic against its null**

| | observed | permutation null |
|---|---|---|
| Count sweep, k 5 → 250 | 0.00 → 0.26–0.40 | 0.00 → **0.32** |
| Size sweep, n 250 → 4000 | 0.04–0.14 throughout | 0.02–0.08 → **0.00–0.03** |

Read without the null, the count sweep looks like recovery. The null rises to meet it. Driving class balance the other way in the size sweep sends the null to zero, which confirms it tracks balance and not the count — and only then is the signal legible (104 discordant replicates to 32, against 249 to 164 in the count sweep).

**What the design could never resolve**

A 500-row subsample carries 6–8% power against these interactions; the GLM arm returned 6.7–10.1%. Theory and experiment agreeing on a null result is what rules out a defect in the pipeline. 80% power would need 17,686 to 48,357 rows. Even the best cell in the study — CHF real arm at `k = 250`, zero failures — has an interquartile width of 0.194 against a reference of 0.0497, nearly four times the quantity being estimated.

## Reflections & Next Steps

The count-not-size mechanism replicates from Bayesian logistic regression into CART-based synthesis, and it reaches further than expected: it holds not just for the output estimate but for the generator's internal behaviour, whether the tree spends a split on the variable at all. That is the part worth carrying forward.

The more transferable result may be the methodological one. Variable-selection diagnostics of the form "does the generator use variable X" cannot be read without a permutation control, and they are only legible where that null is near zero, which on this data means low prevalence. This is a known CART bias placed in a setting where it silently produces the wrong conclusion.

Next steps:
- **A target the design can resolve**: MIMIC-IV's interactions are 0.050–0.071 and need tens of thousands of rows. Repeating on an effect large enough to be estimable at n=500 would separate "the count is binding" from "nothing is estimable here."
- **More than one generator**: the count threshold may be specific to `syn.cart`'s `minbucket = 5`. Sweeping that parameter, or comparing against a parametric synthesizer, would test it.
- **The privacy side**: the report argues that CART's reluctance to isolate a five-patient node is the same reluctance that protects those patients. Nothing here measures the privacy half of that trade-off.

## Repository

```
code/
  ├── mimic-fidelity-analysis.ipynb   # Experiments 1 and 2
  ├── tree-structure-diagnosis.ipynb  # CART probe, permutation null, power analysis
  ├── sql/                            # MIMIC-IV queries, loaded via read_sql()
  └── config.R                        # local DATA_DIR (not tracked)
figures/
  └── *.png                           # report figures
report/
  ├── mimic_fidelity_report.tex       # full writeup
  └── proposal.tex                    # original proposal
```

MIMIC-IV is not redistributable. `config.R` points `DATA_DIR` at a local copy obtained through PhysioNet credentialed access; the notebooks read the `.csv.gz` files directly with DuckDB, so no database server is needed.

## Tools

**Statistical methods**: CART-based synthesis, Monte Carlo replication (M=100), permutation testing, exact paired binomial tests, power analysis
**Language**: R
**Libraries**: synthpop, rpart, duckdb, DBI, dplyr, tidyr, ggplot2
**Data**: MIMIC-IV v2.2 (PhysioNet, credentialed access)

## References

Johnson, A.E.W., Bulgarelli, L., Shen, L., et al. (2023). MIMIC-IV, a freely accessible electronic health record dataset. *Scientific Data*, 10, 1.

Nowok, B., Raab, G.M., & Dibben, C. (2016). synthpop: Bespoke creation of synthetic data in R. *Journal of Statistical Software*, 74(11), 1–26.

Heinze, G., & Schemper, M. (2002). A solution to the problem of separation in logistic regression. *Statistics in Medicine*, 21(16), 2409–2419.

Strobl, C., Boulesteix, A.-L., Zeileis, A., & Hothorn, T. (2007). Bias in random forest variable importance measures. *BMC Bioinformatics*, 8, 25.

Prior work this builds on: [Bayesian Prior Sensitivity](https://github.com/ShengPeiWilliam/bayesian-prior-sensitivity) — the positive-count finding in Bayesian logistic regression that this project tests against a different generator.
