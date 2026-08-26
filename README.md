# Positive-Case Count and Synthetic Data Fidelity: A CART-Based Analysis of MIMIC-IV
[![Full Report](https://img.shields.io/badge/📄_Read_Full_Report-PDF-blue?style=for-the-badge)](report/mimic_fidelity_report.pdf)

CART-based synthetic data generation (`synthpop`) on MIMIC-IV (n=70,954 ICU stays), asking whether fidelity inside a rare stratum is governed by the stratum's positive count or by the size of the table it sits in. Two sweeps hold the generator fixed and vary only the input: one varies the positive count `k` at fixed sample size, the other varies sample size `n` at fixed `k`.

**Every quantity the two sweeps measure follows the count and none follows sample size.** Divergence falls to zero between counts of 40 and 80 in all six series; holding `k = 20` and raising `n` sixteen-fold never brings it near zero. The generator inherits that threshold rather than adding one of its own — what it adds instead is a bias no count removes.

## Motivation

An earlier project, [Bayesian Prior Sensitivity](https://github.com/ShengPeiWilliam/bayesian-prior-sensitivity), found on a 189-row dataset that the point where estimates stop being driven by the prior is set by a predictor's positive-case count, not by sample size. That a logistic fit's precision is set by the event count is not new — it is the events-per-variable literature. What is not established is whether the same threshold governs a **generator**: whether a CART synthesizer preserves a rare-stratum interaction is a question about which variables a tree spends a split on, not about the curvature of a likelihood.

MIMIC-IV is close to the ideal place to ask, since a rare stratum can be built to any size on demand instead of being whatever the data happens to give.

## Key Results

**Full-pool reference** (n=70,954, SEs clustered on patient)

| Comorbidity | Creatinine | Interaction | SE | p |
|---|---|---|---|---|
| CHF | 0.208 | −0.0497 | 0.0141 | 4.4 × 10⁻⁴ |
| Diabetes | 0.234 | −0.0709 | 0.0145 | 1.1 × 10⁻⁶ |
| COPD | 0.188 | +0.0668 | 0.0232 | 4.0 × 10⁻³ |

All three are real and all three are small, roughly three to four times smaller than the creatinine main effect they modify.

**The two levers, matched at roughly eightfold**

| Measurement | count: k 20 → 150 | size: n 500 → 4000 |
|---|---|---|
| Divergence | 6–17% → **zero** (all six series) | −29% to +25%, none near zero |
| Interquartile width | **−42% to −66%** (all six series) | −19% to +10%, no consistent sign |
| CART uses the comorbidity | 26–40% at k=250 | 4–14% at every n |
| GLM detection | — | 2–16%, no trend in n |

Four readings on the same replicates, so not four independent confirmations — but two of them never touch the synthesizer. (Divergence is the interaction running past `|2|`, forty times the full-pool estimate.)

Every series is at zero by `k = 80`, though not at the same count, and a zero over 100 replicates bounds a rate at 3.6% rather than establishing it as nil. What the data supports is a threshold in (40, 80] for each series, not one shared value.

**Sign recovery** (COPD, the only positive interaction, +0.0668)

| | synthetic median |
|---|---|
| k = 20 (count sweep) | −0.177 |
| k = 20, n = 250…4000 (size sweep) | −0.077 to −0.159 |
| k = 40 (count sweep) | **+0.039** |

Raising the count fixes the sign. Raising `n` sixteen-fold does not.

**The diagnostic needs a permutation null**

The permutation null climbs from 0.00 to 0.32 across the count sweep, driven by class balance alone: 1% positive at `k = 5`, 50% at `k = 250`, and a balanced binary variable is structurally easy to split on whether or not it carries signal. Read across the grid the signal survives — a conditional logistic model puts the real column ahead of the shuffled one by an odds ratio of 2.0 at `k = 40` (p = 6 × 10⁻⁶) — but any diagnostic of this form needs permutation calibration before it can be read.

**What synthesis costs when the target is large enough to see** (simulation, known truth)

| True value | Real median | Synthetic median | Real median \|error\| | Synthetic median \|error\| |
|---|---|---|---|---|
| 0.5 | 0.521 | 0.496 | 0.164 | 0.219 |
| 1.5 | 1.73 | 0.891 | 0.292 | 0.713 |
| 3.0 | 2.75 | 1.22 | 0.688 | 1.90 |

The row at 0.5 is **not** evidence that small interactions survive synthesis. An estimator that saturates looks accurate against any target small enough, and this one saturates — the three synthetic medians sit close to `0.70 × √true`. Where the medians agree to within 1%, the per-replicate error is still 34% larger on the synthetic side.

**Attenuation is bought, not incurred.** At a true value of 3 the real arm diverges in 21–29% of replicates while the synthetic arm diverges in 2.5–10.5% — the reverse of the MIMIC-IV pattern. A strong interaction puts the real subsample close to separation; the tree smooths that away, so the synthetic fit converges more readily and converges on the wrong value. At `k = 80` the real median is 2.85 against a truth of 3, and the synthetic median is −0.150.

**What the design could never resolve**

A 500-row subsample carries 6–7% power against these interactions; the GLM arm returned 6.7–10.1%, both within a few points of the 5% a test with no power at all would produce. Theory and experiment landing in the same place rules out a defect in the pipeline and locates the failure in the design.

| | smallest detectable at n=500 | actually present | factor |
|---|---|---|---|
| CHF | 0.472 | 0.0497 | 9.5× |
| Diabetes | 0.485 | 0.0709 | 6.8× |
| COPD | 0.775 | 0.0668 | 11.6× |

Reaching those interactions would need roughly 23,000 to 67,000 rows.

## Method Notes

The choices that are not obvious from the code:

- **The rare stratum is constructed, not found.** Each replicate draws `k` positives and `n − k` negatives from separate pools, so `k` is a controlled variable rather than a random one. One consequence runs through the results: the real subsample can never be inestimable, since it holds exactly `k` positives by construction — only the synthetic table, which generates its own column, can come back with none.
- **The permutation preserves the marginal count.** Shuffling the comorbidity within the subsample keeps the number of positives at `k` and destroys only the row-level pairing with the outcome. Without that, the null would be a differently sized stratum and the comparison would mean nothing.
- **The probe is matched to the generator; the synthetic fit deliberately is not.** The diagnostic tree uses `syn.cart`'s own `rpart.control(minsplit = 20, minbucket = 5, cp = 1e-08)`, and the logistic arm's likelihood-ratio test drops the comorbidity and its interaction together, so both are asked the same question. The synthetic arm, by contrast, is fitted the way a data recipient would fit it — one copy, an ordinary GLM, no correction for the variance synthesis contributes — so the two arms are compared as procedures in use, not as estimators of matched validity.

## Reflections & Next Steps

The count-not-size mechanism carries from Bayesian logistic regression into CART-based synthesis, and holds not just for the output estimate but for whether the tree spends a split on the variable at all. The more transferable result may be the methodological one: a variable-selection diagnostic of the form "does the generator use variable X" cannot be read without a permutation control. A related warning applies to stability metrics — at a large true effect the synthetic arm diverges *less* than the real one, not because it is better behaved but because it has smoothed the effect away. Convergence there is a symptom, not a reassurance.

- **A third sweep**: fix prevalence and raise `k` and `n` together; width should fall as `1/√n`.
- **More than one generator**: the threshold may be specific to `syn.cart`'s `minbucket = 5`; a parametric synthesizer separates "CART" from "synthesis".
- **The privacy side**: CART's reluctance to isolate a five-patient node is the same reluctance that protects those patients, so disclosure risk and fidelity should move together. Nothing here measures it.

## Repository

```
code/
  ├── mimic-fidelity-analysis.ipynb      # count and size sweeps, full-pool reference, missingness
  ├── synthetic-fidelity-analysis.ipynb  # known-truth simulation at three effect sizes
  ├── tree-structure-diagnosis.ipynb     # CART probe, permutation null, conditional logistic
  ├── sql/                               # MIMIC-IV queries, loaded via read_sql()
  └── config.example.R                   # template — copy to config.R for local DATA_DIR
figures/                                 # report figures, written by the notebooks
report/mimic_fidelity_report.tex         # full writeup
```

MIMIC-IV is not redistributable. `config.R` points `DATA_DIR` at a local copy obtained through PhysioNet credentialed access; the notebooks read the `.csv.gz` files directly with DuckDB, so no database server is needed.

## Tools

- **Methods**: CART-based synthesis, Monte Carlo replication (M=100 on MIMIC-IV, M=200 in simulation), permutation testing, conditional logistic regression, cluster-robust standard errors, penalized (Firth) logistic regression, separation detection, power and minimum-detectable-effect analysis
- **Stack**: R — synthpop, rpart, logistf, detectseparation, sandwich, survival, lmtest, duckdb, dplyr, ggplot2
- **Data**: MIMIC-IV v2.2 (PhysioNet, credentialed access)

## References

Full reference list in the [report](report/mimic_fidelity_report.pdf).

Johnson, A.E.W., et al. (2023). MIMIC-IV, a freely accessible electronic health record dataset. *Scientific Data*, 10, 1.

Nowok, B., Raab, G.M., & Dibben, C. (2016). synthpop: Bespoke creation of synthetic data in R. *Journal of Statistical Software*, 74(11), 1–26.

Raab, G.M., Nowok, B., & Dibben, C. (2016). Practical data synthesis for large samples. *Journal of Privacy and Confidentiality*, 7(3), 67–97.

Strobl, C., et al. (2007). Bias in random forest variable importance measures. *BMC Bioinformatics*, 8, 25.

Prior work this builds on: [Bayesian Prior Sensitivity](https://github.com/ShengPeiWilliam/bayesian-prior-sensitivity).
