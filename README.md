# Positive-Case Count and Synthetic Data Fidelity: A CART-Based Analysis of MIMIC-IV
[![Full Report](https://img.shields.io/badge/📄_Read_Full_Report-PDF-blue?style=for-the-badge)](report/mimic_fidelity_report.pdf)

CART-based synthesis (`synthpop`) on MIMIC-IV (n=70,954 ICU stays). Does fidelity inside a rare stratum depend on the stratum's positive count, or on the size of the table it sits in? Two sweeps hold the generator fixed and vary one input each: the positive count `k`, or the sample size `n`.

**The count governs everything the sweeps measure; the sample size governs none of it. But the generator adds no threshold of its own. What it adds is an attenuation no count removes.**

An earlier project, [Bayesian Prior Sensitivity](https://github.com/ShengPeiWilliam/bayesian-prior-sensitivity), found the same count dependence in Bayesian logistic regression. That a logistic fit's precision follows the event count is not new; it is the events-per-variable literature. What is not established is whether the same threshold governs a **generator**, which is a question about which variables a tree spends a split on, not about the curvature of a likelihood.

## What the sweeps found

Matched at roughly eightfold, three comorbidities × two arms:

| | count: k 20 → 150 | size: n 500 → 4000 |
|---|---|---|
| Divergence | 6 to 17% to **zero** | stays at 2 to 24%, never zero |
| Interquartile width | **42 to 66% narrower** | −19% to +10%, no consistent sign |
| CART uses the comorbidity | 26 to 40% at k=250 | 4 to 14% at every n |
| GLM detection | not applicable | 2 to 16%, no trend in n |

Four readings on the same replicates, so not four independent confirmations, but two of them never touch the synthesizer. The threshold is series-dependent, somewhere in (40, 80].

None of this settles whether synthesis recovers the truth: a 500-row subsample carries 6 to 7% power against these interactions, and the smallest effect it could detect is seven to twelve times the one present.

## What synthesis costs

A companion simulation supplies the truth MIMIC-IV lacks, with the same sweep at three known values:

| True value | Real median | Synthetic median | Real \|error\| | Synthetic \|error\| |
|---|---|---|---|---|
| 0.5 | 0.521 | 0.496 | 0.164 | 0.219 |
| 1.5 | 1.73 | 0.891 | 0.292 | 0.713 |
| 3.0 | 2.75 | 1.22 | 0.688 | 1.90 |

The 0.5 row is **not** evidence that small interactions survive. An estimator that pulls toward zero is accurate against any target near zero, and the error columns show what the medians hide.

**Attenuation is bought, not incurred.** At a true value of 3 the real arm diverges in 21 to 29% of replicates and the synthetic arm in 2.5 to 10.5%, the reverse of the MIMIC-IV pattern. The tree smooths away the near-separation, so the synthetic fit converges more readily, and converges on the wrong value.

## Two things that would mislead

- **A split diagnostic without a null.** "Does the generator use variable X" reads as recovery across the count sweep, but the permutation null climbs from 0.00 to 0.32 alongside it, driven by class balance alone.
- **Stability as reassurance.** The synthetic arm diverging *less* is not better behaviour; it is the effect having been smoothed away.

## Notes on the design

- The rare stratum is **constructed, not found**: each replicate draws `k` positives and `n` minus `k` negatives from separate pools. So the real subsample can never be inestimable; only the synthetic table, which generates its own column, can come back empty.
- The permutation **preserves the marginal count**, destroying only the row-level pairing with the outcome.
- The probe is matched to the generator (`syn.cart`'s own control parameters); the synthetic fit deliberately is not, being one copy and an ordinary GLM, as a data recipient would fit it.

Left open: a third sweep at fixed prevalence, a second generator, and the privacy half of the trade-off the report argues for.

## Repository

```
code/
  ├── mimic-fidelity-analysis.ipynb      # count and size sweeps, full-pool reference, missingness
  ├── synthetic-fidelity-analysis.ipynb  # known-truth simulation at three effect sizes
  ├── tree-structure-diagnosis.ipynb     # CART probe, permutation null, conditional logistic
  ├── sql/                               # MIMIC-IV queries, loaded via read_sql()
  └── config.example.R                   # template, copy to config.R for local DATA_DIR
figures/                                 # report figures, written by the notebooks
report/mimic_fidelity_report.tex         # full writeup
```

MIMIC-IV is not included and requires PhysioNet credentialed access. Point `DATA_DIR` in `config.R` at a local copy; the notebooks read the `.csv.gz` files directly with DuckDB.

**Stack**: R with synthpop, rpart, logistf, detectseparation, sandwich, survival, lmtest, duckdb, dplyr, ggplot2

Full reference list in the [report](report/mimic_fidelity_report.pdf).
