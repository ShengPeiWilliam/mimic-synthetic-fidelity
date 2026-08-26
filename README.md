# Positive-Case Count and Synthetic Data Fidelity: A CART-Based Analysis of MIMIC-IV
[![Full Report](https://img.shields.io/badge/📄_Read_Full_Report-PDF-blue?style=for-the-badge)](report/mimic_fidelity_report.pdf)

CART-based synthetic data generation (`synthpop`) on MIMIC-IV (n=70,954 ICU stays), testing whether fidelity inside a rare stratum is governed by the stratum's positive count or by the size of the table it sits in. Two sweeps hold the generator fixed and vary only the input: one varies the positive count `k` at fixed sample size, the other varies sample size `n` at fixed `k`.

Four independent measurements follow the count and none follows sample size. Divergence reaches exactly zero at `k = 80` in all six series; holding `k = 20` and raising `n` sixteen-fold never brings it near zero. The same contrast shows up in estimate spread, in whether CART spends a split on the variable at all, and in a plain logistic fit. But the mechanism does not rescue the analysis: at MIMIC-IV's effect sizes a 500-row subsample carries 6–8% power, and no count reachable inside one makes the estimate usable.

[See key results →](#key-results)

## Motivation

An earlier project of mine, [Bayesian Prior Sensitivity](https://github.com/ShengPeiWilliam/bayesian-prior-sensitivity), found something on a 189-row dataset that I did not expect: the point where estimates stop being driven by the prior is set by a predictor's positive-case count, not by the sample size. Smoking, with 74 positives, settled by n=40. Hypertension, with 12, was still unstable at n=80.

That was one small dataset and one model family, so it was fair to ask whether the pattern was real or a quirk of `birthwt`. This project is the follow-up, and [MIMIC-IV](https://physionet.org/content/mimiciv/2.2/) is close to the ideal place to run it: 70,954 ICU stays, so a rare stratum can be built to any size on demand instead of being whatever the data happens to give.

The generator is different too, which is the point. If the count threshold is a property of Bayesian logistic regression it should not show up in a CART-based synthesizer. If it does, it is a property of small strata rather than of any one method.

## Design Decisions

**Exact positive counts, not random ones**

The rare stratum is constructed, not found. Each replicate draws `k` positives and `n - k` negatives from separate pools:

```r
sub <- bind_rows(
  pos_pool[sample(nrow(pos_pool), k), ],
  neg_pool[sample(nrow(neg_pool), n_fixed - k), ]
)
```

`k` is a controlled variable, not a random one, which is what allows a sweep over it. It has a consequence worth naming: the real subsample can never be inestimable, since it holds exactly `k` positives by construction. Only the synthetic table, which generates its own column, can come back with none.

**The shuffle preserves the marginal count, and the subtraction is what you read**

```r
if (permute) sub[[comorbidity_col]] <- sample(sub[[comorbidity_col]])
```

Permuting within the subsample keeps the number of positives at exactly `k` and destroys only the row-level pairing with the outcome. Without that, the null would be a differently sized stratum and the comparison would mean nothing.

The quantity read off it is `excess = observed − permuted`: how much of the observed rate comes from information rather than from a column that is merely easy to split on. One seed, `k * 1000 + r`, drives the row draw, the synthesis and both tree fits, so the two runs are paired on identical rows and the subtraction has an exact reading:

```
excess = (only_obs − only_perm) / M
```

`only_obs` counts replicates where the real column earned a split and the shuffled one did not, `only_perm` the reverse. Replicates where both did, or neither did, cancel — they say nothing about which column is better. CHF at `k = 40`: 0.21 − 0.07 = 0.14, and (18 − 4) / 100 = 0.14.

A positive excess is not by itself evidence. `binom.test(c(only_obs, only_perm))`, the exact form of McNemar's test, asks whether the imbalance between the two discordant counts beats a coin flip. At `k = 40` that is 18 against 4, p = 0.004. At `k = 250` it is 17 against 20, and the apparent recovery in the raw rate is gone.

**The probe is matched to the generator, and the GLM to the probe**

```r
rpart.control(minsplit = 20, minbucket = 5, cp = 1e-08)
```

These are `synthpop::syn.cart`'s own settings, so the diagnostic tree reflects what the generator does rather than a more permissive stand-in. It is still a proxy: three predictors instead of the full feature set, and no sequential visit order.

The logistic arm is matched the same way. Its likelihood-ratio test drops the comorbidity and its interaction together, so it answers the tree's question — does this variable earn a place at all — instead of a different one. That is what separates "this generator misses the effect" from "nothing could find it at this size."

**Two thresholds, two denominators**

|2| is not a general cutoff, it is 40x the full-pool estimate of ~0.05. Everything is recomputed at 5 and 10 and every comparison points the same way.

Failure rates are shares of all `M` replicates, so converged + inestimable + divergent = 1. Quantiles use converged fits only. Two denominators in the same table, stated explicitly because mixing them is the easiest mistake to make here.

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
