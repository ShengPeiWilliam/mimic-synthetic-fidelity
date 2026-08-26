# Positive-Case Count and Synthetic Data Fidelity: A CART-Based Analysis of MIMIC-IV
[![Full Report](https://img.shields.io/badge/📄_Read_Full_Report-PDF-blue?style=for-the-badge)](report/mimic_fidelity_report.pdf)

CART-based synthetic data generation (`synthpop`) on MIMIC-IV (n=70,954 ICU stays), testing whether fidelity inside a rare stratum is governed by the stratum's positive count or by the size of the table it sits in. Two sweeps hold the generator fixed and vary only the input: one varies the positive count `k` at fixed sample size, the other varies sample size `n` at fixed `k`.

Every quantity measured follows the count and none follows sample size. Divergence falls to zero between counts of 40 and 80 in all six series; holding `k = 20` and raising `n` sixteen-fold never brings it near zero. The same contrast shows up in estimate spread, in whether CART spends a split on the variable at all, and in a plain logistic fit. The generator inherits that threshold rather than adding one of its own — what it adds instead is a bias no count removes. A companion simulation with a known truth puts a number on it: the synthetic median recovers 99% of a true interaction of 0.5 and 41% of one of 3.

The mechanism does not rescue the analysis. At MIMIC-IV's effect sizes a 500-row subsample carries 6–7% power, and the smallest interaction it could detect is seven to twelve times the one actually present.

[See key results →](#key-results)

## Motivation

An earlier project of mine, [Bayesian Prior Sensitivity](https://github.com/ShengPeiWilliam/bayesian-prior-sensitivity), found something on a 189-row dataset that I did not expect: the point where estimates stop being driven by the prior is set by a predictor's positive-case count, not by the sample size. Smoking, with 74 positives, settled by n=40. Hypertension, with 12, was still unstable at n=80.

That a logistic fit's precision is set by the event count is not itself new — it is the events-per-variable literature, and the conditions under which thin event counts produce separation have been characterised since Albert and Anderson (1984). What is not established is whether the same threshold governs a **generator**. Whether a CART-based synthesizer preserves a rare-stratum interaction is a question about which variables a tree spends a split on, not about the curvature of a likelihood, and the two need not share a threshold.

[MIMIC-IV](https://physionet.org/content/mimiciv/2.2/) is close to the ideal place to ask: 70,954 ICU stays, so a rare stratum can be built to any size on demand instead of being whatever the data happens to give.

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

Cell by cell this is underpowered, because each test turns on a handful of discordant pairs. The grid is therefore read with one conditional logistic model instead of twenty-one separate tests, each replicate its own stratum:

```r
clogit(used ~ observed + observed:log2k_c + strata(pair))
```

**The probe is matched to the generator, and the GLM to the probe**

```r
rpart.control(minsplit = 20, minbucket = 5, cp = 1e-08)
```

These are `synthpop::syn.cart`'s own settings, so the diagnostic tree reflects what the generator does rather than a more permissive stand-in. It is still a proxy: three predictors instead of the full feature set, and no sequential visit order.

The logistic arm is matched the same way. Its likelihood-ratio test drops the comorbidity and its interaction together, so it answers the tree's question — does this variable earn a place at all — instead of a different one. That is what separates "this generator misses the effect" from "nothing could find it at this size."

**Two thresholds, two denominators**

|2| is not a general cutoff, it is 40x the full-pool estimate of ~0.05. Everything is recomputed at 5 and 10 and every comparison points the same way. It under-detects at the smallest counts: at `k = 5` a separation test finds 64–76% of real fits separated where |2| flags 18–34%, since `glm` stops after twenty-five iterations rather than running to the boundary. From `k = 20` up the two agree.

Failure rates are shares of all `M` replicates, so converged + inestimable + divergent = 1. Quantiles use converged fits only. Two denominators in the same table, stated explicitly because mixing them is the easiest mistake to make here.

**Rows are not independent, and the pool is not a random sample**

The outcome and the comorbidity flags are recorded per admission, but the pool is indexed on ICU stay: 70,954 stays span 64,308 admissions and 49,563 patients, so 30% of rows repeat a patient. Standard errors on the full-pool fit are clustered on patient, which inflates them by 10–18%.

The 2,227 stays dropped for having no first-day creatinine are not missing at random either. They die in hospital at 28.1% against 10.9%, their median stay is 0.73 days against 1.96, and they carry every comorbidity less often — thinner documentation rather than healthier patients, since a stay ending within hours has neither the time for the draw nor the encounter for the coding. Both sweep arms come from the same pool, so the comparison is unaffected; what the selection limits is the clinical reading of the reference.

**The synthetic arm is fitted the way a recipient would fit it**

One synthetic copy, an ordinary GLM, no correction for the variance synthesis itself contributes. Proper inference would use several copies and a combining rule (Raab et al., 2016), which would widen the synthetic intervals reported here. The naive fit is measured because it is the one a data recipient performs — but the comparison is therefore between procedures as used, not between estimators of matched validity.

## Key Results

**Full-pool reference (n=70,954, standard errors clustered on patient)**

| Comorbidity | Creatinine | Interaction | SE | p |
|---|---|---|---|---|
| CHF | 0.208 | −0.0497 | 0.0141 | 4.4 × 10⁻⁴ |
| Diabetes | 0.234 | −0.0709 | 0.0145 | 1.1 × 10⁻⁶ |
| COPD | 0.188 | +0.0668 | 0.0232 | 4.0 × 10⁻³ |

All three interactions are real and all three are small — roughly three to four times smaller than the creatinine main effect they modify.

**The two levers, matched at roughly eightfold**

| Measurement | count: k 20 → 150 | size: n 500 → 4000 |
|---|---|---|
| Divergence | 6–17% → **zero** (all six series) | −29% to +25%, none near zero |
| Interquartile width | **−42% to −66%** (all six series) | −19% to +10%, no consistent sign |
| CART uses the comorbidity | 26–40% at k=250 | 4–14% at every n |
| GLM detection | — | 2–16%, no trend in n |

Four readings on the same replicates, so not four independent confirmations — but two of them never touch the synthesizer.

**Where the threshold actually sits**

Every series is at zero by `k = 80`, but not at the same count: the CHF and COPD real arms are already there at `k = 40`, where the other four sit at 2–4%. Each zero is 0 of 100 replicates, which bounds a rate at 3.6% rather than establishing it as nil. What the data supports is a threshold in (40, 80] for each series, not one shared value. A penalized fit removes 60–70% of the divergence at `k = 5` but does not move where it stops — the threshold is not an artifact of the estimator.

**The sharpest single comparison**

At `k = 5` and 1% prevalence the synthetic table loses the stratum 4–6% of the time. At `k = 20` and the same 1% prevalence (n = 2000) it never does. Prevalence held equal, count differs, outcome differs.

**Spread against the rate the count predicts**

An interaction's variance draws on both groups, so it scales as `1/k + 1/(n-k)`. Rescaling the `k = 250` width by that factor reproduces the real arm closely — median ratio 0.97, within 0.85 and 1.07. The synthetic arm sits higher at **every one of the nine matched cells** (same comorbidity, same count), median 1.14, sign test p = 0.004. The count accounts for the real spread; the synthetic spread carries a further excess it does not explain.

**Sign recovery (COPD, the only positive interaction, +0.0668)**

| | synthetic median |
|---|---|
| k = 20 (Experiment 1) | −0.177 |
| k = 20, n = 250…4000 (Experiment 2) | −0.077 to −0.159 |
| k = 40 (Experiment 1) | **+0.039** |

Raising the count fixes the sign. Raising `n` sixteen-fold does not.

**Tree diagnostic against its null**

Read cell by cell the diagnostic clears p < 0.01 in two of twenty-one cells. Read across the grid it is unambiguous: the conditional logistic model puts the real column ahead of the shuffled one by an odds ratio of 2.0 at `k = 40` (p = 6 × 10⁻⁶), narrowing by a factor of 0.81 per doubling of the count (p = 0.015) — fitted, from 3.7 at `k = 5` to 1.2 at `k = 250`.

That narrowing is the null rising to meet the observed rate, not the observed rate pulling away. The permutation null climbs from 0.00 to 0.32 across the count sweep, driven by class balance alone: 1% positive at `k = 5`, 50% at `k = 250`, and a balanced binary variable is structurally easy to split on whether or not it carries signal. Any diagnostic of this form needs permutation calibration before it can be read.

**What synthesis costs when the target is large enough to see (simulation, known truth)**

MIMIC-IV cannot answer this: there is no truth to compare against and no power to find one. A companion simulation supplies both, sweeping the same counts at three true values, 200 replicates per cell.

| True value | Real median | Synthetic median | Shortfall | Real median \|error\| | Synthetic median \|error\| |
|---|---|---|---|---|---|
| 0.5 | 0.521 | 0.496 | 1% | 0.164 | 0.219 |
| 1.5 | 1.73 | 0.891 | 41% | 0.292 | 0.713 |
| 3.0 | 2.75 | 1.22 | 59% | 0.688 | 1.90 |

The row at 0.5 is **not** evidence that small interactions survive synthesis. An estimator that saturates looks accurate against any target small enough, and this one saturates — the three synthetic medians sit close to `0.70 × √true`. The error columns say it without the fit: at a true value of 0.5 the medians agree to within 1%, but the per-replicate error is 34% larger on the synthetic side.

**Attenuation is bought, not incurred.** At a true value of 3 the real arm diverges in 21–29% of replicates between `k = 20` and `k = 80` while the synthetic arm diverges in 2.5–10.5% — the reverse of the MIMIC-IV pattern. A strong interaction puts the real subsample close to separation; the tree smooths that away, so the synthetic fit converges more readily and converges on the wrong value. At `k = 80` the real median is 2.85 against a truth of 3, and the synthetic median is −0.150.

Below the threshold the count governs the sign, not just the spread: the real median carries the wrong sign through `k = 5` at a true value of 0.5, through `k = 10` at 1.5 and through `k = 20` at 3. The larger the interaction, the longer the median points the wrong way — the replicates left converging are the ones that did not find it.

**What the design could never resolve**

A 500-row subsample carries 6–7% power against these interactions; the GLM arm returned 6.7–10.1%, both within a few points of the 5% a test with no power at all would produce. Theory and experiment landing in the same place is what rules out a defect in the pipeline and locates the failure in the design.

The gap reads harder as a detectable effect than as a rate:

| | smallest detectable at n=500 | actually present | factor |
|---|---|---|---|
| CHF | 0.472 | 0.0497 | 9.5× |
| Diabetes | 0.485 | 0.0709 | 6.8× |
| COPD | 0.775 | 0.0668 | 11.6× |

Reaching those interactions would need 23,384 to 67,238 rows. Even the best cell in the study — CHF real arm at `k = 250`, zero failures — has an interquartile width of 0.194 against a reference of 0.0497, nearly four times the quantity being estimated.

That floor also explains a choice in the simulation that would otherwise look arbitrary: 0.5, its smallest true value, is of the same order as the smallest effect this design can detect at all. The simulation starts where MIMIC-IV's detectable range ends.

## Reflections & Next Steps

The count-not-size mechanism carries from Bayesian logistic regression into CART-based synthesis, and it reaches further than expected: it holds not just for the output estimate but for the generator's internal behaviour, whether the tree spends a split on the variable at all. The generator does not add a threshold of its own — it inherits the stratum's. What it adds is attenuation, and the count does not remove it.

The more transferable result may be the methodological one. Variable-selection diagnostics of the form "does the generator use variable X" cannot be read without a permutation control, and they are only legible where that null is near zero, which on this data means low prevalence. This is a known CART bias placed in a setting where it silently produces the wrong conclusion.

The second is a warning about stability metrics on synthetic data. At a large true effect the synthetic arm diverges *less* than the real one — not because it is better behaved, but because it has smoothed the effect away. Convergence there is a symptom, not a reassurance.

Next steps, each with what the current results predict:

- **A third sweep**: fix prevalence and raise `k` and `n` together. Width should fall as `1/√n` and sit on the `1/k + 1/(n-k)` curve the real arm already tracks. Neither existing sweep can show this.
- **More than one generator**: the threshold may be specific to `syn.cart`'s `minbucket = 5`. The attenuation mechanism — a tree declining to spend a split on a thinly supported subgroup — predicts that relaxing `minbucket` reduces the shortfall. A parametric synthesizer is the baseline that separates "CART" from "synthesis".
- **The privacy side**: the report argues CART's reluctance to isolate a five-patient node is the same reluctance that protects those patients. If so, a disclosure-risk measure and the fidelity measure should move together. Nothing here measures the privacy half.
- **Proper synthetic inference**: `m > 1` copies with the Raab combining rule, against the naive single-copy fit measured here.
- **External validation**: eICU-CRD has a different admission structure and the same variables.

## Repository

```
code/
  ├── mimic-fidelity-analysis.ipynb      # Experiments 1 and 2, full-pool reference, missingness
  ├── synthetic-fidelity-analysis.ipynb  # Known-truth simulation at three effect sizes
  ├── tree-structure-diagnosis.ipynb     # CART probe, permutation null, conditional logistic
  ├── sql/                               # MIMIC-IV queries, loaded via read_sql()
  ├── config.example.R                   # template — copy to config.R
  └── config.R                           # local DATA_DIR (not tracked)
figures/                                 # report figures, written by the notebooks
report/
  ├── mimic_fidelity_report.tex          # full writeup
  └── proposal.tex                       # original proposal
```

MIMIC-IV is not redistributable. `config.R` points `DATA_DIR` at a local copy obtained through PhysioNet credentialed access; the notebooks read the `.csv.gz` files directly with DuckDB, so no database server is needed.

## Tools

**Statistical methods**: CART-based synthesis, Monte Carlo replication (M=100 on MIMIC-IV, M=200 in simulation), permutation testing, conditional logistic regression, exact paired binomial tests, cluster-robust standard errors, penalized (Firth) logistic regression, separation detection, power and minimum-detectable-effect analysis
**Language**: R
**Libraries**: synthpop, rpart, logistf, detectseparation, sandwich, survival, lmtest, duckdb, DBI, dplyr, tidyr, ggplot2
**Data**: MIMIC-IV v2.2 (PhysioNet, credentialed access)

## References

Johnson, A.E.W., Bulgarelli, L., Shen, L., et al. (2023). MIMIC-IV, a freely accessible electronic health record dataset. *Scientific Data*, 10, 1.

Nowok, B., Raab, G.M., & Dibben, C. (2016). synthpop: Bespoke creation of synthetic data in R. *Journal of Statistical Software*, 74(11), 1–26.

Raab, G.M., Nowok, B., & Dibben, C. (2016). Practical data synthesis for large samples. *Journal of Privacy and Confidentiality*, 7(3), 67–97.

Peduzzi, P., Concato, J., Kemper, E., Holford, T.R., & Feinstein, A.R. (1996). A simulation study of the number of events per variable in logistic regression analysis. *Journal of Clinical Epidemiology*, 49(12), 1373–1379.

van Smeden, M., de Groot, J.A.H., Moons, K.G.M., et al. (2016). No rationale for 1 variable per 10 events criterion for binary logistic regression analysis. *BMC Medical Research Methodology*, 16, 163.

Albert, A., & Anderson, J.A. (1984). On the existence of maximum likelihood estimates in logistic regression models. *Biometrika*, 71(1), 1–10.

Firth, D. (1993). Bias reduction of maximum likelihood estimates. *Biometrika*, 80(1), 27–38.

Heinze, G., & Schemper, M. (2002). A solution to the problem of separation in logistic regression. *Statistics in Medicine*, 21(16), 2409–2419.

Strobl, C., Boulesteix, A.-L., Zeileis, A., & Hothorn, T. (2007). Bias in random forest variable importance measures. *BMC Bioinformatics*, 8, 25.

Prior work this builds on: [Bayesian Prior Sensitivity](https://github.com/ShengPeiWilliam/bayesian-prior-sensitivity) — the positive-count finding in Bayesian logistic regression that this project tests against a different generator.
