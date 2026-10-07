# syphilis_japan

Code and data for:

> Nakajo K, Nishiura H. Estimating the incidence of syphilis from congenital syphilis cases in
> Japan, 2016–2025: a modeling study.

A Bayesian back-calculation model (Stan) reconstructs the quarterly number of newly diagnosed
syphilis cases in women of child-bearing age (15–44 years, six age groups) in Japan from
2016 Q1 to 2025 Q3. Notified syphilis cases and congenital syphilis (CS) cases are the
observed data. The model then estimates gestational syphilis and CS through age-specific
fertility rates, a fertility-rate ratio s(a), the CS risk among viable births and the delay
from diagnosis to delivery.

This repository contains the code for the revised manuscript. The code for the initial
submission is kept in [`archive/initial_submission/`](archive/initial_submission/).

## Repository layout

```
data/      Case.csv  notified syphilis cases by quarter (rows) and age group (age1-age6 = 15-19, ..., 40-44)
           CS.csv    notified CS cases by quarter
stan/      model_primary.stan       primary model
           model_repperiod.stan     sensitivity analysis: reporting fraction by 3 calendar periods
           model_rrprior.stan       sensitivity analysis: lognormal prior on s(a)
           model_delayperiod.stan   sensitivity analysis: delay distribution with 1 or 2 periods
R/         data_prep.R and analysis scripts 01-09 (see below)
results/   saved model output, so that 05-09 can be run without refitting
tables/    tables written by 05-07
figures/   Figures 1-3 of the main text
```

Age-specific fertility rates, annual live births and the CS-risk survey counts (11/76 and
28/375) are entered directly in `R/data_prep.R`.

## Scripts and outputs

Run every script from the repository root, e.g. `Rscript R/01_fit_primary_and_sensitivity.R`.

| Script | What it does | Used for | Run time* |
|---|---|---|---|
| `R/data_prep.R` | Shared data preparation, sourced by the other scripts | – | – |
| `R/01_fit_primary_and_sensitivity.R` | Primary model (1000+1000) and the nine sensitivity analyses for model selection | Table 1, Table 2, Figures 1–3, ST2, ST3, ST4 | 7–9 h |
| `R/02_cs_risk_lhs.R` | CS-risk sensitivity analysis: 10-point Latin hypercube sample of the two CS-risk proportions (500+500) | Table 1 (range column), Figure 3a, ranges quoted in Results | ~6 h |
| `R/03_cs_risk_lhs_full_it.R` | Re-fits the same 10 points, keeping the full quarterly series of diagnosed cases | Figure 2a | ~6 h |
| `R/04_calibration.R` | Posterior predictive interval coverage (90%/95%) of the primary model | ST5 | 1–1.5 h |
| `R/05_posterior_correlation.R` | Posterior correlation matrix of 19 parameters | ST4 | < 1 min |
| `R/06_model_comparison.R` | Convergence diagnostics and LOO comparison | ST2, ST3 | < 1 min |
| `R/07_main_tables.R` | Parameter estimates and fertility-rate ratio scenarios | Table 1, Table 2 | < 1 min |
| `R/08_figure1.R` | Observed vs. predicted notifications | Figure 1 | ~5 min |
| `R/09_figures2_3.R` | Diagnosed cases and gestational syphilis | Figures 2, 3 | ~1 min |

\* 4 chains in parallel on a 4-core machine. MCMC settings are given as warm-up + sampling
iterations per chain; all fits use seed 123.

Supplementary Table ST1 (parameters, priors and their basis) is descriptive and has no
associated script.

Scripts 05–09 use the saved files in `results/` and reproduce the published tables and
figures without refitting. To rerun everything from scratch, run 01–09 in order; 03 reads
the design written by 02, and 05–09 read the output of 01–03.

## Software

The analyses were run with R 4.5.3, CmdStan 2.36.0 and the R packages cmdstanr 0.9.0,
posterior 1.6.1, loo 2.8.0, dplyr 1.1.4, tidyr 1.3.1, ggplot2 4.0.0, patchwork 1.3.2 and
ragg 1.5.0.

```r
install.packages(c("dplyr", "tidyr", "ggplot2", "patchwork", "ragg", "posterior", "loo"))
install.packages("cmdstanr", repos = c("https://stan-dev.r-universe.dev", getOption("repos")))
cmdstanr::install_cmdstan()
```

MCMC results can differ slightly across CmdStan versions, compilers and platforms even with
the same seed.
