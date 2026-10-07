# 01_fit_primary_and_sensitivity.R
#
# Fits the primary model and the nine pre-specified sensitivity analyses used for model
# selection (Methods 2.4; Supplementary Tables ST2-ST3), run sequentially:
#   main            primary model, 1000 warm-up + 1000 sampling
#   repperiod       reporting fraction by 3 calendar periods (model_repperiod.stan), 1000+1000
#   rrprior_sd0.5 / rrprior_sd1.0
#                   s(a) ~ lognormal(0, sd), no lower bound (model_rrprior.stan), 500+500
#   pr_pooled / pr_2016only
#                   no CS-risk breakpoint: p_r = 39/451 or 11/76 for 2016-2025, 500+500
#   delay1 / delay2 delay distribution with one period / two periods (2016-2021, 2022-),
#                   500+500 (the primary model is the three-period case)
#   rrfloor0.8 / rrfloor0.5
#                   lower bound of s(a) relaxed from 1 to 0.8 / 0.5, 500+500
# 4 chains, seed 123. Every fit writes results/revision_summary_<label>.csv (mean, sd,
# 2.5/5/50/95/97.5% quantiles, R-hat, ESS), results/revision_draws_<label>.rds (selected
# variables), results/revision_loo_<label>.rds, and one row of
# results/revision_fit_diagnostics.csv. A failure in one fit does not stop the others.
#
# Run time: roughly 7-9 hours in total on a 4-core machine.
rm(list=ls(all=TRUE))
library(cmdstanr)
library(dplyr)
library(tidyr)
library(posterior)
library(loo)

source("R/data_prep.R")
dir.create("results", showWarnings = FALSE)
SEED <- 123
CHAINS <- 4
DIAG_FILE <- "results/revision_fit_diagnostics.csv"
if (file.exists(DIAG_FILE)) file.remove(DIAG_FILE)

KEEP <- c("i_t", "s_i", "rep", "rep1", "rep2", "rep3", "rr", "rrisk", "p1", "p2", "p3", "p_delay", "i_a",
          "cs_pop", "gst_pop", "cum_cs_w", "log_lik_case", "log_lik_cs")

run_fit <- function(stan_file, label, data_extra, warmup, sampling) {
  cat(sprintf("\n=== %s (%s, %d+%d) ===\n", label, stan_file, warmup, sampling))
  mod <- cmdstan_model(stan_file)
  fit <- mod$sample(data = c(base_stan_data, data_extra), seed = SEED,
                    chains = CHAINS, parallel_chains = CHAINS,
                    iter_warmup = warmup, iter_sampling = sampling, refresh = 100)
  ds <- fit$diagnostic_summary()
  vars <- intersect(KEEP, names(fit$metadata()$stan_variable_sizes))
  draws <- fit$draws(variables = vars)
  saveRDS(draws, sprintf("results/revision_draws_%s.rds", label))
  summ <- summarise_draws(draws, "mean", "sd",
                          ~quantile2(.x, probs = c(0.025, 0.05, 0.5, 0.95, 0.975)),
                          "rhat", "ess_bulk", "ess_tail")
  write.csv(summ, sprintf("results/revision_summary_%s.csv", label), row.names = FALSE)
  ll <- as_draws_matrix(fit$draws(c("log_lik_case", "log_lik_cs")))
  r_eff <- relative_eff(exp(ll), chain_id = rep(seq_len(CHAINS), each = sampling))
  lo <- loo(ll, r_eff = r_eff)
  saveRDS(lo, sprintf("results/revision_loo_%s.rds", label))
  print(lo)
  diag <- data.frame(label = label, warmup = warmup, sampling = sampling,
                     num_divergent = sum(ds$num_divergent),
                     max_rhat = max(summ$rhat, na.rm = TRUE),
                     elpd_loo = lo$estimates["elpd_loo", "Estimate"],
                     elpd_loo_se = lo$estimates["elpd_loo", "SE"],
                     p_loo = lo$estimates["p_loo", "Estimate"],
                     n_bad_pareto_k = sum(lo$diagnostics$pareto_k > 0.7))
  write.table(diag, DIAG_FILE, sep = ",", row.names = FALSE,
              col.names = !file.exists(DIAG_FILE),
              append = file.exists(DIAG_FILE))
  invisible(NULL)
}

jobs <- list(
  list("stan/model_primary.stan", "main",
       list(p_r_fixed = as.array(p_r_primary), rr_floor = 1.0), 1000, 1000),
  list("stan/model_repperiod.stan", "repperiod",
       list(p_r_fixed = as.array(p_r_primary), rr_floor = 1.0), 1000, 1000),
  list("stan/model_rrprior.stan", "rrprior_sd0.5",
       list(p_r_fixed = as.array(p_r_primary), rr_floor = 0.0, rr_prior_sd = 0.5), 500, 500),
  list("stan/model_rrprior.stan", "rrprior_sd1.0",
       list(p_r_fixed = as.array(p_r_primary), rr_floor = 0.0, rr_prior_sd = 1.0), 500, 500),
  list("stan/model_primary.stan", "pr_pooled",
       list(p_r_fixed = as.array(rep(39 / 451, 2)), rr_floor = 1.0), 500, 500),
  list("stan/model_primary.stan", "pr_2016only",
       list(p_r_fixed = as.array(rep(11 / 76, 2)), rr_floor = 1.0), 500, 500),
  # delay-distribution period structure (P=3 is the primary model = "main"):
  list("stan/model_delayperiod.stan", "delay1",
       list(p_r_fixed = as.array(p_r_primary), rr_floor = 1.0, P = 1L,
            delay_period = as.array(rep(1L, T))), 500, 500),
  list("stan/model_delayperiod.stan", "delay2",
       list(p_r_fixed = as.array(p_r_primary), rr_floor = 1.0, P = 2L,
            delay_period = as.array(c(rep(1L, 24), rep(2L, T - 24)))), 500, 500),
  # relaxed lower bound of s(a), no prior:
  list("stan/model_primary.stan", "rrfloor0.8",
       list(p_r_fixed = as.array(p_r_primary), rr_floor = 0.8), 500, 500),
  list("stan/model_primary.stan", "rrfloor0.5",
       list(p_r_fixed = as.array(p_r_primary), rr_floor = 0.5), 500, 500)
)

for (j in jobs) {
  res <- tryCatch(run_fit(j[[1]], j[[2]], j[[3]], j[[4]], j[[5]]),
                  error = function(e) cat(sprintf("\n!!! %s FAILED: %s\n", j[[2]], conditionMessage(e))))
}
cat("\n=== all fits done ===\n")
