# 04_calibration.R
#
# Calibration of the primary model (Supplementary Table ST5): for each observation
# (notified cases Case[c,t], 234 points; CS[t], 39 points), checks whether the observed
# count falls inside the posterior predictive interval (case_prd / cs_prd from the Stan
# generated quantities) at 90% and 95% nominal coverage.
# Primary model, 4 chains x (1000 warm-up + 1000 sampling), seed 123 -- the same settings
# as the "main" fit in 01_fit_primary_and_sensitivity.R.
#
# Outputs: results/calibration_fixedp_primary_{summary,by_category,case_detail,cs_detail}.csv
#
# Run time: roughly 1-1.5 hours on a 4-core machine.

rm(list = ls(all = TRUE))
library(cmdstanr)
library(dplyr)
library(tidyr)
library(posterior)

source("R/data_prep.R")
dir.create("results", showWarnings = FALSE)

CHAINS <- 4
WARMUP <- 1000
SAMPLING <- 1000
SEED <- 123

stan_data <- c(base_stan_data, list(rr_floor = 1.0, p_r_fixed = as.array(p_r_primary)))

cat("=== compiling stan/model_primary.stan ===\n")
mod <- cmdstan_model("stan/model_primary.stan")
cat(sprintf("=== fitting primary model: p_r_fixed = (%.4f, %.4f) ===\n",
            p_r_primary[1], p_r_primary[2]))
fit <- mod$sample(
  data = stan_data,
  seed = SEED,
  chains = CHAINS, parallel_chains = CHAINS,
  iter_warmup = WARMUP, iter_sampling = SAMPLING,
  refresh = 100
)

## ---- calibration: 90%/95% posterior predictive interval coverage ----

coverage_at <- function(draws_mat, observed, probs_lower, probs_upper) {
  lo <- apply(draws_mat, 2, quantile, probs = probs_lower, names = FALSE)
  hi <- apply(draws_mat, 2, quantile, probs = probs_upper, names = FALSE)
  covered <- observed >= lo & observed <= hi
  list(covered = covered, lo = lo, hi = hi)
}

## Case[c,t]: build the (c,t) ordering explicitly and select draws columns by name --
## don't rely on however cmdstanr happens to flatten a [C,T] array, since a silent
## column/row transposition here wouldn't error, just misassign every category label.
case_cat <- rep(seq_len(C), times = T)
case_time <- rep(seq_len(T), each = C)
case_obs <- case_mat[cbind(case_cat, case_time)]   # matches case_mat[c, t] exactly

case_draws_all <- as_draws_matrix(fit$draws("case_prd"))
case_prd_names <- sprintf("case_prd[%d,%d]", case_cat, case_time)
stopifnot(all(case_prd_names %in% colnames(case_draws_all)))
case_draws <- case_draws_all[, case_prd_names]

case_90 <- coverage_at(case_draws, case_obs, 0.05, 0.95)
case_95 <- coverage_at(case_draws, case_obs, 0.025, 0.975)

## CS[t]
cs_draws_all <- as_draws_matrix(fit$draws("cs_prd"))
cs_prd_names <- sprintf("cs_prd[%d]", seq_len(T))
stopifnot(all(cs_prd_names %in% colnames(cs_draws_all)))
cs_draws <- cs_draws_all[, cs_prd_names]
cs_obs <- cs_mat
cs_90 <- coverage_at(cs_draws, cs_obs, 0.05, 0.95)
cs_95 <- coverage_at(cs_draws, cs_obs, 0.025, 0.975)

## ---- overall summary ----
calib_summary <- tibble(
  series = c("Case (overall, all categories x quarters)", "CS (overall, all quarters)"),
  n = c(length(case_obs), length(cs_obs)),
  coverage_90pct = c(mean(case_90$covered), mean(cs_90$covered)),
  coverage_95pct = c(mean(case_95$covered), mean(cs_95$covered))
)

cat("\n=== overall calibration (posterior predictive interval coverage) ===\n")
print(calib_summary)

## ---- per-category breakdown for Case ----
calib_by_category <- tibble(
  category = seq_len(C),
  age_band = c("15-19", "20-24", "25-29", "30-34", "35-39", "40-44"),
  n = as.integer(table(case_cat)),
  coverage_90pct = tapply(case_90$covered, case_cat, mean),
  coverage_95pct = tapply(case_95$covered, case_cat, mean)
)

cat("\n=== Case calibration by age category ===\n")
print(calib_by_category)

## ---- per-observation detail (for follow-up plotting/inspection) ----
case_detail <- tibble(
  category = case_cat, age_band = calib_by_category$age_band[case_cat], time = case_time,
  year = year_vec_T[case_time], observed = case_obs,
  lo90 = case_90$lo, hi90 = case_90$hi, covered_90 = case_90$covered,
  lo95 = case_95$lo, hi95 = case_95$hi, covered_95 = case_95$covered
)

cs_detail <- tibble(
  time = seq_len(T), year = year_vec_T, observed = cs_obs,
  lo90 = cs_90$lo, hi90 = cs_90$hi, covered_90 = cs_90$covered,
  lo95 = cs_95$lo, hi95 = cs_95$hi, covered_95 = cs_95$covered
)

write.csv(calib_summary, "results/calibration_fixedp_primary_summary.csv", row.names = FALSE)
write.csv(calib_by_category, "results/calibration_fixedp_primary_by_category.csv", row.names = FALSE)
write.csv(case_detail, "results/calibration_fixedp_primary_case_detail.csv", row.names = FALSE)
write.csv(cs_detail, "results/calibration_fixedp_primary_cs_detail.csv", row.names = FALSE)

cat("\nSaved: results/calibration_fixedp_primary_summary.csv, ",
    "results/calibration_fixedp_primary_by_category.csv, ",
    "results/calibration_fixedp_primary_case_detail.csv, ",
    "results/calibration_fixedp_primary_cs_detail.csv\n", sep = "")
