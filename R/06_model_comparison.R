# 06_model_comparison.R
#
# Supplementary Table ST2 (convergence diagnostics) and ST3 (leave-one-out cross-validation
# comparison) for the primary model and the nine sensitivity analyses, from the outputs of
# 01_fit_primary_and_sensitivity.R (results/revision_fit_diagnostics.csv and
# results/revision_loo_<label>.rds).
#
# ST3: paired elpd_loo difference (alternative - primary) over the 273 observations
# (234 notified-case and 39 CS data points), with standard error sqrt(n) * sd(pointwise
# differences), as in loo::loo_compare.
#
# Outputs: tables/supplementary_table_ST2_convergence.csv,
#          tables/supplementary_table_ST3_loo_comparison.csv

library(loo)
dir.create("tables", showWarnings = FALSE)

fits <- data.frame(
  label = c("main", "repperiod", "delay1", "delay2", "pr_pooled", "pr_2016only",
            "rrprior_sd0.5", "rrprior_sd1.0", "rrfloor0.8", "rrfloor0.5"),
  name = c("Primary model (main)",
           "Reporting fraction: 3 periods",
           "Delay distribution: 1 period",
           "Delay distribution: 2 periods",
           "CS risk: pooled single value",
           "CS risk: 2016-2021 value only",
           "s(a) prior: lognormal(0, 0.5)",
           "s(a) prior: lognormal(0, 1.0)",
           "s(a) lower bound: 0.8",
           "s(a) lower bound: 0.5"),
  stringsAsFactors = FALSE
)

## ---- ST2: convergence diagnostics ----
diag <- read.csv("results/revision_fit_diagnostics.csv")
st2 <- merge(fits, diag, by = "label", sort = FALSE)
st2 <- st2[match(fits$label, st2$label), ]
st2 <- data.frame(
  Fit = st2$name,
  Setting = sprintf("%d+%d", st2$warmup, st2$sampling),
  Divergent_transitions = st2$num_divergent,
  Max_Rhat = sprintf("%.2f", st2$max_rhat),
  p_loo = sprintf("%.1f", st2$p_loo),
  Pareto_k_gt_0.7_of_273 = st2$n_bad_pareto_k,
  check.names = FALSE
)
print(st2, row.names = FALSE)
write.csv(st2, "tables/supplementary_table_ST2_convergence.csv", row.names = FALSE)

## ---- ST3: LOO comparison against the primary model ----
L <- lapply(fits$label, function(l) readRDS(sprintf("results/revision_loo_%s.rds", l)))
names(L) <- fits$label
cmp <- function(a) {
  d <- L[[a]]$pointwise[, "elpd_loo"] - L[["main"]]$pointwise[, "elpd_loo"]
  se <- sqrt(length(d)) * sd(d)
  c(elpd_diff = sum(d), se = se, ratio = sum(d) / se)
}
alts <- fits$label[-1]
st3 <- t(sapply(alts, cmp))
st3 <- data.frame(
  Alternative = fits$name[-1],
  elpd_loo_diff_alt_minus_primary = sprintf("%.2f", st3[, "elpd_diff"]),
  SE = sprintf("%.2f", st3[, "se"]),
  Ratio_diff_over_SE = sprintf("%.2f", st3[, "ratio"]),
  check.names = FALSE
)
print(st3, row.names = FALSE)
write.csv(st3, "tables/supplementary_table_ST3_loo_comparison.csv", row.names = FALSE)
