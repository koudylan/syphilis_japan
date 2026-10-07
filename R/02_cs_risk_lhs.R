# 02_cs_risk_lhs.R
#
# CS-risk sensitivity analysis (Methods 2.4): the primary model is re-fitted at 10 Latin
# hypercube points for the two CS-risk proportions, sampled from their Beta distributions
#   p_r[1] (2016-2021) ~ Beta(11 + 1, 76 - 11 + 1)   = Beta(12, 66)
#   p_r[2] (2022-2025) ~ Beta(28 + 1, 375 - 28 + 1)  = Beta(29, 348)
# Each margin is stratified into 10 equal-probability strata (one draw per stratum), and
# the strata of p_r[2] are randomly permuted before pairing with p_r[1].
# 4 chains x (500 warm-up + 500 sampling), seed 123.
#
# Outputs:
#   results/lhs_design_fixedp.csv   the 10 (p_r[1], p_r[2]) design points
#   results/summary_fixedp_lhs.csv  posterior summaries at each point (Table 1 range column,
#                                   ranges quoted in Results 3.1, Figure 3a band)
#
# Run time: roughly 6 hours on a 4-core machine.
rm(list=ls(all=TRUE))
library(cmdstanr)
library(dplyr)
library(tidyr)
library(posterior)

source("R/data_prep.R")
dir.create("results", showWarnings = FALSE)

N_LHS <- 10
CHAINS <- 4
WARMUP <- 500
SAMPLING <- 500
SEED <- 123

stan_data_base <- c(base_stan_data, list(rr_floor = 1.0))

cat("=== compiling stan/model_primary.stan ===\n")
mod <- cmdstan_model("stan/model_primary.stan")

run_one <- function(p_r_fixed, label) {
  cat(sprintf("\n=== fitting %s: p_r_fixed = (%.4f, %.4f) ===\n", label, p_r_fixed[1], p_r_fixed[2]))
  stan_data <- c(stan_data_base, list(p_r_fixed = as.array(p_r_fixed)))
  fit <- mod$sample(
    data = stan_data,
    seed = SEED,
    chains = CHAINS, parallel_chains = CHAINS,
    iter_warmup = WARMUP, iter_sampling = SAMPLING,
    refresh = 100
  )
  ds <- fit$diagnostic_summary()
  summ <- fit$summary(c("s_i", "rep", "rr", "rrisk", "i_t[1]", "i_t[39]", "p1", "p2", "p3",
                         "cum_cs_w", "cs_pop", "gst_pop"))
  summ$label <- label
  summ$p_r1_fixed <- p_r_fixed[1]
  summ$p_r2_fixed <- p_r_fixed[2]
  summ$num_divergent <- sum(ds$num_divergent)
  summ$max_rhat <- max(fit$summary()$rhat, na.rm = TRUE)
  summ
}

## ---- Latin hypercube design ----
set.seed(SEED)
a1 <- k_2016 + 1; b1 <- n_2016 - k_2016 + 1   # Beta(12, 66) for p_r[1] (2016-2021)
a2 <- k_2022 + 1; b2 <- n_2022 - k_2022 + 1   # Beta(29, 348) for p_r[2] (2022-2025)

lhs_uniform <- function(n) ((seq_len(n) - 1) + runif(n)) / n
u1 <- lhs_uniform(N_LHS)
u2 <- sample(lhs_uniform(N_LHS))  # independent permutation -> space-filling 2D design

p_r1_lhs <- qbeta(u1, a1, b1)
p_r2_lhs <- qbeta(u2, a2, b2)

lhs_design <- data.frame(point = 1:N_LHS, p_r1 = p_r1_lhs, p_r2 = p_r2_lhs)
print(lhs_design)
write.csv(lhs_design, "results/lhs_design_fixedp.csv", row.names = FALSE)

## ---- fits ----
lhs_summaries <- list()
for (i in 1:N_LHS) {
  lhs_summaries[[i]] <- run_one(c(p_r1_lhs[i], p_r2_lhs[i]), sprintf("lhs_%02d", i))
  lhs_all <- do.call(rbind, lhs_summaries)
  write.csv(lhs_all, "results/summary_fixedp_lhs.csv", row.names = FALSE)
}
cat("\n=== LHS sensitivity analysis done ===\n")
