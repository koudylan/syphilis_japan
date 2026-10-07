# 03_cs_risk_lhs_full_it.R
#
# Re-fits the 10 Latin hypercube points of 02_cs_risk_lhs.R (design read from
# results/lhs_design_fixedp.csv) with the same settings, keeping the full quarterly series
# i_t[1..39] (02 keeps only i_t[1] and i_t[39]). Used for the CS-risk band in Figure 2a.
# 4 chains x (500 warm-up + 500 sampling), seed 123.
#
# Output: results/summary_fixedp_lhs_full_it.csv
#
# Run time: roughly 6 hours on a 4-core machine.
rm(list=ls(all=TRUE))
library(cmdstanr)
library(dplyr)
library(tidyr)
library(posterior)

source("R/data_prep.R")
dir.create("results", showWarnings = FALSE)
SEED <- 123; CHAINS <- 4; WARMUP <- 500; SAMPLING <- 500

stan_data_base <- c(base_stan_data, list(rr_floor = 1.0))

design <- read.csv("results/lhs_design_fixedp.csv")
print(design)

cat("=== compiling stan/model_primary.stan ===\n")
mod <- cmdstan_model("stan/model_primary.stan")

run_one <- function(p_r_fixed, label) {
  cat(sprintf("\n=== fitting %s: p_r_fixed = (%.4f, %.4f) ===\n", label, p_r_fixed[1], p_r_fixed[2]))
  fit <- mod$sample(data = c(stan_data_base, list(p_r_fixed = as.array(p_r_fixed))),
                    seed = SEED, chains = CHAINS, parallel_chains = CHAINS,
                    iter_warmup = WARMUP, iter_sampling = SAMPLING, refresh = 100)
  ds <- fit$diagnostic_summary()
  summ <- fit$summary("i_t")
  summ$label <- label; summ$p_r1_fixed <- p_r_fixed[1]; summ$p_r2_fixed <- p_r_fixed[2]
  summ$num_divergent <- sum(ds$num_divergent); summ$max_rhat <- max(fit$summary()$rhat, na.rm = TRUE)
  summ
}

out <- list()
for (i in seq_len(nrow(design))) {
  lab <- sprintf("lhs_%02d", design$point[i])
  out[[lab]] <- run_one(c(design$p_r1[i], design$p_r2[i]), lab)
  all_so_far <- do.call(rbind, out)
  write.csv(all_so_far, "results/summary_fixedp_lhs_full_it.csv", row.names = FALSE)
  cat(sprintf("--- saved %d / %d fits ---\n", length(out), nrow(design)))
}
cat("\n=== all i_t refits done ===\n")
