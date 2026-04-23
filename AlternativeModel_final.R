# package
rm(list=ls(all=TRUE))
library(cmdstanr)
library(dplyr)
library(tidyr)
library(purrr)
library(ggplot2)
library(patchwork)

case_df <- read.csv("Case.csv")  # number of reported syphilis cases in women
case_mat <- t(as.matrix(case_df[, -1]))  

cs_df <- read.csv("CS.csv")  # number of reported CS cases 
cs_mat <- t(as.matrix(cs_df[, -1]))  
if (is.matrix(cs_mat)) {
  cs_mat <- as.vector(cs_mat)  
}

# bound defined
lower_bound <- c(1, 21, 41, 61, 81, 101)
upper_bound <- c(20, 40, 60, 80, 100, 120)

C <- nrow(case_mat)
T <- ncol(case_mat) # T=2025 Q3
K <- max(upper_bound)

# ensure the length of bound match with rows of case_mat
stopifnot(length(lower_bound) == C, length(upper_bound) == C)

cat_of_k <- integer(K)
for (c in seq_len(C)) {
  cat_of_k[lower_bound[c]:upper_bound[c]] <- c
}
stopifnot(all(cat_of_k >= 1 & cat_of_k <= C))

# quarter to year
t_to_year    <- function(t, start_year = 2016L) start_year + (t - 1L) %/% 4L
year_vec_T   <- t_to_year(1:T)

# list of years
years <- 2016:2025

# tribble for ASFR（/1000/year）
asfr_tbl <- tribble(
  ~year, ~c1,  ~c2,  ~c3,  ~c4,   ~c5,  ~c6,
  2016,  3.8,  28.6, 83.5, 102.7, 57.3, 11.4,
  2017,  3.4,  27.5, 82.1, 102.2, 57.5, 11.4,
  2018,  3.1,  26.6, 81.1, 102.0, 57.4, 11.7,
  2019,  2.8,  24.9, 77.2,  98.5, 55.8, 11.7,
  2020,  2.5,  23.0, 74.7,  97.3, 55.3, 11.8,
  2021,  2.1,  20.8, 72.2,  96.2, 55.5, 12.4,
  2022,  1.7,  18.5, 69.6,  93.9, 53.8, 12.2,
  2023,  1.7,  16.8, 65.0,  90.8, 52.4, 12.5,
  2024,  1.8,  14.1, 56.2,  81.6, 48.4, 11.6,
  2025,  1.8,  14.1, 56.2,  81.6, 48.4, 11.6
)

asfr_tbl <- asfr_tbl %>%
  pivot_longer(-year, names_to = "cat", values_to = "asfr_per_1000") %>%
  mutate(cat = as.integer(sub("c", "", cat))) %>%
  arrange(year, cat)

# scale ASFR to quarter
divide_by4 <- TRUE
asfr_tbl <- asfr_tbl %>%
  mutate(asfr_per_woman = (asfr_per_1000 / 1000) / if (divide_by4) 4 else 1)

# year by category
asfr_mat <- asfr_tbl %>%
  select(year, cat, asfr_per_woman) %>%
  pivot_wider(names_from = cat, values_from = asfr_per_woman) %>%
  arrange(year)

# if year is out of range, clip it to the minimum or maximum value
lookup_asfr <- function(year, cat) {
  yr <- pmax(min(asfr_mat$year), pmin(year, max(asfr_mat$year)))
  row <- match(yr, asfr_mat$year)
  asfr_mat[row, as.character(cat), drop = TRUE]
}

# matrix of age K (quarter) by time T (quarter)
q_base <- array(0, dim = c(K, T))
for (t in 1:T) {
  y <- year_vec_T[t]
  for (k in 1:K) {
    c_idx <- cat_of_k[k]
    q_base[k, t] <- lookup_asfr(y, c_idx)
  }
}

# births_per 100000 (2016-2025)
brt_per100000 <- 100000*c(977242^-1, 946146^-1, 918400^-1, 865239^-1, 840835^-1,
                          811622^-1, 770759^-1, 727288^-1, 686061^-1, 525064^-1)

# vector of proportion of reduced birth rate in women with syphilis in scenario assessment
scen_mult <- c(0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0)　

# list of data
stan_data <- list( 
  T = T,
  C = C,
  K = K,
  Y = length(years),
  W = length(scen_mult),
  lower_bound = as.array(lower_bound),
  upper_bound = as.array(upper_bound),
  cat_of_k = as.array(cat_of_k),
  Case = case_mat,
  CS = cs_mat,
  q_base = q_base,
  brt_per100000 = brt_per100000,
  scen_mult = scen_mult
)

mod <- cmdstan_model("syp_AltModel.stan") # alternative model

fit2 <- mod$sample( # for alternative model
  data = stan_data,
  seed = 123,
  chains = 4, parallel_chains = 4,
  iter_warmup = 2000, iter_sampling = 1000,
  adapt_delta = 0.99
)

# summary of parameters
summ_one <- function(fit, vars, probs = c(0.025, 0.5, 0.975)) {
  summarise_draws(
    fit2$draws(variables = vars),
    rhat,
    ess_bulk,
    ess_tail,
    q = function(x) quantile2(x, probs = probs)  
  )
}

tab2 <- summ_one(fit2, vars = c("s_i", "rr", "rep1", "rep2","rep3","rrisk", "p1", "p2", "p3", "cum_cs_w"))
tab2 <- readRDS("tab2_alternative_model.rds")
print(tab2, n = Inf) # alternative model


#　save summary for each CS scenario
## alternative model
# base (cs risk 15%)
sum_case_prd_baseA <- summarise_draws( # for predicted age-dependent syphilis cases
  fit$draws(variables = "case_prd"),
  q2.5  = ~quantile(.x, 0.025),
  q50   = ~quantile(.x, 0.5),
  q97.5 = ~quantile(.x, 0.975)
)
saveRDS(sum_case_prd_baseA, "sum_case_prd_baseA.rds")

sum_cs_prd_baseA <- summarise_draws( # for predicted congenital syphilis cases
  fit$draws(variables = "cs_prd"),
  q2.5  = ~quantile(.x, 0.025),
  q50   = ~quantile(.x, 0.5),
  q97.5 = ~quantile(.x, 0.975)
)
saveRDS(sum_cs_prd_baseA, "sum_cs_prd_baseA.rds")

sum_i_t_baseA <- summarise_draws( # for diagnosed syphilis
  fit$draws(variables = "i_t"),
  q2.5  = ~quantile(.x, 0.025),
  q50   = ~quantile(.x, 0.5),
  q97.5 = ~quantile(.x, 0.975)
)
saveRDS(sum_i_t_baseA, "sum_i_t_baseA.rds")

sum_i_c_baseA <- summarise_draws( # for diagnosed syphilis by age category
  fit$draws(variables = "i_c"),
  q2.5  = ~quantile(.x, 0.025),
  q50   = ~quantile(.x, 0.5),
  q97.5 = ~quantile(.x, 0.975)
)
saveRDS(sum_i_c_baseA, "sum_i_c_baseA.rds")

sum_ip_c_baseA <- summarise_draws( # for diagnosed gestational syphilis by age category
  fit$draws(variables = "ip_c"),
  q2.5  = ~quantile(.x, 0.025),
  q50   = ~quantile(.x, 0.5),
  q97.5 = ~quantile(.x, 0.975)
)
saveRDS(sum_ip_c_baseA, "sum_ip_c_baseA.rds")

sum_gst_pop_baseA <- summarise_draws( # for gestational syphilis per 100,000 births
  fit$draws(variables = "gst_pop"),
  q2.5  = ~quantile(.x, 0.025),
  q50   = ~quantile(.x, 0.5),
  q97.5 = ~quantile(.x, 0.975)
)
saveRDS(sum_gst_pop_baseA, "sum_gst_pop_baseA.rds")

# sensitivity analysis 1 (cs risk 10%)
sum_case_prd_sa1A <- summarise_draws( # for predicted age-dependent syphilis cases
  fit$draws(variables = "case_prd"),
  q2.5  = ~quantile(.x, 0.025),
  q50   = ~quantile(.x, 0.5),
  q97.5 = ~quantile(.x, 0.975)
)
saveRDS(sum_case_prd_sa1A, "sum_case_prd_sa1A.rds")

sum_cs_prd_sa1A <- summarise_draws( # for predicted congenital syphilis cases
  fit$draws(variables = "cs_prd"),
  q2.5  = ~quantile(.x, 0.025),
  q50   = ~quantile(.x, 0.5),
  q97.5 = ~quantile(.x, 0.975)
)
saveRDS(sum_cs_prd_sa1A, "sum_cs_prd_sa1A.rds")

sum_i_t_sa1A <- summarise_draws( # for diagnosed syphilis
  fit$draws(variables = "i_t"),
  q2.5  = ~quantile(.x, 0.025),
  q50   = ~quantile(.x, 0.5),
  q97.5 = ~quantile(.x, 0.975)
)
saveRDS(sum_i_t_sa1A, "sum_i_t_sa1A.rds")

sum_i_c_sa1A <- summarise_draws( # for diagnosed syphilis by age category
  fit$draws(variables = "i_c"),
  q2.5  = ~quantile(.x, 0.025),
  q50   = ~quantile(.x, 0.5),
  q97.5 = ~quantile(.x, 0.975)
)
saveRDS(sum_i_c_sa1A, "sum_i_c_sa1A.rds")

sum_ip_c_sa1A <- summarise_draws( # for diagnosed gestational syphilis by age category
  fit$draws(variables = "ip_c"),
  q2.5  = ~quantile(.x, 0.025),
  q50   = ~quantile(.x, 0.5),
  q97.5 = ~quantile(.x, 0.975)
)
saveRDS(sum_ip_c_sa1A, "sum_ip_c_sa1A.rds")

sum_gst_pop_sa1A <- summarise_draws( # for gestational syphilis per 100,000 births
  fit$draws(variables = "gst_pop"),
  q2.5  = ~quantile(.x, 0.025),
  q50   = ~quantile(.x, 0.5),
  q97.5 = ~quantile(.x, 0.975)
)
saveRDS(sum_gst_pop_sa1A, "sum_gst_pop_sa1A.rds")

# sensitivity analysis 2 (cs risk 20%)
sum_case_prd_sa2A <- summarise_draws( # for predicted age-dependent syphilis cases
  fit$draws(variables = "case_prd"),
  q2.5  = ~quantile(.x, 0.025),
  q50   = ~quantile(.x, 0.5),
  q97.5 = ~quantile(.x, 0.975)
)
saveRDS(sum_case_prd_sa2A, "sum_case_prd_sa2A.rds")

sum_cs_prd_sa2A <- summarise_draws( # for predicted congenital syphilis cases
  fit$draws(variables = "cs_prd"),
  q2.5  = ~quantile(.x, 0.025),
  q50   = ~quantile(.x, 0.5),
  q97.5 = ~quantile(.x, 0.975)
)
saveRDS(sum_cs_prd_sa2A, "sum_cs_prd_sa2A.rds")

sum_i_t_sa2A <- summarise_draws( # for diagnosed syphilis
  fit$draws(variables = "i_t"),
  q2.5  = ~quantile(.x, 0.025),
  q50   = ~quantile(.x, 0.5),
  q97.5 = ~quantile(.x, 0.975)
)
saveRDS(sum_i_t_sa2A, "sum_i_t_sa2A.rds")

sum_i_c_sa2A <- summarise_draws( # for diagnosed syphilis by age category
  fit$draws(variables = "i_c"),
  q2.5  = ~quantile(.x, 0.025),
  q50   = ~quantile(.x, 0.5),
  q97.5 = ~quantile(.x, 0.975)
)
saveRDS(sum_i_c_sa2A, "sum_i_c_sa2A.rds")

sum_ip_c_sa2A <- summarise_draws( # for diagnosed gestational syphilis by age category
  fit$draws(variables = "ip_c"),
  q2.5  = ~quantile(.x, 0.025),
  q50   = ~quantile(.x, 0.5),
  q97.5 = ~quantile(.x, 0.975)
)
saveRDS(sum_ip_c_sa2A, "sum_ip_c_sa2A.rds")

sum_gst_pop_sa2A <- summarise_draws( # for gestational syphilis per 100,000 births
  fit$draws(variables = "gst_pop"),
  q2.5  = ~quantile(.x, 0.025),
  q50   = ~quantile(.x, 0.5),
  q97.5 = ~quantile(.x, 0.975)
)
saveRDS(sum_gst_pop_sa2A, "sum_gst_pop_sa2A.rds")


## WAIC and LOO calculation
library(loo)

# log_lik
draws_array <- fit$draws(c("log_lik_case", "log_lik_cs"))

# change to draws_matrix
draws_mat <- as_draws_matrix(draws_array)

# get name of column
case_names <- grep("^log_lik_case\\[", colnames(draws_mat), value = TRUE)
cs_names   <- grep("^log_lik_cs\\[", colnames(draws_mat), value = TRUE)

# format table
log_lik_case_mat <- draws_mat[, case_names]
log_lik_cs_mat   <- draws_mat[, cs_names]

# bind by column
log_lik_all <- cbind(log_lik_case_mat, log_lik_cs_mat)

# Convert to the input format for loo::waic
log_lik_all_pointwise <- as.matrix(log_lik_all)

# WAIC
waic_result <- waic(log_lik_all_pointwise)
print(waic_result)

# LOO
loo_result <- loo(log_lik_all_pointwise)
print(loo_result)

# visualization of Pareto k diagnosis 
plot(loo_result) # k>0.7 is highlighted by red

# Pareto k
k_values <- loo_result$diagnostics$pareto_k

# identify data with k > 0.7
problematic_points <- which(k_values > 0.7)

# check
print(problematic_points)
print(k_values[problematic_points])

# identify index of data with k > 0.7
T_case <- 39 # number of time per one category 

# vector of index (below is example)
i <- c(50, 176, 200, 206, 241)  

c <- (i - 1) %/% T_case + 1
t <- (i - 1) %% T_case + 1

cat("問題のインデックス", i, "は、カテゴリー", c, "、時点", t, "に対応します\n")



