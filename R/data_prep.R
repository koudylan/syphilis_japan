# data_prep.R
#
# Shared data preparation, sourced by every script that fits a model or reconstructs model
# quantities. Run all scripts from the repository root.
#
# Creates:
#   case_mat       notified syphilis cases, women aged 15-44 years [6 age categories x 39 quarters]
#   cs_mat         notified congenital syphilis (CS) cases [39 quarters], 2016 Q1 - 2025 Q3
#   q_base         quarterly baseline fertility per woman [120 fine age indices x 39 quarters],
#                  from age-specific fertility rates (2025 uses the 2024 value)
#   brt_per100000  100,000 / annual live births, 2016-2025 (2025: January-September, provisional)
#   scen_mult      scenario proportions for the fertility-rate ratio scenario analysis (Table 2)
#   p_r_primary    CS risk used in the primary analysis: 11/76 (2016-2021), 28/375 (2022-2025)
#   base_stan_data data list shared by all Stan models (model-specific inputs such as
#                  p_r_fixed and rr_floor are added by each script)

library(dplyr)
library(tidyr)

case_df <- read.csv("data/Case.csv")
case_mat <- t(as.matrix(case_df[, -1]))

cs_df <- read.csv("data/CS.csv")
cs_mat <- t(as.matrix(cs_df[, -1]))
if (is.matrix(cs_mat)) {
  cs_mat <- as.vector(cs_mat)
}

# six 5-year age categories (15-19, ..., 40-44), each split into 20 fine age indices
lower_bound <- c(1, 21, 41, 61, 81, 101)
upper_bound <- c(20, 40, 60, 80, 100, 120)

C <- nrow(case_mat)
T <- ncol(case_mat)
K <- max(upper_bound)

stopifnot(length(lower_bound) == C, length(upper_bound) == C)

cat_of_k <- integer(K)
for (c in seq_len(C)) {
  cat_of_k[lower_bound[c]:upper_bound[c]] <- c
}
stopifnot(all(cat_of_k >= 1 & cat_of_k <= C))

t_to_year  <- function(t, start_year = 2016L) start_year + (t - 1L) %/% 4L
year_vec_T <- t_to_year(1:T)

years <- 2016:2025

# age-specific fertility rate per 1,000 women (vital statistics; 2025 = 2024 value)
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

# annual rate per 1,000 -> quarterly rate per woman
asfr_tbl <- asfr_tbl %>%
  mutate(asfr_per_woman = (asfr_per_1000 / 1000) / 4)

asfr_mat <- asfr_tbl %>%
  select(year, cat, asfr_per_woman) %>%
  pivot_wider(names_from = cat, values_from = asfr_per_woman) %>%
  arrange(year)

lookup_asfr <- function(year, cat) {
  yr <- pmax(min(asfr_mat$year), pmin(year, max(asfr_mat$year)))
  row <- match(yr, asfr_mat$year)
  asfr_mat[row, as.character(cat), drop = TRUE]
}

q_base <- array(0, dim = c(K, T))
for (t in 1:T) {
  y <- year_vec_T[t]
  for (k in 1:K) {
    c_idx <- cat_of_k[k]
    q_base[k, t] <- lookup_asfr(y, c_idx)
  }
}

# annual live births, 2016-2025 (2025: January-September)
brt_per100000 <- 100000*c(977242^-1, 946146^-1, 918400^-1, 865239^-1, 840835^-1,
                          811622^-1, 770759^-1, 727288^-1, 686061^-1, 525064^-1)

scen_mult <- c(0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0)

# CS risk surveys (Japan Society of Obstetrics and Gynecology): CS cases / pregnancies
# with syphilis at >= 22 weeks
k_2016 <- 11; n_2016 <- 76     # Suzuki et al., applied to 2016-2021
k_2022 <- 28; n_2022 <- 375    # Hayata et al., applied to 2022-2025
p_r_primary <- c(k_2016 / n_2016, k_2022 / n_2022)

base_stan_data <- list(
  T = T, C = C, K = K, Y = length(years), W = length(scen_mult),
  lower_bound = as.array(lower_bound),
  upper_bound = as.array(upper_bound),
  cat_of_k = as.array(cat_of_k),
  Case = case_mat, CS = cs_mat, q_base = q_base,
  brt_per100000 = brt_per100000, scen_mult = scen_mult
)
