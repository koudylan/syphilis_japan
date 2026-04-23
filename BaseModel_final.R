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

mod <- cmdstan_model("syp_BaseModel_ver1.stan")  

fit <- mod$sample(
  data = stan_data,
  seed = 123,
  chains = 4, parallel_chains = 4,
  iter_warmup = 1000, iter_sampling = 1000
)

# summary of parameters for each CS assumption
summ_one <- function(fit, vars, probs = c(0.025, 0.5, 0.975)) {
  summarise_draws(
    fit$draws(variables = vars),
    rhat,
    ess_bulk,
    ess_tail,
    q = function(x) quantile2(x, probs = probs)  
  )
}
tab1 <- summ_one(fit, vars = c("s_i", "rr", "rep","rrisk", "p1", "p2", "p3", "cum_cs_w"))
saveRDS(tab1, "tab1_base.rds")
tab1 <- readRDS("tab1_base.rds")
print(tab1, n = Inf)
# tab1 <- summ_one(fit, vars = c("s_i", "rr", "rep","rrisk", "p1", "p2", "p3", "cum_cs_w"))
# saveRDS(tab1, "tab1_sa1.rds")
# tab1 <- readRDS("tab1_sa1.rds")
# print(tab1, n = Inf)
# tab1 <- summ_one(fit, vars = c("s_i", "rr", "rep","rrisk", "p1", "p2", "p3", "cum_cs_w"))
# saveRDS(tab1, "tab1_sa2.rds")
# tab1 <- readRDS("tab1_sa2.rds")
# print(tab1, n = Inf)

# save summary for each CS scenario
# base (cs risk 15%)
sum_case_prd_base <- summarise_draws( # for predicted age-dependent syphilis cases
  fit$draws(variables = "case_prd"),
  q2.5  = ~quantile(.x, 0.025),
  q50   = ~quantile(.x, 0.5),
  q97.5 = ~quantile(.x, 0.975)
)
saveRDS(sum_case_prd_base, "sum_case_prd_base.rds")

sum_cs_prd_base <- summarise_draws( # for predicted congenital syphilis cases
  fit$draws(variables = "cs_prd"),
  q2.5  = ~quantile(.x, 0.025),
  q50   = ~quantile(.x, 0.5),
  q97.5 = ~quantile(.x, 0.975)
)
saveRDS(sum_cs_prd_base, "sum_cs_prd_base.rds")

sum_i_t_base <- summarise_draws( # for diagnosed syphilis
  fit$draws(variables = "i_t"),
  q2.5  = ~quantile(.x, 0.025),
  q50   = ~quantile(.x, 0.5),
  q97.5 = ~quantile(.x, 0.975)
)
saveRDS(sum_i_t_base, "sum_i_t_base.rds")

sum_i_c_base <- summarise_draws( # for diagnosed syphilis by age category
  fit$draws(variables = "i_c"),
  q2.5  = ~quantile(.x, 0.025),
  q50   = ~quantile(.x, 0.5),
  q97.5 = ~quantile(.x, 0.975)
)
saveRDS(sum_i_c_base, "sum_i_c_base.rds")

sum_ip_c_base <- summarise_draws( # for diagnosed gestational syphilis by age category
  fit$draws(variables = "ip_c"),
  q2.5  = ~quantile(.x, 0.025),
  q50   = ~quantile(.x, 0.5),
  q97.5 = ~quantile(.x, 0.975)
)
saveRDS(sum_ip_c_base, "sum_ip_c_base.rds")

sum_gst_pop_base <- summarise_draws( # for gestational syphilis per 100,000 births
  fit$draws(variables = "gst_pop"),
  q2.5  = ~quantile(.x, 0.025),
  q50   = ~quantile(.x, 0.5),
  q97.5 = ~quantile(.x, 0.975)
)
saveRDS(sum_gst_pop_base, "sum_gst_pop_base.rds")

sum_cum_cs_w_base <- summarise_draws( # 
  fit$draws(variables = "cum_cs_w"),
  q2.5  = ~quantile(.x, 0.025),
  q50   = ~quantile(.x, 0.5),
  q97.5 = ~quantile(.x, 0.975)
)
saveRDS(sum_cum_cs_w_base, "sum_cum_cs_w_base.rds")

# sensitivity analysis 1 (cs risk 10%)
sum_case_prd_sa1 <- summarise_draws( # for predicted age-dependent syphilis cases
  fit$draws(variables = "case_prd"),
  q2.5  = ~quantile(.x, 0.025),
  q50   = ~quantile(.x, 0.5),
  q97.5 = ~quantile(.x, 0.975)
)
saveRDS(sum_case_prd_sa1, "sum_case_prd_sa1.rds")

sum_cs_prd_sa1 <- summarise_draws( # for predicted congenital syphilis cases
  fit$draws(variables = "cs_prd"),
  q2.5  = ~quantile(.x, 0.025),
  q50   = ~quantile(.x, 0.5),
  q97.5 = ~quantile(.x, 0.975)
)
saveRDS(sum_cs_prd_sa1, "sum_cs_prd_sa1.rds")

sum_i_t_sa1 <- summarise_draws( # for diagnosed syphilis
  fit$draws(variables = "i_t"),
  q2.5  = ~quantile(.x, 0.025),
  q50   = ~quantile(.x, 0.5),
  q97.5 = ~quantile(.x, 0.975)
)
saveRDS(sum_i_t_sa1, "sum_i_t_sa1.rds")

sum_i_c_sa1 <- summarise_draws( # for diagnosed syphilis by age category
  fit$draws(variables = "i_c"),
  q2.5  = ~quantile(.x, 0.025),
  q50   = ~quantile(.x, 0.5),
  q97.5 = ~quantile(.x, 0.975)
)
saveRDS(sum_i_c_sa1, "sum_i_c_sa1.rds")

sum_ip_c_sa1 <- summarise_draws( # for diagnosed gestational syphilis by age category
  fit$draws(variables = "ip_c"),
  q2.5  = ~quantile(.x, 0.025),
  q50   = ~quantile(.x, 0.5),
  q97.5 = ~quantile(.x, 0.975)
)
saveRDS(sum_ip_c_sa1, "sum_ip_c_sa1.rds")

sum_gst_pop_sa1 <- summarise_draws( # for gestational syphilis per 100,000 births
  fit$draws(variables = "gst_pop"),
  q2.5  = ~quantile(.x, 0.025),
  q50   = ~quantile(.x, 0.5),
  q97.5 = ~quantile(.x, 0.975)
)
saveRDS(sum_gst_pop_sa1, "sum_gst_pop_sa1.rds")

# sensitivity analysis 2 (cs risk 20%)
sum_case_prd_sa2 <- summarise_draws( # for predicted age-dependent syphilis cases
  fit$draws(variables = "case_prd"),
  q2.5  = ~quantile(.x, 0.025),
  q50   = ~quantile(.x, 0.5),
  q97.5 = ~quantile(.x, 0.975)
)
saveRDS(sum_case_prd_sa2, "sum_case_prd_sa2.rds")

sum_cs_prd_sa2 <- summarise_draws( # for predicted congenital syphilis cases
  fit$draws(variables = "cs_prd"),
  q2.5  = ~quantile(.x, 0.025),
  q50   = ~quantile(.x, 0.5),
  q97.5 = ~quantile(.x, 0.975)
)
saveRDS(sum_cs_prd_sa2, "sum_cs_prd_sa2.rds")

sum_i_t_sa2 <- summarise_draws( # for diagnosed syphilis
  fit$draws(variables = "i_t"),
  q2.5  = ~quantile(.x, 0.025),
  q50   = ~quantile(.x, 0.5),
  q97.5 = ~quantile(.x, 0.975)
)
saveRDS(sum_i_t_sa2, "sum_i_t_sa2.rds")

sum_i_c_sa2 <- summarise_draws( # for diagnosed syphilis by age category
  fit$draws(variables = "i_c"),
  q2.5  = ~quantile(.x, 0.025),
  q50   = ~quantile(.x, 0.5),
  q97.5 = ~quantile(.x, 0.975)
)
saveRDS(sum_i_c_sa2, "sum_i_c_sa2.rds")

sum_ip_c_sa2 <- summarise_draws( # for diagnosed gestational syphilis by age category
  fit$draws(variables = "ip_c"),
  q2.5  = ~quantile(.x, 0.025),
  q50   = ~quantile(.x, 0.5),
  q97.5 = ~quantile(.x, 0.975)
)
saveRDS(sum_ip_c_sa2, "sum_ip_c_sa2.rds")

sum_gst_pop_sa2 <- summarise_draws( # for gestational syphilis per 100,000 births
  fit$draws(variables = "gst_pop"),
  q2.5  = ~quantile(.x, 0.025),
  q50   = ~quantile(.x, 0.5),
  q97.5 = ~quantile(.x, 0.975)
)
saveRDS(sum_gst_pop_sa2, "sum_gst_pop_sa2.rds")


## Figure 1A: observed and predicted notified cases
# ---- observed data (Case.csv -> tidy) ----
obs_case_df <- as.data.frame(case_mat) %>%
  mutate(c = row_number()) %>%
  pivot_longer(-c, names_to = "t_raw", values_to = "observed") %>%
  group_by(c) %>%
  mutate(
    t = row_number(),
    category = paste0("Cat_", c)
  ) %>%
  ungroup() %>%
  select(c, t, category, observed)

# "case_prd[c,t]" -> tidy DF (c, t, lower, median, upper)
readRDS("sum_case_prd_base.rds")
idx <- regmatches(
  sum_case_prd_base$variable,
  regexec("^case_prd\\[(\\d+),(\\d+)\\]$", sum_case_prd_base$variable)
)
idx <- do.call(rbind, lapply(idx, `[`, 2:3))
idx <- apply(idx, 2, as.integer)

df_case_prd <- tibble(
  c = idx[, 1],
  t = idx[, 2],
  lower  = sum_case_prd_base$'2.5%',
  median = sum_case_prd_base$'50%',
  upper  = sum_case_prd_base$'97.5%'
) %>%
  left_join(obs_case_df, by = c("c", "t")) %>%
  mutate(category = paste0("Cat_", c)) %>%
  arrange(c, t)

# ---- labels ----
category_labels <- c(
  Cat_1 = "15-19",
  Cat_2 = "20-24",
  Cat_3 = "25-29",
  Cat_4 = "30-34",
  Cat_5 = "35-39",
  Cat_6 = "40-44"
)

start_year  <- 2016
year_breaks <- seq(1, T, by = 4)
year_labels <- start_year + (seq_along(year_breaks) - 1)

# ---- plot ----
cb_palette_6 <- c(
  Cat_1 = "#0072B2", # blue
  Cat_2 = "#E69F00", # orange
  Cat_3 = "#009E73", # bluish green
  Cat_4 = "#CC79A7", # magenta (red replacement)
  Cat_5 = "#56B4E9", # sky blue
  Cat_6 = "#000000"  # black
)

p1 <- ggplot(df_case_prd, aes(x = t, group = category)) +
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = category), alpha = 0.2, colour = NA) +
  geom_line(aes(y = median, colour = category), linewidth = 1) +
  geom_point(aes(y = observed, colour = category), size = 1.5, na.rm = TRUE) +
  scale_x_continuous(breaks = year_breaks, labels = year_labels) +
  scale_colour_manual(
    values = cb_palette_6,
    breaks = names(category_labels),
    labels = category_labels
  ) +
  scale_fill_manual(
    values = cb_palette_6,
    breaks = names(category_labels),
    labels = category_labels
  ) +
  labs(
    x = "Year", y = "Persons",
    colour = "Age group", fill = "Age group"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black"),
    axis.ticks = element_line(colour = "black")
  )+
  labs(tag = "a") +
  theme(
    plot.tag = element_text(face = "bold", size = 14),
    plot.tag.position = c(0.01, 0.99)
  )

## Figure 1B:observed and predicted congenital syphilis 
# "cs_prd[t]" -> tidy DF + observed
df_cs <- tibble(
  time = as.integer(sub("^cs_prd\\[(\\d+)\\]$", "\\1", sum_cs_prd_base$variable)),
  lower  = sum_cs_prd_base$'2.5%',
  median = sum_cs_prd_base$'50%',
  upper  = sum_cs_prd_base$'97.5%'
) %>%
  arrange(time) %>%
  mutate(observed = cs_mat)

# ---- plot ----
p2 <- ggplot(df_cs, aes(x = time)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), fill = "#56B4E9", alpha = 0.25) +  # sky blue
  geom_line(aes(y = median), color = "#0072B2", linewidth = 1) +                   # blue
  geom_point(aes(y = observed), color = "#CC79A7", size = 1.5, na.rm = TRUE) +     # magenta (red replacement)
  scale_x_continuous(breaks = year_breaks, labels = year_labels) +
  labs(x = "Year", y = "Persons") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black"),
    axis.ticks = element_line(colour = "black")
  )+
  labs(tag = "b") +
  theme(
    plot.tag = element_text(face = "bold", size = 14),
    plot.tag.position = c(0.01, 0.99)
  )


figure1 <- p1 / p2

# save as tiff
library(ragg)
ggsave(
  filename = "figure1.tiff",
  plot = figure1,
  device = ragg::agg_tiff,
  width = 107,
  height = 140,
  units = "mm",
  dpi = 300,
  compression = "lzw"
)

## Figure 2a: estimated diagnosed syphilis for all CS risk assumptions
sum_base <- readRDS("sum_i_t_base.rds") 
sum_sa1 <- readRDS("sum_i_t_sa1.rds") 
sum_sa2 <- readRDS("sum_i_t_sa2.rds") 
to_i_t_df <- function(sum_df, model_label) { 
  data.frame( 
    time = as.integer(sub("^i_t\\[(\\d+)\\]$", "\\1", sum_df$variable)), 
    lower = sum_df[["2.5%"]], 
    median = sum_df[["50%"]], 
    upper = sum_df[["97.5%"]], model = model_label ) } 
df_i <- rbind( 
  to_i_t_df(sum_base, "15%"), 
  to_i_t_df(sum_sa1, "10%"), 
  to_i_t_df(sum_sa2, "20%") )

# Color-blind safe + high contrast + red→magenta
csrisk_palette <- c(
  "10%" = "#0072B2",  # blue
  "15%" = "#CC79A7",  # magenta (red replacement)
  "20%" = "#009E73"   # bluish green
)

p3 <- ggplot(df_i, aes(x = time, group = model)) +
  geom_ribbon(
    aes(ymin = lower, ymax = upper, fill = model),
    alpha = 0.15, colour = NA
  ) +
  geom_line(
    aes(y = median, colour = model),
    linewidth = 1
  ) +
  scale_colour_manual(
    values = csrisk_palette,
    breaks = c("10%", "15%", "20%")
  ) +
  scale_fill_manual(
    values = csrisk_palette,
    breaks = c("10%", "15%", "20%")
  ) +
  scale_x_continuous(breaks = year_breaks, labels = year_labels) +
  labs(
    x = "Year", y = "Quarterly number of \nnew diagnosis",
    colour = "CS risk", fill = "CS risk"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black"),
    axis.ticks = element_line(colour = "black")
  ) +
  labs(tag = "a") +
  theme(
    plot.tag = element_text(face = "bold", size = 14),
    plot.tag.position = c(0.001, 0.99)
  )
p3

## Figure 2b and supplement: estimated diagnosed syphilis by age category 
# "i_c[c,t]" -> tidy DF (c, time, lower, median, upper)
sum_i_c_base <- readRDS("sum_i_c_base.rds") 
sum_i_c_sa1 <- readRDS("sum_i_c_sa1.rds") 
sum_i_c_sa2 <- readRDS("sum_i_c_sa2.rds") 
idx <- regmatches(
  sum_i_c_base$variable,
  regexec("^i_c\\[(\\d+),(\\d+)\\]$", sum_i_c_base$variable)
)
# idx <- regmatches(
#   sum_i_c_sa1$variable,
#   regexec("^i_c\\[(\\d+),(\\d+)\\]$", sum_i_c_sa1$variable)
# )
# idx <- regmatches(
#   sum_i_c_sa2$variable,
#   regexec("^i_c\\[(\\d+),(\\d+)\\]$", sum_i_c_sa2$variable)
# )
idx <- do.call(rbind, lapply(idx, `[`, 2:3))
idx <- apply(idx, 2, as.integer)

df_i_c <- tibble(
  c = idx[, 1],
  time = idx[, 2],
  lower  = sum_i_c_base$'2.5%',
  median = sum_i_c_base$'50%',
  upper  = sum_i_c_base$'97.5%'
) %>%
  mutate(category = paste0("Cat_", c)) %>%
  arrange(c, time)
# df_i_c <- tibble(
#   c = idx[, 1],
#   time = idx[, 2],
#   lower  = sum_i_c_sa1$'2.5%',
#   median = sum_i_c_sa1$'50%',
#   upper  = sum_i_c_sa1$'97.5%'
# ) %>%
#   mutate(category = paste0("Cat_", c)) %>%
#   arrange(c, time)
# df_i_c <- tibble(
#   c = idx[, 1],
#   time = idx[, 2],
#   lower  = sum_i_c_sa2$'2.5%',
#   median = sum_i_c_sa2$'50%',
#   upper  = sum_i_c_sa2$'97.5%'
# ) %>%
#   mutate(category = paste0("Cat_", c)) %>%
#   arrange(c, time)

# ---- labels ----
# ---- color-blind safe palette (6 categories) ----
cb_palette_6 <- c(
  Cat_1 = "#0072B2", # blue
  Cat_2 = "#E69F00", # orange
  Cat_3 = "#009E73", # bluish green
  Cat_4 = "#CC79A7", # magenta (red replacement)
  Cat_5 = "#56B4E9", # sky blue
  Cat_6 = "#000000"  # black
)

# Fix the legend order and plotting order
df_i_c <- df_i_c %>%
  mutate(category = factor(category, levels = names(category_labels)))

# ---- plot ----
p4sa2 <- ggplot(df_i_c, aes(x = time, group = category)) +
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = category), alpha = 0.2, colour = NA) +
  geom_line(aes(y = median, colour = category), linewidth = 1) +
  scale_x_continuous(breaks = year_breaks, labels = year_labels) +
  scale_colour_manual(
    values = cb_palette_6,
    breaks = names(category_labels),
    labels = category_labels
  ) +
  scale_fill_manual(
    values = cb_palette_6,
    breaks = names(category_labels),
    labels = category_labels
  ) +
  labs(
    x = "Year", y = "Quarterly number of \nnew diagnosis",
    colour = "Age group", fill = "Age group"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black"),
    axis.ticks = element_line(colour = "black")
  ) +
  labs(tag = "b") +
  theme(
    plot.tag = element_text(face = "bold", size = 14),
    plot.tag.position = c(0.001, 0.99)
  )
p4
figure2 <- p3 / p4 # Main
figure2bSup <- p4sa1 / p4sa2 # Supplement (same plot for all CS risk assumptions)

ggsave(
  filename = "figure2.tiff",
  plot = figure2,
  device = ragg::agg_tiff,
  width = 107,
  height = 140,
  units = "mm",
  dpi = 300,
  compression = "lzw"
)
ggsave(
  filename = "figure2bSup.tiff",
  plot = figure2bSup,
  device = ragg::agg_tiff,
  width = 107,
  height =140,
  units = "mm",
  dpi = 300,
  compression = "lzw"
)


## Figure 3a and supplement: diagnosed gestational syphilis by age category 
# "ip_c[c,t]" -> tidy DF (c, time, lower, median, upper)
sum_ip_c_base <-readRDS("sum_ip_c_base.rds") 
sum_ip_c_sa1 <- readRDS("sum_ip_c_sa1.rds") 
sum_ip_c_sa2 <- readRDS("sum_ip_c_sa2.rds") 
idx <- regmatches(
  sum_ip_c_base$variable,
  regexec("^ip_c\\[(\\d+),(\\d+)\\]$", sum_ip_c_base$variable)
)
# idx <- regmatches(
#   sum_ip_c_sa1$variable,
#   regexec("^ip_c\\[(\\d+),(\\d+)\\]$", sum_ip_c_sa1$variable)
# )
# idx <- regmatches(
#   sum_ip_c_sa2$variable,
#   regexec("^ip_c\\[(\\d+),(\\d+)\\]$", sum_ip_c_sa2$variable)
# )
idx <- do.call(rbind, lapply(idx, `[`, 2:3))
idx <- apply(idx, 2, as.integer)

df_ip_c <- tibble(
  c = idx[, 1],
  time = idx[, 2],
  lower  = sum_ip_c_base$'2.5%',
  median = sum_ip_c_base$'50%',
  upper  = sum_ip_c_base$'97.5%'
) %>%
  mutate(category = paste0("Cat_", c)) %>%
  arrange(c, time)
# df_ip_c <- tibble(
#   c = idx[, 1],
#   time = idx[, 2],
#   lower  = sum_ip_c_sa1$'2.5%',
#   median = sum_ip_c_sa1$'50%',
#   upper  = sum_ip_c_sa1$'97.5%'
# ) %>%
#   mutate(category = paste0("Cat_", c)) %>%
#   arrange(c, time)
# df_ip_c <- tibble(
#   c = idx[, 1],
#   time = idx[, 2],
#   lower  = sum_ip_c_sa2$'2.5%',
#   median = sum_ip_c_sa2$'50%',
#   upper  = sum_ip_c_sa2$'97.5%'
# ) %>%
#   mutate(category = paste0("Cat_", c)) %>%
#   arrange(c, time)

# ---- labels ----
# ---- color-blind safe palette (6 categories) ----
cb_palette_6 <- c(
  Cat_1 = "#0072B2", # blue
  Cat_2 = "#E69F00", # orange
  Cat_3 = "#009E73", # bluish green
  Cat_4 = "#CC79A7", # magenta (red replacement)
  Cat_5 = "#56B4E9", # sky blue
  Cat_6 = "#000000"  # black
)

# Fix the legend order and plotting order
df_ip_c <- df_ip_c %>%
  mutate(category = factor(category, levels = names(category_labels)))

# ---- plot ----
p5sa2 <- ggplot(df_ip_c, aes(x = time, group = category)) +
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = category), alpha = 0.2, colour = NA) +
  geom_line(aes(y = median, colour = category), linewidth = 1) +
  scale_x_continuous(breaks = year_breaks, labels = year_labels) +
  scale_colour_manual(
    values = cb_palette_6,
    breaks = names(category_labels),
    labels = category_labels
  ) +
  scale_fill_manual(
    values = cb_palette_6,
    breaks = names(category_labels),
    labels = category_labels
  ) +
  labs(
    x = "Year", y = "Quarterly number of \ngestational syphilis",
    colour = "Age group", fill = "Age group"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black"),
    axis.ticks = element_line(colour = "black")
  ) +
  labs(tag = "b") +
  theme(
    plot.tag = element_text(face = "bold", size = 14),
    plot.tag.position = c(0.001, 0.99)
  )
p5
figure3aSup <- p5sa1 / p5sa2 # Supplement (same plot for all CS risk assumptions)

ggsave(
  filename = "figure3aSup.tiff",
  plot = figure3aSup,
  device = ragg::agg_tiff,
  width = 107,
  height = 140,
  units = "mm",
  dpi = 300,
  compression = "lzw"
)

## Figure 3b: gestational syphilis per 100,000 births for all CS risk assumptions
sum_base <- readRDS("sum_gst_pop_base.rds")
sum_sa1  <- readRDS("sum_gst_pop_sa1.rds")
sum_sa2  <- readRDS("sum_gst_pop_sa2.rds")

to_gst_pop_df <- function(sum_df, model_label) { 
  data.frame( 
    time  = 2015 + as.integer(sub("^gst_pop\\[(\\d+)\\]$", "\\1", sum_df$variable)),  # 1->2016, ..., 10->2025
    lower = sum_df[["2.5%"]], 
    median = sum_df[["50%"]], 
    upper = sum_df[["97.5%"]], 
    model = model_label
  ) 
}

df_gst_pop <- rbind( 
  to_gst_pop_df(sum_base, "15%"), 
  to_gst_pop_df(sum_sa1, "10%"), 
  to_gst_pop_df(sum_sa2, "20%")
)

# Fix the legend order 
df_gst_pop$model <- factor(df_gst_pop$model, levels = c("10%", "15%", "20%"))

p6 <- ggplot(df_gst_pop, aes(x = time, group = model)) +
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = model), alpha = 0.15, colour = NA) +
  geom_line(aes(y = median, colour = model), linewidth = 1) +
  scale_colour_manual(values = csrisk_palette, breaks = c("10%", "15%", "20%")) +
  scale_fill_manual(values = csrisk_palette, breaks = c("10%", "15%", "20%")) +
  scale_x_continuous(
    breaks = 2016:2025,
    labels = 2016:2025
  ) +
  labs(
    x = "Year", y = "Yearly number of \ngestational syphilis per 100,000 birth",
    colour = "CS risk", fill = "CS risk"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black"),
    axis.ticks = element_line(colour = "black")
  ) +
  labs(tag = "b") +
  theme(
    plot.tag = element_text(face = "bold", size = 14),
    plot.tag.position = c(0.001, 0.99)
  )

p6

figure3 <- p5 / p6

library(ragg)
ggsave(
  filename = "figure3.tiff",
  plot = figure3,
  device = ragg::agg_tiff,
  width = 107,
  height = 140,
  units = "mm",
  dpi = 300,
  compression = "lzw"
)


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



