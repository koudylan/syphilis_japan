# 08_figure1.R
#
# Figure 1: observed and model-predicted quarterly notifications of syphilis in women of
# child-bearing age by age group (a) and of CS (b), primary model (median and 95%
# posterior predictive interval).
#
# The posterior predictive draws case_prd / cs_prd are reconstructed from the primary fit's
# saved draws (results/revision_draws_main.rds: i_t, i_a, rep, rr, p1, p2, p3) with the same
# transforms as the generated quantities block of stan/model_primary.stan, followed by
# Poisson sampling (seed 123). As a check, the reconstructed cs_pop is compared with the
# Stan-computed values in results/revision_summary_main.csv.
#
# Outputs: results/fig1_predictive_draws.rds, figures/figure1.tiff

suppressMessages({library(posterior); library(dplyr); library(tidyr)})
library(ggplot2)
library(patchwork)
library(ragg)

source("R/data_prep.R")
dir.create("figures", showWarnings = FALSE)

## ---- posterior predictive draws ----
r_t <- ifelse(1:T <= 24, p_r_primary[1], p_r_primary[2])
period_of_t <- ifelse(1:T <= 16, 1L, ifelse(1:T <= 24, 2L, 3L))

dr <- readRDS("results/revision_draws_main.rds")
i_t_d  <- as_draws_matrix(subset_draws(dr, variable = "i_t"))     # ndraws x T
i_a_d  <- as_draws_matrix(subset_draws(dr, variable = "i_a"))     # ndraws x K
rep_d  <- as.numeric(as_draws_matrix(subset_draws(dr, variable = "rep")))
rr_d   <- as_draws_matrix(subset_draws(dr, variable = "rr"))      # ndraws x C
p1_d   <- as_draws_matrix(subset_draws(dr, variable = "p1"))      # ndraws x 3
p2_d   <- as_draws_matrix(subset_draws(dr, variable = "p2"))
p3_d   <- as_draws_matrix(subset_draws(dr, variable = "p3"))
ord_it <- order(as.integer(sub("i_t\\[(\\d+)\\]", "\\1", colnames(i_t_d)))); i_t_d <- i_t_d[, ord_it]
ord_ia <- order(as.integer(sub("i_a\\[(\\d+)\\]", "\\1", colnames(i_a_d)))); i_a_d <- i_a_d[, ord_ia]
ord_rr <- order(as.integer(sub("rr\\[(\\d+)\\]", "\\1", colnames(rr_d)))); rr_d <- rr_d[, ord_rr]
nd <- nrow(i_t_d)
cat("ndraws =", nd, "\n")

# g[s, t] for s=1,2,3, built per-draw from p1/p2/p3 according to period_of_t
# e_case_c[c,t] = i_t[t] * sum(i_a[bounds_c]) * rep
sum_ia_c <- sapply(seq_len(C), function(c) rowSums(i_a_d[, lower_bound[c]:upper_bound[c], drop = FALSE]))  # nd x C

set.seed(123)
cs_pop_check <- matrix(NA_real_, nd, length(years))

case_prd_draws <- array(NA_integer_, dim = c(nd, C, T))
cs_prd_draws   <- matrix(NA_integer_, nd, T)

for (i in seq_len(nd)) {
  it <- i_t_d[i, ]; ia <- i_a_d[i, ]; rp <- rep_d[i]; rr <- rr_d[i, ]
  g <- rbind(p1_d[i,], p2_d[i,], p3_d[i,])[period_of_t, , drop = FALSE]  # T x 3 (row=t, col=s)

  # e_case_c
  e_case_c <- sweep(matrix(it, nrow = C, ncol = T, byrow = TRUE), 1, sum_ia_c[i, ], "*") * rp

  # e_cs via 3-term convolution
  e_cs_a <- matrix(0, K, T)
  for (s in 1:3) {
    if (K - (s - 1) < 1 || T - (s - 1) < 1) next
    k_idx <- s:K; kk <- k_idx - (s - 1)
    t_idx <- s:T; tt <- t_idx - (s - 1)
    qb <- q_base[cbind(rep(kk, times = length(tt)), rep(tt, each = length(kk)))]
    qb <- matrix(qb, nrow = length(kk), ncol = length(tt))
    add <- (ia[kk] * rr[cat_of_k[kk]]) * qb
    add <- sweep(add, 2, it[tt] * r_t[tt] * g[tt, s], "*")
    e_cs_a[k_idx, t_idx] <- e_cs_a[k_idx, t_idx] + add
  }
  e_cs <- colSums(e_cs_a)

  case_prd_draws[i, , ] <- matrix(rpois(C * T, as.vector(e_case_c)), C, T)
  cs_prd_draws[i, ] <- rpois(T, e_cs)

  for (y in seq_along(years)) {
    idx <- if (y <= 9) (4*y-3):(4*y) else (4*y-3):(4*y-1)
    cs_pop_check[i, y] <- brt_per100000[y] * sum(e_cs[idx])
  }
  if (i %% 1000 == 0) cat("draw", i, "done\n")
}

# validation vs saved cs_pop summary
main_summ <- read.csv("results/revision_summary_main.csv")
cs_pop_saved <- main_summ %>% filter(grepl("^cs_pop\\[", variable)) %>%
  mutate(idx = as.integer(sub("cs_pop\\[(\\d+)\\]", "\\1", variable))) %>% arrange(idx)
cat("\nValidation: reconstructed vs saved cs_pop (mean):\n")
print(data.frame(year = years, reconstructed = colMeans(cs_pop_check), saved = cs_pop_saved$mean))

saveRDS(list(case_prd = case_prd_draws, cs_prd = cs_prd_draws, case_mat = case_mat, cs_mat = cs_mat,
             years = years, C = C, T = T), "results/fig1_predictive_draws.rds")
cat("\nSaved results/fig1_predictive_draws.rds\n")

## ---- plot ----
pred <- readRDS("results/fig1_predictive_draws.rds")
C <- pred$C; T <- pred$T

## ---- panel a: Case ----
obs_case_df <- as.data.frame(pred$case_mat) %>%
  mutate(c = row_number()) %>%
  pivot_longer(-c, names_to = "t_raw", values_to = "observed") %>%
  group_by(c) %>% mutate(t = row_number(), category = paste0("Cat_", c)) %>%
  ungroup() %>% select(c, t, category, observed)

q_case <- apply(pred$case_prd, c(2, 3), quantile, probs = c(0.025, 0.5, 0.975), na.rm = TRUE)  # 3 x C x T
grid <- expand.grid(c = seq_len(C), t = seq_len(T))
df_case_prd <- grid %>%
  mutate(lower = q_case[1, , ][cbind(c, t)], median = q_case[2, , ][cbind(c, t)], upper = q_case[3, , ][cbind(c, t)]) %>%
  left_join(obs_case_df, by = c("c", "t")) %>% mutate(category = paste0("Cat_", c)) %>% arrange(c, t)

category_labels <- c(Cat_1 = "15-19", Cat_2 = "20-24", Cat_3 = "25-29",
                     Cat_4 = "30-34", Cat_5 = "35-39", Cat_6 = "40-44")
start_year <- 2016; year_breaks <- seq(1, T, by = 4); year_labels <- start_year + (seq_along(year_breaks) - 1)
cb_palette_6 <- c(Cat_1 = "#0072B2", Cat_2 = "#E69F00", Cat_3 = "#009E73",
                  Cat_4 = "#CC79A7", Cat_5 = "#56B4E9", Cat_6 = "#000000")

p1 <- ggplot(df_case_prd, aes(x = t, group = category)) +
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = category), alpha = 0.2, colour = NA) +
  geom_line(aes(y = median, colour = category), linewidth = 1) +
  geom_point(aes(y = observed, colour = category), size = 1.5, na.rm = TRUE) +
  scale_x_continuous(breaks = year_breaks, labels = year_labels) +
  scale_colour_manual(values = cb_palette_6, breaks = names(category_labels), labels = category_labels) +
  scale_fill_manual(values = cb_palette_6, breaks = names(category_labels), labels = category_labels) +
  labs(x = "Year", y = "Persons", colour = "Age group", fill = "Age group") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
        axis.ticks = element_line(colour = "black")) +
  labs(tag = "a") + theme(plot.tag = element_text(face = "bold", size = 14), plot.tag.position = c(0.01, 0.99))

## ---- panel b: CS ----
q_cs <- apply(pred$cs_prd, 2, quantile, probs = c(0.025, 0.5, 0.975), na.rm = TRUE)
df_cs <- tibble(time = seq_len(T), lower = q_cs[1, ], median = q_cs[2, ], upper = q_cs[3, ], observed = pred$cs_mat)

p2 <- ggplot(df_cs, aes(x = time)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), fill = "#56B4E9", alpha = 0.25) +
  geom_line(aes(y = median), color = "#0072B2", linewidth = 1) +
  geom_point(aes(y = observed), color = "#CC79A7", size = 1.5, na.rm = TRUE) +
  scale_x_continuous(breaks = year_breaks, labels = year_labels) +
  labs(x = "Year", y = "Persons") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
        axis.ticks = element_line(colour = "black")) +
  labs(tag = "b") + theme(plot.tag = element_text(face = "bold", size = 14), plot.tag.position = c(0.01, 0.99))

figure1 <- p1 / p2
ggsave(filename = "figures/figure1.tiff", plot = figure1, device = ragg::agg_tiff,
       width = 107, height = 140, units = "mm", dpi = 300, compression = "lzw")
cat("Saved figures/figure1.tiff
")
