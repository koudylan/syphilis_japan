# 09_figures2_3.R
#
# Figure 2: (a) quarterly number of newly diagnosed syphilis cases in women of child-bearing
#           age, primary model (median and 95% CrI) with the range of posterior medians
#           across the 10-point CS-risk sensitivity analysis
#           (results/summary_fixedp_lhs_full_it.csv); (b) the same by age group.
# Figure 3: (a) yearly gestational syphilis per 100,000 live births, primary model with the
#           CS-risk sensitivity range (results/summary_fixedp_lhs.csv); (b) quarterly
#           gestational syphilis by age group.
# Panels b are reconstructed from the primary fit's saved draws
# (results/revision_draws_main.rds: i_t, i_a, rr) with the same transforms as the generated
# quantities block of stan/model_primary.stan (i_c, ip_c).
#
# Outputs: figures/figure2.tiff, figures/figure3.tiff

library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)
library(ragg)
library(posterior)

source("R/data_prep.R")
dir.create("figures", showWarnings = FALSE)

category_labels <- c(Cat_1 = "15-19", Cat_2 = "20-24", Cat_3 = "25-29",
                     Cat_4 = "30-34", Cat_5 = "35-39", Cat_6 = "40-44")
start_year <- 2016; year_breaks <- seq(1, T, by = 4); year_labels <- start_year + (seq_along(year_breaks) - 1)
cb_palette_6 <- c(Cat_1 = "#0072B2", Cat_2 = "#E69F00", Cat_3 = "#009E73",
                  Cat_4 = "#CC79A7", Cat_5 = "#56B4E9", Cat_6 = "#000000")

dr <- readRDS("results/revision_draws_main.rds")
i_t <- as_draws_matrix(subset_draws(dr, variable="i_t")); i_t <- i_t[, order(as.integer(sub("i_t\\[(\\d+)\\]","\\1",colnames(i_t))))]
i_a <- as_draws_matrix(subset_draws(dr, variable="i_a")); i_a <- i_a[, order(as.integer(sub("i_a\\[(\\d+)\\]","\\1",colnames(i_a))))]
rr  <- as_draws_matrix(subset_draws(dr, variable="rr"));  rr  <- rr[, order(as.integer(sub("rr\\[(\\d+)\\]","\\1",colnames(rr))))]

## ---- Figure 2a: total i_t, primary line + p_r LHS band ----
i_t_q <- apply(i_t, 2, quantile, probs = c(0.025, 0.5, 0.975))
df_it <- tibble(t = 1:T, lower = i_t_q[1,], median = i_t_q[2,], upper = i_t_q[3,])

lhs_it <- read.csv("results/summary_fixedp_lhs_full_it.csv") %>% filter(grepl("^i_t\\[", variable)) %>%
  mutate(t = as.integer(sub("i_t\\[(\\d+)\\]", "\\1", variable)))
band_it <- lhs_it %>% group_by(t) %>% summarise(band_lo = min(median), band_hi = max(median))
df_it <- df_it %>% left_join(band_it, by = "t")

p3 <- ggplot(df_it, aes(x = t)) +
  geom_ribbon(aes(ymin = band_lo, ymax = band_hi), fill = "#D55E00", alpha = 0.15) +
  geom_ribbon(aes(ymin = lower, ymax = upper), fill = "#0072B2", alpha = 0.2) +
  geom_line(aes(y = median), colour = "#0072B2", linewidth = 1) +
  scale_x_continuous(breaks = year_breaks, labels = year_labels) +
  labs(x = "Year", y = "Quarterly number of \nnew diagnosis") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
        axis.ticks = element_line(colour = "black")) +
  labs(tag = "a") + theme(plot.tag = element_text(face = "bold", size = 14), plot.tag.position = c(0.001, 0.99))

## ---- Figure 2b: i_c by age category, primary model ----
sum_ia_c <- sapply(1:C, function(c) rowSums(i_a[, lower_bound[c]:upper_bound[c], drop = FALSE]))  # nd x C
i_c_q <- array(NA_real_, dim = c(3, C, T))
for (t in 1:T) i_c_q[, , t] <- apply(sweep(sum_ia_c, 1, i_t[, t], "*"), 2, quantile, probs = c(0.025, 0.5, 0.975))
df_i_c <- expand.grid(c = 1:C, t = 1:T) %>%
  mutate(lower = i_c_q[1, , ][cbind(c, t)], median = i_c_q[2, , ][cbind(c, t)], upper = i_c_q[3, , ][cbind(c, t)],
        category = factor(paste0("Cat_", c), levels = names(category_labels)))

p4 <- ggplot(df_i_c, aes(x = t, group = category)) +
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = category), alpha = 0.2, colour = NA) +
  geom_line(aes(y = median, colour = category), linewidth = 1) +
  scale_x_continuous(breaks = year_breaks, labels = year_labels) +
  scale_colour_manual(values = cb_palette_6, breaks = names(category_labels), labels = category_labels) +
  scale_fill_manual(values = cb_palette_6, breaks = names(category_labels), labels = category_labels) +
  labs(x = "Year", y = "Quarterly number of \nnew diagnosis", colour = "Age group", fill = "Age group") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
        axis.ticks = element_line(colour = "black")) +
  labs(tag = "b") + theme(plot.tag = element_text(face = "bold", size = 14), plot.tag.position = c(0.001, 0.99))

figure2 <- p3 / p4
ggsave("figures/figure2.tiff", figure2, device = ragg::agg_tiff, width = 107, height = 140, units = "mm", dpi = 300, compression = "lzw")

## ---- Figure 3b: ip_c (gestational syphilis) by age category, primary model ----
ip_c_q <- array(NA_real_, dim = c(3, C, T))
for (t in 1:T) {
  wsum <- sapply(1:C, function(c) rowSums(sweep(i_a[, lower_bound[c]:upper_bound[c], drop = FALSE], 2, q_base[lower_bound[c]:upper_bound[c], t], "*")))
  ipc <- sweep(wsum, 1, i_t[, t], "*") * rr
  ip_c_q[, , t] <- apply(ipc, 2, quantile, probs = c(0.025, 0.5, 0.975))
}
df_ip_c <- expand.grid(c = 1:C, t = 1:T) %>%
  mutate(lower = ip_c_q[1, , ][cbind(c, t)], median = ip_c_q[2, , ][cbind(c, t)], upper = ip_c_q[3, , ][cbind(c, t)],
        category = factor(paste0("Cat_", c), levels = names(category_labels)))

p5 <- ggplot(df_ip_c, aes(x = t, group = category)) +
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = category), alpha = 0.2, colour = NA) +
  geom_line(aes(y = median, colour = category), linewidth = 1) +
  scale_x_continuous(breaks = year_breaks, labels = year_labels) +
  scale_colour_manual(values = cb_palette_6, breaks = names(category_labels), labels = category_labels) +
  scale_fill_manual(values = cb_palette_6, breaks = names(category_labels), labels = category_labels) +
  labs(x = "Year", y = "Quarterly number of \ngestational syphilis", colour = "Age group", fill = "Age group") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
        axis.ticks = element_line(colour = "black")) +
  labs(tag = "b") + theme(plot.tag = element_text(face = "bold", size = 14), plot.tag.position = c(0.001, 0.99))

## ---- Figure 3a: gst_pop per 100,000 births, primary line + p_r LHS band ----
gst_draws <- as_draws_matrix(subset_draws(dr, variable = "gst_pop"))
gst_draws <- gst_draws[, order(as.integer(sub("gst_pop\\[(\\d+)\\]", "\\1", colnames(gst_draws))))]
gst_q <- apply(gst_draws, 2, quantile, probs = c(0.025, 0.5, 0.975))
df_gst <- tibble(year = years, lower = gst_q[1,], median = gst_q[2,], upper = gst_q[3,])

lhs_gst <- read.csv("results/summary_fixedp_lhs.csv") %>% filter(grepl("^gst_pop\\[", variable)) %>%
  mutate(idx = as.integer(sub("gst_pop\\[(\\d+)\\]", "\\1", variable)), year = 2015 + idx)
band_gst <- lhs_gst %>% group_by(year) %>% summarise(band_lo = min(median), band_hi = max(median))
df_gst <- df_gst %>% left_join(band_gst, by = "year")

p6 <- ggplot(df_gst, aes(x = year)) +
  geom_ribbon(aes(ymin = band_lo, ymax = band_hi), fill = "#D55E00", alpha = 0.15) +
  geom_ribbon(aes(ymin = lower, ymax = upper), fill = "#0072B2", alpha = 0.2) +
  geom_line(aes(y = median), colour = "#0072B2", linewidth = 1) +
  scale_x_continuous(breaks = 2016:2025, labels = 2016:2025) +
  labs(x = "Year", y = "Yearly number of\ngestational syphilis\nper 100,000 live births") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
        axis.ticks = element_line(colour = "black")) +
  labs(tag = "a") + theme(plot.tag = element_text(face = "bold", size = 14), plot.tag.position = c(0.001, 0.99))

figure3 <- p6 / p5
ggsave("figures/figure3.tiff", figure3, device = ragg::agg_tiff, width = 107, height = 140, units = "mm", dpi = 300, compression = "lzw")

cat("Saved figures/figure2.tiff, figures/figure3.tiff
")
