# 07_main_tables.R
#
# Table 1: parameter estimates of the primary model (posterior median and 95% credible
#          interval, results/revision_summary_main.csv) and the range of posterior medians
#          across the 10-point CS-risk Latin hypercube sensitivity analysis
#          (results/summary_fixedp_lhs.csv).
# Table 2: cumulative CS cases, 2016-2025, when a proportion (0%, 10%, ..., 100%) of women
#          with syphilis is given the fertility rate of all women in the same age group
#          (cum_cs_w, primary model; median and 95% credible interval).
# Also writes the CS-risk sensitivity ranges for every summarised quantity, including the
# diagnosed-case incidence in 2016 Q1 and 2025 Q3 quoted in Results.
#
# Outputs: tables/table1.csv, tables/table2.csv, tables/cs_risk_lhs_ranges.csv

library(dplyr)
dir.create("tables", showWarnings = FALSE)

main <- read.csv("results/revision_summary_main.csv")
lhs  <- read.csv("results/summary_fixedp_lhs.csv")

lhs_range <- lhs %>%
  group_by(variable) %>%
  summarise(lhs_median_min = min(median), lhs_median_max = max(median), .groups = "drop")
write.csv(lhs_range, "tables/cs_risk_lhs_ranges.csv", row.names = FALSE)

## ---- Table 1 ----
t1_rows <- data.frame(
  variable = c(sprintf("rr[%d]", 1:6), "rep",
               sprintf("p1[%d]", 1:3), sprintf("p2[%d]", 1:3), sprintf("p3[%d]", 1:3)),
  parameter = c(paste("Fertility-rate ratio,", c("15-19", "20-24", "25-29", "30-34", "35-39", "40-44")),
                "Reporting fraction",
                paste("Waiting interval from Dx to CS delivery, 2016-2019,", c("<1Q", ">1Q <2Q", ">2Q <3Q")),
                paste("Waiting interval from Dx to CS delivery, 2020-2021,", c("<1Q", ">1Q <2Q", ">2Q <3Q")),
                paste("Waiting interval from Dx to CS delivery, 2022-2025,", c("<1Q", ">1Q <2Q", ">2Q <3Q"))),
  stringsAsFactors = FALSE
)
table1 <- t1_rows %>%
  left_join(main %>% select(variable, median = q50, lower = q2.5, upper = q97.5), by = "variable") %>%
  left_join(lhs_range, by = "variable") %>%
  mutate(across(c(median, lower, upper, lhs_median_min, lhs_median_max), ~ signif(.x, 3)))
print(table1, row.names = FALSE)
write.csv(table1, "tables/table1.csv", row.names = FALSE)

## ---- Table 2 ----
table2 <- main %>%
  filter(grepl("^cum_cs_w\\[", variable)) %>%
  mutate(w = as.integer(sub("cum_cs_w\\[(\\d+)\\]", "\\1", variable))) %>%
  arrange(w) %>%
  transmute(proportion_pct = (w - 1) * 10,
            cumulative_cs_median = q50, lower = q2.5, upper = q97.5)
print(table2, row.names = FALSE)
write.csv(table2, "tables/table2.csv", row.names = FALSE)
