# 05_posterior_correlation.R
#
# Supplementary Table ST4: full pairwise posterior correlation matrix among 19 parameters
# of the primary model (reporting fraction, j(t) at 2016 Q1 and 2025 Q3, s(a) for the six
# age categories, random-walk SD, and the three delay-distribution simplexes), computed
# from the primary fit's saved draws (results/revision_draws_main.rds, from
# 01_fit_primary_and_sensitivity.R).
#
# Outputs:
#   results/posterior_correlation_matrix_main.csv   19 x 19 matrix (3 decimals)
#   results/posterior_correlation_pairs_main.csv    all 171 pairs, sorted by |r|
#   tables/supplementary_table_ST4_posterior_correlations.md

library(posterior)
dir.create("tables", showWarnings = FALSE)

d <- readRDS("results/revision_draws_main.rds")
vn <- variables(d)

target_vars <- c(
  "rep",
  "i_t[1]", "i_t[39]",
  "rr[1]", "rr[2]", "rr[3]", "rr[4]", "rr[5]", "rr[6]",
  "s_i",
  "p1[1]", "p1[2]", "p1[3]",
  "p2[1]", "p2[2]", "p2[3]",
  "p3[1]", "p3[2]", "p3[3]"
)

missing_vars <- setdiff(target_vars, vn)
if (length(missing_vars) > 0) {
  stop("Missing variables in draws object: ", paste(missing_vars, collapse = ", "))
}

draws_df <- as_draws_df(subset_draws(d, variable = target_vars))
mat <- as.matrix(draws_df[, target_vars])

cor_mat <- cor(mat)          # dimnames = target_vars (Stan names), used for indexing below
cor_mat_raw <- cor_mat

# Nice display labels matching the manuscript's notation
labels <- c(
  "rep" = "rho (reporting fraction)",
  "i_t[1]" = "j(t) 2016 Q1",
  "i_t[39]" = "j(t) 2025 Q3",
  "rr[1]" = "s(a) 15-19y",
  "rr[2]" = "s(a) 20-24y",
  "rr[3]" = "s(a) 25-29y",
  "rr[4]" = "s(a) 30-34y",
  "rr[5]" = "s(a) 35-39y",
  "rr[6]" = "s(a) 40-44y",
  "s_i" = "s_i (random-walk SD)",
  "p1[1]" = "f(u) 2016-2019, u=1",
  "p1[2]" = "f(u) 2016-2019, u=2",
  "p1[3]" = "f(u) 2016-2019, u=3",
  "p2[1]" = "f(u) 2020-2021, u=1",
  "p2[2]" = "f(u) 2020-2021, u=2",
  "p2[3]" = "f(u) 2020-2021, u=3",
  "p3[1]" = "f(u) 2022-2025, u=1",
  "p3[2]" = "f(u) 2022-2025, u=2",
  "p3[3]" = "f(u) 2022-2025, u=3"
)

rownames(cor_mat) <- labels[target_vars]
colnames(cor_mat) <- labels[target_vars]

out_mat <- round(cor_mat, 3)
write.csv(out_mat, "results/posterior_correlation_matrix_main.csv", row.names = TRUE)

# Long format (one row per unique pair, excluding diagonal), sorted by |r| descending
pairs <- t(combn(target_vars, 2))
pair_df <- data.frame(
  var1 = labels[pairs[, 1]],
  var2 = labels[pairs[, 2]],
  r = mapply(function(a, b) cor_mat_raw[a, b], pairs[, 1], pairs[, 2])
)
pair_df$r <- round(pair_df$r, 3)
pair_df <- pair_df[order(-abs(pair_df$r)), ]
write.csv(pair_df, "results/posterior_correlation_pairs_main.csv", row.names = FALSE)

cat("Wrote results/posterior_correlation_matrix_main.csv (", nrow(out_mat), "x", ncol(out_mat), ")\n")
cat("Wrote results/posterior_correlation_pairs_main.csv (", nrow(pair_df), "pairs)\n")
cat("\nTop 15 pairs by |r|:\n")
print(head(pair_df, 15))
cat("\nNumber of pairs with |r| >= 0.2:", sum(abs(pair_df$r) >= 0.2), "out of", nrow(pair_df), "\n")

## ---- Supplementary Table ST4 (markdown, lower triangle, 2 decimals) ----
mat <- out_mat   # the 3-decimal matrix written to the CSV above
short <- c(
  "rho", "j.1", "j.39",
  "s.1", "s.2", "s.3", "s.4", "s.5", "s.6",
  "s_i",
  "f1.1", "f1.2", "f1.3",
  "f2.1", "f2.2", "f2.3",
  "f3.1", "f3.2", "f3.3"
)
rownames(mat) <- short
colnames(mat) <- short

# Lower triangle only (incl. diagonal), 2 decimals, blank upper triangle
m2 <- round(mat, 2)
disp <- matrix("", nrow(m2), ncol(m2), dimnames = dimnames(m2))
for (i in seq_len(nrow(m2))) {
  for (j in seq_len(i)) {
    disp[i, j] <- sprintf("%.2f", m2[i, j])
  }
}

lines <- c(
  "# Supplementary Table ST4. Full posterior correlation matrix (primary model)",
  "",
  "Pairwise posterior correlations among the 19 parameters used throughout the sensitivity",
  "analyses (Methods 2.4), computed from the primary model's saved draws (4 chains x",
  "(1000 warm-up + 1000 sampling)). Lower triangle shown (matrix is symmetric); diagonal = 1.",
  "The two correlations reported in Results (rho-j(t)) are highlighted in bold.",
  "",
  "**Abbreviations.** rho = reporting fraction; j.1/j.39 = diagnosed-case incidence j(t) at",
  "2016 Q1 / 2025 Q3; s.1-s.6 = fertility-rate ratio s(a) for the six age categories (15-19,",
  "20-24, 25-29, 30-34, 35-39, 40-44 years); s_i = SD of the random-walk prior on log j(t);",
  "f1.1-f1.3 / f2.1-f2.3 / f3.1-f3.3 = the three-point delay-distribution simplex f(u) for",
  "2016-2019 / 2020-2021 / 2022-2025.",
  ""
)

header <- paste0("| | ", paste(short, collapse = " | "), " |")
sep <- paste0("|---|", paste(rep("---", length(short)), collapse = "|"), "|")
lines <- c(lines, header, sep)

for (i in seq_len(nrow(disp))) {
  row_vals <- disp[i, ]
  # bold the two pairs reported in Results: rho-j.1, rho-j.39
  if (short[i] %in% c("j.1", "j.39")) {
    idx_rho <- which(short == "rho")
    if (row_vals[idx_rho] != "") row_vals[idx_rho] <- paste0("**", row_vals[idx_rho], "**")
  }
  lines <- c(lines, paste0("| **", short[i], "** | ", paste(row_vals, collapse = " | "), " |"))
}

writeLines(lines, "tables/supplementary_table_ST4_posterior_correlations.md")
cat("Wrote tables/supplementary_table_ST4_posterior_correlations.md\n")
