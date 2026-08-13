# data-raw/subtype_power.R
#
# Can an ~85%-accurate SOX6/CALB1 classifier detect a shift in subtype
# proportions between groups?
#
# The operating characteristics are recomputed from your own annotation
# rather than hardcoded, so this is an independent check rather than a
# rerun of numbers read off a screen.

library(Matrix)
devtools::load_all()

dat <- readRDS("data-raw/derived/kamath_controls_expr.rds")

# ---- recompute the confusion from scratch ---------------------------
da    <- grepl("^(SOX6|CALB1)_", dat$cell_type)
axis  <- ifelse(grepl("^SOX6", dat$cell_type[da]), "SOX6", "CALB1")

cand_axis <- derive_markers(dat$expr[, da], axis, donor = dat$donor[da])
mk_axis   <- as_marker_table(cand_axis, n = 12, min_gap = 0.1,
                             source = "SOX6/CALB1 axis, Kamath control donors")

pred_axis <- annotate_cells(dat$expr[, da], markers = mk_axis, seed = 42)$cell_type

cm <- table(true = axis, predicted = pred_axis)
print(cm)

sens <- cm["SOX6", "SOX6"]   / sum(cm["SOX6", ])    # correctly call SOX6
fpr  <- cm["CALB1", "SOX6"]  / sum(cm["CALB1", ])   # wrongly call SOX6
p_ctrl <- mean(axis == "SOX6")

cat(sprintf("\nsensitivity %.3f | false positive rate %.3f | attenuation %.3f\n",
            sens, fpr, sens - fpr))
cat(sprintf("baseline SOX6 share in these controls: %.3f\n\n", p_ctrl))

# A measured proportion is a distorted version of the true one.
observed <- function(p) p * sens + (1 - p) * fpr

for (p in c(0.1, 0.3, 0.6, 0.9)) {
  cat(sprintf("true share %.2f measures as %.2f\n", p, observed(p)))
}

# ---- power ----------------------------------------------------------
# Donors are the unit of analysis, not cells: cells within a donor are
# not independent, and treating them as such inflates significance
# dramatically. `donor_sd` is between-donor variability in the true
# proportion -- set it from your own data if you can.

simulate <- function(eff, k = 8, n = 500, donor_sd = 0.06,
                     reps = 2000, perfect = FALSE, seed = 1) {
  set.seed(seed)
  f <- if (perfect) identity else observed
  mean(replicate(reps, {
    p_a <- pmin(pmax(stats::rnorm(k, p_ctrl, donor_sd), 0.01), 0.99)
    p_b <- pmin(pmax(stats::rnorm(k, p_ctrl - eff, donor_sd), 0.01), 0.99)
    a <- stats::rbinom(k, n, f(p_a)) / n
    b <- stats::rbinom(k, n, f(p_b)) / n
    stats::t.test(a, b)$p.value < 0.05
  }))
}

cat("\npower, 8 donors per group, 500 DA cells each\n")
cat(sprintf("%-12s %-12s %-10s %-10s\n", "true shift", "measured", "power", "if perfect"))
for (eff in c(0.02, 0.05, 0.10, 0.15, 0.20)) {
  cat(sprintf("%8.0f pp %9.1f pp %10.2f %10.2f\n",
              100 * eff,
              100 * (observed(p_ctrl) - observed(p_ctrl - eff)),
              simulate(eff), simulate(eff, perfect = TRUE)))
}

cat(sprintf("\nfalse positives when there is no true difference: %.3f (nominal 0.05)\n",
            simulate(0)))

# ---- how many donors would you actually need? -----------------------
cat("\ndonors per group needed for 0.80 power\n")
for (eff in c(0.05, 0.10, 0.15)) {
  k <- 3
  while (k < 60 && simulate(eff, k = k) < 0.80) k <- k + 1
  cat(sprintf("  %2.0f pp shift: %s donors\n", 100 * eff,
              if (k >= 60) "60+" else as.character(k)))
}

# ---- the assumption this all rests on -------------------------------
# The bias cancels between groups only if it is the SAME in both. If a
# condition pushes cells toward weakly-marked subtypes, sensitivity drops
# in that group alone and the comparison becomes biased rather than just
# attenuated. That cannot be tested here -- it needs labelled cells from
# the condition of interest.
cat("\nNOTE: this assumes classifier accuracy is equal in both groups.\n")
cat("If disease shifts cells toward weakly-marked subtypes, it is not.\n")
