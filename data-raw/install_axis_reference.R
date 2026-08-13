# data-raw/install_axis_reference.R
#
# Add the SOX6/CALB1 axis as a third bundled reference.
#
# Rationale: the nine-way subtype call runs at ~63% accuracy, the
# two-class axis at ~85%. For questions about the major dopaminergic
# division -- including the SOX6/AGTR1 population depleted in Parkinson's
# -- the coarser call is the more honest instrument. Simulation says an
# 85% classifier costs almost nothing in power against a perfect one
# (0.78 vs 0.82 for a 10-point shift): donor count is the constraint.
#
# Run after install_reference.R.

library(Matrix)
devtools::load_all()

out_dir <- "data-raw/derived"
dat <- readRDS(file.path(out_dir, "kamath_controls_expr.rds"))

da   <- grepl("^(SOX6|CALB1)_", dat$cell_type)
axis <- ifelse(grepl("^SOX6", dat$cell_type[da]), "SOX6", "CALB1")

cand_axis <- derive_markers(dat$expr[, da], axis, donor = dat$donor[da])

midbrain_da_axis <- as_marker_table(
  cand_axis, n = 12, min_gap = 0.10,
  source = paste("Control donors, human substantia nigra pars compacta;",
                 "Kamath et al. 2022, Broad SCP1768 / GEO GSE178265;",
                 "SOX6 vs CALB1 axis derived within dopaminergic cells")
)

print(table(midbrain_da_axis$cell_type))
cat("\ntop markers per class:\n")
for (c in unique(midbrain_da_axis$cell_type)) {
  cat(" ", c, ":", paste(head(midbrain_da_axis$gene[midbrain_da_axis$cell_type == c], 6),
                         collapse = ", "), "\n")
}

save(midbrain_da_axis, file = "data/midbrain_da_axis.rda", compress = "bzip2")
saveRDS(cand_axis, file.path(out_dir, "candidates_kamath_da_axis.rds"))

# ---------------------------------------------------------------------
# Record the operating characteristics -- these are what a user needs in
# order to interpret a measured proportion, and they belong with the data
# rather than in a chat log.
# ---------------------------------------------------------------------

pred <- annotate_cells(dat$expr[, da], markers = midbrain_da_axis, seed = 42)$cell_type
cm <- table(true = axis, predicted = pred)
print(cm)

sens <- cm["SOX6", "SOX6"]  / sum(cm["SOX6", ])
fpr  <- cm["CALB1", "SOX6"] / sum(cm["CALB1", ])

cat(sprintf("\nsensitivity (SOX6) %.3f | false positive rate %.3f\n", sens, fpr))
cat(sprintf("attenuation factor %.3f -- a true difference measures ~%.0f%% of its size\n",
            sens - fpr, 100 * (sens - fpr)))
cat(sprintf("fixed point %.3f -- proportions above this are understated, below overstated\n",
            fpr / (1 - sens + fpr)))

cat("\nThese hold only if accuracy is equal across the groups being\n")
cat("compared. If a condition shifts cells toward weakly-marked states,\n")
cat("the error stops cancelling and becomes bias.\n")

message("\nnext: devtools::document(), then test() and check()")

# ---------------------------------------------------------------------
# Usage
# ---------------------------------------------------------------------
#   labels <- annotate_cells(expr)
#   labels <- annotate_subtypes(expr, labels,
#                               parent  = "Dopaminergic neuron",
#                               markers = midbrain_da_axis)
#   table(labels$subtype)
