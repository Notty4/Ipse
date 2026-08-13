# data-raw/kamath_two_level.R
#
# The first Kamath run derived at SUBTYPE resolution across all lineages
# at once. That is the wrong granularity for lineage markers: every Astro
# subtype's strongest competitor is another Astro subtype, so exclusivity
# collapses and you get subtype-discriminating genes instead of markers
# that say "this is an astrocyte".
#
# So derive twice:
#   LEVEL 1  broad lineage classes -> the identity layer of the reference
#   LEVEL 2  DA subtypes, within DA cells only -> the resolution the
#            thesis question needs, and what Agarwal could not provide
#
# Run after kamath_pipeline.R has cached the expression subset.

library(Matrix)
devtools::load_all()

out_dir <- "data-raw/derived"
dat <- readRDS(file.path(out_dir, "kamath_controls_expr.rds"))

# =====================================================================
# Clean up the label set
# =====================================================================

# Clusters to drop, on evidence from the derivation run:
#
#   Olig_PLXDC2   every top marker has donor_consistency of exactly 0.000
#                 -- nothing about it reproduces across donors. Its earlier
#                 neuronal contamination (NEFL, SNAP25) disappeared once
#                 the other oligodendrocyte subtypes were present to
#                 compete, but the donor irreproducibility did not.
#
# MG_TSPO_VIM was dropped in an earlier version and is now KEPT: with the
# full lineage set present its markers include AIF1 (consistency 1.0),
# TSPO, PYCARD and FCGR1A -- genuine microglial genes, not just ribosomal
# ones. At broad level it merges into Microglia regardless.
#
# `da` and `nonda` no longer exist now that the coarse nurr files are
# excluded from label assembly, but are left here defensively.
drop_clusters <- c("Olig_PLXDC2", "da", "nonda")

usable <- !(dat$cell_type %in% drop_clusters)
message(sprintf("dropping %d cells in %d suspect clusters",
                sum(!usable), length(drop_clusters)))

expr    <- dat$expr[, usable]
subtype <- dat$cell_type[usable]
donor   <- dat$donor[usable]

# =====================================================================
# LEVEL 1 -- broad lineage classes
# =====================================================================

broad <- character(length(subtype))
broad[grepl("^Astro_",   subtype)] <- "Astrocyte"
broad[grepl("^MG_",      subtype)] <- "Microglia"
broad[grepl("^Macro_",   subtype)] <- "Macrophage"
broad[grepl("^Olig_",    subtype)] <- "Oligodendrocyte"
broad[grepl("^OPC",      subtype)] <- "OPC"
broad[grepl("^Endo",     subtype)] <- "Endothelial"
broad[grepl("^Ependyma", subtype)] <- "Ependymal"
broad[grepl("^Ex_",      subtype)] <- "Excitatory neuron"
broad[grepl("^Inh_",     subtype)] <- "Inhibitory neuron"
broad[grepl("^(SOX6|CALB1)_", subtype)] <- "Dopaminergic neuron"

unmapped <- unique(subtype[broad == ""])
if (length(unmapped) > 0) {
  warning("unmapped subtypes (excluded): ", paste(unmapped, collapse = ", "),
          call. = FALSE)
}
keep_broad <- broad != ""
print(table(broad[keep_broad]))

cand_broad <- derive_markers(
  expr[, keep_broad], broad[keep_broad], donor = donor[keep_broad]
)

cat("\n=== LEVEL 1: broad lineage ceilings ===\n")
print(round(tapply(cand_broad$detect_gap, cand_broad$cell_type, max), 3))

for (t in sort(unique(cand_broad$cell_type))) {
  cat("\n---", t, "---\n")
  print(utils::head(cand_broad[cand_broad$cell_type == t,
                               c("gene", "detect_gap", "donor_consistency",
                                 "n_donors", "competitor")], 10))
}
saveRDS(cand_broad, file.path(out_dir, "candidates_kamath_broad.rds"))

# =====================================================================
# LEVEL 2 -- DA subtypes, within DA cells only
# =====================================================================
# Restricting to DA cells is the point: each subtype's competitor is now
# another DA subtype, which is exactly the discrimination wanted. Run at
# this level, a small detect_gap is informative rather than a failure --
# subtypes are genuinely similar.

is_da <- broad == "Dopaminergic neuron"
message(sprintf("\n%d DA cells across %d subtypes", sum(is_da),
                length(unique(subtype[is_da]))))
print(table(subtype[is_da]))

cand_da <- derive_markers(
  expr[, is_da], subtype[is_da], donor = donor[is_da]
)

cat("\n=== LEVEL 2: DA subtype ceilings ===\n")
print(round(tapply(cand_da$detect_gap, cand_da$cell_type, max), 3))

for (t in sort(unique(cand_da$cell_type))) {
  cat("\n---", t, "---\n")
  print(utils::head(cand_da[cand_da$cell_type == t,
                            c("gene", "detect_gap", "donor_consistency",
                              "n_donors", "competitor")], 10))
}
saveRDS(cand_da, file.path(out_dir, "candidates_kamath_da_subtypes.rds"))

# =====================================================================
# NEXT
# =====================================================================
# 1. Build the broad reference and test Agarwal's markers against it --
#    the first non-circular evaluation:
#
#      mk_agarwal <- readRDS(file.path(out_dir, "markers_agarwal.rds"))
#      benchmark_methods(expr[, keep_broad], broad[keep_broad],
#                        markers = mk_agarwal)
#
#    Agarwal has no Macrophage or Ependymal, so those cells have no
#    matching truth; exclude them or expect them to count as errors.
#
# 2. Cross-cohort consensus for the shared classes:
#
#      cons <- consensus_markers(
#        list(agarwal = readRDS(file.path(out_dir, "candidates_GSE140231.rds")),
#             kamath  = cand_broad),
#        top_n = 50, min_cohorts = 2)
#
#    Agarwal's labels must match Kamath's broad names exactly, or nothing
#    will be shared -- rename before combining.
