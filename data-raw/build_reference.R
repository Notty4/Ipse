# data-raw/build_reference.R
#
# Build the Ipse marker reference from the Kamath broad-lineage
# candidates, then test it -- including the first non-circular
# evaluation available: Agarwal's markers applied to Kamath's cells.
#
# Run after kamath_two_level.R.

library(Matrix)
devtools::load_all()

out_dir <- "data-raw/derived"
cand_broad <- readRDS(file.path(out_dir, "candidates_kamath_broad.rds"))
cand_da    <- readRDS(file.path(out_dir, "candidates_kamath_da_subtypes.rds"))

# =====================================================================
# STEP 1 -- inspect DA candidates before assigning layers
# =====================================================================
# Only LMX1B and EBF3 in the top 10 are lineage-determining; the rest
# (TH, SLC6A3, SLC18A2, DDC, DRD2) are effector genes. Taking the top n
# blindly would leave the identity layer with two markers and put the
# functional programme where it can assign labels -- which is exactly
# what the layer split exists to prevent.

cat("=== DA candidates, top 30 -- pick identity genes from here ===\n")
print(head(cand_broad[cand_broad$cell_type == "Dopaminergic neuron",
                      c("gene", "detect_gap", "donor_consistency")], 30))

# =====================================================================
# STEP 2 -- curation decisions
# =====================================================================

# Effector genes: functional programme, scored but never assign a label.
# A DA neuron that has lost these is still a DA neuron.
effector_genes <- c(
  "TH", "SLC6A3", "SLC18A2", "DDC", "DRD2",       # dopaminergic programme
  "MBP", "CNP", "MAG", "ERMN", "CARNS1",          # myelination
  "GAD2", "SLC17A6", "SLC17A7"                    # neurotransmitter identity
)

# Drop: housekeeping, ribosomal, and MHC genes that track cell state or
# library depth rather than lineage.
junk <- c("RPS27", "TMSB10", "SELM", "JUND", "HLA-A", "ITM2C", "VEGFB",
          "EEF1B2", "RPL26", "RPL31", "RPL7", "RPL4", "UBA52", "RPS25",
          "RPS3A", "RPL27", "RPL9", "RPL10A", "RPL23A")

markers_kamath <- as_marker_table(
  cand_broad,
  n        = 10,
  min_gap  = 0.15,
  exclude  = junk,
  effector = effector_genes,
  source   = "Control donors, human SNpc, Kamath et al. 2022 (SCP1768/GSE178265)"
)

print(table(markers_kamath$cell_type, markers_kamath$layer))

# Every cell type must retain identity markers, or it cannot be labelled.
no_identity <- setdiff(unique(markers_kamath$cell_type),
                       unique(markers_kamath$cell_type[markers_kamath$layer == "identity"]))
if (length(no_identity) > 0) {
  warning("no identity markers for: ", paste(no_identity, collapse = ", "),
          call. = FALSE)
}

saveRDS(markers_kamath, file.path(out_dir, "markers_kamath.rds"))

# =====================================================================
# STEP 3 -- benchmark on the cohort it came from (self-annotation)
# =====================================================================
# This is the ceiling, not a validation: the markers were derived from
# these labels. A poor score here means something is broken; a good one
# proves little.

dat <- readRDS(file.path(out_dir, "kamath_controls_expr.rds"))

broad <- character(length(dat$cell_type))
broad[grepl("^Astro_",   dat$cell_type)] <- "Astrocyte"
broad[grepl("^MG_",      dat$cell_type)] <- "Microglia"
broad[grepl("^Macro_",   dat$cell_type)] <- "Macrophage"
broad[grepl("^Olig_",    dat$cell_type)] <- "Oligodendrocyte"
broad[grepl("^OPC",      dat$cell_type)] <- "OPC"
broad[grepl("^Endo",     dat$cell_type)] <- "Endothelial"
broad[grepl("^Ependyma", dat$cell_type)] <- "Ependymal"
broad[grepl("^Ex_",      dat$cell_type)] <- "Excitatory neuron"
broad[grepl("^Inh_",     dat$cell_type)] <- "Inhibitory neuron"
broad[grepl("^(SOX6|CALB1)_", dat$cell_type)] <- "Dopaminergic neuron"
ok <- broad != "" & !(dat$cell_type %in% c("Olig_PLXDC2"))

cat("\n=== self-annotation on Kamath (ceiling, not validation) ===\n")
print(benchmark_methods(dat$expr[, ok], broad[ok], markers = markers_kamath))

# =====================================================================
# STEP 4 -- CROSS-COHORT: Agarwal's markers on Kamath's cells
# =====================================================================
# The first genuinely non-circular number in this project. Expect it to
# be worse than either self-annotation figure; that gap is the honest
# measure of how far a marker set transfers.

ag_path <- file.path(out_dir, "markers_agarwal.rds")
if (file.exists(ag_path)) {
  mk_ag <- readRDS(ag_path)

  # Agarwal has no Macrophage, Ependymal, Excitatory or Inhibitory class,
  # so cells of those types have no reachable truth label. Restrict to the
  # classes both cohorts define, or the score punishes Agarwal for
  # categories it never had.
  shared <- intersect(unique(mk_ag$cell_type), unique(broad))
  cat("\nclasses shared with Agarwal:", paste(shared, collapse = ", "), "\n")

  cross <- ok & broad %in% shared
  cat("\n=== CROSS-COHORT: Agarwal markers -> Kamath cells ===\n")
  print(benchmark_methods(dat$expr[, cross], broad[cross], markers = mk_ag))

  cat("\n=== same cells, Kamath's own markers (for comparison) ===\n")
  print(benchmark_methods(dat$expr[, cross], broad[cross], markers = markers_kamath))
} else {
  message("markers_agarwal.rds not found -- skipping the cross-cohort test")
}

# =====================================================================
# NEXT
# =====================================================================
# 1. If the cross-cohort number holds up, build the consensus reference:
#      cons <- consensus_markers(
#        list(agarwal = readRDS(file.path(out_dir, "candidates_GSE140231.rds")),
#             kamath  = cand_broad),
#        top_n = 50, min_cohorts = 2)
#    Agarwal's class names must match Kamath's broad names exactly.
#
# 2. Add the DA subtypes as a second reference from cand_da, for use on
#    cells already labelled dopaminergic.
#
# 3. Replace data/midbrain_markers.rda via data-raw/prepare_markers.R.
