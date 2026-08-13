# data-raw/install_reference.R
#
# Replace the placeholder bundled markers with the validated reference.
#
# The originals were canonical genes typed from memory, never checked
# against data. These are derived from control human substantia nigra
# (Kamath et al. 2022, SCP1768/GSE178265, 8 control donors), ranked by
# exclusivity against the strongest competing lineage, and validated by
# applying an independently derived marker set (Agarwal et al. 2020,
# GSE140231, 5 control donors) to these cells -- 0.93 balanced accuracy
# across the six shared classes.
#
# Run from the Ipse project root, after build_reference.R.

library(Matrix)
devtools::load_all()

out_dir <- "data-raw/derived"
cand_broad <- readRDS(file.path(out_dir, "candidates_kamath_broad.rds"))
cand_da    <- readRDS(file.path(out_dir, "candidates_kamath_da_subtypes.rds"))

kamath_source <- paste(
  "Control donors, human substantia nigra pars compacta;",
  "Kamath et al. 2022, Broad SCP1768 / GEO GSE178265"
)

# ---------------------------------------------------------------------
# Layer assignment -- the one judgement the derivation cannot make
# ---------------------------------------------------------------------
# Effector genes describe what a cell is doing, not what lineage it
# belongs to. A dopaminergic neuron that has lost TH is still a
# dopaminergic neuron, and putting TH in the identity layer would make it
# invisible. Identity is left to lineage-determining transcription
# factors and stable surface markers.

effector_genes <- c(
  "TH", "SLC6A3", "SLC18A2", "DDC", "DRD2",     # dopaminergic programme
  "MBP", "CNP", "MAG", "ERMN", "CARNS1",        # myelination
  "GAD2", "SLC17A6", "SLC17A7"                  # neurotransmitter identity
)

# Housekeeping, ribosomal and MHC genes track library depth or activation
# state rather than lineage.
junk <- c("RPS27", "TMSB10", "SELM", "JUND", "HLA-A", "ITM2C", "VEGFB",
          "EEF1B2", "RPL26", "RPL31", "RPL7", "RPL4", "UBA52", "RPS25",
          "RPS3A", "RPL27", "RPL9", "RPL10A", "RPL23A")

midbrain_markers <- as_marker_table(
  cand_broad, n = 16, min_gap = 0.10,
  exclude = junk, effector = effector_genes, source = kamath_source
)

# Keep the reactive-astrocyte section from the previous reference. It is
# canonical rather than cohort-derived, and is opt-in, so it never
# influences a label.
midbrain_markers <- rbind(
  midbrain_markers,
  data.frame(
    cell_type = "Astrocyte", gene = c("GFAP", "VIM"),
    layer = "state", context = "reactive",
    source = "Canonical reactive-astrocyte markers; no single cohort",
    weight = 1, stringsAsFactors = FALSE
  )
)
rownames(midbrain_markers) <- NULL

cat("=== bundled reference ===\n")
print(table(midbrain_markers$cell_type, midbrain_markers$layer))

# ---------------------------------------------------------------------
# DA subtypes -- a separate reference, for cells already called DA
# ---------------------------------------------------------------------
# Derived within DA cells only, so each subtype competes against the
# others rather than against glia. Gaps here are smaller than in the
# lineage reference by nature: subtypes are genuinely similar, and a gap
# of 0.2 between two DA populations means more than 0.2 between a neuron
# and an oligodendrocyte.
#
# CALB1_RBP4 is excluded: its top candidates were housekeeping and
# ribosomal genes, which is what a subtype looks like when it has no
# distinguishing expression of its own.

midbrain_da_subtypes <- as_marker_table(
  cand_da[cand_da$cell_type != "CALB1_RBP4", ],
  n = 12, min_gap = 0.05, exclude = junk,
  source = paste(kamath_source, "-- DA subtypes, derived within DA cells only")
)

cat("\n=== DA subtype reference ===\n")
print(table(midbrain_da_subtypes$cell_type))

# ---------------------------------------------------------------------
# Write into the package
# ---------------------------------------------------------------------

save(midbrain_markers, file = "data/midbrain_markers.rda", compress = "bzip2")
save(midbrain_da_subtypes, file = "data/midbrain_da_subtypes.rda", compress = "bzip2")

message("bundled references written to data/")
message("next: document midbrain_da_subtypes in R/data.R, then devtools::document()")

# ---------------------------------------------------------------------
# Sanity check: the bundled reference should annotate its own cohort well
# ---------------------------------------------------------------------

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

set.seed(42)
idx <- which(broad != "" & dat$cell_type != "Olig_PLXDC2")
idx <- unlist(lapply(split(idx, broad[idx]),
                     function(i) if (length(i) <= 2000) i else sample(i, 2000)))

cat("\n=== bundled reference, balanced sample of its own cohort ===\n")
print(benchmark_methods(dat$expr[, idx], broad[idx], markers = midbrain_markers))
