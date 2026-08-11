# data-raw/prepare_markers.R
#
# Builds the `midbrain_markers` reference table that ships inside the
# package (data/midbrain_markers.rda), so annotate_midbrain_cells() never
# needs an external atlas download.
#
# Run this script, then:
#   usethis::use_data(midbrain_markers, overwrite = TRUE)
# whenever the marker list changes.

midbrain_markers <- data.frame(
  cell_type = c(
    "Dopaminergic neuron", "Dopaminergic neuron", "Dopaminergic neuron", "Dopaminergic neuron",
    "Astrocyte", "Astrocyte", "Astrocyte",
    "Microglia", "Microglia", "Microglia",
    "Oligodendrocyte", "Oligodendrocyte", "Oligodendrocyte",
    "OPC", "OPC",
    "Endothelial", "Endothelial"
  ),
  gene = c(
    "TH", "SLC6A3", "NR4A2", "SLC18A2",
    "GFAP", "AQP4", "SLC1A3",
    "P2RY12", "C1QA", "CX3CR1",
    "MBP", "PLP1", "MOG",
    "PDGFRA", "CSPG4",
    "CLDN5", "PECAM1"
  ),
  weight = 1,
  stringsAsFactors = FALSE
)

# --- TODO --------------------------------------------------------------
# This is a placeholder starter set (canonical markers, not curated against
# your thesis panel or a differential-expression check). Before relying on
# this for real annotation, replace/expand it with a curated list — e.g.
# cross-checked against the Agarwal midbrain reference (GSE140231) markers
# used in the thesis pipeline, or literature-validated marker sets per
# subtype (SNpc vs VTA dopaminergic neurons, etc).
# -------------------------------------------------------------------------

usethis::use_data(midbrain_markers, overwrite = TRUE)
