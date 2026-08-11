# data-raw/prepare_markers.R
#
# Builds the `midbrain_markers` reference table that ships inside the
# package (data/midbrain_markers.rda).
#
# Run this script, then:
#   usethis::use_data(midbrain_markers, overwrite = TRUE)
# whenever the marker list changes.
#
# ---------------------------------------------------------------------
# LAYERS
#
# Markers are split into layers so that cell *identity* is scored
# separately from cell *state*:
#
#   identity  Lineage-determining genes (mostly transcription factors).
#             These, and only these, decide the assigned label. Kept
#             separate so that a label never encodes disease status.
#
#   effector  Functional-programme genes for that lineage. Reflect what
#             the cell is currently doing, and can be lost in disease
#             without the cell ceasing to be that lineage. Scored and
#             reported, but never used to assign a label.
#
#   state     Condition-associated signatures (reactive, degenerative,
#             disease-associated). These form opt-in *sections*: they are
#             never scored unless the user names them. The `context`
#             column names the section; the `source` column records where
#             the signature came from, so a user can judge for themselves
#             whether a given analysis is independent of it.
#
# The layer split is what lets the package express "retains dopaminergic
# lineage identity but has lost the effector programme" rather than
# forcing a single label onto a cell that is mid-transition.
# ---------------------------------------------------------------------

m <- function(cell_type, gene, layer, context = "canonical",
              source = NA_character_, weight = 1) {
  data.frame(cell_type = cell_type, gene = gene, layer = layer,
             context = context, source = source, weight = weight,
             stringsAsFactors = FALSE)
}

midbrain_markers <- rbind(
  # --- Dopaminergic neuron ---------------------------------------------
  m("Dopaminergic neuron", c("NR4A2", "PITX3", "LMX1B"), "identity"),
  m("Dopaminergic neuron", c("TH", "SLC6A3", "SLC18A2", "DDC"), "effector"),

  # --- Astrocyte --------------------------------------------------------
  m("Astrocyte", c("AQP4", "SLC1A3", "ALDH1L1"), "identity"),
  # GFAP is strongly reactive-state dependent, so it is deliberately NOT
  # an identity marker: using it as one conflates reactive astrocytes with
  # astrocyte identity.
  m("Astrocyte", c("GFAP", "VIM"), "state", context = "reactive",
    source = "Canonical reactive-astrocyte markers; no single cohort"),

  # --- Microglia --------------------------------------------------------
  m("Microglia", c("P2RY12", "CX3CR1", "CSF1R"), "identity"),
  m("Microglia", c("C1QA", "C1QB"), "effector"),

  # --- Oligodendrocyte --------------------------------------------------
  m("Oligodendrocyte", c("ST18", "PLP1", "MOG"), "identity"),
  m("Oligodendrocyte", c("MBP", "MAG"), "effector"),

  # --- OPC --------------------------------------------------------------
  m("OPC", c("PDGFRA", "CSPG4", "OLIG1"), "identity"),

  # --- Endothelial ------------------------------------------------------
  m("Endothelial", c("CLDN5", "FLT1"), "identity"),
  m("Endothelial", c("PECAM1"), "effector")
)

rownames(midbrain_markers) <- NULL

# --- TODO --------------------------------------------------------------
# 1. This remains a starter set of canonical markers, not curated against
#    a reference dataset. Derive replacements from control samples of a
#    neurotypical-donor atlas, and cross-check candidates across at least
#    two independent cohorts before trusting them.
#
# 2. No `state` markers with a disease `context` are bundled. Add them in
#    the form below, but note the circularity risk: a state signature
#    derived from the same data you then apply it to will trivially
#    "detect" itself. Derive from independent data.
#
#      m("Astrocyte", c("GENE1", "GENE2"), "state", context = "SCZ",
#        source = "Author et al. 2024, GSE123456, control-vs-SCZ DE"),
#
#    Always fill in `source`. A named, reference-anchored group is the
#    point of a section -- it is what makes a result comparable across
#    datasets rather than an opaque cluster index.
#
#    Use cell_type = "Any" for a signature that applies across cell types
#    rather than to one lineage.
#
# 3. Dopaminergic subtypes (SNpc vs VTA; the SOX6 / CALB1 axis) are not
#    yet represented. That is the resolution most likely to be worth
#    having, and needs subtype-level markers from a subtype-resolved
#    reference.
# -------------------------------------------------------------------------

usethis::use_data(midbrain_markers, overwrite = TRUE)
