#' Midbrain marker gene reference table
#'
#' A layered marker gene reference for major midbrain cell types, bundled
#' with the package so that [annotate_cells()] can run without an
#' external reference download.
#'
#' @section Layers:
#' Markers are separated by what they tell you about a cell:
#' \describe{
#'   \item{`identity`}{Lineage-determining genes, mostly transcription
#'     factors. These alone decide the assigned label.}
#'   \item{`effector`}{Functional-programme genes for that lineage. These
#'     can be lost in disease without the cell ceasing to belong to the
#'     lineage, so they are scored and reported but never used to assign a
#'     label.}
#'   \item{`state`}{Condition-associated signatures. The `context` column
#'     names the condition (e.g. `"reactive"`). Use `cell_type = "Any"` for
#'     a signature that applies across all cell types.}
#' }
#'
#' The split exists so that a label can never encode disease status: if
#' disease-associated genes helped assign labels, downstream comparisons of
#' cell type composition between conditions would be partly circular.
#'
#' @format A data frame with one row per (cell type, marker gene) pair:
#' \describe{
#'   \item{cell_type}{Character. Cell type label, or `"Any"`.}
#'   \item{gene}{Character. Marker gene symbol.}
#'   \item{layer}{Character. `"identity"`, `"effector"`, or `"state"`.}
#'   \item{context}{Character. `"canonical"`, or a named condition for
#'     state markers. This is the name passed to the `sections` argument
#'     of [annotate_cells()].}
#'   \item{source}{Character. Where the signature came from (cohort,
#'     accession, publication), or `NA` for canonical markers. Recorded so
#'     that a user can judge whether a given analysis is independent of
#'     the data a section was derived from.}
#'   \item{weight}{Numeric. Relative weight within its marker set.}
#' }
#'
#' @source Starter set of canonical markers, not curated against a
#'   reference dataset. See `data-raw/prepare_markers.R`, which also notes
#'   the circularity risk when adding disease state signatures.
#'
#' @examples
#' data(midbrain_markers)
#' table(midbrain_markers$layer)
#' subset(midbrain_markers, layer == "identity")
"midbrain_markers"

#' Dopaminergic neuron subtype marker reference
#'
#' Markers distinguishing dopaminergic neuron subtypes from one another,
#' for use on cells already labelled `"Dopaminergic neuron"` by
#' [annotate_cells()] with [midbrain_markers].
#'
#' @section Why this is separate:
#' These markers were derived within dopaminergic cells only, so each
#' subtype's strongest competitor is another subtype rather than a glial
#' lineage. Deriving subtypes alongside lineages collapses their scores:
#' every astrocyte subtype's nearest competitor is another astrocyte, and
#' the resulting markers separate subtypes rather than identifying
#' lineages. Judge the gaps here on a different scale to those in
#' [midbrain_markers] — a detection gap of 0.2 between two dopaminergic
#' populations is a stronger signal than 0.2 between a neuron and an
#' oligodendrocyte.
#'
#' @format A data frame with the same columns as [midbrain_markers].
#'   All markers are in the `identity` layer.
#'
#' @source Control donors, human substantia nigra pars compacta; Kamath
#'   et al. 2022, Broad Single Cell Portal SCP1768 / GEO GSE178265. The
#'   `CALB1_RBP4` subtype is absent: its candidates were housekeeping and
#'   ribosomal genes, which is what a subtype looks like when it has no
#'   distinguishing expression of its own.
#'
#' @examples
#' data(midbrain_da_subtypes)
#' unique(midbrain_da_subtypes$cell_type)
"midbrain_da_subtypes"
