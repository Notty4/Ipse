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
