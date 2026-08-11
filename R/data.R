#' Midbrain marker gene reference table
#'
#' A marker gene reference for major midbrain cell types, bundled with the
#' package so that [annotate_midbrain_cells()] can run without an external
#' reference download. This is a starter set — see `data-raw/prepare_markers.R`
#' for the build script and a note on curating it further.
#'
#' @format A data frame with one row per (cell type, marker gene) pair and
#'   the following columns:
#' \describe{
#'   \item{cell_type}{Character. Cell type label, e.g. `"Dopaminergic neuron"`.}
#'   \item{gene}{Character. Marker gene symbol.}
#'   \item{weight}{Numeric. Relative weight of this marker within its cell
#'     type's score (default 1 for all markers).}
#' }
#'
#' @source Starter set of canonical midbrain markers; not yet curated
#'   against a specific reference dataset. See `data-raw/prepare_markers.R`.
#'
#' @examples
#' data(midbrain_markers)
#' head(midbrain_markers)
#' unique(midbrain_markers$cell_type)
"midbrain_markers"
