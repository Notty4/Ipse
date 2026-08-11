#' Annotate midbrain cell types from expression data
#'
#' Assigns each cell (or spot) the cell type it scores highest against,
#' using the bundled midbrain marker reference by default. No external
#' atlas or reference download is required.
#'
#' Scores are z-scored across cells by default, which makes them comparable
#' across cell types whose markers differ in absolute abundance. See
#' [score_markers()] for the available methods and their trade-offs — the
#' choice matters, and the right one depends on your data.
#'
#' @param expr A genes x cells numeric matrix (dense or sparse `Matrix`) of
#'   normalised expression values. Row names must be gene symbols, column
#'   names cell/spot IDs. Pass the full matrix rather than a marker-only
#'   subset, so that control genes can be drawn.
#' @param markers A marker reference data frame (see [midbrain_markers]).
#'   Defaults to the bundled midbrain reference.
#' @param method Scoring method passed to [score_markers()].
#' @param min_score Minimum winning score required to accept a call; cells
#'   scoring below this for every cell type are labelled `"Unassigned"`.
#'   Default `-Inf` (always call the top-scoring type). With the
#'   background-corrected methods, `0` is a natural starting threshold: it
#'   requires a cell to exceed its expression-matched background.
#' @param min_margin Minimum gap between the best and second-best score
#'   required to accept a call. Cells where the top two types are within
#'   `min_margin` are labelled `"Ambiguous"`. Default `0` (no check).
#' @param seed Optional integer seed, passed to [score_markers()]. Control
#'   gene sampling is random, so set this for reproducible labels.
#' @param return_scores If `TRUE`, also return the full score matrix.
#'
#' @return A data frame with one row per cell: `cell_id`, `cell_type`,
#'   `score` (the winning score), and `margin` (gap to the runner-up). If
#'   `return_scores = TRUE`, returns `list(labels = <data.frame>,
#'   scores = <matrix>)` instead.
#'
#' @examples
#' \dontrun{
#' labels <- annotate_midbrain_cells(expr_matrix, seed = 42)
#' table(labels$cell_type)
#'
#' # Flag weak and ambiguous calls rather than forcing every cell
#' labels <- annotate_midbrain_cells(
#'   expr_matrix, min_score = 0, min_margin = 0.1, seed = 42
#' )
#' }
#' @export
annotate_midbrain_cells <- function(expr,
                                    markers = midbrain_markers,
                                    method = c("zscore", "ucell", "control", "mean_weighted"),
                                    min_score = -Inf,
                                    min_margin = 0,
                                    seed = NULL,
                                    return_scores = FALSE) {
  method <- match.arg(method)

  scores <- score_markers(
    expr,
    markers = markers,
    method = method,
    seed = seed
  )

  top_type <- apply(scores, 2, function(col) {
    if (all(is.na(col))) return(NA_character_)
    rownames(scores)[which.max(col)]
  })
  top_score <- apply(scores, 2, function(col) {
    if (all(is.na(col))) return(NA_real_)
    max(col, na.rm = TRUE)
  })
  margin <- apply(scores, 2, function(col) {
    col <- sort(col[!is.na(col)], decreasing = TRUE)
    if (length(col) < 2) return(NA_real_)
    col[1] - col[2]
  })

  # Order matters: an ambiguous-but-strong call is more informative to the
  # user than a bare "Unassigned", so margin is applied after min_score.
  top_type[!is.na(top_score) & top_score < min_score] <- "Unassigned"
  if (min_margin > 0) {
    top_type[!is.na(margin) & margin < min_margin &
               top_type != "Unassigned"] <- "Ambiguous"
  }

  labels <- data.frame(
    cell_id = colnames(expr),
    cell_type = top_type,
    score = top_score,
    margin = margin,
    row.names = NULL,
    stringsAsFactors = FALSE
  )

  if (return_scores) {
    return(list(labels = labels, scores = scores))
  }

  labels
}

#' List the cell types covered by a marker reference
#'
#' @param markers A marker reference data frame. Defaults to the bundled
#'   midbrain reference.
#' @return A character vector of unique cell type labels.
#' @examples
#' list_midbrain_cell_types()
#' @export
list_midbrain_cell_types <- function(markers = midbrain_markers) {
  unique(markers$cell_type)
}
