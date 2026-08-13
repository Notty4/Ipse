#' Sub-label cells of one type using a subtype reference
#'
#' Takes labels from [annotate_cells()] and resolves one broad class into
#' its subtypes — dopaminergic neurons into the SOX6 and CALB1 populations,
#' for instance. Cells of other types are left alone.
#'
#' @section Why this is a separate pass:
#' Subtype scores are computed **within the parent class only**. That is
#' the whole point: a dopaminergic subtype's meaningful competitor is
#' another dopaminergic subtype, not an oligodendrocyte. Scoring subtypes
#' alongside broad lineages collapses their scores, because every subtype's
#' nearest neighbour is a sibling and exclusivity vanishes — the same
#' reason [midbrain_da_subtypes] is derived separately from
#' [midbrain_markers].
#'
#' A consequence worth knowing: with composition-dependent methods such as
#' `"zscore"`, the scores here depend on which cells were assigned to the
#' parent class. Mislabelled cells entering this pass will distort it.
#'
#' @section Reading the output:
#' Expect smaller margins than at the lineage level. Subtypes are
#' genuinely similar, so a `subtype_margin` of 0.05 is not necessarily a
#' failure — but a population of low-margin cells concentrated in one
#' condition is a result, not noise. Setting `min_margin` marks those
#' `"Ambiguous"` rather than forcing a call.
#'
#' @param expr A genes x cells matrix, as passed to [annotate_cells()].
#' @param labels The data frame returned by [annotate_cells()], with one
#'   row per column of `expr` in the same order.
#' @param parent The `cell_type` value to resolve, e.g.
#'   `"Dopaminergic neuron"`.
#' @param markers A subtype marker reference. Defaults to
#'   [midbrain_da_subtypes].
#' @param method Scoring method, passed to [score_markers()].
#' @param min_score Minimum winning score to accept a subtype call; below
#'   it, cells are `"Unassigned"`.
#' @param min_margin Minimum gap between best and second-best subtype
#'   score; below it, cells are `"Ambiguous"`.
#' @param seed Optional integer seed, passed to [score_markers()].
#'
#' @return `labels` with three columns added: `subtype`, `subtype_score`
#'   and `subtype_margin`. All three are `NA` for cells whose `cell_type`
#'   is not `parent`.
#'
#' @examples
#' \dontrun{
#' labels <- annotate_cells(expr, seed = 42)
#' labels <- annotate_subtypes(expr, labels, parent = "Dopaminergic neuron",
#'                             seed = 42)
#' table(labels$subtype, useNA = "ifany")
#'
#' # SOX6_AGTR1 is the population depleted in Parkinson's disease
#' sum(labels$subtype == "SOX6_AGTR1", na.rm = TRUE)
#' }
#' @export
annotate_subtypes <- function(expr,
                              labels,
                              parent,
                              markers = midbrain_da_subtypes,
                              method = c("zscore", "ucell", "control", "mean_weighted"),
                              min_score = -Inf,
                              min_margin = 0,
                              seed = NULL) {
  method <- match.arg(method)
  markers <- .validate_markers(markers)

  if (!is.data.frame(labels) || !all(c("cell_id", "cell_type") %in% names(labels))) {
    stop("`labels` must be the data frame returned by annotate_cells().",
         call. = FALSE)
  }
  if (nrow(labels) != ncol(expr)) {
    stop("`labels` must have one row per column of `expr`.", call. = FALSE)
  }
  if (length(parent) != 1L || is.na(parent)) {
    stop("`parent` must be a single cell type name.", call. = FALSE)
  }

  labels$subtype <- NA_character_
  labels$subtype_score <- NA_real_
  labels$subtype_margin <- NA_real_

  target <- which(labels$cell_type == parent)
  if (length(target) == 0) {
    warning(sprintf("No cells labelled '%s'; nothing to sub-label. Present: %s",
                    parent, paste(sort(unique(labels$cell_type)), collapse = ", ")),
            call. = FALSE)
    return(labels)
  }

  scores <- score_markers(
    expr[, target, drop = FALSE],
    markers = markers, layer = "identity",
    method = method, seed = seed
  )

  top <- apply(scores, 2, function(col) {
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

  top[!is.na(top_score) & top_score < min_score] <- "Unassigned"
  if (min_margin > 0) {
    top[!is.na(margin) & margin < min_margin & top != "Unassigned"] <- "Ambiguous"
  }

  labels$subtype[target] <- top
  labels$subtype_score[target] <- top_score
  labels$subtype_margin[target] <- margin

  labels
}
