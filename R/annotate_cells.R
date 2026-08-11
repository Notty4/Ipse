#' Annotate midbrain cell types from expression data
#'
#' Assigns each cell (or spot) the cell type it scores highest against, and
#' — separately — reports how that cell's effector programme and any
#' condition-associated state signatures are scoring.
#'
#' @section Why identity and state are kept separate:
#' The label is assigned from `identity` markers **only**. Effector and
#' state markers are scored and returned, but never contribute to the
#' label. This is deliberate: if disease-associated genes helped decide the
#' label, then the annotation would partly encode disease status, and any
#' downstream comparison of cell type composition between conditions would
#' be partly circular.
#'
#' Keeping them apart lets a cell be described as retaining its lineage
#' identity while having lost its effector programme — which a single
#' categorical label cannot express. For cells mid-transition, the `margin`
#' column and the full score matrix are more informative than the label,
#' and a cluster of low-margin cells concentrated in one condition is a
#' result rather than a QC failure.
#'
#' @param expr A genes x cells numeric matrix (dense or sparse `Matrix`) of
#'   normalised expression values. Row names must be gene symbols, column
#'   names cell/spot IDs. Pass the full matrix rather than a marker-only
#'   subset.
#' @param markers A marker reference data frame (see [midbrain_markers]).
#' @param method Scoring method passed to [score_markers()].
#' @param min_score Minimum winning identity score required to accept a
#'   call; cells below it for every cell type are labelled `"Unassigned"`.
#' @param min_margin Minimum gap between best and second-best identity
#'   score; cells below it are labelled `"Ambiguous"`. Default `0`.
#' @param sections Which state sections to run, by name — e.g.
#'   `"SCZ"`, or `c("SCZ", "reactive")`. Sections are **opt-in**: the
#'   default `NULL` runs none, so no condition-associated signature is
#'   scored unless you ask for it. Pass `"all"` for every section in
#'   `markers`. See [list_state_sections()] for what is available and
#'   where each signature came from.
#' @param score_effector If `TRUE` (default), also score the `effector`
#'   layer and return it as a column.
#' @param seed Optional integer seed, passed to [score_markers()].
#' @param return_scores If `TRUE`, also return the underlying score
#'   matrices.
#'
#' @return A data frame with one row per cell: `cell_id`, `cell_type`,
#'   `score`, `margin`, `effector_score`, and one `state_<context>` column
#'   per section named in `sections`. When any section ran, the result
#'   carries a `"sections"` attribute recording which ones and their
#'   sources.
#'   Effector and state scores are taken for the cell's *assigned* cell
#'   type, and are `NA` where that type has no markers in the layer.
#'   If `return_scores = TRUE`, returns a list of `labels` and the score
#'   matrices instead.
#'
#' @examples
#' \dontrun{
#' labels <- annotate_midbrain_cells(expr_matrix, seed = 42)
#' table(labels$cell_type)
#'
#' # Cells that kept lineage identity but lost the effector programme
#' subset(labels, cell_type == "Dopaminergic neuron" & effector_score < 0)
#'
#' # Opt in to a named section, then inspect its provenance
#' labels <- annotate_midbrain_cells(expr_matrix, sections = "reactive")
#' attr(labels, "sections")
#' }
#' @export
annotate_midbrain_cells <- function(expr,
                                    markers = midbrain_markers,
                                    method = c("zscore", "ucell", "control", "mean_weighted"),
                                    min_score = -Inf,
                                    min_margin = 0,
                                    sections = NULL,
                                    score_effector = TRUE,
                                    seed = NULL,
                                    return_scores = FALSE) {
  method <- match.arg(method)
  markers <- .validate_markers(markers)

  identity_scores <- score_markers(
    expr, markers = markers, layer = "identity",
    method = method, seed = seed
  )

  top_type <- apply(identity_scores, 2, function(col) {
    if (all(is.na(col))) return(NA_character_)
    rownames(identity_scores)[which.max(col)]
  })
  top_score <- apply(identity_scores, 2, function(col) {
    if (all(is.na(col))) return(NA_real_)
    max(col, na.rm = TRUE)
  })
  margin <- apply(identity_scores, 2, function(col) {
    col <- sort(col[!is.na(col)], decreasing = TRUE)
    if (length(col) < 2) return(NA_real_)
    col[1] - col[2]
  })

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

  extra_scores <- list()

  if (score_effector && any(markers$layer == "effector")) {
    eff <- score_markers(expr, markers = markers, layer = "effector",
                         method = method, seed = seed)
    labels$effector_score <- .score_for_assigned(eff, top_type)
    extra_scores$effector <- eff
  }

  run_sections <- .resolve_sections(markers, sections)
  state_markers <- markers[markers$layer == "state", , drop = FALSE]

  for (ctx in run_sections) {
    ctx_markers <- state_markers[state_markers$context == ctx, , drop = FALSE]
    st <- score_markers(expr, markers = ctx_markers, layer = "state",
                        method = method, seed = seed)
    labels[[paste0("state_", ctx)]] <- .score_for_assigned(st, top_type)
    extra_scores[[paste0("state_", ctx)]] <- st
  }

  # Record which sections ran and where their signatures came from, so the
  # provenance travels with the result rather than living only in the call.
  attr(labels, "sections") <- if (length(run_sections) == 0) {
    NULL
  } else {
    sect <- list_state_sections(markers)
    sect[sect$context %in% run_sections, , drop = FALSE]
  }

  if (return_scores) {
    return(c(list(labels = labels, identity = identity_scores), extra_scores))
  }

  labels
}

#' Pull each cell's score for the cell type it was assigned
#'
#' A signature marked `"Any"` applies to every cell type; otherwise a cell
#' whose assigned type is absent from the score matrix (including
#' `"Unassigned"` and `"Ambiguous"`) gets `NA`.
#' @keywords internal
#' @noRd
.score_for_assigned <- function(score_matrix, assigned) {
  if (identical(rownames(score_matrix), "Any")) {
    return(as.numeric(score_matrix["Any", ]))
  }
  vapply(seq_along(assigned), function(i) {
    ct <- assigned[i]
    if (is.na(ct) || !(ct %in% rownames(score_matrix))) return(NA_real_)
    score_matrix[ct, i]
  }, numeric(1))
}

#' List the cell types covered by a marker reference
#'
#' @param markers A marker reference data frame.
#' @param layer Restrict to cell types having markers in this layer, or
#'   `NULL` for any layer. Defaults to `"identity"` — the layer that
#'   actually determines labels.
#' @return A character vector of unique cell type labels.
#' @examples
#' list_midbrain_cell_types()
#' @export
list_midbrain_cell_types <- function(markers = midbrain_markers,
                                     layer = "identity") {
  markers <- .validate_markers(markers)
  if (!is.null(layer)) markers <- markers[markers$layer %in% layer, , drop = FALSE]
  unique(markers$cell_type)
}
