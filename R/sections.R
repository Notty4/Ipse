#' List the state sections available in a marker reference
#'
#' State sections are condition-associated signatures (disease states,
#' reactive states) that are **never scored unless you ask for them** by
#' name via the `sections` argument of [annotate_midbrain_cells()].
#'
#' The point of a section is interpretability: it gives a named,
#' reference-anchored cell group that can be compared against published
#' work, rather than an opaque cluster index that means nothing outside the
#' run that produced it. The `source` column records where each signature
#' came from, so you can judge whether a given analysis is independent of
#' the data the signature was derived from.
#'
#' @param markers A marker reference data frame. Defaults to the bundled
#'   midbrain reference.
#'
#' @return A data frame with one row per section: `context` (the name you
#'   pass to `sections`), `cell_types` it applies to, `n_genes`, and
#'   `source`. Zero rows if the reference defines no sections.
#'
#' @examples
#' list_state_sections()
#' @export
list_state_sections <- function(markers = midbrain_markers) {
  markers <- .validate_markers(markers)
  state <- markers[markers$layer == "state", , drop = FALSE]

  if (nrow(state) == 0) {
    return(data.frame(context = character(0), cell_types = character(0),
                      n_genes = integer(0), source = character(0),
                      stringsAsFactors = FALSE))
  }

  contexts <- unique(state$context)
  do.call(rbind, lapply(contexts, function(ctx) {
    sub <- state[state$context == ctx, , drop = FALSE]
    src <- unique(stats::na.omit(sub$source))
    data.frame(
      context = ctx,
      cell_types = paste(unique(sub$cell_type), collapse = ", "),
      n_genes = length(unique(sub$gene)),
      source = if (length(src) == 0) NA_character_ else paste(src, collapse = "; "),
      row.names = NULL,
      stringsAsFactors = FALSE
    )
  }))
}

#' Resolve and validate a requested set of sections
#' @keywords internal
#' @noRd
.resolve_sections <- function(markers, sections) {
  available <- unique(markers$context[markers$layer == "state"])

  if (is.null(sections)) return(character(0))

  if (identical(sections, "all")) return(available)

  unknown <- setdiff(sections, available)
  if (length(unknown) > 0) {
    stop(sprintf(
      "Unknown section(s): %s. Available: %s. See list_state_sections().",
      paste(unknown, collapse = ", "),
      if (length(available) == 0) "none in this reference"
      else paste(available, collapse = ", ")
    ), call. = FALSE)
  }
  sections
}
