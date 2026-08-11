#' Combine marker candidates across independent cohorts
#'
#' A gene that ranks highly in one dataset may be reflecting that dataset's
#' donors, dissection, or protocol. Requiring a candidate to replicate
#' across independent cohorts is a stronger filter than any single-dataset
#' statistic, and is usually worth more than moving to a larger dataset.
#'
#' @param candidate_list A named list of data frames as returned by
#'   [derive_markers()]. Names are used as cohort labels.
#' @param top_n Consider only each cohort's top `top_n` candidates per cell
#'   type. Set `NULL` to use all candidates that passed `min_detect`.
#' @param min_cohorts Minimum number of cohorts a gene must appear in to be
#'   retained. Defaults to all of them.
#'
#' @return A data frame with one row per (cell type, gene) surviving the
#'   filter: `cell_type`, `gene`, `n_cohorts`, `cohorts` (comma-separated),
#'   `mean_specificity`, `min_specificity`, and `worst_rank` — the poorest
#'   rank it achieved in any cohort, which is the conservative view.
#'   Sorted by cell type then `worst_rank`.
#'
#' @examples
#' \dontrun{
#' consensus <- consensus_markers(
#'   list(kamath = cand_kamath, agarwal = cand_agarwal, siletti = cand_siletti),
#'   top_n = 50,
#'   min_cohorts = 2
#' )
#' }
#' @export
consensus_markers <- function(candidate_list, top_n = 50, min_cohorts = NULL) {
  if (!is.list(candidate_list) || length(candidate_list) < 2) {
    stop("`candidate_list` must be a list of at least two candidate tables.",
         call. = FALSE)
  }
  if (is.null(names(candidate_list)) || any(names(candidate_list) == "")) {
    names(candidate_list) <- paste0("cohort", seq_along(candidate_list))
  }
  if (is.null(min_cohorts)) min_cohorts <- length(candidate_list)

  required <- c("cell_type", "gene", "specificity", "rank")
  trimmed <- lapply(names(candidate_list), function(nm) {
    df <- candidate_list[[nm]]
    missing <- setdiff(required, colnames(df))
    if (length(missing) > 0) {
      stop(sprintf("Candidate table '%s' is missing column(s): %s",
                   nm, paste(missing, collapse = ", ")), call. = FALSE)
    }
    if (!is.null(top_n)) df <- df[df$rank <= top_n, , drop = FALSE]
    df$cohort <- nm
    df[, c("cell_type", "gene", "specificity", "rank", "cohort")]
  })

  all_cand <- do.call(rbind, trimmed)
  key <- paste(all_cand$cell_type, all_cand$gene, sep = "\r")
  split_cand <- split(all_cand, key)

  rows <- lapply(split_cand, function(df) {
    data.frame(
      cell_type = df$cell_type[1],
      gene = df$gene[1],
      n_cohorts = length(unique(df$cohort)),
      cohorts = paste(sort(unique(df$cohort)), collapse = ", "),
      mean_specificity = mean(df$specificity),
      min_specificity = min(df$specificity),
      worst_rank = max(df$rank),
      row.names = NULL,
      stringsAsFactors = FALSE
    )
  })

  out <- do.call(rbind, rows)
  out <- out[out$n_cohorts >= min_cohorts, , drop = FALSE]
  out <- out[order(out$cell_type, out$worst_rank), , drop = FALSE]
  rownames(out) <- NULL
  out
}

#' Convert marker candidates into a reference table
#'
#' Takes the output of [derive_markers()] or [consensus_markers()] and
#' produces a table in the format [annotate_midbrain_cells()] expects.
#'
#' Layer assignment is deliberately **not** automated. Whether a gene is
#' lineage-determining (`identity`) or part of a functional programme
#' (`effector`) is a biological judgement, and getting it wrong defeats the
#' purpose of separating them: an effector gene mislabelled as identity
#' will make cells that have lost that programme unrecognisable. Everything
#' defaults to `identity`; name the effector genes explicitly.
#'
#' @param candidates A data frame from [derive_markers()] or
#'   [consensus_markers()].
#' @param n Maximum number of genes to keep per cell type.
#' @param min_gap Quality floor applied before taking the top `n`:
#'   candidates scoring below this on `detect_gap` (or `specificity` if
#'   detection is unavailable) are discarded. Prevents a cell type with
#'   only two real markers from being padded out with noise, which
#'   degrades annotation rather than improving it. Set `NULL` to disable.
#'
#'   The right value is dataset-dependent, not universal: with few cells
#'   per type, background genes reach a detection gap of 0.15 or so by
#'   chance alone, so the default of 0.1 will not clear them. Inspect the
#'   distribution of `detect_gap` in your candidates and set the floor
#'   above the noise, rather than trusting the default.
#' @param effector Character vector of gene symbols to assign to the
#'   `effector` layer instead of `identity`.
#' @param source A string recording where these markers came from — cohort,
#'   accession, or publication. Strongly recommended.
#'
#' @return A data frame with `cell_type`, `gene`, `layer`, `context`,
#'   `source`, and `weight`, ready to pass as `markers`.
#'
#' @examples
#' \dontrun{
#' markers <- as_marker_table(
#'   consensus, n = 6, min_gap = 0.2,
#'   effector = c("TH", "SLC6A3", "SLC18A2", "DDC"),
#'   source = "Control samples, GSE178265 + GSE140231"
#' )
#' annotate_midbrain_cells(expr, markers = markers)
#' }
#' @export
as_marker_table <- function(candidates, n = 6, min_gap = 0.1,
                            effector = character(0),
                            source = NA_character_) {
  if (!all(c("cell_type", "gene") %in% colnames(candidates))) {
    stop("`candidates` must have `cell_type` and `gene` columns.", call. = FALSE)
  }

  order_col <- if ("worst_rank" %in% colnames(candidates)) "worst_rank" else "rank"
  if (!order_col %in% colnames(candidates)) {
    stop("`candidates` must have a `rank` or `worst_rank` column.", call. = FALSE)
  }

  # Taking the top n regardless of quality quietly admits noise when a cell
  # type has fewer than n genuine markers -- which is common, and which
  # degrades annotation rather than improving it. Apply a floor first.
  gap_col <- if ("detect_gap" %in% colnames(candidates)) "detect_gap" else
             if ("min_specificity" %in% colnames(candidates)) "min_specificity" else
             if ("specificity" %in% colnames(candidates)) "specificity" else NULL

  if (!is.null(min_gap) && !is.null(gap_col)) {
    before <- unique(candidates$cell_type)
    candidates <- candidates[candidates[[gap_col]] >= min_gap, , drop = FALSE]
    lost <- setdiff(before, unique(candidates$cell_type))
    if (length(lost) > 0) {
      warning(sprintf(
        "No candidate reached min_gap = %s for cell type(s): %s. They are absent from the marker table.",
        format(min_gap), paste(lost, collapse = ", ")), call. = FALSE)
    }
    if (nrow(candidates) == 0) {
      stop("No candidates passed `min_gap`. Lower it, or check the reference data.",
           call. = FALSE)
    }
  }

  candidates <- candidates[order(candidates$cell_type, candidates[[order_col]]), ]
  keep <- unlist(lapply(split(seq_len(nrow(candidates)), candidates$cell_type),
                        function(i) utils::head(i, n)), use.names = FALSE)
  sel <- candidates[sort(keep), , drop = FALSE]

  out <- data.frame(
    cell_type = sel$cell_type,
    gene = sel$gene,
    layer = ifelse(sel$gene %in% effector, "effector", "identity"),
    context = "canonical",
    source = source,
    weight = 1,
    row.names = NULL,
    stringsAsFactors = FALSE
  )

  no_identity <- setdiff(unique(out$cell_type),
                         unique(out$cell_type[out$layer == "identity"]))
  if (length(no_identity) > 0) {
    warning(sprintf(
      "Cell type(s) left with no identity markers after assigning effectors: %s. These cannot be labelled.",
      paste(no_identity, collapse = ", ")), call. = FALSE)
  }
  out
}
