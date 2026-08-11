#' Derive candidate marker genes from an annotated reference dataset
#'
#' Ranks genes by how *exclusively* they mark each cell type, given an
#' annotated expression matrix. Intended to be run on control (or
#' neurotypical) samples only: identity markers derived from disease tissue
#' risk encoding disease status into the labels they later assign.
#'
#' @section Why exclusivity rather than fold change:
#' The usual differential-expression ranking asks whether a gene is higher
#' in this cell type than in the average of the rest. That tolerates a gene
#' being equally high in one *other* type, which is exactly the failure
#' that defeats marker scoring in practice (see
#' `vignette("scoring-methods")`). This function instead compares each cell
#' type against its strongest competitor: `specificity` is the gap between
#' mean expression in the type and mean expression in whichever *other*
#' type expresses it most. A gene shared by two lineages scores near zero
#' however abundant it is.
#'
#' `detect_gap` applies the same idea to detection rate, and is the
#' default ranking. On simulated data with a gene made abundant in every
#' cell type but slightly higher in one — the ambient/spillover pattern
#' that defeats marker scoring — `specificity` still ranks it second,
#' because a large raw gap survives log transformation. `detect_gap` gives
#' it a score of zero, since it is detected everywhere. Ranking on
#' detection is therefore the safer default, particularly for spatial
#' panels; `specificity` is available when expression level matters more
#' than presence/absence.
#'
#' @section Donor consistency:
#' When `donor` is supplied, `donor_consistency` reports the fraction of
#' donors in which the gene remains positively specific to the type. A
#' marker driven by one donor will rank highly on `specificity` but poorly
#' here, which is usually the more informative signal.
#'
#' @param expr A genes x cells matrix (dense or sparse) of normalised,
#'   log-transformed expression. Row names must be gene symbols.
#' @param labels Character vector of cell type labels, one per column of
#'   `expr`.
#' @param donor Optional vector of donor/sample IDs, one per column, used
#'   to compute `donor_consistency`.
#' @param panel Optional character vector of gene symbols to restrict to —
#'   e.g. the genes on your spatial panel. Applied before ranking, so the
#'   returned ranks reflect the genes you can actually measure.
#' @param rank_by Which statistic to rank on: `"detect_gap"` (default) or
#'   `"specificity"`. See the section above for why detection is the safer
#'   default.
#' @param min_cells Cell types with fewer cells than this are dropped, with
#'   a warning.
#' @param min_detect Minimum detection rate within the type for a gene to
#'   be considered a candidate at all.
#'
#' @return A data frame, one row per (cell type, gene) candidate, sorted by
#'   cell type then descending `specificity`, with columns: `cell_type`,
#'   `gene`, `mean_in`, `detect_in`, `mean_out_max`, `competitor` (the
#'   other cell type that expresses it most), `specificity`, `detect_gap`,
#'   `donor_consistency`, and `rank` within the cell type.
#'
#' @seealso [consensus_markers()] to combine candidates across cohorts,
#'   [as_marker_table()] to convert them into a reference table.
#'
#' @examples
#' \dontrun{
#' candidates <- derive_markers(
#'   expr   = counts[, meta$condition == "control"],
#'   labels = meta$cell_type[meta$condition == "control"],
#'   donor  = meta$donor[meta$condition == "control"],
#'   panel  = xenium_panel_genes
#' )
#' subset(candidates, cell_type == "Dopaminergic neuron" & rank <= 10)
#' }
#' @export
derive_markers <- function(expr,
                           labels,
                           donor = NULL,
                           panel = NULL,
                           rank_by = c("detect_gap", "specificity"),
                           min_cells = 20,
                           min_detect = 0.1) {
  rank_by <- match.arg(rank_by)
  .validate_expr(expr)

  if (length(labels) != ncol(expr)) {
    stop("`labels` must have one entry per column of `expr`.", call. = FALSE)
  }
  if (!is.null(donor) && length(donor) != ncol(expr)) {
    stop("`donor` must have one entry per column of `expr`.", call. = FALSE)
  }

  if (!is.null(panel)) {
    keep <- intersect(rownames(expr), panel)
    if (length(keep) == 0) {
      stop("None of the genes in `panel` are present in `expr`.", call. = FALSE)
    }
    if (length(keep) < length(panel)) {
      warning(sprintf("%d of %d panel genes not found in `expr`; ignoring those.",
                      length(panel) - length(keep), length(panel)), call. = FALSE)
    }
    expr <- expr[keep, , drop = FALSE]
  }

  labels <- as.character(labels)
  types <- sort(unique(labels[!is.na(labels)]))

  small <- types[vapply(types, function(t) sum(labels == t, na.rm = TRUE), integer(1)) < min_cells]
  if (length(small) > 0) {
    warning(sprintf("Dropping cell type(s) with fewer than %d cells: %s",
                    min_cells, paste(small, collapse = ", ")), call. = FALSE)
    types <- setdiff(types, small)
  }
  if (length(types) < 2) {
    stop("At least two cell types are needed to assess specificity.", call. = FALSE)
  }

  stats_by_type <- .type_summaries(expr, labels, types)
  out <- .rank_specificity(stats_by_type, types, min_detect)

  if (!is.null(donor)) {
    out$donor_consistency <- .donor_consistency(expr, labels, donor, types, out)
  } else {
    out$donor_consistency <- NA_real_
  }

  out <- out[order(out$cell_type, -out[[rank_by]]), , drop = FALSE]
  out$rank <- stats::ave(out[[rank_by]], out$cell_type,
                         FUN = function(x) rank(-x, ties.method = "first"))
  attr(out, "rank_by") <- rank_by
  rownames(out) <- NULL
  out
}

#' Mean expression and detection rate per gene per cell type
#' @keywords internal
#' @noRd
.type_summaries <- function(expr, labels, types) {
  means <- matrix(NA_real_, nrow(expr), length(types),
                  dimnames = list(rownames(expr), types))
  detect <- means

  for (t in types) {
    idx <- which(labels == t)
    sub <- expr[, idx, drop = FALSE]
    means[, t] <- Matrix::rowMeans(sub)
    detect[, t] <- Matrix::rowMeans(sub > 0)
  }
  list(means = means, detect = detect)
}

#' Compare each cell type against its single strongest competitor
#' @keywords internal
#' @noRd
.rank_specificity <- function(s, types, min_detect) {
  rows <- lapply(types, function(t) {
    others <- setdiff(types, t)
    other_means <- s$means[, others, drop = FALSE]

    competitor_idx <- max.col(other_means, ties.method = "first")
    mean_out_max <- other_means[cbind(seq_len(nrow(other_means)), competitor_idx)]
    competitor <- others[competitor_idx]

    other_detect <- s$detect[, others, drop = FALSE]
    detect_out_max <- apply(other_detect, 1, max)

    keep <- s$detect[, t] >= min_detect
    data.frame(
      cell_type = t,
      gene = rownames(s$means)[keep],
      mean_in = s$means[keep, t],
      detect_in = s$detect[keep, t],
      mean_out_max = mean_out_max[keep],
      competitor = competitor[keep],
      specificity = s$means[keep, t] - mean_out_max[keep],
      detect_gap = s$detect[keep, t] - detect_out_max[keep],
      row.names = NULL,
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

#' Fraction of donors in which a candidate stays positively specific
#' @keywords internal
#' @noRd
.donor_consistency <- function(expr, labels, donor, types, candidates) {
  donors <- unique(donor[!is.na(donor)])
  key <- paste(candidates$cell_type, candidates$gene, sep = "\r")
  hits <- stats::setNames(numeric(length(key)), key)
  usable <- 0L

  for (d in donors) {
    in_d <- donor == d & !is.na(donor)
    d_types <- types[vapply(types, function(t) sum(labels[in_d] == t, na.rm = TRUE),
                            integer(1)) >= 3]
    if (length(d_types) < 2) next
    usable <- usable + 1L

    s <- .type_summaries(expr[, in_d, drop = FALSE], labels[in_d], d_types)
    for (t in d_types) {
      others <- setdiff(d_types, t)
      spec <- s$means[, t] - apply(s$means[, others, drop = FALSE], 1, max)
      pos <- names(spec)[spec > 0]
      k <- paste(t, pos, sep = "\r")
      k <- k[k %in% key]
      hits[k] <- hits[k] + 1
    }
  }

  if (usable == 0) {
    warning("No donor had enough cells in two or more cell types; ",
            "donor_consistency is NA.", call. = FALSE)
    return(rep(NA_real_, length(key)))
  }
  as.numeric(hits[key]) / usable
}
