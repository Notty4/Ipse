#' Score cells against marker gene sets
#'
#' Computes a per-cell-type score for every cell in an expression matrix,
#' from a marker gene reference table. This is the scoring engine behind
#' [annotate_midbrain_cells()]; call it directly if you want raw scores
#' rather than hard-assigned labels.
#'
#' @section Scoring methods:
#' Raw mean expression is a poor basis for comparing scores *across* cell
#' types: highly-expressed markers (e.g. `MBP`) inflate their own set
#' relative to lowly-expressed ones (e.g. `SLC6A3`), and per-cell library
#' size leaks into the result. The available corrections differ in *what*
#' they normalise against, which matters more than it might appear — see
#' `vignette("scoring-methods")` for a benchmark.
#'
#' \describe{
#'   \item{`"zscore"` (default)}{Each gene is standardised across cells, then
#'     the weighted mean of the standardised markers is taken. This is the
#'     only method here that compares a gene against *its own distribution
#'     across cells*, which is what detects enrichment relative to a gene's
#'     own baseline. Consequently it is the most robust to markers that
#'     differ in absolute abundance. Its weakness is the flip side: scores
#'     are relative to the cells present, so it degrades on datasets
#'     dominated by a single cell type, and scores are not comparable
#'     between separately-processed batches.}
#'   \item{`"ucell"`}{Mann-Whitney U statistic on within-cell gene ranks,
#'     following the UCell formulation. Scores lie in `[0, 1]`, are
#'     independent of the other cells in the matrix (so batches and tiles
#'     are comparable, and cells can be scored in chunks), and are robust to
#'     sparsity. The trade-off: because it ranks *within* a cell, a
#'     constitutively abundant marker set ranks highly in every cell, so it
#'     does not correct for abundance imbalance or ambient/spillover
#'     contamination.}
#'   \item{`"control"`}{Background-matched scoring in the style of Seurat's
#'     `AddModuleScore`: genes are binned by mean expression, control genes
#'     are drawn from matching bins, and the control mean is subtracted.
#'     Note this is an *additive* correction, so it removes an offset but
#'     not a difference in dynamic range; it helps less than expected when
#'     marker sets differ in abundance. Requires enough non-marker genes in
#'     `expr` to build the control pool.}
#'   \item{`"mean_weighted"`}{Weighted mean of raw marker expression. Kept
#'     as an uncorrected reference point; not recommended for real use.}
#' }
#'
#' @param expr A genes x cells numeric matrix (dense or sparse `Matrix`) of
#'   normalised (typically log-normalised) expression values. Row names must
#'   be gene symbols.
#' @param markers A marker reference data frame with columns `cell_type`,
#'   `gene`, `weight`, and optionally `layer` and `context` (see
#'   [midbrain_markers]). Defaults to the bundled midbrain reference. A
#'   table without `layer`/`context` is treated as all-identity.
#' @param layer Which marker layer(s) to score: `"identity"` (the default),
#'   `"effector"`, `"state"`, a vector of these, or `NULL` for all. See
#'   [midbrain_markers] for what the layers mean.
#' @param method Scoring method: `"zscore"`, `"ucell"`, `"control"`, or
#'   `"mean_weighted"`. See Scoring methods.
#' @param n_control Number of control genes sampled per marker gene
#'   (`method = "control"` only).
#' @param n_bins Number of expression bins used to match control genes to
#'   markers (`method = "control"` only).
#' @param max_rank Rank cutoff below which genes are treated as tied
#'   (`method = "ucell"` only).
#' @param seed Optional integer seed. Control gene sampling is random, so
#'   set this for reproducible scores when `method = "control"`.
#'
#' @return A cell-types x cells numeric matrix of scores. A row is all `NA`
#'   (with a warning) if none of that cell type's markers were found in
#'   `expr`.
#'
#' @examples
#' \dontrun{
#' scores <- score_markers(expr_matrix, seed = 42)
#' scores[, 1:5]
#' }
#' @export
score_markers <- function(expr,
                          markers = midbrain_markers,
                          layer = "identity",
                          method = c("zscore", "ucell", "control", "mean_weighted"),
                          n_control = 100,
                          n_bins = 25,
                          max_rank = 1500,
                          seed = NULL) {
  method <- match.arg(method)

  .validate_expr(expr)
  markers <- .validate_markers(markers)
  markers <- .filter_layer(markers, layer)

  if (!is.null(seed)) {
    old_seed <- .capture_seed()
    on.exit(.restore_seed(old_seed), add = TRUE)
    set.seed(seed)
  }

  all_marker_genes <- unique(markers$gene)

  # Method-specific one-off preparation, done once rather than per cell type.
  prep <- switch(
    method,
    zscore  = list(expr = .zscore_rows(expr)),
    control = .prepare_control_bins(expr, all_marker_genes, n_bins),
    ucell   = list(ranks = .rank_cells(expr, max_rank)),
    list()
  )
  if (method == "zscore") expr <- prep$expr

  cell_types <- unique(markers$cell_type)
  scores <- matrix(
    NA_real_,
    nrow = length(cell_types),
    ncol = ncol(expr),
    dimnames = list(cell_types, colnames(expr))
  )

  for (ct in cell_types) {
    ct_markers <- markers[markers$cell_type == ct, ]
    present <- ct_markers$gene %in% rownames(expr)

    if (!any(present)) {
      warning(sprintf(
        "None of the markers for '%s' were found in `expr`; scores set to NA.",
        ct
      ))
      next
    }

    genes <- ct_markers$gene[present]
    weights <- ct_markers$weight[present]
    weights <- weights / sum(weights)

    if (method == "ucell") {
      scores[ct, ] <- .ucell_score(prep$ranks, genes, max_rank)
      next
    }

    marker_score <- as.numeric(crossprod(expr[genes, , drop = FALSE], weights))

    if (method == "control" && !is.null(prep$bins)) {
      control_genes <- .sample_control_genes(genes, prep$bins, prep$pool, n_control)
      if (length(control_genes) > 0) {
        control_score <- Matrix::colMeans(expr[control_genes, , drop = FALSE])
        marker_score <- marker_score - as.numeric(control_score)
      }
    }

    scores[ct, ] <- marker_score
  }

  scores
}

#' Standardise each gene (row) across cells
#'
#' Zero-variance genes are returned as all-zero rather than `NaN`.
#' @keywords internal
#' @noRd
.zscore_rows <- function(expr) {
  mu <- Matrix::rowMeans(expr)
  # var = E[x^2] - E[x]^2, computed this way to stay sparse-friendly
  sq <- Matrix::rowMeans(expr * expr)
  sdev <- sqrt(pmax(sq - mu^2, 0)) * sqrt(ncol(expr) / max(ncol(expr) - 1, 1))

  zero_var <- sdev < .Machine$double.eps
  sdev[zero_var] <- 1

  out <- as.matrix(expr)
  out <- (out - mu) / sdev
  out[zero_var, ] <- 0
  dimnames(out) <- dimnames(expr)
  out
}

#' Bin genes by mean expression to build a control gene pool
#'
#' Returns `NULL` bins (with a warning) when there are too few non-marker
#' genes to form a meaningful background.
#' @keywords internal
#' @noRd
.prepare_control_bins <- function(expr, marker_genes, n_bins) {
  pool <- setdiff(rownames(expr), marker_genes)

  if (length(pool) < 10) {
    warning(
      "Too few non-marker genes in `expr` to build a control pool ",
      sprintf("(%d found). ", length(pool)),
      "Falling back to uncorrected weighted mean scores; ",
      "pass the full expression matrix, not just marker genes, ",
      "or use method = \"zscore\".",
      call. = FALSE
    )
    return(list(bins = NULL, pool = pool))
  }

  gene_means <- Matrix::rowMeans(expr)
  n_bins <- min(n_bins, max(length(gene_means) %/% 5, 2))

  breaks <- stats::quantile(gene_means, probs = seq(0, 1, length.out = n_bins + 1))
  breaks <- unique(breaks)
  if (length(breaks) < 3) {
    bins <- stats::setNames(rep(1L, length(gene_means)), names(gene_means))
  } else {
    bins <- cut(gene_means, breaks = breaks, include.lowest = TRUE, labels = FALSE)
    bins <- stats::setNames(as.integer(bins), names(gene_means))
  }

  list(bins = bins, pool = pool)
}

#' Sample expression-matched control genes for a marker set
#' @keywords internal
#' @noRd
.sample_control_genes <- function(genes, bins, pool, n_control) {
  pool_bins <- bins[pool]
  controls <- character(0)

  for (g in genes) {
    target_bin <- bins[[g]]
    candidates <- pool[pool_bins == target_bin]
    if (length(candidates) == 0) next
    n_take <- min(n_control, length(candidates))
    controls <- c(controls, sample(candidates, n_take, replace = FALSE))
  }

  unique(controls)
}

#' @keywords internal
#' @noRd
.capture_seed <- function() {
  if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) {
    get(".Random.seed", envir = globalenv(), inherits = FALSE)
  } else {
    NULL
  }
}

#' @keywords internal
#' @noRd
.restore_seed <- function(old_seed) {
  if (is.null(old_seed)) {
    if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) {
      rm(".Random.seed", envir = globalenv())
    }
  } else {
    assign(".Random.seed", old_seed, envir = globalenv())
  }
}

#' @keywords internal
#' @noRd
.validate_expr <- function(expr) {
  if (!is.matrix(expr) && !methods::is(expr, "Matrix")) {
    stop("`expr` must be a matrix or Matrix (genes x cells).", call. = FALSE)
  }
  if (is.null(rownames(expr))) {
    stop("`expr` must have gene symbols as row names.", call. = FALSE)
  }
  if (is.null(colnames(expr))) {
    stop("`expr` must have cell IDs as column names.", call. = FALSE)
  }
  if (anyDuplicated(rownames(expr))) {
    stop("`expr` has duplicated gene symbols in its row names.", call. = FALSE)
  }
}

#' Validate a marker table, backfilling optional columns
#'
#' `layer` and `context` are optional: a table without them is treated as
#' all-identity, canonical, which keeps pre-layer marker tables working.
#' @keywords internal
#' @noRd
.validate_markers <- function(markers) {
  required_cols <- c("cell_type", "gene", "weight")
  missing_cols <- setdiff(required_cols, colnames(markers))
  if (length(missing_cols) > 0) {
    stop(
      sprintf("`markers` is missing column(s): %s", paste(missing_cols, collapse = ", ")),
      call. = FALSE
    )
  }
  if (is.null(markers$layer))   markers$layer <- "identity"
  if (is.null(markers$context)) markers$context <- "canonical"
  if (is.null(markers$source))  markers$source  <- NA_character_

  valid <- c("identity", "effector", "state")
  bad <- setdiff(unique(markers$layer), valid)
  if (length(bad) > 0) {
    stop(sprintf("Unknown marker layer(s): %s. Must be one of: %s",
                 paste(bad, collapse = ", "), paste(valid, collapse = ", ")),
         call. = FALSE)
  }
  markers
}

#' Subset a marker table to one or more layers
#' @keywords internal
#' @noRd
.filter_layer <- function(markers, layer) {
  if (is.null(layer)) return(markers)
  out <- markers[markers$layer %in% layer, , drop = FALSE]
  if (nrow(out) == 0) {
    stop(sprintf("No markers found in layer(s): %s", paste(layer, collapse = ", ")),
         call. = FALSE)
  }
  out
}

#' Rank genes within each cell, descending by expression
#'
#' Genes falling below `max_rank` are collapsed to `max_rank + 1`, which
#' both bounds the score and makes the computation robust to the long tail
#' of lowly/zero-expressed genes typical of sparse spatial data.
#' @keywords internal
#' @noRd
.rank_cells <- function(expr, max_rank) {
  dense <- as.matrix(expr)
  ranks <- apply(dense, 2, function(x) {
    r <- rank(-x, ties.method = "average")
    r[r > max_rank] <- max_rank + 1
    r
  })
  rownames(ranks) <- rownames(dense)
  colnames(ranks) <- colnames(dense)
  ranks
}

#' Mann-Whitney U based enrichment score for one marker set
#'
#' Follows the UCell formulation: the U statistic of the marker ranks is
#' normalised by its theoretical maximum, so scores lie in `[0, 1]` and are
#' comparable across marker sets of different size *and* across cell types
#' whose markers differ in absolute abundance.
#' @keywords internal
#' @noRd
.ucell_score <- function(ranks, genes, max_rank) {
  n <- length(genes)
  sub <- ranks[genes, , drop = FALSE]
  rank_sum <- colSums(sub)
  u <- rank_sum - n * (n + 1) / 2
  1 - u / (n * max_rank)
}
