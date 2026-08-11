#' Compare scoring methods against known labels
#'
#' Runs every scoring method over the same data and marker set and reports
#' how well each recovers labels you already trust. Use it whenever you
#' have annotated data, because which method wins is a property of your
#' dataset rather than a fact about the methods.
#'
#' @section Why there is no universal default:
#' The methods fail in different circumstances, and those circumstances are
#' properties of your data:
#'
#' \itemize{
#'   \item `"zscore"` standardises each gene across cells, so it is robust
#'     when marker sets differ in absolute abundance — but it degrades
#'     badly under **class imbalance**. A cell type comprising most of the
#'     dataset pulls up its own markers' means, compressing its members'
#'     scores, while a 1% type's markers sit near zero for everyone. Weakly
#'     expressing members of the abundant type then fall below the rare
#'     type's baseline and are misassigned to it.
#'   \item `"ucell"` scores each cell independently, so composition cannot
#'     distort it — but it cannot correct for a marker set that is
#'     constitutively abundant across all cell types.
#'   \item `"control"` and `"mean_weighted"` are unaffected by composition
#'     but do not correct abundance imbalance either.
#' }
#'
#' Real data is nearly always class-imbalanced, and often also carries
#' ambient contamination, so both failure modes can be present at once.
#' Measure rather than assume.
#'
#' @param expr A genes x cells matrix, as for [annotate_cells()].
#' @param labels Trusted cell type labels, one per column of `expr`.
#' @param markers A marker reference data frame.
#' @param methods Which methods to compare.
#' @param calibrate Which calibration settings to compare; every
#'   combination of `methods` and `calibrate` is evaluated. See
#'   [annotate_cells()].
#' @param seed Passed to [score_markers()] for reproducibility.
#'
#' @return A data frame with one row per method and calibration setting:
#'   `method`, `calibrate`, `accuracy`
#'   overall, `balanced_accuracy` (mean of per-class recall, which does not
#'   let one large class mask failures on small ones), and `worst_class`
#'   with its recall. Sorted by balanced accuracy.
#'
#' @examples
#' \dontrun{
#' benchmark_methods(expr, obj$cell_type, markers = my_markers)
#' }
#' @export
benchmark_methods <- function(expr,
                              labels,
                              markers = midbrain_markers,
                              methods = c("zscore", "ucell", "control", "mean_weighted"),
                              calibrate = "none",
                              seed = 42) {
  .validate_expr(expr)
  if (length(labels) != ncol(expr)) {
    stop("`labels` must have one entry per column of `expr`.", call. = FALSE)
  }

  labels <- as.character(labels)
  classes <- sort(unique(labels[!is.na(labels)]))

  # A class with no identity markers can never be predicted, so it scores
  # zero recall for every method and drags balanced accuracy down
  # uniformly -- which hides the actual differences between methods.
  markers <- .validate_markers(markers)
  labelable <- unique(markers$cell_type[markers$layer == "identity"])
  impossible <- setdiff(classes, labelable)
  if (length(impossible) > 0) {
    warning(sprintf(
      "No identity markers for: %s. These can never be predicted, so their recall is 0 for every method and balanced accuracy is not comparable. Consider excluding these cells before benchmarking.",
      paste(impossible, collapse = ", ")), call. = FALSE)
  }

  grid <- expand.grid(method = methods, calibrate = calibrate,
                      stringsAsFactors = FALSE)

  rows <- lapply(seq_len(nrow(grid)), function(i) {
    m <- grid$method[i]
    cal <- grid$calibrate[i]
    pred <- annotate_cells(expr, markers = markers, method = m,
                           calibrate = cal, seed = seed)$cell_type

    recall <- vapply(classes, function(cl) {
      idx <- which(labels == cl)
      if (length(idx) == 0) return(NA_real_)
      mean(pred[idx] == cl)
    }, numeric(1))

    worst <- names(recall)[which.min(recall)]
    data.frame(
      method = m,
      calibrate = cal,
      accuracy = mean(pred == labels, na.rm = TRUE),
      balanced_accuracy = mean(recall, na.rm = TRUE),
      worst_class = worst,
      worst_recall = min(recall, na.rm = TRUE),
      row.names = NULL,
      stringsAsFactors = FALSE
    )
  })

  out <- do.call(rbind, rows)
  out <- out[order(-out$balanced_accuracy), , drop = FALSE]
  rownames(out) <- NULL
  out
}
