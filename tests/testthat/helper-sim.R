# Shared simulated datasets for the test suite.
# testthat sources helper-*.R before running any test file.

#' Well-separated data: markers are near-silent outside their own cell type.
make_expr <- function(n_cells = 60, n_bg = 300, seed = 1) {
  set.seed(seed)
  mk <- midbrain_markers
  bg <- if (n_bg > 0) paste0("BG", seq_len(n_bg)) else character(0)
  all_genes <- c(bg, unique(mk$gene))
  truth <- rep(unique(mk$cell_type), length.out = n_cells)

  base <- stats::setNames(rep(0.05, length(all_genes)), all_genes)
  if (n_bg > 0) base[bg] <- stats::runif(n_bg, 0.1, 2)

  expr <- matrix(0, length(all_genes), n_cells,
                 dimnames = list(all_genes, paste0("cell", seq_len(n_cells))))
  for (j in seq_len(n_cells)) {
    lam <- base
    lam[mk$gene[mk$cell_type == truth[j]]] <- 8
    expr[, j] <- stats::rpois(length(all_genes), lam)
  }
  list(expr = log1p(expr), truth = truth)
}

#' Abundance-imbalanced data: oligodendrocyte markers are constitutively
#' abundant in every cell, mimicking ambient RNA / transcript spillover.
make_spillover_expr <- function(n_cells = 120, n_bg = 300, seed = 7) {
  set.seed(seed)
  mk <- midbrain_markers
  bg <- paste0("BG", seq_len(n_bg))
  all_genes <- c(bg, unique(mk$gene))
  truth <- rep(unique(mk$cell_type), length.out = n_cells)

  base <- stats::setNames(rep(0.5, length(all_genes)), all_genes)
  base[bg] <- stats::runif(n_bg, 0.1, 2)
  base[mk$gene[mk$cell_type == "Oligodendrocyte"]] <- 6

  expr <- matrix(0, length(all_genes), n_cells,
                 dimnames = list(all_genes, paste0("cell", seq_len(n_cells))))
  for (j in seq_len(n_cells)) {
    lam <- base
    own <- mk$gene[mk$cell_type == truth[j]]
    lam[own] <- lam[own] * 4
    expr[, j] <- stats::rpois(length(all_genes), lam)
  }
  list(expr = log1p(expr), truth = truth)
}

all_methods <- c("zscore", "ucell", "control", "mean_weighted")

#' Cells that retain lineage identity but have lost their effector
#' programme -- the pattern the identity/effector split exists to capture.
make_lineage_shift_expr <- function(n_cells = 240, n_bg = 800, seed = 11) {
  set.seed(seed)
  mk <- midbrain_markers
  bg <- paste0("BG", seq_len(n_bg))
  all_genes <- c(bg, unique(mk$gene))
  truth <- rep(unique(mk$cell_type[mk$layer == "identity"]), length.out = n_cells)

  da <- which(truth == "Dopaminergic neuron")
  affected <- da[seq(1, length(da), by = 2)]

  base <- stats::setNames(rep(0.05, length(all_genes)), all_genes)
  base[bg] <- stats::runif(n_bg, 0.1, 2)

  expr <- matrix(0, length(all_genes), n_cells,
                 dimnames = list(all_genes, paste0("cell", seq_len(n_cells))))
  for (j in seq_len(n_cells)) {
    lam <- base
    own <- mk[mk$cell_type == truth[j], ]
    lam[own$gene[own$layer %in% c("identity", "effector")]] <- 8
    if (j %in% affected) {
      lam[own$gene[own$layer == "effector"]] <- 0.05
      lam[mk$gene[mk$layer == "state" & mk$context == "reactive"]] <- 6
    }
    expr[, j] <- stats::rpois(length(all_genes), lam)
  }
  list(expr = log1p(expr), truth = truth, affected = affected)
}
