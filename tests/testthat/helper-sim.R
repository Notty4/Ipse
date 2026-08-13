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

#' A small annotated reference with three kinds of planted gene:
#' exclusive (should win), shared between two types (fold change likes
#' these), and abundant everywhere but higher in one (the spillover trap).
make_reference <- function(n_per = 60, seed = 5) {
  set.seed(seed)
  types <- c("DA", "Astro", "Oligo", "Micro")
  labels <- rep(types, each = n_per)
  donors <- rep(paste0("D", 1:4), length.out = length(labels))

  genes <- c(paste0("TRUE_", types), paste0("SHARED_", types),
             paste0("ABUND_", types), paste0("BG", 1:60))
  lam <- matrix(0.3, length(genes), length(labels), dimnames = list(genes, NULL))

  for (i in seq_along(types)) {
    t <- types[i]
    cols <- labels == t
    lam[paste0("TRUE_", t), cols] <- 12
    lam[paste0("SHARED_", t), cols] <- 12
    lam[paste0("SHARED_", t), labels == types[(i %% 4) + 1]] <- 12
    lam[paste0("ABUND_", t), ] <- 20
    lam[paste0("ABUND_", t), cols] <- 34
  }

  expr <- log1p(matrix(stats::rpois(length(lam), lam), nrow = nrow(lam),
                       dimnames = list(genes, paste0("c", seq_along(labels)))))
  list(expr = expr, labels = labels, donors = donors, types = types)
}

#' A dataset with a parent class that splits into subtypes, plus two
#' unrelated lineages. Returns matching broad and subtype marker tables.
make_subtype_expr <- function(n_per = 120, seed = 21) {
  set.seed(seed)
  subs <- c("SubA", "SubB", "SubC")
  others <- c("Glia1", "Glia2")
  truth_sub <- c(rep(subs, each = n_per), rep(NA, length(others) * n_per))
  truth_broad <- c(rep("Parent", length(subs) * n_per),
                   rep(others, each = n_per))
  n <- length(truth_broad)

  parent_mk <- paste0("PARENT_", 1:6)
  sub_mk <- unlist(lapply(subs, function(s) paste0(s, "_", 1:5)))
  other_mk <- unlist(lapply(others, function(o) paste0(o, "_", 1:6)))
  bg <- paste0("BG", 1:300)
  all_genes <- c(parent_mk, sub_mk, other_mk, bg)

  lam <- matrix(0.2, length(all_genes), n, dimnames = list(all_genes, paste0("c", 1:n)))
  lam[bg, ] <- stats::runif(length(bg), 0.1, 2)
  lam[parent_mk, truth_broad == "Parent"] <- 10
  for (o in others) lam[paste0(o, "_", 1:6), truth_broad == o] <- 10
  for (s in subs) lam[paste0(s, "_", 1:5), which(truth_sub == s)] <- 8

  expr <- log1p(matrix(stats::rpois(length(lam), lam), nrow = nrow(lam),
                       dimnames = dimnames(lam)))

  mk_broad <- data.frame(
    cell_type = c(rep("Parent", 6), rep(others, each = 6)),
    gene = c(parent_mk, other_mk),
    layer = "identity", context = "canonical",
    source = NA_character_, weight = 1, stringsAsFactors = FALSE)

  mk_sub <- data.frame(
    cell_type = rep(subs, each = 5), gene = sub_mk,
    layer = "identity", context = "canonical",
    source = NA_character_, weight = 1, stringsAsFactors = FALSE)

  list(expr = expr, truth_broad = truth_broad, truth_sub = truth_sub,
       markers = mk_broad, subtype_markers = mk_sub)
}
