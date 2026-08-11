test_that("every method returns a correctly shaped score matrix", {
  d <- make_expr()
  for (m in all_methods) {
    s <- score_markers(d$expr, method = m, seed = 42)
    expect_equal(dim(s), c(length(unique(midbrain_markers$cell_type)), ncol(d$expr)),
                 info = m)
    expect_false(anyNA(s), info = m)
  }
})

test_that("every method recovers cell identity on well-separated data", {
  d <- make_expr()
  for (m in all_methods) {
    lab <- annotate_midbrain_cells(d$expr, method = m, seed = 42)
    expect_gt(mean(lab$cell_type == d$truth), 0.95)
  }
})

test_that("zscore tolerates marker abundance imbalance better than raw means", {
  d <- make_spillover_expr()
  acc <- function(m) {
    mean(annotate_midbrain_cells(d$expr, method = m, seed = 42)$cell_type == d$truth)
  }
  expect_gt(acc("zscore"), acc("mean_weighted"))
})

test_that("ucell scores are bounded in [0, 1]", {
  d <- make_expr()
  s <- score_markers(d$expr, method = "ucell")
  expect_true(all(s >= 0 & s <= 1))
})

test_that("ucell scores do not depend on which other cells are present", {
  d <- make_expr()
  full <- score_markers(d$expr, method = "ucell")
  subset <- score_markers(d$expr[, 1:10, drop = FALSE], method = "ucell")
  expect_equal(unname(full[, 1:10]), unname(subset), tolerance = 1e-8)
})

test_that("zscore scores DO depend on the other cells present", {
  # Documents the known trade-off, so a future refactor can't silently
  # change it without a test failing.
  d <- make_expr()
  full <- score_markers(d$expr, method = "zscore")
  subset <- score_markers(d$expr[, 1:10, drop = FALSE], method = "zscore")
  expect_false(isTRUE(all.equal(unname(full[, 1:10]), unname(subset))))
})

test_that("control method warns and degrades gracefully without a background pool", {
  # A matrix containing nothing but the identity markers themselves leaves
  # no genes to draw expression-matched controls from.
  genes <- unique(midbrain_markers$gene[midbrain_markers$layer == "identity"])
  set.seed(3)
  expr <- matrix(stats::rpois(length(genes) * 20, 4), nrow = length(genes),
                 dimnames = list(genes, paste0("cell", 1:20)))
  expr <- log1p(expr)

  expect_warning(score_markers(expr, method = "control", seed = 42), "control pool")
})

test_that("seed makes control scores reproducible", {
  d <- make_expr()
  a <- score_markers(d$expr, method = "control", seed = 99)
  b <- score_markers(d$expr, method = "control", seed = 99)
  expect_equal(a, b)
})

test_that("scoring does not disturb the caller's RNG stream", {
  d <- make_expr()
  set.seed(123)
  before <- stats::runif(1)
  set.seed(123)
  invisible(score_markers(d$expr, method = "control", seed = 7))
  after <- stats::runif(1)
  expect_equal(before, after)
})

test_that("zscore handles zero-variance genes without producing NaN", {
  d <- make_expr()
  d$expr["GFAP", ] <- 1  # constant across all cells
  s <- score_markers(d$expr, method = "zscore")
  expect_false(anyNA(s))
})

test_that("duplicated gene symbols are rejected", {
  d <- make_expr()
  e <- rbind(d$expr, d$expr["TH", , drop = FALSE])
  rownames(e)[nrow(e)] <- "TH"
  expect_error(score_markers(e), "duplicated")
})

test_that("min_margin flags ambiguous cells", {
  d <- make_expr()
  lab <- annotate_midbrain_cells(d$expr, method = "zscore", min_margin = Inf)
  expect_true(all(lab$cell_type == "Ambiguous"))
})
