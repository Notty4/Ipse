make_imbalanced <- function(n = 2000, seed = 3) {
  set.seed(seed)
  types <- c("Oligo", "Astro", "Micro", "Endo")
  props <- c(0.70, 0.14, 0.14, 0.02)
  truth <- sample(types, n, replace = TRUE, prob = props)

  bg <- paste0("BG", 1:400)
  mk <- unlist(lapply(types, function(t) paste0(t, "_", 1:5)))
  all_genes <- c(bg, mk)

  lam <- matrix(0.2, length(all_genes), n, dimnames = list(all_genes, paste0("c", 1:n)))
  lam[bg, ] <- stats::runif(length(bg), 0.1, 2)
  for (t in types) {
    g <- paste0(t, "_", 1:5)
    lam[g, truth == t] <- 10
    if (t == "Oligo") lam[g, truth != t] <- 3.5   # ambient bleed
  }
  expr <- log1p(matrix(stats::rpois(length(lam), lam), nrow = nrow(lam),
                       dimnames = dimnames(lam)))

  markers <- data.frame(
    cell_type = rep(types, each = 5), gene = mk, layer = "identity",
    context = "canonical", source = NA_character_, weight = 1,
    stringsAsFactors = FALSE)

  list(expr = expr, truth = truth, markers = markers)
}

test_that("benchmark_methods returns one row per method with expected columns", {
  d <- make_imbalanced()
  b <- benchmark_methods(d$expr, d$truth, markers = d$markers)

  expect_equal(nrow(b), 4)
  expect_true(all(c("method", "accuracy", "balanced_accuracy",
                    "worst_class", "worst_recall") %in% names(b)))
  expect_true(all(b$accuracy >= 0 & b$accuracy <= 1))
})

test_that("results are sorted by balanced accuracy", {
  d <- make_imbalanced()
  b <- benchmark_methods(d$expr, d$truth, markers = d$markers)
  expect_equal(b$balanced_accuracy, sort(b$balanced_accuracy, decreasing = TRUE))
})

test_that("zscore is penalised under class imbalance with ambient bleed", {
  # Documents the failure mode: the abundant class leaks into the rare one
  # because z-scoring compresses the scores of a class that dominates the
  # dataset. If a future change fixes this, the test should be updated,
  # not deleted.
  d <- make_imbalanced()
  b <- benchmark_methods(d$expr, d$truth, markers = d$markers)

  z <- b$accuracy[b$method == "zscore"]
  u <- b$accuracy[b$method == "ucell"]
  expect_lt(z, u)
})

test_that("mismatched labels are rejected", {
  d <- make_imbalanced()
  expect_error(benchmark_methods(d$expr, d$truth[-1], markers = d$markers),
               "one entry per column")
})

test_that("a restricted method set is honoured", {
  d <- make_imbalanced()
  b <- benchmark_methods(d$expr, d$truth, markers = d$markers,
                         methods = c("ucell", "control"))
  expect_setequal(b$method, c("ucell", "control"))
})
