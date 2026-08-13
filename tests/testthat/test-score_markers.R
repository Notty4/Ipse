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
    lab <- annotate_cells(d$expr, method = m, seed = 42)
    expect_gt(mean(lab$cell_type == d$truth), 0.95)
  }
})

test_that("zscore tolerates marker abundance imbalance better than raw means", {
  d <- make_spillover_expr()
  acc <- function(m) {
    mean(annotate_cells(d$expr, method = m, seed = 42)$cell_type == d$truth)
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
  lab <- annotate_cells(d$expr, method = "zscore", min_margin = Inf)
  expect_true(all(lab$cell_type == "Ambiguous"))
})

test_that("scoring keeps the expression matrix sparse", {
  # ucell used to call as.matrix() on the whole matrix, which is tens of
  # GB on a real dataset. Assert the mechanism rather than a memory
  # figure: gc() accounting is not portable, but a densified matrix
  # provably is not sparse.
  set.seed(4)
  g <- 3000; n <- 800
  m <- methods::as(Matrix::rsparsematrix(g, n, density = 0.03,
                                         rand.x = function(k) stats::rpois(k, 3) + 1),
                   "CsparseMatrix")
  rownames(m) <- c(unique(midbrain_markers$gene),
                   paste0("BG", seq_len(g - length(unique(midbrain_markers$gene)))))
  colnames(m) <- paste0("c", seq_len(n))

  expect_true(methods::is(.rank_cells(m, 1500), "sparseMatrix"))

  for (meth in c("zscore", "ucell", "control", "mean_weighted")) {
    s <- score_markers(m, method = meth, seed = 42)
    expect_false(anyNA(s), info = meth)
    # the input must be untouched -- no method may densify it in place
    expect_true(methods::is(m, "sparseMatrix"), info = meth)
  }
})

test_that("ucell ranking matches a dense reference implementation", {
  set.seed(1)
  g <- 400; n <- 60
  m <- methods::as(Matrix::rsparsematrix(g, n, density = 0.25,
                                         rand.x = function(k) stats::rpois(k, 5) + 1),
                   "CsparseMatrix")
  rownames(m) <- c(unique(midbrain_markers$gene),
                   paste0("BG", seq_len(g - length(unique(midbrain_markers$gene)))))
  colnames(m) <- paste0("c", seq_len(n))

  max_rank <- 150
  dense_ranks <- apply(as.matrix(m), 2, function(x) {
    r <- rank(-x, ties.method = "average")
    r[r > max_rank] <- max_rank + 1
    r
  })
  rownames(dense_ranks) <- rownames(m)

  genes <- intersect(midbrain_markers$gene[midbrain_markers$cell_type == "Astrocyte"],
                     rownames(m))
  nn <- length(genes)
  expected <- 1 - (colSums(dense_ranks[genes, , drop = FALSE]) - nn * (nn + 1) / 2) /
    (nn * max_rank)

  actual <- .ucell_score(.rank_cells(m, max_rank), genes, max_rank)
  expect_equal(unname(actual), unname(expected), tolerance = 1e-10)
})

test_that("gene moments are chunked without changing the result", {
  # rowMeans(expr * expr) would materialise a whole second matrix; on real
  # data that is several GB and fails. Chunking must be exact, not
  # approximate, and independent of chunk size.
  set.seed(2)
  g <- 600; n <- 400
  m <- methods::as(Matrix::rsparsematrix(g, n, density = 0.2,
                                         rand.x = function(k) stats::rpois(k, 4) + 1),
                   "CsparseMatrix")
  rownames(m) <- paste0("g", seq_len(g))
  colnames(m) <- paste0("c", seq_len(n))

  mu_ref <- Matrix::rowMeans(m)
  sq_ref <- Matrix::rowMeans(m * m)
  sd_ref <- sqrt(pmax(sq_ref - mu_ref^2, 0)) * sqrt(n / (n - 1))
  zv <- sd_ref < .Machine$double.eps
  sd_ref[zv] <- 1; mu_ref[zv] <- 0

  got <- .gene_moments(m, chunk = 50)
  expect_equal(unname(got$mu), unname(mu_ref))
  expect_equal(unname(got$sdev), unname(sd_ref))
  expect_equal(.gene_moments(m, chunk = 37), .gene_moments(m, chunk = 9999))
})

test_that("scoring works without Matrix attached to the search path", {
  # Unqualified calls to Matrix generics (crossprod, rowSums...) resolve
  # to the base versions when the package is installed rather than
  # sourced, and base::crossprod cannot handle a sparse matrix. This only
  # shows up under R CMD check, never in an interactive session where the
  # user has run library(Matrix). Every such call must be Matrix:: qualified.
  set.seed(6)
  g <- 500; n <- 200
  m <- methods::as(Matrix::rsparsematrix(g, n, density = 0.1,
                                         rand.x = function(k) stats::rpois(k, 3) + 1),
                   "CsparseMatrix")
  rownames(m) <- c(unique(midbrain_markers$gene),
                   paste0("BG", seq_len(g - length(unique(midbrain_markers$gene)))))
  colnames(m) <- paste0("c", seq_len(n))

  # Evaluate in an environment where the Matrix package is not attached,
  # mimicking how the installed package sees the world.
  was_attached <- "package:Matrix" %in% search()
  if (was_attached) {
    detach("package:Matrix", unload = FALSE, character.only = TRUE)
    on.exit(suppressMessages(library(Matrix)), add = TRUE)
  }

  for (meth in c("zscore", "ucell", "control", "mean_weighted")) {
    expect_silent(s <- score_markers(m, method = meth, seed = 42))
    expect_false(anyNA(s), info = meth)
  }
})
