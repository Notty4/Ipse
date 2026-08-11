test_that("annotate_midbrain_cells assigns one row per cell", {
  d <- make_expr()
  labels <- annotate_midbrain_cells(d$expr, seed = 42)

  expect_equal(nrow(labels), ncol(d$expr))
  expect_true(all(c("cell_id", "cell_type", "score", "margin") %in% colnames(labels)))
  expect_equal(labels$cell_id, colnames(d$expr))
  expect_true(all(labels$cell_type %in%
                    c(unique(midbrain_markers$cell_type), "Unassigned", "Ambiguous")))
})

test_that("known cell identities are recovered on well-separated data", {
  d <- make_expr()
  labels <- annotate_midbrain_cells(d$expr, seed = 42)
  expect_gt(mean(labels$cell_type == d$truth), 0.95)
})

test_that("min_score labels low-confidence cells as Unassigned", {
  d <- make_expr()
  labels <- annotate_midbrain_cells(d$expr, min_score = Inf, seed = 42)
  expect_true(all(labels$cell_type == "Unassigned"))
})

test_that("return_scores returns both labels and the score matrix", {
  d <- make_expr()
  out <- annotate_midbrain_cells(d$expr, seed = 42, return_scores = TRUE)

  expect_named(out, c("labels", "scores"))
  expect_equal(nrow(out$labels), ncol(d$expr))
  expect_equal(ncol(out$scores), ncol(d$expr))
  # the reported score must match the winning row of the score matrix
  expect_equal(out$labels$score, unname(apply(out$scores, 2, max)))
})

test_that("margin is the gap between the best and second-best score", {
  d <- make_expr(n_cells = 12)
  out <- annotate_midbrain_cells(d$expr, seed = 42, return_scores = TRUE)
  expected <- apply(out$scores, 2, function(col) {
    s <- sort(col, decreasing = TRUE)
    s[1] - s[2]
  })
  expect_equal(out$labels$margin, unname(expected))
})

test_that("score_markers errors informatively on missing row names", {
  expr <- matrix(1:10, nrow = 2)
  expect_error(score_markers(expr), "row names")
})

test_that("score_markers errors informatively on a malformed marker table", {
  d <- make_expr()
  bad_markers <- data.frame(cell_type = "X", gene = "TH")
  expect_error(score_markers(d$expr, markers = bad_markers), "weight")
})

test_that("a cell type with no markers present yields NA with a warning", {
  d <- make_expr()
  drop <- midbrain_markers$gene[midbrain_markers$cell_type == "Microglia"]
  reduced <- d$expr[setdiff(rownames(d$expr), drop), , drop = FALSE]

  expect_warning(s <- score_markers(reduced, method = "zscore"), "Microglia")
  expect_true(all(is.na(s["Microglia", ])))
})

test_that("list_midbrain_cell_types returns the expected labels", {
  expect_setequal(list_midbrain_cell_types(), unique(midbrain_markers$cell_type))
})
