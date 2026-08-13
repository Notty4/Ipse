test_that("subtypes are resolved within the parent class only", {
  d <- make_subtype_expr()
  lab <- annotate_cells(d$expr, markers = d$markers, seed = 42)
  out <- annotate_subtypes(d$expr, lab, parent = "Parent",
                           markers = d$subtype_markers, seed = 42)

  expect_true(all(c("subtype", "subtype_score", "subtype_margin") %in% names(out)))
  # cells outside the parent class are untouched
  expect_true(all(is.na(out$subtype[out$cell_type != "Parent"])))
  expect_false(anyNA(out$subtype[out$cell_type == "Parent"]))
})

test_that("known subtypes are recovered", {
  d <- make_subtype_expr()
  lab <- annotate_cells(d$expr, markers = d$markers, seed = 42)
  out <- annotate_subtypes(d$expr, lab, parent = "Parent",
                           markers = d$subtype_markers, seed = 42)

  is_parent <- out$cell_type == "Parent"
  expect_gt(mean(out$subtype[is_parent] == d$truth_sub[is_parent], na.rm = TRUE), 0.9)
})

test_that("scoring within the parent beats scoring across everything", {
  # The reason this is a separate pass: subtype markers compared against
  # glia rather than sibling subtypes lose their discriminating power.
  d <- make_subtype_expr()
  lab <- annotate_cells(d$expr, markers = d$markers, seed = 42)
  is_parent <- lab$cell_type == "Parent"

  within <- annotate_subtypes(d$expr, lab, parent = "Parent",
                              markers = d$subtype_markers, seed = 42)
  acc_within <- mean(within$subtype[is_parent] == d$truth_sub[is_parent], na.rm = TRUE)

  # naive alternative: score every cell at once, then keep the parent ones
  all_scores <- score_markers(d$expr, markers = d$subtype_markers, seed = 42)
  naive <- rownames(all_scores)[apply(all_scores, 2, which.max)]
  acc_naive <- mean(naive[is_parent] == d$truth_sub[is_parent], na.rm = TRUE)

  expect_gte(acc_within, acc_naive)
})

test_that("min_margin flags ambiguous subtype calls", {
  d <- make_subtype_expr()
  lab <- annotate_cells(d$expr, markers = d$markers, seed = 42)
  out <- annotate_subtypes(d$expr, lab, parent = "Parent",
                           markers = d$subtype_markers, min_margin = Inf, seed = 42)

  expect_true(all(out$subtype[out$cell_type == "Parent"] == "Ambiguous"))
})

test_that("min_score flags weak subtype calls", {
  d <- make_subtype_expr()
  lab <- annotate_cells(d$expr, markers = d$markers, seed = 42)
  out <- annotate_subtypes(d$expr, lab, parent = "Parent",
                           markers = d$subtype_markers, min_score = Inf, seed = 42)

  expect_true(all(out$subtype[out$cell_type == "Parent"] == "Unassigned"))
})

test_that("an absent parent class warns and leaves labels unchanged", {
  d <- make_subtype_expr()
  lab <- annotate_cells(d$expr, markers = d$markers, seed = 42)

  expect_warning(
    out <- annotate_subtypes(d$expr, lab, parent = "NOT_A_TYPE",
                             markers = d$subtype_markers, seed = 42),
    "No cells labelled"
  )
  expect_true(all(is.na(out$subtype)))
  expect_equal(out$cell_type, lab$cell_type)
})

test_that("mismatched labels are rejected", {
  d <- make_subtype_expr()
  lab <- annotate_cells(d$expr, markers = d$markers, seed = 42)

  expect_error(annotate_subtypes(d$expr, lab[-1, ], parent = "Parent",
                                 markers = d$subtype_markers),
               "one row per column")
  expect_error(annotate_subtypes(d$expr, data.frame(x = 1), parent = "Parent",
                                 markers = d$subtype_markers),
               "annotate_cells")
})

test_that("the broad labels are never altered by sub-labelling", {
  d <- make_subtype_expr()
  lab <- annotate_cells(d$expr, markers = d$markers, seed = 42)
  out <- annotate_subtypes(d$expr, lab, parent = "Parent",
                           markers = d$subtype_markers, seed = 42)

  expect_equal(out$cell_type, lab$cell_type)
  expect_equal(out$score, lab$score)
})
