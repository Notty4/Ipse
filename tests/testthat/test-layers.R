test_that("marker tables without layer/context columns still work", {
  legacy <- data.frame(
    cell_type = c("A", "A", "B", "B"),
    gene = c("TH", "NR4A2", "GFAP", "AQP4"),
    weight = 1,
    stringsAsFactors = FALSE
  )
  d <- make_expr()
  s <- score_markers(d$expr, markers = legacy, method = "zscore")
  expect_equal(sort(rownames(s)), c("A", "B"))
})

test_that("invalid layer names are rejected", {
  bad <- midbrain_markers
  bad$layer[1] <- "nonsense"
  expect_error(score_markers(make_expr()$expr, markers = bad), "Unknown marker layer")
})

test_that("requesting an absent layer errors rather than silently returning nothing", {
  identity_only <- midbrain_markers[midbrain_markers$layer == "identity", ]
  expect_error(
    score_markers(make_expr()$expr, markers = identity_only, layer = "state"),
    "No markers found in layer"
  )
})

test_that("labels are decided by identity markers alone", {
  d <- make_expr()
  # Wildly inflating an effector marker of another lineage must not change
  # any label, since effector markers never contribute to assignment.
  altered <- d$expr
  altered["MBP", ] <- altered["MBP", ] + 50

  a <- annotate_cells(d$expr, seed = 42)
  b <- annotate_cells(altered, seed = 42)
  expect_equal(a$cell_type, b$cell_type)
})

test_that("state markers never contribute to the assigned label", {
  d <- make_expr()
  altered <- d$expr
  altered[c("GFAP", "VIM"), ] <- altered[c("GFAP", "VIM"), ] + 50

  a <- annotate_cells(d$expr, seed = 42)
  b <- annotate_cells(altered, seed = 42)
  expect_equal(a$cell_type, b$cell_type)
})

test_that("effector score is returned by default", {
  d <- make_expr()
  labels <- annotate_cells(d$expr, seed = 42)
  expect_true("effector_score" %in% names(labels))
})

test_that("score_effector = FALSE omits the effector column", {
  d <- make_expr()
  labels <- annotate_cells(d$expr, score_effector = FALSE, seed = 42)
  expect_false("effector_score" %in% names(labels))
})

test_that("lost effector programme is detected while identity is retained", {
  d <- make_lineage_shift_expr()
  labels <- annotate_cells(d$expr, seed = 42)

  # affected cells keep their dopaminergic label ...
  expect_true(all(labels$cell_type[d$affected] == "Dopaminergic neuron"))
  # ... but score markedly lower on the effector layer
  intact <- setdiff(which(d$truth == "Dopaminergic neuron"), d$affected)
  expect_lt(mean(labels$effector_score[d$affected]),
            mean(labels$effector_score[intact]))
})

test_that("effector score is NA for cell types lacking effector markers", {
  d <- make_expr()
  labels <- annotate_cells(d$expr, seed = 42)
  # OPC has identity markers only in the bundled reference
  opc <- labels$cell_type == "OPC"
  if (any(opc)) expect_true(all(is.na(labels$effector_score[opc])))
})

test_that("a shared 'Any' state signature applies to every cell type", {
  d <- make_expr()
  shared <- rbind(
    midbrain_markers,
    data.frame(cell_type = "Any", gene = c("BG1", "BG2"), layer = "state",
               context = "shared", source = "test fixture", weight = 1,
               stringsAsFactors = FALSE)
  )
  labels <- annotate_cells(d$expr, markers = shared,
                                    sections = "shared", seed = 42)
  expect_true("state_shared" %in% names(labels))
  expect_false(anyNA(labels$state_shared))
})

test_that("return_scores exposes identity, effector and state matrices", {
  d <- make_expr()
  out <- annotate_cells(d$expr, seed = 42, sections = "all",
                                 return_scores = TRUE)
  expect_true(all(c("labels", "identity", "effector", "state_reactive") %in% names(out)))
  expect_equal(ncol(out$identity), ncol(d$expr))
})

test_that("list_cell_types reflects the identity layer", {
  identity_types <- unique(
    midbrain_markers$cell_type[midbrain_markers$layer == "identity"]
  )
  expect_setequal(list_cell_types(), identity_types)
})
