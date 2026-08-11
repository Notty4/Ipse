test_that("sections are opt-in: none run by default", {
  d <- make_expr()
  labels <- annotate_midbrain_cells(d$expr, seed = 42)

  expect_false(any(grepl("^state_", names(labels))))
  expect_null(attr(labels, "sections"))
})

test_that("a named section is scored when requested", {
  d <- make_expr()
  labels <- annotate_midbrain_cells(d$expr, sections = "reactive", seed = 42)

  expect_true("state_reactive" %in% names(labels))
})

test_that("sections = 'all' runs every section in the reference", {
  d <- make_expr()
  labels <- annotate_midbrain_cells(d$expr, sections = "all", seed = 42)
  available <- unique(midbrain_markers$context[midbrain_markers$layer == "state"])

  expect_true(all(paste0("state_", available) %in% names(labels)))
})

test_that("an unknown section name errors and lists what is available", {
  d <- make_expr()
  expect_error(
    annotate_midbrain_cells(d$expr, sections = "NOT_A_SECTION", seed = 42),
    "Unknown section"
  )
  expect_error(
    annotate_midbrain_cells(d$expr, sections = "NOT_A_SECTION", seed = 42),
    "reactive"   # names the available section
  )
})

test_that("provenance travels with the result", {
  d <- make_expr()
  labels <- annotate_midbrain_cells(d$expr, sections = "reactive", seed = 42)
  prov <- attr(labels, "sections")

  expect_s3_class(prov, "data.frame")
  expect_equal(prov$context, "reactive")
  expect_false(is.na(prov$source))
})

test_that("running a section never changes the assigned labels", {
  d <- make_expr()
  without <- annotate_midbrain_cells(d$expr, seed = 42)
  with    <- annotate_midbrain_cells(d$expr, sections = "all", seed = 42)

  expect_equal(without$cell_type, with$cell_type)
  expect_equal(without$score, with$score)
})

test_that("list_state_sections reports contexts, cell types, size and source", {
  s <- list_state_sections()
  expect_true(all(c("context", "cell_types", "n_genes", "source") %in% names(s)))
  expect_true("reactive" %in% s$context)
  expect_gt(s$n_genes[s$context == "reactive"], 0)
})

test_that("list_state_sections returns zero rows when none are defined", {
  identity_only <- midbrain_markers[midbrain_markers$layer != "state", ]
  s <- list_state_sections(identity_only)

  expect_equal(nrow(s), 0)
  expect_true(all(c("context", "source") %in% names(s)))
})

test_that("requesting a section from a reference that has none errors clearly", {
  d <- make_expr()
  identity_only <- midbrain_markers[midbrain_markers$layer != "state", ]
  expect_error(
    annotate_midbrain_cells(d$expr, markers = identity_only,
                            sections = "reactive", seed = 42),
    "none in this reference"
  )
})

test_that("a user-supplied section is picked up without touching the package", {
  d <- make_expr()
  custom <- rbind(
    midbrain_markers,
    data.frame(cell_type = "Dopaminergic neuron", gene = c("BG1", "BG2"),
               layer = "state", context = "MY_COHORT",
               source = "Someone et al. 2026, GSE999999", weight = 1,
               stringsAsFactors = FALSE)
  )
  labels <- annotate_midbrain_cells(d$expr, markers = custom,
                                    sections = "MY_COHORT", seed = 42)

  expect_true("state_MY_COHORT" %in% names(labels))
  expect_equal(attr(labels, "sections")$source, "Someone et al. 2026, GSE999999")
})
