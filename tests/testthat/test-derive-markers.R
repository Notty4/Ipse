test_that("derive_markers returns one row per candidate with expected columns", {
  r <- make_reference()
  cand <- derive_markers(r$expr, r$labels)

  expect_true(all(c("cell_type", "gene", "mean_in", "detect_in", "mean_out_max",
                    "competitor", "specificity", "detect_gap",
                    "donor_consistency", "rank") %in% names(cand)))
  expect_setequal(unique(cand$cell_type), r$types)
})

test_that("the exclusive marker ranks first for each cell type", {
  r <- make_reference()
  cand <- derive_markers(r$expr, r$labels)

  for (t in r$types) {
    top <- cand$gene[cand$cell_type == t & cand$rank == 1]
    expect_equal(top, paste0("TRUE_", t))
  }
})

test_that("genes shared between two types are penalised, unlike fold change", {
  r <- make_reference()
  cand <- derive_markers(r$expr, r$labels)

  # A gene shared by DA and Astro should not be a top-3 candidate for either
  shared <- cand[cand$cell_type == "DA" & cand$gene == "SHARED_DA", ]
  expect_gt(shared$rank, 3)
  # and its specificity should be near zero rather than large
  expect_lt(abs(shared$specificity), 0.5)
})

test_that("detect_gap rejects the abundant-everywhere trap", {
  r <- make_reference()
  cand <- derive_markers(r$expr, r$labels, rank_by = "detect_gap")

  abund <- cand[cand$cell_type == "DA" & cand$gene == "ABUND_DA", ]
  true_marker <- cand[cand$cell_type == "DA" & cand$gene == "TRUE_DA", ]

  expect_lt(abund$detect_gap, 0.05)
  expect_gt(true_marker$detect_gap, 0.5)
  expect_gt(abund$rank, true_marker$rank)
})

test_that("competitor names the cell type that most contests the gene", {
  r <- make_reference()
  cand <- derive_markers(r$expr, r$labels)
  # SHARED_DA is planted high in DA and Astro
  row <- cand[cand$cell_type == "DA" & cand$gene == "SHARED_DA", ]
  expect_equal(row$competitor, "Astro")
})

test_that("donor_consistency is high for true markers and NA without donors", {
  r <- make_reference()
  with_donor <- derive_markers(r$expr, r$labels, donor = r$donors)
  no_donor <- derive_markers(r$expr, r$labels)

  true_da <- with_donor[with_donor$cell_type == "DA" & with_donor$gene == "TRUE_DA", ]
  expect_equal(true_da$donor_consistency, 1)
  expect_true(all(is.na(no_donor$donor_consistency)))
})

test_that("panel restriction limits candidates to measurable genes", {
  r <- make_reference()
  panel <- c("TRUE_DA", "TRUE_Astro", "TRUE_Oligo", "TRUE_Micro", "BG1", "BG2")
  cand <- derive_markers(r$expr, r$labels, panel = panel)

  expect_true(all(cand$gene %in% panel))
})

test_that("genes absent from expr are reported when restricting to a panel", {
  r <- make_reference()
  expect_warning(
    derive_markers(r$expr, r$labels, panel = c("TRUE_DA", "TRUE_Astro", "NOT_A_GENE")),
    "not found"
  )
})

test_that("mismatched labels are rejected", {
  r <- make_reference()
  expect_error(derive_markers(r$expr, r$labels[-1]), "one entry per column")
})

test_that("cell types below min_cells are dropped with a warning", {
  r <- make_reference()
  labels <- r$labels
  labels[labels == "Micro"][1:55] <- "DA"   # leave Micro very small
  expect_warning(derive_markers(r$expr, labels, min_cells = 20), "fewer than")
})

test_that("consensus_markers keeps only genes replicating across cohorts", {
  r1 <- make_reference(seed = 1)
  r2 <- make_reference(seed = 2)
  c1 <- derive_markers(r1$expr, r1$labels)
  c2 <- derive_markers(r2$expr, r2$labels)

  cons <- consensus_markers(list(a = c1, b = c2), top_n = 5, min_cohorts = 2)

  expect_true(all(cons$n_cohorts == 2))
  # the planted exclusive markers should survive replication
  expect_true("TRUE_DA" %in% cons$gene[cons$cell_type == "DA"])
})

test_that("consensus_markers reports the worst rank achieved, not the best", {
  r1 <- make_reference(seed = 1)
  c1 <- derive_markers(r1$expr, r1$labels)
  c2 <- c1
  c2$rank[c2$gene == "TRUE_DA" & c2$cell_type == "DA"] <- 9

  cons <- consensus_markers(list(a = c1, b = c2), top_n = 20, min_cohorts = 2)
  expect_equal(cons$worst_rank[cons$gene == "TRUE_DA" & cons$cell_type == "DA"], 9)
})

test_that("consensus_markers rejects a malformed candidate table", {
  r <- make_reference()
  cand <- derive_markers(r$expr, r$labels)
  expect_error(
    consensus_markers(list(a = cand, b = cand[, c("cell_type", "gene")])),
    "missing column"
  )
})

test_that("as_marker_table produces a table annotate_cells accepts", {
  r <- make_reference()
  cand <- derive_markers(r$expr, r$labels)
  mk <- as_marker_table(cand, n = 3, min_gap = 0.3, source = "test fixture")

  expect_true(all(c("cell_type", "gene", "layer", "context", "source", "weight")
                  %in% names(mk)))
  expect_true(all(table(mk$cell_type) <= 3))

  labels <- annotate_cells(r$expr, markers = mk, seed = 42)
  expect_equal(nrow(labels), ncol(r$expr))
  # the derived markers should recover the cell types they were derived from
  expect_gt(mean(labels$cell_type == r$labels), 0.9)
})

test_that("as_marker_table assigns only the named genes to the effector layer", {
  r <- make_reference()
  cand <- derive_markers(r$expr, r$labels)
  mk <- as_marker_table(cand, n = 3, min_gap = NULL, effector = "TRUE_DA")

  expect_equal(mk$layer[mk$gene == "TRUE_DA"], "effector")
  expect_true(all(mk$layer[mk$gene != "TRUE_DA"] == "identity"))
})

test_that("as_marker_table warns if a cell type is left with no identity markers", {
  r <- make_reference()
  cand <- derive_markers(r$expr, r$labels)
  top_da <- cand$gene[cand$cell_type == "DA"][1:3]

  expect_warning(as_marker_table(cand, n = 3, min_gap = NULL, effector = top_da),
                 "cannot be labelled")
})

test_that("min_gap keeps noise out of the marker table", {
  r <- make_reference()
  cand <- derive_markers(r$expr, r$labels)

  padded <- as_marker_table(cand, n = 3, min_gap = NULL)
  filtered <- as_marker_table(cand, n = 3, min_gap = 0.3)

  expect_lt(nrow(filtered), nrow(padded))
  # only the genuinely exclusive genes clear a realistic floor
  expect_true(all(grepl("^TRUE_", filtered$gene)))
})

test_that("min_gap warns when some cell types lose all their candidates", {
  r <- make_reference()
  cand <- derive_markers(r$expr, r$labels)
  # Cripple one cell type's candidates so it alone falls below the floor
  cand$detect_gap[cand$cell_type == "Micro"] <- 0.01

  expect_warning(as_marker_table(cand, n = 3, min_gap = 0.3), "Micro")
})

test_that("min_gap errors rather than returning an empty table", {
  r <- make_reference()
  cand <- derive_markers(r$expr, r$labels)
  expect_error(
    suppressWarnings(as_marker_table(cand, n = 3, min_gap = 0.99)),
    "No candidates passed"
  )
})

test_that("depth correction stops a high-RNA cell type manufacturing markers", {
  # A cell type with no genuine markers, but 3x the RNA content of the
  # others, will otherwise have background genes promoted as markers --
  # which then attract cells away from every other type at annotation time.
  set.seed(9)
  types <- c("Deep", "A", "B", "C")
  n <- 1200
  truth <- rep(types, each = n / 4)

  mk <- c(paste0("A_", 1:5), paste0("B_", 1:5), paste0("C_", 1:5))
  bg <- paste0("BG", 1:200)
  all_genes <- c(mk, bg)

  lam <- matrix(0.3, length(all_genes), n, dimnames = list(all_genes, paste0("c", 1:n)))
  for (t in c("A", "B", "C")) lam[paste0(t, "_", 1:5), truth == t] <- 10
  lam <- sweep(lam, 2, ifelse(truth == "Deep", 3.0, 1.0), "*")
  expr <- log1p(matrix(stats::rpois(length(lam), lam), nrow = nrow(lam),
                       dimnames = dimnames(lam)))

  raw <- derive_markers(expr, truth, depth_correct = FALSE)
  fixed <- derive_markers(expr, truth, depth_correct = TRUE)

  spurious_raw <- max(raw$detect_gap[raw$cell_type == "Deep"])
  spurious_fixed <- max(fixed$detect_gap[fixed$cell_type == "Deep"])
  genuine_fixed <- max(fixed$detect_gap[fixed$gene %in% paste0("A_", 1:5) &
                                          fixed$cell_type == "A"])

  expect_lt(spurious_fixed, spurious_raw / 2)
  # and a real marker must remain clearly above the spurious ceiling
  expect_gt(genuine_fixed, spurious_fixed * 3)
})

test_that("depth correction leaves genuine markers ranked first", {
  r <- make_reference()
  cand <- derive_markers(r$expr, r$labels, depth_correct = TRUE)
  for (t in r$types) {
    expect_equal(cand$gene[cand$cell_type == t & cand$rank == 1], paste0("TRUE_", t))
  }
})
