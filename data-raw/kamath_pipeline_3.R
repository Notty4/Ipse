# data-raw/kamath_pipeline.R
#
# Derive markers from Kamath et al. 2022 (Broad SCP1768).
#
# Two things about this dataset shape the script:
#
#  1. METADATA_PD.tsv.gz has NO cell type column. Cell types live in the
#     per-lineage *_UMAP.tsv files, which we concatenate to build the
#     label vector. da_UMAP.tsv carries the DA SUBTYPES (SOX6_AGTR1 etc),
#     which is the resolution Agarwal could not provide.
#
#  2. The metadata covers 1.49M rows across four species, two brain
#     regions and three conditions. Human + substantia nigra + Ctrl is a
#     hard requirement before deriving anything.
#
# The 4.74 GB matrix is never loaded whole. We decide which cells we want
# first, then stream the file once keeping only those columns.
#
#   source("data-raw/kamath_pipeline.R")

library(Matrix)
devtools::load_all()

data_dir <- "C:/Users/User/Downloads/SCP1768"
out_dir  <- "data-raw/derived"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

cache_expr <- file.path(out_dir, "kamath_controls_expr.rds")
cap <- 2000        # cells per type; marker ranks stabilise well below this

# =====================================================================
# STAGE A -- which cells do we want?
# =====================================================================

read_scp <- function(path) {
  x <- read.delim(path, stringsAsFactors = FALSE, check.names = FALSE)
  if (tolower(as.character(x[1, 1])) %in% c("type", "group")) x <- x[-1, , drop = FALSE]
  x
}

## --- cell types, assembled from the per-lineage UMAP files -----------
# Keep ONLY the lineage-specific files. The nurr files must be excluded:
# nurrneg_liger_UMAP.tsv labels every NURR-negative nucleus "nonda", and
# because it sorts before olig_UMAP and opc_UMAP alphabetically, the
# first-occurrence dedup below would stamp every oligodendrocyte and OPC
# as "nonda" and discard their real labels. The cross-species and macaque
# files are excluded for the same reason -- they re-label DA cells.
umap_files <- list.files(data_dir, pattern = "_UMAP\\.tsv$", full.names = TRUE)
umap_files <- umap_files[!grepl("^(liger_species|new_macaque|nurrneg|nurrpos)",
                                basename(umap_files))]
message("building labels from:")
print(basename(umap_files))

lab <- do.call(rbind, lapply(umap_files, function(f) {
  x <- read_scp(f)
  if (!"Cell_Type" %in% names(x)) {
    warning("no Cell_Type column in ", basename(f), "; skipping", call. = FALSE)
    return(NULL)
  }
  data.frame(NAME = x$NAME, cell_type = x$Cell_Type,
             origin = sub("_UMAP\\.tsv$", "", basename(f)),
             stringsAsFactors = FALSE)
}))

# A cell appearing in two files would otherwise be counted twice.
lab <- lab[!duplicated(lab$NAME), , drop = FALSE]
message(sprintf("%d labelled cells across %d types",
                nrow(lab), length(unique(lab$cell_type))))
print(table(lab$cell_type))

## --- metadata filters -------------------------------------------------
meta <- read_scp(file.path(data_dir, "METADATA_PD.tsv.gz"))

eligible <-
  meta$species__ontology_label == "Homo sapiens" &
  meta$organ__ontology_label   == "substantia nigra pars compacta" &
  meta$Status                  == "Ctrl"

message(sprintf("%d of %d metadata rows are human SNpc controls",
                sum(eligible), nrow(meta)))

meta_ok <- meta[eligible, c("NAME", "donor_id"), drop = FALSE]

## --- combine and downsample ------------------------------------------
sel <- merge(lab, meta_ok, by = "NAME")
message(sprintf("%d cells are both labelled and eligible", nrow(sel)))
print(table(sel$cell_type))
stopifnot(nrow(sel) > 0)

set.seed(42)
idx <- unlist(lapply(split(seq_len(nrow(sel)), sel$cell_type), function(i) {
  if (length(i) <= cap) i else sample(i, cap)
}), use.names = FALSE)
sel <- sel[sort(idx), , drop = FALSE]
message(sprintf("downsampled to %d cells", nrow(sel)))
print(table(sel$cell_type))

# =====================================================================
# STAGE B -- stream the wanted columns out of the 4.74 GB matrix
# =====================================================================

if (file.exists(cache_expr)) {
  message("loading cached expression subset")
  dat <- readRDS(cache_expr)
} else {

  barcodes <- readLines(gzfile(file.path(data_dir, "Homo_bcd.tsv.gz")))
  features <- read.delim(gzfile(file.path(data_dir, "Homo_features.tsv.gz")),
                         header = FALSE, stringsAsFactors = FALSE)
  genes <- make.unique(features[[1]])
  message(sprintf("expecting %d genes x %d cells", length(genes), length(barcodes)))

  want_cols <- match(sel$NAME, barcodes)
  if (anyNA(want_cols)) {
    warning(sprintf("%d selected cells absent from the barcode file; dropping",
                    sum(is.na(want_cols))), call. = FALSE)
    sel <- sel[!is.na(want_cols), , drop = FALSE]
    want_cols <- want_cols[!is.na(want_cols)]
  }
  ord <- order(want_cols)
  sel <- sel[ord, , drop = FALSE]
  want_cols <- want_cols[ord]

  message("streaming matrix -- this will take several minutes")
  con <- gzfile(file.path(data_dir, "Homo_matrix.mtx.gz"), "r")
  repeat { line <- readLines(con, n = 1); if (!startsWith(line, "%")) break }
  dims <- as.numeric(strsplit(trimws(line), "\\s+")[[1]])
  message(sprintf("header: %d rows x %d cols, %d non-zero",
                  dims[1], dims[2], dims[3]))
  stopifnot(dims[2] == length(barcodes))

  keep   <- logical(dims[2]); keep[want_cols] <- TRUE
  newidx <- integer(dims[2]); newidx[want_cols] <- seq_along(want_cols)

  ii <- jj <- xx <- list(); k <- 0L; seen <- 0
  repeat {
    d <- scan(con, what = list(integer(), integer(), double()),
              nmax = 5e6, quiet = TRUE)
    if (length(d[[1]]) == 0) break
    seen <- seen + length(d[[1]])
    s <- keep[d[[2]]]
    if (any(s)) {
      k <- k + 1L
      ii[[k]] <- d[[1]][s]; jj[[k]] <- newidx[d[[2]][s]]; xx[[k]] <- d[[3]][s]
    }
    cat(sprintf("\r  %.0f%% scanned", 100 * seen / dims[3]))
  }
  close(con); cat("\n")

  counts <- sparseMatrix(i = unlist(ii), j = unlist(jj), x = unlist(xx),
                         dims = c(dims[1], length(want_cols)),
                         dimnames = list(genes, sel$NAME))
  rm(ii, jj, xx); invisible(gc())

  libsize <- Matrix::colSums(counts); libsize[libsize == 0] <- 1
  expr <- counts
  expr@x <- log1p(expr@x / rep.int(libsize, diff(expr@p)) * 1e4)

  dat <- list(expr = expr, cell_type = sel$cell_type, donor = sel$donor_id)
  saveRDS(dat, cache_expr)
  message("cached to ", cache_expr)
}

message(sprintf("%d genes x %d cells, %d donors",
                nrow(dat$expr), ncol(dat$expr), length(unique(dat$donor))))

# =====================================================================
# STAGE C -- derive
# =====================================================================

cand_kamath <- derive_markers(dat$expr, dat$cell_type, donor = dat$donor)

print(round(tapply(cand_kamath$detect_gap, cand_kamath$cell_type, max), 3))
for (t in sort(unique(cand_kamath$cell_type))) {
  cat("\n---", t, "---\n")
  print(utils::head(cand_kamath[cand_kamath$cell_type == t,
                                c("gene", "detect_gap", "donor_consistency",
                                  "n_donors", "competitor")], 12))
}

saveRDS(cand_kamath, file.path(out_dir, "candidates_kamath.rds"))

# =====================================================================
# NEXT -- the first non-circular test available
# =====================================================================
#   mk_agarwal <- readRDS("data-raw/derived/markers_agarwal.rds")
#   benchmark_methods(dat$expr, dat$cell_type, markers = mk_agarwal)
#
# Kamath's labels are finer than Agarwal's, so map them to Agarwal's broad
# classes first or most cells will have no matching truth label.
