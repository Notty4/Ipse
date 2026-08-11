# data-raw/GSE140231_pipeline.R
#
# Full pipeline for Agarwal et al. 2020 (GSE140231): load -> QC -> cluster
# -> label -> derive candidate markers.
#
# Written to be re-run safely. The clustered object is cached, so a second
# run skips straight to the labelling step rather than redoing twenty
# minutes of work -- which is what makes losing the session survivable.
#
# Run from the Ipse project root.

library(Seurat)
library(Matrix)
devtools::load_all()          # or library(Ipse) once installed

data_dir <- "C:/Users/User/Downloads/GSE140231_RAW"
out_dir  <- "data-raw/derived"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

cache <- file.path(out_dir, "GSE140231_clustered.rds")

# =====================================================================
# STAGES A + B -- load, QC, cluster   (skipped if the cache exists)
# =====================================================================

if (file.exists(cache)) {
  message("loading cached clustered object")
  obj <- readRDS(cache)
} else {

  ## --- A: read the per-sample tar.gz archives ------------------------
  archives <- list.files(data_dir, pattern = "\\.tar\\.gz$", full.names = TRUE)
  stopifnot(length(archives) > 0)

  read_sample <- function(archive) {
    dest <- file.path(tempdir(), "GSE140231",
                      tools::file_path_sans_ext(
                        tools::file_path_sans_ext(basename(archive))))
    dir.create(dest, recursive = TRUE, showWarnings = FALSE)
    utils::untar(archive, exdir = dest)

    inner <- list.dirs(dest, recursive = FALSE)
    inner <- if (length(inner) > 0) inner[1] else dest

    genes    <- utils::read.delim(file.path(inner, "genes.tsv"),
                                  header = FALSE, stringsAsFactors = FALSE)
    barcodes <- readLines(file.path(inner, "barcodes.tsv"))
    m <- Matrix::readMM(file.path(inner, "matrix.mtx"))

    rownames(m) <- make.unique(genes[[2]])
    colnames(m) <- barcodes
    methods::as(m, "CsparseMatrix")
  }

  mats <- lapply(archives, function(a) { message("reading ", basename(a)); read_sample(a) })
  names(mats) <- sub("^(GSM[0-9]+)_(.*)\\.tar\\.gz$", "\\2", basename(archives))

  # C = cortex, N = nigra; trailing digit = donor. Cell counts match the
  # paper's split (~5,943 nigra / ~10,706 cortex), which is what confirms
  # the letters are the right way round.
  region <- ifelse(grepl("_C[0-9]", names(mats)), "cortex", "nigra")
  donor  <- sub(".*_[CN]([0-9]+).*", "\\1", names(mats))
  print(data.frame(sample = names(mats), region, donor,
                   n_cells = vapply(mats, ncol, integer(1)), row.names = NULL))

  # Nigra only. The matched cortex samples are worth keeping for the
  # region-identification work later, but including them here would put
  # cortical cell types into the specificity ranking.
  keep_s <- which(region == "nigra")
  mats <- mats[keep_s]; donor <- donor[keep_s]

  common <- Reduce(intersect, lapply(mats, rownames))
  counts <- do.call(cbind, lapply(seq_along(mats), function(i) {
    m <- mats[[i]][common, , drop = FALSE]
    colnames(m) <- paste(names(mats)[i], colnames(m), sep = "_")
    m
  }))
  sample_id <- rep(names(mats), vapply(mats, ncol, integer(1)))
  donor_id  <- rep(donor,       vapply(mats, ncol, integer(1)))

  obj <- CreateSeuratObject(counts, project = "GSE140231",
                            min.cells = 3, min.features = 200)
  obj$sample <- sample_id[match(colnames(obj), colnames(counts))]
  obj$donor  <- donor_id[match(colnames(obj), colnames(counts))]

  ## --- B: QC, normalise, cluster -------------------------------------
  obj[["percent_mt"]] <- PercentageFeatureSet(obj, pattern = "^MT-")

  # 5% is a whole-cell figure and permissive for nuclei; it retained 5,945
  # of ~6,100 last time, matching the paper's 5,943.
  obj <- subset(obj, subset = nFeature_RNA > 200 &
                              nFeature_RNA < 6000 &
                              percent_mt < 5)
  message(sprintf("%d nuclei after QC", ncol(obj)))

  obj <- NormalizeData(obj)
  obj <- FindVariableFeatures(obj, nfeatures = 2000)
  obj <- ScaleData(obj)
  obj <- RunPCA(obj, npcs = 30)
  obj <- FindNeighbors(obj, dims = 1:20)
  obj <- FindClusters(obj, resolution = 0.5)
  obj <- RunUMAP(obj, dims = 1:20)

  saveRDS(obj, cache)
  message("cached to ", cache)
}

# =====================================================================
# Check the clustering matches what the labels below were derived from
# =====================================================================

message(sprintf("%d clusters", length(levels(obj$seurat_clusters))))
print(table(obj$seurat_clusters))

e <- GetAssayData(obj, layer = "data")
cat(sprintf("TH+ in cluster 12: %.2f (expect ~0.76)\n",
            mean(e["TH", obj$seurat_clusters == "12"] > 0)))
cat(sprintf("SLC32A1+ in cluster 11: %.2f (expect ~0.29)\n",
            mean(e["SLC32A1", obj$seurat_clusters == "11"] > 0)))

# If those two are badly off, the cluster numbering has shifted and the
# labels below would attach to the wrong clusters. Stop and re-derive the
# mapping from FindAllMarkers rather than trusting it.

# =====================================================================
# Labels
# =====================================================================
#   0,1,4,6  oligodendrocyte (OPALIN, QDPR, CA2; 4 is a cholesterol-
#                             biosynthesis state, 1 splits partly on XIST)
#   3        astrocyte       (AQP4, GJA1, ETNPPL, SLC14A1)
#   5,8      microglia       (C1QB/C1QC, CD14, HLA-DR; 8 activated)
#   7        OPC             (PDGFRA, VCAN, TNR, GPR17)
#   10       endothelial     (CLDN5, FLT1 -- also carries pericyte markers
#                             RGS5/NDUFA4L2, so really "vascular")
#   11       GABAergic       (SLC32A1, GAD2, SOX14, OTX2, GATA3)
#   12       dopaminergic    (TH, SLC6A3, ALDH1A1, SLC10A4)
#   2        NA -- S100B/LGALS1/TUBB2B, no clean canonical match
#   9        NA -- clearly neuronal, subtype not resolvable here
#   13       NA -- ependymal (ciliary genes), but donor 1 only

cluster_labels <- c(
  "0"  = "Oligodendrocyte",
  "1"  = "Oligodendrocyte",
  "2"  = NA_character_,
  "3"  = "Astrocyte",
  "4"  = "Oligodendrocyte",
  "5"  = "Microglia",
  "6"  = "Oligodendrocyte",
  "7"  = "OPC",
  "8"  = "Microglia",
  "9"  = NA_character_,
  "10" = "Endothelial",
  "11" = "GABAergic neuron",
  "12" = "Dopaminergic neuron",
  "13" = NA_character_
)

obj$cell_type <- unname(cluster_labels[as.character(obj$seurat_clusters)])
print(table(obj$cell_type, useNA = "ifany"))

keep <- !is.na(obj$cell_type)
stopifnot(sum(keep) > 0)
message(sprintf("%d of %d cells labelled", sum(keep), length(keep)))

# =====================================================================
# Derive candidate markers
# =====================================================================

cand <- derive_markers(
  expr   = GetAssayData(obj, layer = "data")[, keep],
  labels = obj$cell_type[keep],
  donor  = obj$donor[keep]
)

print(quantile(cand$detect_gap, c(0.5, 0.9, 0.99, 1)))

# Per-type ceilings, not a global floor: oligodendrocytes top out near
# 0.30 because ambient myelin transcripts are detected in every cell type,
# so a global min_gap of 0.3 would delete them entirely.
print(round(tapply(cand$detect_gap, cand$cell_type, max), 3))

for (t in sort(unique(cand$cell_type))) {
  cat("\n---", t, "---\n")
  print(utils::head(cand[cand$cell_type == t,
                         c("gene", "detect_gap", "specificity",
                           "donor_consistency", "n_donors", "competitor")], 15))
}

saveRDS(cand, file.path(out_dir, "candidates_GSE140231.rds"))
message("saved candidates")
