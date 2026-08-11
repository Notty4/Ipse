# data-raw/01_derive_GSE140231.R
#
# Derive candidate identity markers from Agarwal et al. 2020 (GSE140231):
# ~17,000 nuclei from matched substantia nigra and cortex of 5 control
# post-mortem donors. All donors are controls, so no disease subsetting is
# needed -- which is why this is the first cohort to run.
#
# STAGES
#   A. Inspect what GEO actually gave you and load it
#   B. QC, normalise, cluster
#   C. Label the clusters (the step that needs your judgement)
#   D. Derive markers and save candidates
#
# Stage C is not automatable and should not be rushed: every marker you
# derive downstream inherits whatever mistakes you make naming clusters.
#
# NOTE: the file-reading in Stage A is written defensively because GEO
# deposits vary in layout. Run Stage A's inspection block first and adjust
# the reader to match what you actually see, rather than assuming.

library(Seurat)
library(Matrix)
library(Ipse)

data_dir <- "~/data/GSE140231"     # <- where you unpacked the GEO download
out_dir  <- "data-raw/derived"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ---------------------------------------------------------------------
# STAGE A -- inspect, then load
# ---------------------------------------------------------------------

files <- list.files(data_dir, recursive = TRUE, full.names = TRUE)
print(basename(files))
print(utils::head(readLines(files[1], n = 3)))   # peek at the format

# GSE140231 supplementary files are per-sample UMI matrices. Depending on
# the download you will have either (a) one dense .txt/.csv per sample, or
# (b) 10x-style barcodes/features/matrix triplets per sample. Use whichever
# branch matches what printed above; delete the other.

## --- branch (a): one delimited matrix per sample -------------------
read_sample_dense <- function(path) {
  m <- utils::read.delim(path, row.names = 1, check.names = FALSE)
  as(as.matrix(m), "dgCMatrix")
}

sample_files <- grep("\\.(txt|csv|tsv)(\\.gz)?$", files, value = TRUE)
mats <- lapply(sample_files, read_sample_dense)
names(mats) <- sub("_.*$", "", basename(sample_files))   # GSM id as sample name

## --- branch (b): 10x triplets, one directory per sample ------------
# sample_dirs <- list.dirs(data_dir, recursive = FALSE)
# mats <- lapply(sample_dirs, Seurat::Read10X)
# names(mats) <- basename(sample_dirs)

# Build one object, keeping sample identity so it can be used as `donor`
# later. Genes are intersected across samples so the matrices bind.
common <- Reduce(intersect, lapply(mats, rownames))
message(sprintf("%d genes shared across %d samples", length(common), length(mats)))

counts <- do.call(cbind, lapply(names(mats), function(s) {
  m <- mats[[s]][common, , drop = FALSE]
  colnames(m) <- paste(s, colnames(m), sep = "_")
  m
}))
sample_id <- sub("_.*$", "", colnames(counts))

obj <- CreateSeuratObject(counts, project = "GSE140231", min.cells = 3, min.features = 200)
obj$sample <- sample_id

# ---------------------------------------------------------------------
# STAGE B -- QC, normalise, cluster
# ---------------------------------------------------------------------

obj[["percent_mt"]] <- PercentageFeatureSet(obj, pattern = "^MT-")
VlnPlot(obj, c("nFeature_RNA", "nCount_RNA", "percent_mt"), ncol = 3)

# Set these from the plot above rather than trusting the defaults --
# single-nucleus data has lower mitochondrial content than whole cells,
# so a 5% cutoff borrowed from a scRNA-seq tutorial is too permissive.
obj <- subset(obj, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent_mt < 5)

obj <- NormalizeData(obj)
obj <- FindVariableFeatures(obj, nfeatures = 2000)
obj <- ScaleData(obj)
obj <- RunPCA(obj, npcs = 30)
ElbowPlot(obj, ndims = 30)

obj <- FindNeighbors(obj, dims = 1:20)
obj <- FindClusters(obj, resolution = 0.5)
obj <- RunUMAP(obj, dims = 1:20)
DimPlot(obj, group.by = "sample")     # check donors mix rather than split

# ---------------------------------------------------------------------
# STAGE C -- label the clusters   (YOUR JUDGEMENT, NOT THE SCRIPT'S)
# ---------------------------------------------------------------------
#
# The bundled markers give a first pass, but treat them as a suggestion.
# They are canonical genes, not a validated reference -- that is the whole
# reason you are doing this.

expr <- GetAssayData(obj, layer = "data")
first_pass <- annotate_cells(expr, seed = 42)

# Assign per CLUSTER, not per cell: cluster-level calls are far more stable
# than single-cell calls, and clusters are what you are naming.
cluster_call <- tapply(first_pass$cell_type, obj$seurat_clusters,
                       function(x) names(sort(table(x), decreasing = TRUE))[1])
print(cluster_call)

# Now check it independently. Unsupervised cluster markers should agree
# with the call; where they disagree, believe these, not the first pass.
cluster_markers <- FindAllMarkers(obj, only.pos = TRUE,
                                  min.pct = 0.25, logfc.threshold = 0.25)
top10 <- do.call(rbind, lapply(split(cluster_markers, cluster_markers$cluster),
                               function(d) utils::head(d[order(-d$avg_log2FC), ], 10)))
print(top10[, c("cluster", "gene", "avg_log2FC", "pct.1", "pct.2")])

# Fill this in yourself once you have looked at the above. Any cluster you
# are unsure about should be labelled NA and excluded -- a wrong label is
# far more damaging than a missing one, because every marker derived from
# it will be wrong too.
cluster_labels <- c(
  "0" = NA_character_,
  "1" = NA_character_
  # ... one entry per cluster
)

obj$cell_type <- unname(cluster_labels[as.character(obj$seurat_clusters)])
table(obj$cell_type, useNA = "ifany")

# ---------------------------------------------------------------------
# STAGE D -- derive
# ---------------------------------------------------------------------

keep <- !is.na(obj$cell_type)
cand <- derive_markers(
  expr   = GetAssayData(obj, layer = "data")[, keep],
  labels = obj$cell_type[keep],
  donor  = obj$sample[keep]
)

# Look at the detection-gap distribution BEFORE choosing a threshold. The
# default min_gap of 0.1 is a placeholder, not a recommendation: set the
# floor above wherever the noise sits in your data.
hist(cand$detect_gap, breaks = 60,
     main = "detect_gap, GSE140231", xlab = "detection gap")

for (t in unique(cand$cell_type)) {
  cat("\n---", t, "---\n")
  print(utils::head(cand[cand$cell_type == t,
                         c("gene", "detect_gap", "specificity",
                           "donor_consistency", "competitor")], 15))
}

saveRDS(cand, file.path(out_dir, "candidates_GSE140231.rds"))

# Next: repeat for a second cohort (GSE178265 -- controls only, annotations
# from Broad Single Cell Portal SCP1768), then combine:
#
#   cons <- consensus_markers(
#     list(agarwal = cand_agarwal, kamath = cand_kamath),
#     top_n = 50, min_cohorts = 2
#   )
#   markers <- as_marker_table(cons, n = 6, min_gap = <from the histogram>,
#                             effector = c("TH", "SLC6A3", "SLC18A2", "DDC"),
#                             source = "Control donors, GSE140231 + GSE178265")
