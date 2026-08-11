# midbrainType

Automated cell type annotation for midbrain single-cell and spatial
transcriptomic data — from a raw expression matrix straight to per-cell
labels, with the marker gene reference bundled inside the package. No
external atlas or reference download required.

This is the first region covered under a planned region-by-region rollout
(midbrain now; other regions later, as separate reference tables or
companion packages).

## Status

Early skeleton (v0.0.0.9000). The bundled marker set in
`data-raw/prepare_markers.R` is a starter list of canonical markers — not
yet curated against a specific reference (e.g. the Agarwal midbrain atlas,
GSE140231) or validated against real Xenium/spatial data. Treat calls as
provisional until that curation pass happens.

## Installation

```r
# install.packages("devtools")
devtools::install_github("Notty4/midbrainType")
```

## Usage

```r
library(midbrainType)

# expr: a genes x cells matrix of normalised expression, with gene symbols
# as row names and cell/spot IDs as column names.
labels <- annotate_midbrain_cells(expr)
table(labels$cell_type)

# Raw per-cell-type scores, if you want them instead of hard labels:
scores <- score_markers(expr)

# What cell types does the bundled reference cover?
list_midbrain_cell_types()
```

## How it works

Each cell type in the marker reference gets a score per cell, and each cell
is assigned the highest-scoring type. Four scoring methods are available —
`zscore` (default), `ucell`, `control`, and `mean_weighted` — which differ
in what they normalise against. See `vignette("scoring-methods")` for a
benchmark and guidance on choosing.

Cells can be flagged rather than forced into a call:

```r
labels <- annotate_midbrain_cells(
  expr,
  min_score  = 0,    # must beat background, else "Unassigned"
  min_margin = 0.1   # top two types must separate, else "Ambiguous"
)
```

### A caveat worth reading

Benchmarking on simulated data showed that when marker sets differ in
absolute abundance — which happens routinely via ambient RNA and transcript
spillover in spatial assays — *every* scoring method degrades badly, and the
best of them still lands under 50% accuracy. Only `zscore` degrades
gracefully, which is why it is the default, but the honest conclusion is
that the scoring formula is not where this problem gets solved. Marker
curation matters more.

## Package layout

- `R/` — `annotate_midbrain_cells()`, `score_markers()`,
  `list_midbrain_cell_types()`
- `data/midbrain_markers.rda` — bundled marker reference
- `data-raw/prepare_markers.R` — script that builds the bundled reference;
  edit this and re-run to update the marker list
- `tests/testthat/` — unit tests
- `vignettes/` — worked example

## Roadmap

- Curate the marker set against a real reference (thesis / literature)
  rather than the current canonical starter list — currently the highest-value
  work, per the benchmark above
- Validate on real Xenium/spatial data, not just simulated counts
- Handle ambient RNA / transcript spillover explicitly, rather than hoping
  the scoring method absorbs it
- Extend beyond the midbrain to further regions
