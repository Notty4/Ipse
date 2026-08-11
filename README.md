# Ipse

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
devtools::install_github("Notty4/Ipse")
```

## Usage

```r
library(Ipse)

# expr: a genes x cells matrix of normalised expression, with gene symbols
# as row names and cell/spot IDs as column names.
labels <- annotate_cells(expr)
table(labels$cell_type)

# Raw per-cell-type scores, if you want them instead of hard labels:
scores <- score_markers(expr)

# What cell types does the bundled reference cover?
list_cell_types()
```

## How it works

Each cell type in the marker reference gets a score per cell, and each cell
is assigned the highest-scoring type. Four scoring methods are available —
`zscore` (default), `ucell`, `control`, and `mean_weighted` — which differ
in what they normalise against. See `vignette("scoring-methods")` for a
benchmark and guidance on choosing.

### Identity, effector, and opt-in sections

Markers are layered. **Identity** markers (lineage transcription factors)
are the only ones that assign a label. **Effector** markers score the
functional programme separately, so a cell that retains its lineage but
has lost that programme is described rather than mislabelled. **State**
markers form named *sections* — condition-associated signatures that are
never scored unless you ask for them:

```r
list_state_sections()          # what is available, and where it came from

labels <- annotate_cells(expr, sections = "reactive")
attr(labels, "sections")       # provenance travels with the result
```

A section gives a named, reference-anchored cell group you can compare
against published work, instead of an opaque cluster index. Sections never
influence the assigned label — that stays derived from identity markers
alone, so annotation cannot encode condition status.

Cells can be flagged rather than forced into a call:

```r
labels <- annotate_cells(
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

## Deriving your own markers

The bundled markers are a placeholder. To build a real reference from an
annotated dataset — run this on **control samples only**, so that labels
can never encode disease status:

```r
candidates <- derive_markers(
  expr   = counts,          # genes x cells, log-normalised
  labels = meta$cell_type,
  donor  = meta$donor,      # optional: adds donor_consistency
  panel  = xenium_genes     # optional: restrict to what you can measure
)

# Require replication across independent cohorts
consensus <- consensus_markers(
  list(cohortA = cand_a, cohortB = cand_b),
  top_n = 50, min_cohorts = 2
)

markers <- as_marker_table(
  consensus, n = 6, min_gap = 0.3,
  effector = c("TH", "SLC6A3", "SLC18A2"),
  source = "Control samples, <accessions>"
)
```

Candidates are ranked by **exclusivity**, not fold change: each cell type
is compared against its single strongest competitor, so a gene shared by
two lineages scores near zero however abundant it is. The default ranking
uses detection rate (`detect_gap`) rather than expression level, because a
gene that is abundant everywhere but slightly higher in one type still
scores well on expression and poorly on detection — and that is precisely
the pattern that defeats annotation in spatial data.

Layer assignment is not automated: name your `effector` genes explicitly.

## Package layout

- `R/` — `annotate_cells()`, `score_markers()`,
  `list_cell_types()`
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
