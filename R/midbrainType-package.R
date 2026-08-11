#' @keywords internal
"_PACKAGE"

# `midbrain_markers` is a lazy-loaded package dataset used as a default
# argument value. R CMD check's static analysis can't see that binding,
# so declare it here to silence a spurious note.
utils::globalVariables("midbrain_markers")

devtools::document()
devtools::test()     # expect 21 tests now, not 9
devtools::check()

writeLines(c(
  "#' @keywords internal",
  "#' @importFrom utils globalVariables",
  '"_PACKAGE"',
  "",
  "# `midbrain_markers` is a lazy-loaded package dataset used as a default",
  "# argument value. R CMD check's static analysis can't see that binding,",
  "# so declare it here to silence a spurious note.",
  'globalVariables("midbrain_markers")'
), "R/midbrainType-package.R")

devtools::document()   # regenerates NAMESPACE with the importFrom line
devtools::check()