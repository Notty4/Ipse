#' @keywords internal
#' @importFrom utils globalVariables
"_PACKAGE"

# `midbrain_markers` is a lazy-loaded package dataset used as a default
# argument value. R CMD check's static analysis can't see that binding,
# so declare it here to silence a spurious note.
globalVariables("midbrain_markers")
