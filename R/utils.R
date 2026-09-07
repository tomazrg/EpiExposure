# EpiExposure namespace imports
#
# Package functions now use explicit namespace qualification (`pkg::fun()`) for
# nearly all non-base functions. This keeps the package namespace small and
# makes dependencies visible at the call site.
#
# `setNames()` remains imported from `stats` because it is intentionally called
# unqualified by current internal EpiExposure code. If those calls are changed
# to `stats::setNames()` in a future cleanup, this import can also be removed.

#' @importFrom stats setNames
NULL
