# EpiExposure global-variable declarations
#
# No `utils::globalVariables()` declarations are currently required.
#
# Earlier EpiExposure versions used non-standard evaluation in data-manipulation
# code and registered column names here to suppress "no visible binding for
# global variable" notes during `R CMD check`.
#
# The current audited implementation accesses columns explicitly (for example,
# with `$`, `[[`, character column names, or ordinary local variables) and does
# not rely on unquoted data-masked symbols in package code. Keeping obsolete
# declarations would therefore provide no benefit and could mask future
# namespace/programming mistakes that `R CMD check` should report.
#
# If non-standard evaluation is intentionally reintroduced in a future version,
# prefer explicit pronouns/selectors (for example `.data[[name]]`) or
# namespace-qualified tidy-evaluation helpers. Add `utils::globalVariables()`
# only when an unavoidable check note is understood and documented.

NULL
