# EpiExposure - internal spaMM spatial specification helpers (not exported)

# Spatial effects are distinct from the conventional random intercept. These
# helpers construct only the spatial term and never modify epidemic histories.

.epix_spatial_quote_name <- function(name) {
  if (!is.character(name) || length(name) != 1L || is.na(name) ||
      !nzchar(trimws(name)) || grepl("`", name, fixed = TRUE)) {
    stop("Spatial column names must be non-empty and cannot contain backticks.",
         call. = FALSE)
  }
  if (identical(make.names(name), name)) name else paste0("`", name, "`")
}

.build_spamm_spatial_term <- function(spatial_effect = NULL,
                                      spatial_structure = "matern",
                                      spatial_group = NULL) {
  if (is.null(spatial_effect)) return(NULL)
  if (!is.character(spatial_effect) || length(spatial_effect) != 2L ||
      anyNA(spatial_effect) || any(!nzchar(trimws(spatial_effect))) ||
      anyDuplicated(spatial_effect)) {
    stop("`spatial_effect` must contain exactly two distinct, non-empty coordinate column names.",
         call. = FALSE)
  }
  if (!is.character(spatial_structure) || length(spatial_structure) != 1L ||
      is.na(spatial_structure) || !nzchar(trimws(spatial_structure)) ||
      tolower(trimws(spatial_structure)) != "matern") {
    stop("Only `spatial_structure = 'matern'` is currently supported.",
         call. = FALSE)
  }
  if (!is.null(spatial_group) &&
      (!is.character(spatial_group) || length(spatial_group) != 1L ||
       is.na(spatial_group) || !nzchar(trimws(spatial_group)))) {
    stop("`spatial_group` must be NULL or one non-empty column name.", call. = FALSE)
  }
  coords <- vapply(spatial_effect, .epix_spatial_quote_name, character(1))
  term <- paste0("Matern(1 | ", coords[[1L]], " + ", coords[[2L]])
  if (!is.null(spatial_group)) {
    term <- paste0(term, " %in% ", .epix_spatial_quote_name(spatial_group))
  }
  paste0(term, ")")
}

.epix_validate_spatial_spec <- function(data, model_engine,
                                        spatial_effect = NULL,
                                        spatial_structure = "matern",
                                        spatial_group = NULL,
                                        forbidden = character(0)) {
  if (is.null(spatial_effect)) {
    if (!is.null(spatial_group)) {
      stop("`spatial_group` can only be supplied when `spatial_effect` is supplied.",
           call. = FALSE)
    }
    # `spatial_structure` is deliberately inert for legacy non-spatial fits.
    return(list(effect = NULL, structure = NULL, group = NULL, term = NULL,
                has_spatial_effect = FALSE))
  }
  if (!identical(model_engine, "spamm")) {
    stop("Spatial effects are supported only with `model_engine = 'spamm'`. ",
         "Use 'spamm' or remove `spatial_effect`.", call. = FALSE)
  }
  # Validate names before accessing data or forming R formula expressions.
  term <- .build_spamm_spatial_term(spatial_effect, spatial_structure, spatial_group)
  if (length(intersect(spatial_effect, forbidden)) ||
      any(startsWith(spatial_effect, "cb_"))) {
    stop("Coordinate names in `spatial_effect` cannot be a response, time, ",
         "exposure, or reserved `cb_*` column.", call. = FALSE)
  }
  if (!is.null(spatial_group) &&
      (spatial_group %in% forbidden || startsWith(spatial_group, "cb_") ||
       spatial_group %in% spatial_effect)) {
    stop("`spatial_group` must be distinct from response, time, exposure, ",
         "coordinate, and reserved `cb_*` columns.", call. = FALSE)
  }
  spatial_columns <- unique(c(spatial_effect, spatial_group))
  missing <- setdiff(spatial_columns, names(data))
  if (length(missing)) {
    stop("Spatial column(s) missing from `data`: ", paste(missing, collapse = ", "),
         ".", call. = FALSE)
  }
  for (coord in spatial_effect) {
    value <- data[[coord]]
    if (!is.numeric(value)) {
      stop("Spatial coordinate column '", coord, "' must be numeric.", call. = FALSE)
    }
    if (anyNA(value) || any(!is.finite(value))) {
      stop("Spatial coordinate column '", coord,
           "' must contain only finite non-missing values.", call. = FALSE)
    }
  }
  if (nrow(unique(data[, spatial_effect, drop = FALSE])) < 2L) {
    stop("`spatial_effect` requires at least two distinct coordinate pairs.",
         call. = FALSE)
  }
  if (!is.null(spatial_group)) {
    value <- data[[spatial_group]]
    if (anyNA(value)) {
      stop("`spatial_group` ('", spatial_group,
           "') cannot contain missing values.", call. = FALSE)
    }
    # Whole-number numeric years often arrive as double after data import.
    integer_like <- is.numeric(value) && all(is.finite(value)) &&
      all(value == floor(value))
    if (!is.factor(value) && !is.character(value) && !integer_like) {
      stop("`spatial_group` must be a factor, character, or integer-valued column.",
           call. = FALSE)
    }
  }
  list(effect = spatial_effect, structure = "matern", group = spatial_group,
       term = term, has_spatial_effect = TRUE)
}

# Validate group-level metadata *before* costly grouped CV, while leaving the
# epidemic definition independent of coordinate pairs or spatial group levels.
.epix_validate_spatial_constancy <- function(data, spatial_effect = NULL,
                                             spatial_group = NULL,
                                             epi_id = "epi_id") {
  if (is.null(spatial_effect)) return(invisible(TRUE))
  for (column in unique(c(spatial_effect, spatial_group))) {
    unique_pairs <- unique(data[, c(epi_id, column), drop = FALSE])
    bad <- duplicated(as.character(unique_pairs[[epi_id]]))
    if (any(bad)) {
      stop("Spatial metadata column '", column,
           "' is not constant within epi_id '",
           as.character(unique_pairs[[epi_id]][which(bad)[1L]]), "'.",
           call. = FALSE)
    }
  }
  invisible(TRUE)
}
