# Image loading, the per-image processing wrapper, specimen records, and the
# CSV form of the specimen store.
# Specification: docs/specs/2026-11-17_app_rewrite_spec.md, sections 6.1, 6.5.

SCHEMA_VERSION <- 1L
CSV_COLUMNS <- c("schema_version", "id", "specimen", "group", "px_per_mm", "base_row",
                 "threshold", "polarity", "speck", "direction", "manual_landmarks", "cropped",
                 "row", "x", "y")

# Load an image for processing, resized so its longer side is at most
# `max_side` pixels. Returns list(gray (matrix [x, y], values 0..1), raster,
# width, height, factor), where factor is working pixels per original pixel.
load_working_image <- function(path, max_side = 1200L) {
  im <- imager::load.image(path)
  if (imager::depth(im) > 1L) im <- imager::frame(im, 1)
  if (imager::spectrum(im) == 4L) im <- imager::rm.alpha(im)
  if (imager::spectrum(im) == 2L) im <- imager::channel(im, 1)
  if (max(im) > 1) im <- im / 255
  w <- imager::width(im); h <- imager::height(im)
  factor <- 1
  if (max(w, h) > max_side) {
    factor <- max_side / max(w, h)
    im <- imager::resize(im, round(w * factor), round(h * factor), interpolation_type = 3)
    im[im < 0] <- 0; im[im > 1] <- 1
  }
  gray <- if (imager::spectrum(im) == 3L) imager::grayscale(im) else im
  list(gray = as.matrix(gray[, , 1, 1]), raster = grDevices::as.raster(im),
       width = imager::width(im), height = imager::height(im), factor = factor)
}

# Process one image region.
#   crop:      c(x1, x2, y1, y2) in plot coordinates (pixel i spans i-1..i), or NULL
#   tip, base: optional manual landmarks in plot coordinates
# Returns a list with ok and warnings, and, when ok, the outline in plot
# coordinates, the standardized result, the landmarks in plot coordinates, and
# the settings used.
process_image <- function(img, crop = NULL, threshold = NULL, polarity = "auto", speck_side = 3L,
                          direction = "up", tip = NULL, base = NULL) {
  w <- img$width; h <- img$height
  whole <- c(0, w, 0, h)
  if (is.null(crop)) crop <- whole
  crop <- as.numeric(crop)
  # Part of the crop that lies inside the image, as pixel index ranges.
  x1 <- max(1, floor(crop[1]) + 1); x2 <- min(w, ceiling(crop[2]))
  y1 <- max(1, floor(crop[3]) + 1); y2 <- min(h, ceiling(crop[4]))
  if (length(crop) != 4L || !all(is.finite(crop)) || x2 - x1 + 1 < 20 || y2 - y1 + 1 < 20)
    return(list(ok = FALSE, warnings = "The crop is too small or lies outside the image. Drag a larger rectangle around the artifact."))
  xi <- seq.int(x1, x2); yi <- seq.int(y1, y2)
  seg <- segment_image(img$gray[xi, yi, drop = FALSE], threshold, polarity, speck_side)
  if (is.null(seg$mask)) return(list(ok = FALSE, warnings = seg$warnings, threshold = seg$threshold))
  # The mask and the plot coordinates of its lower-index corner, for display.
  shown <- list(mask = seg$mask, mask_origin = c(xi[1] - 1, yi[1] - 1))
  o <- mask_outline(seg$mask)
  if (is.null(o) || nrow(o) < 100L)
    return(list(ok = FALSE, threshold = seg$threshold, mask = shown$mask, mask_origin = shown$mask_origin,
                warnings = c(seg$warnings, "The detected outline has fewer than 100 points, which is too small to analyse. Use a larger image of the artifact or adjust the threshold.")))
  # Pixel index i has its centre at plot coordinate i - 0.5.
  o <- cbind(x = o[, 1] + xi[1] - 1.5, y = o[, 2] + yi[1] - 1.5)
  flip <- function(p) if (is.null(p)) NULL else c(p[1], -p[2])
  std <- tryCatch(standardize_outline(prepare_outline(o), direction, flip(tip), flip(base)),
                  error = function(e) e)
  if (inherits(std, "error"))
    return(list(ok = FALSE, threshold = seg$threshold, outline_plot = o,
                mask = shown$mask, mask_origin = shown$mask_origin,
                warnings = c(seg$warnings, conditionMessage(std))))
  list(ok = TRUE, warnings = c(seg$warnings, std$warnings), outline_plot = o, std = std,
       mask = shown$mask, mask_origin = shown$mask_origin,
       threshold = seg$threshold, polarity = seg$polarity,
       tip_plot = c(std$tip_input[1], -std$tip_input[2]),
       base_plot = c(std$base_input[1], -std$base_input[2]),
       settings = list(threshold = seg$threshold, polarity = seg$polarity, speck = as.integer(speck_side),
                       direction = direction, manual_landmarks = !is.null(tip) && !is.null(base),
                       cropped = !isTRUE(all.equal(as.numeric(crop), whole))))
}

# Scale conversion between working pixels (after resizing) and original pixels
# of an image whose resize factor is `factor` (working / original). A shared
# scale is held in original pixels per millimetre so that it transfers between
# images resized by different factors.
ppm_to_original <- function(ppm_working, factor) ppm_working / factor
ppm_to_working <- function(ppm_original, factor) ppm_original * factor

# A specimen record. `outline` is the oriented contour (tip up, tip first) in
# working pixels of its image; `base_row` is the row of the base landmark. The
# corresponding points are derived from the contour.
new_specimen <- function(id, name, group, outline, base_row, px_per_mm = NA_real_, settings = list()) {
  outline <- cbind(x = outline[, 1], y = outline[, 2])
  list(id = id, name = name, group = group, outline = outline, base_row = as.integer(base_row),
       points = points_from_outline(outline, as.integer(base_row)),
       px_per_mm = if (is.null(px_per_mm) || !is.finite(px_per_mm)) NA_real_ else px_per_mm,
       settings = settings)
}

# Reason why an oriented contour cannot be analysed, or NULL when it can. A
# usable contour has finite coordinates, no repeated consecutive vertex,
# positive (counterclockwise) area, the tip (row 1) above the base (row
# base_row), and yields finite measurements, 200 distinct corresponding points,
# and Fourier coefficients.
MAX_COORDINATE <- 1e6
outline_problem <- function(outline, base_row) {
  n <- nrow(outline)
  if (!all(is.finite(outline))) return("it contains values that are not finite numbers")
  if (max(abs(outline)) > MAX_COORDINATE) return("it contains coordinates larger than 1,000,000")
  if (is.na(base_row) || base_row < 2L || base_row > n) return("its base row is outside the outline")
  seg <- sqrt(rowSums((outline[c(2:n, 1L), , drop = FALSE] - outline)^2))
  if (any(seg <= 0)) return("it repeats a vertex")
  ext <- c(diff(range(outline[, 1])), diff(range(outline[, 2])))
  if (any(ext <= 0)) return("it has no width or no length")
  area <- polygon_moments(outline)$area
  if (!is.finite(area) || area <= 1e-9 * prod(ext)) return("it does not enclose an area in counterclockwise order")
  if (outline[1, 2] <= outline[base_row, 2]) return("its tip is not above its base")
  derived <- tryCatch({
    pts <- points_from_outline(outline, base_row)
    step <- sqrt(rowSums((pts[c(2:nrow(pts), 1L), , drop = FALSE] - pts)^2))
    m <- measure_outline(outline)
    f <- efa_vector(efa_forward(center_scale(pts), 30L))
    a <- asymmetry_score(pts, mirror_index())
    all(is.finite(pts)) && all(step > 0) && all(is.finite(f)) && is.finite(a) &&
      all(is.finite(unlist(m[c("length_px", "width_px", "area_px2", "length_width_ratio", "widest_position")])))
  }, error = function(e) FALSE)
  if (!isTRUE(derived)) return("its corresponding points or measurements cannot be computed")
  NULL
}

# TRUE when two contours have the same base row, the same number of vertices,
# and coordinates that differ by at most 1e-9 pixels. (A contour written to
# CSV and read back differs from the original by rounding of that order.)
same_contour <- function(a, b) {
  identical(dim(a$outline), dim(b$outline)) && identical(a$base_row, b$base_row) &&
    max(abs(a$outline - b$outline)) <= 1e-9
}

# One row of measurements per specimen.
specimen_measurements <- function(specimens) {
  if (length(specimens) == 0L) return(NULL)
  pair <- mirror_index()
  out <- do.call(rbind, lapply(specimens, function(s) {
    cbind(data.frame(specimen = s$name, group = s$group, stringsAsFactors = FALSE),
          measure_outline(s$outline, s$px_per_mm),
          px_per_mm = s$px_per_mm,
          asymmetry = asymmetry_score(s$points, pair))
  }))
  rownames(out) <- NULL
  out
}

# Data frame in schema version 1: one row per contour vertex.
specimens_to_csv <- function(specimens) {
  val <- function(x, default) if (is.null(x)) default else x
  out <- do.call(rbind, lapply(specimens, function(s) {
    st <- s$settings
    data.frame(schema_version = SCHEMA_VERSION, id = s$id, specimen = s$name, group = s$group,
               px_per_mm = s$px_per_mm, base_row = s$base_row,
               threshold = val(st$threshold, NA_real_), polarity = val(st$polarity, NA_character_),
               speck = val(st$speck, NA_integer_), direction = val(st$direction, NA_character_),
               manual_landmarks = val(st$manual_landmarks, NA), cropped = val(st$cropped, NA),
               row = seq_len(nrow(s$outline)), x = s$outline[, 1], y = s$outline[, 2],
               stringsAsFactors = FALSE)
  }))
  rownames(out) <- NULL
  out
}

# Read specimens from a schema version 1 CSV. Stops with a plain message when
# the file is not in that form. Returns a list of specimen records named by id.
specimens_from_csv <- function(path, max_specimens = 500L) {
  # Read everything as text first so that a file with other columns is
  # rejected by the column check and not by a type conversion.
  df <- tryCatch(utils::read.csv(path, stringsAsFactors = FALSE, colClasses = "character", na.strings = character(0)),
                 error = function(e) NULL)
  if (!is.null(df) && all(CSV_COLUMNS %in% names(df))) {
    text_cols <- c("id", "specimen", "group")
    for (col in setdiff(names(df), text_cols)) {
      v <- df[[col]]
      v[v %in% c("", "NA")] <- NA
      df[[col]] <- utils::type.convert(v, as.is = TRUE)
    }
  }
  if (is.null(df) || !all(CSV_COLUMNS %in% names(df)))
    stop("This file is not an outline file saved by this app. It must contain the columns: ",
         paste(CSV_COLUMNS, collapse = ", "), ".")
  if (nrow(df) == 0L) stop("The outline file contains no rows.")
  if (any(is.na(df$schema_version)) || any(df$schema_version != SCHEMA_VERSION))
    stop("This outline file has a schema version this app cannot read. Expected version ", SCHEMA_VERSION, ".")
  if (!is.numeric(df$x) || !is.numeric(df$y) || anyNA(df$x) || anyNA(df$y))
    stop("The x and y columns of the outline file must contain only numbers.")
  if (!all(is.finite(df$x)) || !all(is.finite(df$y)) || max(abs(c(df$x, df$y))) > MAX_COORDINATE)
    stop("The x and y columns of the outline file must contain finite numbers no larger than 1,000,000.")
  if (!is.numeric(df$row) || anyNA(df$row))
    stop("The row column of the outline file must contain only numbers.")
  if (anyNA(df$id) || any(!nzchar(trimws(df$id))))
    stop("Every row of the outline file must have an id.")
  df$group[is.na(df$group)] <- ""
  df$group <- trimws(df$group)    # labels are trimmed wherever they enter the app
  df$specimen[is.na(df$specimen)] <- ""
  ids <- unique(df$id)
  if (length(ids) > max_specimens) stop("The outline file holds more than ", max_specimens, " specimens.")
  na_null <- function(v) if (length(v) == 0L || is.na(v)) NULL else v
  out <- lapply(ids, function(id) {
    d <- df[df$id == id, , drop = FALSE]
    d <- d[order(d$row), , drop = FALSE]
    base_row <- suppressWarnings(as.integer(d$base_row[1]))
    label <- if (nzchar(d$specimen[1])) d$specimen[1] else id
    if (nrow(d) < 10L || is.na(base_row) || base_row < 2L || base_row > nrow(d) ||
        !isTRUE(all.equal(as.numeric(d$row), as.numeric(seq_len(nrow(d))))))
      stop("Specimen '", label, "' in the outline file has an incomplete outline.")
    one_value <- function(col) length(unique(d[[col]])) == 1L
    if (!all(vapply(c("specimen", "group", "px_per_mm", "base_row"), one_value, logical(1))))
      stop("Specimen '", label, "' in the outline file has rows that disagree about its name, group, scale, or base row.")
    ppm <- suppressWarnings(as.numeric(d$px_per_mm[1]))
    if (!is.na(ppm) && (!is.finite(ppm) || ppm <= 0))
      stop("Specimen '", label, "' in the outline file has a scale that is not a positive number.")
    problem <- outline_problem(cbind(d$x, d$y), base_row)
    if (!is.null(problem))
      stop("Specimen '", label, "' in the outline file has an outline that cannot be analysed: ", problem, ".")
    settings <- list(threshold = na_null(suppressWarnings(as.numeric(d$threshold[1]))),
                     polarity = na_null(as.character(d$polarity[1])),
                     speck = na_null(suppressWarnings(as.integer(d$speck[1]))),
                     direction = na_null(as.character(d$direction[1])),
                     manual_landmarks = na_null(as.logical(d$manual_landmarks[1])),
                     cropped = na_null(as.logical(d$cropped[1])))
    new_specimen(id, d$specimen[1], d$group[1], cbind(d$x, d$y), base_row,
                 suppressWarnings(as.numeric(d$px_per_mm[1])), settings[!vapply(settings, is.null, logical(1))])
  })
  stats::setNames(out, ids)
}

# Merge loaded specimens into a store (both lists named by id). A loaded
# specimen whose id exists replaces the stored one, including its name, group,
# scale, and settings, when their contours are the same (see same_contour),
# and otherwise is added under a new id. Display names are made unique.
merge_specimens <- function(store, loaded) {
  for (s in loaded) {
    id <- s$id
    if (!is.null(store[[id]]) && !same_contour(store[[id]], s)) {
      k <- 2L
      while (!is.null(store[[paste0(id, " #", k)]])) k <- k + 1L
      id <- paste0(id, " #", k)
      s$id <- id
    }
    s$name <- unique_specimen_name(s$name, store[setdiff(names(store), id)])
    store[[id]] <- s
  }
  store
}

# `name`, or `name (2)`, `name (3)`, ... if other specimens already use it.
unique_specimen_name <- function(name, others) {
  taken <- vapply(others, function(s) s$name, character(1))
  name <- trimws(name)
  if (!nzchar(name)) name <- "specimen"
  stem <- name; i <- 2L
  while (name %in% taken) { name <- sprintf("%s (%d)", stem, i); i <- i + 1L }
  name
}
