# Image segmentation: grayscale image -> binary mask -> closed outline.
# Specification: docs/specs/2026-11-17_app_rewrite_spec.md, section 6.1.

# Otsu's threshold on a 256-bin histogram of values in [0, 1].
otsu_threshold <- function(values, bins = 256L) {
  values <- as.numeric(values)
  values <- values[is.finite(values)]
  if (length(values) == 0L) stop("No pixel values supplied.")
  rng <- range(values)
  if (rng[1] == rng[2]) return(rng[1])
  breaks <- seq(0, 1, length.out = bins + 1L)
  values <- pmin(pmax(values, 0), 1)
  counts <- tabulate(findInterval(values, breaks, rightmost.closed = TRUE), nbins = bins)
  mids <- (breaks[-1] + breaks[-(bins + 1L)]) / 2
  p <- counts / sum(counts)
  w0 <- cumsum(p)
  m0 <- cumsum(p * mids)
  mt <- m0[bins]
  w1 <- 1 - w0
  between <- ifelse(w0 > 0 & w1 > 0, (mt * w0 - m0)^2 / (w0 * w1), NA_real_)
  best <- which(between == max(between, na.rm = TRUE))
  # Upper edge of the last bin assigned to the lower class. When several bins
  # tie (an empty gap between two modes), take the middle of the gap.
  breaks[floor(stats::median(best)) + 1L]
}

# Values on the border of a matrix.
border_values <- function(m) {
  nx <- nrow(m); ny <- ncol(m)
  c(m[1, ], m[nx, ], m[-c(1, nx), 1], m[-c(1, nx), ny])
}

# Shift a logical matrix by k cells along dimension `dim`, filling with FALSE.
shift_matrix <- function(m, k, dim) {
  if (k == 0L) return(m)
  n <- if (dim == 1L) nrow(m) else ncol(m)
  out <- matrix(FALSE, nrow(m), ncol(m))
  if (abs(k) >= n) return(out)
  src <- if (k > 0L) seq_len(n - k) else seq.int(1L - k, n)
  dst <- src + k
  if (dim == 1L) out[dst, ] <- m[src, ] else out[, dst] <- m[, src]
  out
}

# Morphological opening with a square structuring element of side `side`
# pixels. Cells outside the matrix count as background. Side 0 or 1 does nothing.
open_square <- function(fg, side) {
  side <- as.integer(side)
  if (side <= 1L) return(fg)
  off <- seq.int(-((side - 1L) %/% 2L), side %/% 2L)
  pass <- function(m, offsets, op) {
    for (d in 1:2) m <- Reduce(op, lapply(offsets, function(k) shift_matrix(m, k, d)))
    m
  }
  pass(pass(fg, off, `&`), -off, `|`)
}

# 4-connected components of the TRUE cells of a logical matrix. Returns
# list(lab, table); table has one row per component: label, area, touches_border.
component_table <- function(fg) {
  lab <- as.matrix(imager::label(imager::as.cimg(fg * 1), high_connectivity = FALSE))
  ids <- sort(unique(lab[fg]))
  if (length(ids) == 0L) {
    return(list(lab = lab, table = data.frame(label = integer(0), area = integer(0),
                                              touches_border = logical(0))))
  }
  areas <- tabulate(match(lab[fg], ids), nbins = length(ids))
  border_ids <- unique(border_values(lab)[border_values(fg)])
  list(lab = lab, table = data.frame(label = ids, area = areas, touches_border = ids %in% border_ids))
}

# Segment the artifact from a grayscale matrix [x, y] with values in [0, 1].
# The artifact is the largest 4-connected region that does not touch the border
# of the matrix. With polarity "auto" both classes (lighter than the threshold,
# and at or below it) are searched. Specks are removed first by an opening with
# a square of side `speck_side` pixels. If every region touches the border, the
# largest region overall is used and a warning is returned. Holes in the
# selected region are filled.
# Returns list(mask, threshold, polarity, warnings); mask is NULL if nothing
# was found.
segment_image <- function(gray, threshold = NULL, polarity = c("auto", "light", "dark"),
                          speck_side = 3L) {
  polarity <- match.arg(polarity)
  g <- as.matrix(gray)
  if (is.null(threshold)) threshold <- otsu_threshold(g)
  above <- g > threshold
  classes <- switch(polarity, auto = c("light", "dark"), polarity)

  candidates <- lapply(classes, function(cl) {
    fg <- open_square(if (cl == "light") above else !above, speck_side)
    ct <- component_table(fg)
    list(class = cl, fg = fg, lab = ct$lab, table = ct$table)
  })

  best <- NULL
  for (interior_only in c(TRUE, FALSE)) {
    for (cand in candidates) {
      tab <- cand$table
      if (interior_only) tab <- tab[!tab$touches_border, , drop = FALSE]
      if (nrow(tab) == 0L) next
      i <- which.max(tab$area)
      if (is.null(best) || tab$area[i] > best$area) {
        best <- list(area = tab$area[i], label = tab$label[i], cand = cand, interior = interior_only)
      }
    }
    if (!is.null(best)) break
  }
  if (is.null(best)) {
    return(list(mask = NULL, threshold = threshold, polarity = polarity,
                warnings = "No artifact was found. Adjust the threshold or the artifact brightness setting."))
  }

  mask <- best$cand$lab == best$label & best$cand$fg
  # Fill holes: regions outside the mask that do not connect to the border.
  outside <- component_table(!mask)
  open_ids <- outside$table$label[outside$table$touches_border]
  mask <- mask | (!mask & !(outside$lab %in% open_ids))

  # Competing objects: regions of every searched class, under the same border
  # rule as the selection, with at least one quarter of the selected area.
  # Regions that lie inside the filled mask are holes of the artifact and are
  # not counted.
  n_other <- 0L
  for (cand in candidates) {
    tab <- cand$table
    if (best$interior) tab <- tab[!tab$touches_border, , drop = FALSE]
    for (l in tab$label[tab$area >= 0.25 * best$area]) {
      if (identical(cand$class, best$cand$class) && l == best$label) next
      if (!all(mask[cand$lab == l & cand$fg])) n_other <- n_other + 1L
    }
  }

  warnings <- character(0)
  frac <- mean(mask)
  if (!best$interior)
    warnings <- c(warnings, "Every detected region touches the edge of the image or crop, so the artifact may be cut off or merged with the background. Redraw the crop so that background surrounds the artifact.")
  if (frac < 0.01)
    warnings <- c(warnings, "The detected artifact covers less than 1 percent of the image or crop. Crop closer to the artifact.")
  if (frac > 0.90)
    warnings <- c(warnings, "The detected region covers more than 90 percent of the image or crop. Check the artifact brightness setting.")
  if (n_other > 0)
    warnings <- c(warnings, sprintf("%d other object(s) at least one quarter the size of the selected one were detected and ignored. Crop to the artifact if the wrong object was selected.", n_other))
  list(mask = mask, threshold = threshold, polarity = best$cand$class, warnings = warnings)
}

# Closed boundary of a mask in pixel-index coordinates (x right, y down): the
# 0.5 iso-contour of the mask padded with one background pixel on every side,
# which guarantees closure. The longest contour is returned without a repeated
# end point.
mask_outline <- function(mask) {
  nx <- nrow(mask); ny <- ncol(mask)
  padded <- matrix(0, nx + 2L, ny + 2L)
  padded[2:(nx + 1L), 2:(ny + 1L)] <- mask * 1
  cs <- imager::contours(imager::as.cimg(padded), levels = 0.5)
  if (length(cs) == 0L) return(NULL)
  n <- vapply(cs, function(z) length(z$x), integer(1))
  z <- cs[[which.max(n)]]
  xy <- cbind(x = z$x - 1, y = z$y - 1)
  if (all(xy[1, ] == xy[nrow(xy), ])) xy <- xy[-nrow(xy), , drop = FALSE]
  xy
}
