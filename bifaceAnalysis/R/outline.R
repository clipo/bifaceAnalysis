# Outline standardization: orientation, tip and base landmarks, corresponding
# points, and linear measurements. All functions take and return plain matrices.
# Specification: docs/specs/2026-11-17_app_rewrite_spec.md, sections 6.2 to 6.4.

N_EDGE <- 100L                 # points per edge
N_POINTS <- 2L * N_EDGE        # points per specimen
TIP_INDEX <- 1L
BASE_INDEX <- N_EDGE + 1L

direction_vector <- function(direction) {
  switch(direction,
         up = c(0, 1), down = c(0, -1), left = c(-1, 0), right = c(1, 0),
         stop("Unknown tip direction: ", direction))
}

# Signed area, centroid, and central second moments of a closed polygon.
polygon_moments <- function(xy) {
  x <- xy[, 1]; y <- xy[, 2]
  xn <- c(x[-1], x[1]); yn <- c(y[-1], y[1])
  cr <- x * yn - xn * y
  A <- sum(cr) / 2
  cx <- sum((x + xn) * cr) / (6 * A)
  cy <- sum((y + yn) * cr) / (6 * A)
  ixx <- sum((x^2 + x * xn + xn^2) * cr) / 12
  iyy <- sum((y^2 + y * yn + yn^2) * cr) / 12
  ixy <- sum((x * yn + 2 * x * y + 2 * xn * yn + xn * y) * cr) / 24
  cov <- matrix(c(ixx / A - cx^2, ixy / A - cx * cy,
                  ixy / A - cx * cy, iyy / A - cy^2), 2, 2)
  list(area = A, centroid = c(cx, cy), cov = cov)
}

# Convert an image-coordinate outline (y down) to y-up, counterclockwise, with
# consecutive duplicate points removed.
prepare_outline <- function(xy) {
  xy <- cbind(x = xy[, 1], y = -xy[, 2])
  n <- nrow(xy)
  if (n > 1L && all(abs(xy[1, ] - xy[n, ]) < 1e-12)) xy <- xy[-n, , drop = FALSE]
  prev <- xy[c(nrow(xy), seq_len(nrow(xy) - 1L)), , drop = FALSE]
  xy <- xy[sqrt(rowSums((xy - prev)^2)) > 1e-12, , drop = FALSE]
  if (nrow(xy) < 3L) stop("The outline has fewer than three distinct points.")
  if (polygon_moments(xy)$area < 0) xy <- xy[nrow(xy):1, , drop = FALSE]
  xy
}

# Unit vector along the tip-base axis, pointing toward the tip: the principal
# axis of the polygon's area with the larger absolute cosine to the declared
# direction (the major axis on an exact tie). When the principal moments differ
# by less than `min_ratio` the declared direction itself is returned.
# Returns list(u, isotropic).
axis_direction <- function(xy, direction, min_ratio = 1.1) {
  d <- direction_vector(direction)
  e <- eigen(polygon_moments(xy)$cov, symmetric = TRUE)
  if (e$values[1] / e$values[2] < min_ratio) return(list(u = d, isotropic = TRUE))
  dots <- as.numeric(crossprod(e$vectors, d))
  # eigen() lists the major axis first and which.max takes the first maximum,
  # so an exact 45 degree tie resolves to the major axis.
  k <- which.max(round(abs(dots), 9))
  u <- e$vectors[, k] * sign(dots[k])
  list(u = u / sqrt(sum(u^2)), isotropic = FALSE)
}

# Rotation matrix R such that R %*% u = (0, 1).
axis_rotation <- function(u) matrix(c(u[2], u[1], -u[1], u[2]), 2, 2)

# Rotate points so that unit vector u maps to (0, 1).
rotate_axis_up <- function(xy, u) {
  out <- xy %*% t(axis_rotation(u))
  colnames(out) <- c("x", "y")
  out
}

# Crossings of the vertical line x = x0 with a closed polygon. An edge from
# vertex i to i+1 crosses when x_i <= x0 < x_{i+1} or x_{i+1} <= x0 < x_i, so a
# vertex on the line is counted once.
# Returns a data frame: edge (index of the edge's first vertex), y.
vertical_crossings <- function(xy, x0) {
  n <- nrow(xy)
  x <- xy[, 1]; y <- xy[, 2]
  j <- c(2:n, 1L)
  i <- which((x <= x0 & x[j] > x0) | (x > x0 & x[j] <= x0))
  t <- (x0 - x[i]) / (x[j[i]] - x[i])
  data.frame(edge = i, y = y[i] + t * (y[j[i]] - y[i]))
}

# Points at given arc-length fractions along an open polyline.
resample_path <- function(path, fractions) {
  s <- c(0, cumsum(sqrt(rowSums(diff(path)^2))))
  target <- fractions * s[length(s)]
  cbind(x = stats::approx(s, path[, 1], target, ties = "ordered")$y,
        y = stats::approx(s, path[, 2], target, ties = "ordered")$y)
}

# Position of each of the N_POINTS points from tip (0) to base (1).
point_positions <- function() {
  c((0:(N_EDGE - 1L)) / N_EDGE, 1, ((N_EDGE - 1L):1) / N_EDGE)
}

# Index of the point at the same tip-to-base position on the opposite edge.
mirror_index <- function() {
  out <- N_POINTS + 2L - seq_len(N_POINTS)
  out[TIP_INDEX] <- TIP_INDEX
  out
}

# Corresponding points from an oriented outline whose first row is the tip and
# whose row `base_row` is the base: each directed edge (tip to base, base to
# tip) sampled at arc-length fractions j / N_EDGE, j = 0 .. N_EDGE - 1.
points_from_outline <- function(outline, base_row) {
  n <- nrow(outline)
  if (base_row < 2L || base_row > n) stop("The base row is outside the outline.")
  left <- outline[1:base_row, , drop = FALSE]
  right <- rbind(outline[base_row:n, , drop = FALSE], outline[1, , drop = FALSE])
  fr <- (0:(N_EDGE - 1L)) / N_EDGE
  pts <- rbind(resample_path(left, fr), resample_path(right, fr))
  rownames(pts) <- NULL
  pts
}

# Cross-section width of an oriented outline at height y: the distance between
# the outermost crossings of the horizontal line, spanning any notch.
cross_section_width <- function(outline, y) {
  cr <- vertical_crossings(outline[, 2:1, drop = FALSE], y)
  if (nrow(cr) < 2L) 0 else diff(range(cr$y))
}

# TRUE when the cross-section 5 percent of the length below the tip is wider
# than 35 percent of the maximum width.
blunt_tip <- function(outline, depth = 0.05, share = 0.35) {
  yr <- range(outline[, 2])
  cross_section_width(outline, yr[2] - depth * diff(yr)) > share * diff(range(outline[, 1]))
}

# Standardize a prepared outline.
#   xy:        prepared outline (y up, counterclockwise)
#   direction: declared tip direction ("up", "down", "left", "right")
#   tip, base: optional manual landmarks as c(x, y) in the same coordinates;
#              each is snapped to the nearest outline vertex
# Returns list(points (N_POINTS x 2, tip first), outline (oriented contour, tip
# up, tip first, base at row base_row), base_row, tip, base, tip_input and
# base_input (landmarks in the input coordinates), axis, warnings).
standardize_outline <- function(xy, direction = "up", tip = NULL, base = NULL) {
  warnings <- character(0)
  manual <- !is.null(tip) && !is.null(base)
  n <- nrow(xy)

  if (manual) {
    nearest <- function(p) which.min((xy[, 1] - p[1])^2 + (xy[, 2] - p[2])^2)
    ti <- nearest(as.numeric(tip)); bi <- nearest(as.numeric(base))
    if (ti == bi) stop("The tip and base landmarks fall on the same outline point. Click two different places.")
    u <- xy[ti, ] - xy[bi, ]
    u <- u / sqrt(sum(u^2))
    rot <- rotate_axis_up(xy, u)
    rot <- rot[c(ti:n, seq_len(ti - 1L)), , drop = FALSE]
    bi <- ((bi - ti) %% n) + 1L
    left <- rot[1:bi, , drop = FALSE]
    right <- rbind(rot[bi:n, , drop = FALSE], rot[1, , drop = FALSE])
  } else {
    ax <- axis_direction(xy, direction)
    u <- ax$u
    if (ax$isotropic)
      warnings <- c(warnings, "This outline is nearly as wide as it is long, so its long axis could not be determined. The declared tip direction was used as the axis. Check the tip and base markers, and set them by hand if they are misplaced.")
    rot <- rotate_axis_up(xy, u)
    cx <- polygon_moments(rot)$centroid[1]
    cand <- which(rot[, 2] >= max(rot[, 2]) - 1e-9)
    if (length(cand) > 1L) {
      # Several vertices share the maximum: nearest the midline, then smaller x.
      off <- round(abs(rot[cand, 1] - cx), 9)
      cand <- cand[off == min(off)]
      cand <- cand[which.min(rot[cand, 1])]
    }
    rot <- rot[c(cand:n, seq_len(cand - 1L)), , drop = FALSE]
    cr <- vertical_crossings(rot, cx)
    if (nrow(cr) < 2L) stop("The midline does not cross the outline.")
    k <- which.min(cr$y)
    e <- cr$edge[k]
    bpt <- c(x = cx, y = cr$y[k])
    if (nrow(cr) > 2L)
      warnings <- c(warnings, "The midline crosses the outline more than twice (for example through a basal notch or a barb). The base marker was placed at the crossing farthest from the tip. Check it, and set the tip and base by hand if it is misplaced.")
    # The half-open crossing rule means the base point can coincide with the
    # edge's first vertex but never with its second.
    on_vertex <- sqrt(sum((rot[e, ] - bpt)^2)) < 1e-12
    left <- if (on_vertex) rot[1:e, , drop = FALSE] else rbind(rot[1:e, , drop = FALSE], bpt)
    rest <- if (e < n) rot[(e + 1L):n, , drop = FALSE] else rot[0, , drop = FALSE]
    right <- rbind(bpt, rest, rot[1, , drop = FALSE])
  }

  outline <- rbind(left, right[-c(1L, nrow(right)), , drop = FALSE])
  rownames(outline) <- NULL
  base_row <- nrow(left)
  if (!manual && blunt_tip(outline))
    warnings <- c(warnings, "The tip end of this outline is blunt or broken. The tip marker is at the point farthest along the axis, which may not be where the original tip was. You can set the tip and base by hand.")

  pts <- points_from_outline(outline, base_row)
  back <- function(p) as.numeric(matrix(p, 1) %*% axis_rotation(u))
  list(points = pts, outline = outline, base_row = base_row,
       tip = pts[TIP_INDEX, ], base = pts[BASE_INDEX, ],
       tip_input = back(pts[TIP_INDEX, ]), base_input = back(pts[BASE_INDEX, ]),
       axis = u, warnings = warnings)
}

# Cross-section width at `n_levels` heights: the midpoints of equal slices of
# the length. `height` is the fraction of length above the lowest point.
width_profile <- function(outline, n_levels = 200L) {
  yr <- range(outline[, 2])
  levels <- yr[1] + (seq_len(n_levels) - 0.5) / n_levels * diff(yr)
  data.frame(height = (levels - yr[1]) / diff(yr),
             width = vapply(levels, function(l) cross_section_width(outline, l), numeric(1)))
}

# Linear measurements of an oriented outline (tip up), in the outline's pixel
# units and in millimetres when px_per_mm is given.
#   length: extent along the axis; width: extent at right angles to it;
#   widest_position: height of the widest cross-section as a fraction of length
#   above the lowest point (ties go to the level nearest the base).
measure_outline <- function(outline, px_per_mm = NA_real_) {
  len <- diff(range(outline[, 2]))
  wid <- diff(range(outline[, 1]))
  prof <- width_profile(outline)
  area <- abs(polygon_moments(outline)$area)
  s <- if (length(px_per_mm) == 1L && is.finite(px_per_mm) && px_per_mm > 0) px_per_mm else NA_real_
  data.frame(
    length_px = len, width_px = wid, area_px2 = area,
    length_mm = len / s, width_mm = wid / s, area_mm2 = area / s^2,
    length_width_ratio = len / wid,
    widest_position = prof$height[which.max(round(prof$width, 9))]
  )
}
