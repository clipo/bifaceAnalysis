# Elliptical Fourier analysis of closed outlines (Kuhl and Giardina formulas).
# The curve parameter t is cumulative chord length, with t = 0 at the first
# point. Specification: docs/specs/2026-11-17_app_rewrite_spec.md, section 7.2.

# Coefficients for harmonics 1..n_harmonics of a closed polygon (rows = points,
# first point not repeated). Returns list(a, b, c, d, x0, y0, perimeter), where
# x0 and y0 are the arc-length means of x and y.
efa_forward <- function(xy, n_harmonics) {
  K <- nrow(xy)
  if (n_harmonics < 1L) stop("At least one harmonic is required.")
  nxt <- c(2:K, 1L)
  dx <- xy[nxt, 1] - xy[, 1]
  dy <- xy[nxt, 2] - xy[, 2]
  dt <- sqrt(dx^2 + dy^2)
  if (any(dt < 1e-12)) stop("The outline contains repeated consecutive points.")
  t1 <- cumsum(dt)            # parameter at the end of each segment
  t0 <- t1 - dt               # parameter at the start of each segment
  per <- t1[K]
  h <- seq_len(n_harmonics)
  w <- 2 * pi * h / per
  dcos <- cos(outer(t1, w)) - cos(outer(t0, w))     # K x H
  dsin <- sin(outer(t1, w)) - sin(outer(t0, w))
  scale <- per / (2 * pi^2 * h^2)
  list(
    a = scale * colSums((dx / dt) * dcos),
    b = scale * colSums((dx / dt) * dsin),
    c = scale * colSums((dy / dt) * dcos),
    d = scale * colSums((dy / dt) * dsin),
    x0 = sum(dt * (xy[, 1] + xy[nxt, 1]) / 2) / per,
    y0 = sum(dt * (xy[, 2] + xy[nxt, 2]) / 2) / per,
    perimeter = per
  )
}

# Coefficient vector ordered a1..aH, b1..bH, c1..cH, d1..dH.
efa_vector <- function(coef) {
  H <- length(coef$a)
  stats::setNames(c(coef$a, coef$b, coef$c, coef$d),
                  paste0(rep(c("a", "b", "c", "d"), each = H), seq_len(H)))
}

# Specimens-by-4H coefficient matrix from a named list of point matrices.
efa_matrix <- function(shapes, n_harmonics) {
  m <- t(vapply(shapes, function(s) efa_vector(efa_forward(s, n_harmonics)),
                numeric(4L * n_harmonics)))
  rownames(m) <- names(shapes)
  m
}

# Parameter fraction (t / T) of each vertex of a closed polygon, 0 at vertex 1.
vertex_fractions <- function(xy) {
  dt <- sqrt(rowSums((xy[c(2:nrow(xy), 1L), , drop = FALSE] - xy)^2))
  c(0, cumsum(dt)[-nrow(xy)]) / sum(dt)
}

# Outline from a coefficient vector (a1..aH, b1.., c1.., d1..) using the first
# n_harmonics harmonics. Evaluated at n_points equally spaced parameter values,
# or at the parameter fractions `at` (t / T) when given.
efa_inverse <- function(vec, n_points = 200L, n_harmonics = length(vec) / 4, x0 = 0, y0 = 0, at = NULL) {
  H <- length(vec) / 4
  use <- seq_len(min(n_harmonics, H))
  a <- vec[use]; b <- vec[H + use]; cc <- vec[2 * H + use]; d <- vec[3 * H + use]
  theta <- if (is.null(at)) 2 * pi * (seq_len(n_points) - 1) / n_points else 2 * pi * at
  ct <- cos(outer(theta, use)); st <- sin(outer(theta, use))
  cbind(x = x0 + as.numeric(ct %*% a + st %*% b),
        y = y0 + as.numeric(ct %*% cc + st %*% d))
}

# Power of each harmonic, (a^2 + b^2 + c^2 + d^2) / 2, averaged over specimens,
# with cumulative shares of total power including and excluding the first
# harmonic. With no power beyond the first harmonic those shares are zero.
harmonic_power <- function(coef_matrix) {
  H <- ncol(coef_matrix) / 4
  idx <- function(h) c(h, H + h, 2 * H + h, 3 * H + h)
  power <- vapply(seq_len(H), function(h) mean(rowSums(coef_matrix[, idx(h), drop = FALSE]^2) / 2),
                  numeric(1))
  rest <- sum(power[-1])
  share_rest <- if (rest > 0) power[-1] / rest else rep(0, H - 1L)
  data.frame(harmonic = seq_len(H), power = power,
             cumulative_all = cumsum(power / sum(power)),
             cumulative_excluding_first = c(NA_real_, cumsum(share_rest)))
}

# Smallest number of harmonics whose cumulative power, excluding the first
# harmonic, reaches `target`; bounded to [lower, upper]. When no harmonic
# reaches the target (no power beyond the first), `lower` is returned.
choose_harmonics <- function(power_table, target = 0.99, lower = 6L, upper = 30L) {
  h <- which(power_table$cumulative_excluding_first >= target)[1]
  if (is.na(h)) h <- lower
  as.integer(min(max(h, lower), upper, nrow(power_table)))
}
