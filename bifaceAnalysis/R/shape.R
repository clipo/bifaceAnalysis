# Shape analysis on corresponding points: Procrustes alignment, symmetry, PCA,
# clustering, and region statistics. Configurations are k x 2 matrices.
# Specification: docs/specs/2026-11-17_app_rewrite_spec.md, section 7.

centroid_size <- function(m) sqrt(sum(sweep(m, 2, colMeans(m))^2))

# Centre a configuration and scale it to unit centroid size.
center_scale <- function(m) {
  m <- sweep(m, 2, colMeans(m))
  m / sqrt(sum(m^2))
}

# Rotation (no reflection) of x that minimizes the summed squared distance to
# y. Both must be centred. Returns the rotated x.
rotate_onto <- function(x, y) {
  s <- svd(crossprod(x, y))
  d <- sign(det(s$u %*% t(s$v)))
  x %*% (s$u %*% diag(c(1, d)) %*% t(s$v))
}

# Partial Procrustes distance between two configurations.
procrustes_distance <- function(x, y) {
  x <- center_scale(x); y <- center_scale(y)
  sqrt(sum((rotate_onto(x, y) - y)^2))
}

# Generalized Procrustes alignment (translation, unit centroid size, rotation)
# of a named list of k x 2 matrices. The aligned set is rotated as a whole so
# that the mean shape's base-to-tip direction points up.
# Returns list(aligned, mean, iterations).
gpa <- function(shapes, tip = 1L, base = nrow(shapes[[1]]) %/% 2L + 1L,
                tol = 1e-10, max_iter = 100L) {
  x <- lapply(shapes, center_scale)
  ref <- x[[1]]
  ss_old <- Inf
  for (iter in seq_len(max_iter)) {
    x <- lapply(x, rotate_onto, y = ref)
    m <- Reduce(`+`, x) / length(x)
    m <- m / sqrt(sum(m^2))
    ss <- sum(vapply(x, function(s) sum((s - m)^2), numeric(1)))
    ref <- m
    if (abs(ss_old - ss) < tol) break
    ss_old <- ss
  }
  u <- m[tip, ] - m[base, ]
  u <- u / sqrt(sum(u^2))
  R <- matrix(c(u[2], u[1], -u[1], u[2]), 2, 2)
  spin <- function(s) { o <- s %*% t(R); colnames(o) <- c("x", "y"); o }
  list(aligned = lapply(x, spin), mean = spin(m), iterations = iter)
}

# Mirror image of a configuration about the vertical axis, with points
# relabelled by `pair` so that each keeps its tip-to-base position.
mirror_shape <- function(m, pair) {
  r <- cbind(-m[, 1], m[, 2])[pair, , drop = FALSE]
  colnames(r) <- colnames(m)
  r
}

# Partial Procrustes distance between a configuration and its relabelled
# mirror image.
asymmetry_score <- function(m, pair) procrustes_distance(mirror_shape(m, pair), m)

# Average of the centred unit-size configuration and its relabelled mirror
# image rotated onto it, rescaled to unit centroid size.
symmetrize_shape <- function(m, pair) {
  x <- center_scale(m)
  r <- rotate_onto(center_scale(mirror_shape(m, pair)), x)
  s <- (x + r) / 2
  colnames(s) <- c("x", "y")
  s / sqrt(sum(s^2))
}

# PCA of a specimens-by-variables matrix (centred, unscaled). Only components
# whose variance exceeds 1e-12 of the total are kept. `rank` is their number
# and is 0 when all rows are identical.
shape_pca <- function(x) {
  p <- stats::prcomp(x, center = TRUE, scale. = FALSE)
  var <- p$sdev^2
  total <- sum(var)
  keep <- if (total > 0) which(var > 1e-12 * total) else integer(0)
  keep <- keep[keep <= nrow(x) - 1L]
  list(scores = p$x[, keep, drop = FALSE], loadings = p$rotation[, keep, drop = FALSE],
       center = p$center, sdev = p$sdev[keep],
       explained = if (total > 0) var[keep] / total else numeric(0), rank = length(keep))
}

# Leave-one-out stability of the first two components. Each row of x is left
# out in turn and the PCA refitted on the remaining rows; x is not recomputed.
# Returns NULL with fewer than 4 rows or fewer than 2 components. Otherwise:
#   plane:     product of the cosines of the principal angles between the
#              full-data and refitted PC1-PC2 planes. It is unchanged when the
#              two components exchange order or rotate within the plane.
#   axis:      |cosine| between full-data and refitted loading vectors of the
#              same rank (columns PC1, PC2).
#   separable: whether each component's variance is at least `min_ratio` times
#              the next one's and, for PC2, at most 1 / min_ratio of PC1's.
#              Direction agreement is meaningful only for separable components.
pca_stability <- function(x, min_ratio = 1.5) {
  n <- nrow(x)
  full <- shape_pca(x)
  if (n < 4L || full$rank < 2L) return(NULL)
  V <- full$loadings[, 1:2, drop = FALSE]
  axis <- matrix(NA_real_, n, 2L, dimnames = list(rownames(x), c("PC1", "PC2")))
  plane <- stats::setNames(rep(NA_real_, n), rownames(x))
  for (i in seq_len(n)) {
    fit <- shape_pca(x[-i, , drop = FALSE])
    if (fit$rank < 2L) next
    W <- fit$loadings[, 1:2, drop = FALSE]
    axis[i, ] <- abs(colSums(V * W))
    plane[i] <- prod(pmin(svd(crossprod(V, W))$d, 1))
  }
  v <- full$sdev^2
  above_next <- function(j) j == length(v) || v[j] >= min_ratio * v[j + 1L]
  list(axis = axis, plane = plane,
       separable = c(PC1 = above_next(1L), PC2 = v[1] >= min_ratio * v[2] && above_next(2L)))
}

# Ward hierarchical clustering on Euclidean distances.
shape_clusters <- function(x) stats::hclust(stats::dist(x), method = "ward.D2")

# Mean silhouette width of a partition, from a distance object or matrix. For
# specimen i, with a the mean distance to the other members of its cluster and
# b the smallest mean distance to another cluster, the width is
# (b - a) / max(a, b). It is 0 for a specimen alone in its cluster and 0 when
# a = b = 0.
mean_silhouette <- function(d, cl) {
  d <- as.matrix(d)
  ids <- sort(unique(cl))
  if (length(ids) < 2L) return(NA_real_)
  s <- vapply(seq_along(cl), function(i) {
    own <- cl == cl[i]
    if (sum(own) == 1L) return(0)
    a <- sum(d[i, own]) / (sum(own) - 1L)
    b <- min(vapply(setdiff(ids, cl[i]), function(k) mean(d[i, cl == k]), numeric(1)))
    if (max(a, b) == 0) 0 else (b - a) / max(a, b)
  }, numeric(1))
  mean(s)
}

# Mean silhouette width for k = 2 .. n - 1 clusters cut from `tree`.
silhouette_by_k <- function(x, tree = shape_clusters(x)) {
  n <- nrow(x)
  if (n < 4L) return(data.frame(k = integer(0), silhouette = numeric(0)))
  d <- stats::dist(x)
  ks <- 2:(n - 1L)
  data.frame(k = ks, silhouette = vapply(ks, function(k) mean_silhouette(d, stats::cutree(tree, k)), numeric(1)))
}

# Region of each point from its tip-to-base position and three boundaries:
# tip [0, b1), blade [b1, b2), shoulder [b2, b3), base [b3, 1].
REGION_NAMES <- c("Tip", "Blade", "Shoulder", "Base")
assign_regions <- function(positions, boundaries = c(0.15, 0.60, 0.80)) {
  factor(REGION_NAMES[findInterval(positions, sort(boundaries)) + 1L], levels = REGION_NAMES)
}

# Variance of aligned coordinates at each point (x variance + y variance).
point_variance <- function(aligned) {
  arr <- simplify2array(aligned)                # k x 2 x n
  apply(arr, 1, function(p) sum(apply(p, 1, stats::var)))
}

# Specimens-by-(2 * points) matrix of the coordinates of the points in `idx`.
region_block <- function(aligned, idx) {
  t(vapply(aligned, function(s) as.numeric(s[idx, , drop = FALSE]), numeric(2L * length(idx))))
}

# RV coefficient between two specimens-by-variables blocks:
# tr(XX'YY') / sqrt(tr(XX'XX') tr(YY'YY')) for column-centred X and Y.
# NA when either block has no variance.
rv_coefficient <- function(x, y) {
  x <- scale(x, scale = FALSE); y <- scale(y, scale = FALSE)
  kx <- tcrossprod(x); ky <- tcrossprod(y)
  den <- sqrt(sum(kx^2) * sum(ky^2))
  if (!is.finite(den) || den < 1e-24) return(NA_real_)
  sum(kx * ky) / den
}

# Observed RV for each pair of regions with a reference: the mean and central
# 95 percent range of RV when the rows of one block are re-paired at random
# with the rows of the other, holding the aligned coordinates fixed. The
# caller's random number state is restored.
region_rv <- function(aligned, regions, n_perm = 999L, seed = 1L) {
  present <- levels(droplevels(regions))
  pairs <- utils::combn(present, 2, simplify = FALSE)
  n <- length(aligned)
  had_seed <- exists(".Random.seed", envir = globalenv())
  old <- if (had_seed) get(".Random.seed", envir = globalenv()) else NULL
  on.exit(if (had_seed) assign(".Random.seed", old, envir = globalenv())
          else if (exists(".Random.seed", envir = globalenv())) rm(".Random.seed", envir = globalenv()))
  set.seed(seed)
  do.call(rbind, lapply(pairs, function(p) {
    x <- region_block(aligned, which(regions == p[1]))
    y <- region_block(aligned, which(regions == p[2]))
    obs <- rv_coefficient(x, y)
    perm <- if (is.na(obs)) NA_real_ else
      vapply(seq_len(n_perm), function(i) rv_coefficient(x, y[sample.int(n), , drop = FALSE]), numeric(1))
    data.frame(region_a = p[1], region_b = p[2], rv = obs,
               reference_mean = mean(perm),
               reference_low = unname(stats::quantile(perm, 0.025, na.rm = TRUE)),
               reference_high = unname(stats::quantile(perm, 0.975, na.rm = TRUE)),
               stringsAsFactors = FALSE)
  }))
}
