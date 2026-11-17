test_that("alignment matches geomorph::gpagen reference coordinates up to one rotation", {
  g <- gpa(REF$points)
  mine <- simplify2array(g$aligned)
  ref_mean <- apply(REF$gpagen, c(1, 2), mean)
  s <- svd(crossprod(g$mean, ref_mean))
  R <- s$u %*% diag(c(1, sign(det(s$u %*% t(s$v))))) %*% t(s$v)
  for (i in seq_along(g$aligned))
    expect_lt(max(abs(g$aligned[[i]] %*% R - REF$gpagen[, , i])), 1e-5)
  expect_equal(unname(vapply(g$aligned, centroid_size, numeric(1))), rep(1, 4), tolerance = 1e-10)
  expect_equal(unname(g$mean[1, 1] - g$mean[101, 1]), 0, tolerance = 1e-10)   # mean tip-base axis is vertical
  expect_gt(g$mean[1, 2], g$mean[101, 2])
  expect_equal(g$aligned, REF$aligned, tolerance = 1e-10)
})

test_that("PCA scores equal prcomp of the Momocs coefficients up to sign", {
  p <- shape_pca(efa_matrix(REF$aligned, 20))
  q <- prcomp(REF$momocs_coe)
  expect_identical(p$rank, 3L)
  expect_equal(abs(unname(p$scores)), abs(unname(q$x[, 1:3])), tolerance = 1e-8, scale = 1)
  expect_equal(sum(p$explained), 1)
})

test_that("components beyond the data's rank are not reported", {
  set.seed(1)
  expect_identical(shape_pca(matrix(rnorm(30), 3))$rank, 2L)
  same <- matrix(rep(1:5, each = 4), 4)
  expect_identical(shape_pca(same)$rank, 0L)
  line <- outer(c(-1, 0, 1, 2), 1:6)
  expect_identical(shape_pca(line)$rank, 1L)
  expect_null(pca_stability(matrix(rnorm(30), 3)))
  expect_null(pca_stability(line))
})

test_that("symmetric shapes have no asymmetry and unit size", {
  pair <- mirror_index()
  for (p in REF$points) {
    s <- symmetrize_shape(p, pair)
    expect_lt(asymmetry_score(s, pair), 1e-10)
    expect_equal(centroid_size(s), 1, tolerance = 1e-12)
    expect_gt(asymmetry_score(p, pair), 0.01)
  }
})

test_that("plane agreement is unaffected by an exchange of component order", {
  x <- rbind(c(1.01, 0), c(-1.01, 0), c(0, 1), c(0, -1))
  rownames(x) <- letters[1:4]
  st <- pca_stability(x)
  expect_equal(unname(st$plane), rep(1, 4), tolerance = 1e-8)
  expect_identical(unname(st$separable), c(FALSE, FALSE))
  expect_equal(unname(st$axis["a", ]), c(0, 0), tolerance = 1e-8)   # same-rank cosines mislead here
  y <- cbind(c(-3, -1, 1, 3, 0, 0), c(0, 0.2, -0.2, 0, 0.5, -0.5), c(0.01, -0.01, 0.01, -0.01, 0.02, -0.02))
  expect_true(pca_stability(y)$separable[["PC1"]])
})

test_that("silhouette widths equal cluster::silhouette", {
  skip_if_not_installed("cluster")
  set.seed(3)
  x <- matrix(rnorm(60), 12); d <- dist(x); tree <- shape_clusters(x)
  ref <- function(cl, d) mean(cluster::silhouette(cl, d)[, "sil_width"])
  for (k in 2:5) expect_equal(mean_silhouette(d, cutree(tree, k)), ref(cutree(tree, k), d), tolerance = 1e-10)
  single <- c(1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 3)                       # cluster 3 is a singleton
  expect_equal(mean_silhouette(d, single), ref(single, d), tolerance = 1e-10)
  dd <- dist(c(0, 0, 0, 3, 3, 7)); cl <- c(1, 1, 1, 2, 2, 3)            # duplicates
  expect_equal(mean_silhouette(dd, cl), ref(cl, dd), tolerance = 1e-10)
  expect_equal(mean_silhouette(dist(c(0, 0, 0, 0)), c(1, 1, 2, 2)), 0)  # a = b = 0
  s <- silhouette_by_k(x, tree)
  expect_identical(s$k, 2:11)
  expect_equal(nrow(silhouette_by_k(x[1:3, ])), 0L)
})

test_that("regions partition the points", {
  reg <- assign_regions(point_positions(), c(0.80, 0.15, 0.60))        # unsorted input is sorted
  expect_identical(as.integer(table(reg)), c(29L, 90L, 40L, 41L))
  expect_identical(as.character(reg[c(1, 15, 16, 101, 187, 188)]), c("Tip", "Tip", "Blade", "Base", "Tip", "Tip"))
  expect_identical(as.character(assign_regions(c(0.15, 0.60, 0.80, 1))), c("Blade", "Shoulder", "Base", "Base"))
})

test_that("the RV coefficient behaves as specified", {
  set.seed(2)
  x <- matrix(rnorm(40), 8); y <- matrix(rnorm(24), 8)
  expect_equal(rv_coefficient(x, x), 1, tolerance = 1e-12)
  xc <- scale(x, scale = FALSE); yc <- scale(y, scale = FALSE)
  sxy <- crossprod(xc, yc); sxx <- crossprod(xc); syy <- crossprod(yc)
  expect_equal(rv_coefficient(x, y), sum(diag(sxy %*% t(sxy))) / sqrt(sum(diag(sxx %*% sxx)) * sum(diag(syy %*% syy))), tolerance = 1e-12)
  expect_true(is.na(rv_coefficient(x, matrix(1, 8, 3))))
  expect_equal(rv_coefficient(x, x %*% diag(c(2, 2, 2, 2, 2))), 1, tolerance = 1e-12)
})

test_that("the region reference is reproducible and leaves the random state alone", {
  al <- REF$aligned; reg <- assign_regions(point_positions())
  set.seed(99); before <- .Random.seed
  a <- region_rv(al, reg, n_perm = 99L)
  expect_identical(.Random.seed, before)
  expect_identical(a, region_rv(al, reg, n_perm = 99L))
  expect_equal(nrow(a), 6L)
  expect_true(all(a$reference_low <= a$reference_high))
  expect_true(all(a$rv >= 0 & a$rv <= 1))
  pv <- point_variance(al)
  expect_length(pv, 200L); expect_true(all(pv >= 0))
  expect_equal(pv[5], var(sapply(al, function(s) s[5, 1])) + var(sapply(al, function(s) s[5, 2])))
})
