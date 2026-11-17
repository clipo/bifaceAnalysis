test_that("(i) forward coefficients match Momocs 1.4.1 reference values", {
  mine <- efa_matrix(REF$aligned, 20)
  expect_equal(unname(mine), unname(REF$momocs_coe), tolerance = 1e-8, scale = 1)
  expect_identical(colnames(mine), c(paste0("a", 1:20), paste0("b", 1:20), paste0("c", 1:20), paste0("d", 1:20)))
  expect_identical(colnames(REF$momocs_coe), c(paste0("A", 1:20), paste0("B", 1:20), paste0("C", 1:20), paste0("D", 1:20)))
})

test_that("(ii) a circle is a single harmonic", {
  th <- 2 * pi * (0:1999) / 2000
  v <- efa_vector(efa_forward(cbind(2 * cos(th), 2 * sin(th)), 10))
  expect_equal(unname(v[c("a1", "d1")]), c(2, 2), tolerance = 1e-4, scale = 1)
  expect_lt(max(abs(v[!names(v) %in% c("a1", "d1")])), 1e-8)
})

test_that("(iii) an ellipse under arc-length parameterization matches Momocs, not its semi-axes", {
  v <- efa_vector(efa_forward(REF$ellipse, 5))
  expect_equal(unname(v), unname(REF$momocs_ellipse), tolerance = 1e-8, scale = 1)
  expect_equal(unname(round(v[c("a1", "a3")], 4)), c(2.6245, 0.2418))
})

test_that("(iv) the inverse reproduces a configuration at its own vertices", {
  p <- REF$aligned[[1]]
  f <- efa_forward(p, 100)
  rec <- efa_inverse(efa_vector(f), x0 = f$x0, y0 = f$y0, at = vertex_fractions(p))
  err <- max(sqrt(rowSums((rec - p)^2)))
  expect_lt(err, 0.005 * centroid_size(p))
  expect_equal(vertex_fractions(p)[1], 0)
  expect_equal(nrow(efa_inverse(efa_vector(f), n_points = 50L)), 50L)
})

test_that("(v) harmonic power uses the a, b, c, d of the same harmonic", {
  m <- REF$momocs_coe
  H <- 20
  direct <- vapply(1:H, function(h) mean((m[, h]^2 + m[, H + h]^2 + m[, 2 * H + h]^2 + m[, 3 * H + h]^2) / 2), numeric(1))
  pw <- harmonic_power(efa_matrix(REF$aligned, H))
  expect_equal(pw$power, direct, tolerance = 1e-10)
  expect_equal(pw$cumulative_all[H], 1); expect_equal(pw$cumulative_excluding_first[H], 1)
  expect_true(is.na(pw$cumulative_excluding_first[1]))
})

test_that("the automatic number of harmonics is bounded and defined without residual power", {
  pw <- harmonic_power(efa_matrix(REF$aligned, 30))
  h <- choose_harmonics(pw)
  expect_gte(h, 6L); expect_lte(h, 30L)
  expect_gte(pw$cumulative_excluding_first[h], 0.99)
  if (h > 6L) expect_lt(pw$cumulative_excluding_first[h - 1L], 0.99)
  zero <- matrix(0, 3, 120); zero[, 1] <- 1; zero[, 91] <- 1      # a1 = d1 = 1 only
  pz <- harmonic_power(zero)
  expect_false(any(is.nan(pz$cumulative_excluding_first[-1])))
  expect_identical(choose_harmonics(pz), 6L)
  smooth <- pw; smooth$cumulative_excluding_first <- c(NA, rep(1, 29))
  expect_identical(choose_harmonics(smooth), 6L)                   # lower bound applies below 6
})

test_that("the constant terms are arc-length means and need not be zero", {
  tri <- rbind(c(0, 0), c(4, 0), c(0, 3)); tri <- sweep(tri, 2, colMeans(tri))
  f <- efa_forward(tri, 3)
  expect_equal(c(f$x0, f$y0), c(1 / 6, 0), tolerance = 1e-10)
  expect_error(efa_forward(rbind(c(0, 0), c(0, 0), c(1, 1)), 2), "repeated")
})
