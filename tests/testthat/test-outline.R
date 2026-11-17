test_that("the result does not depend on rotation or on the contour's start vertex", {
  d <- densify(KITE_ASYM, 12L)                       # 84 vertices
  base <- standardize_yup(d, "up")
  cases <- list(c(0, "up"), c(37, "up"), c(90, "left"), c(180, "down"), c(265, "right"))
  for (cs in cases) for (k in c(0, 7, 23)) {
    r <- standardize_yup(shift_start(rotate_deg(d, as.numeric(cs[1])), k), cs[2])
    expect_equal(r$points, base$points, tolerance = 1e-8, info = paste(cs[1], k))
  }
  expect_equal(standardize_yup(KITE_ASYM, "up")$points, base$points, tolerance = 1e-8)   # added collinear vertices change nothing
})

test_that("tip and base are placed where the geometry puts them", {
  r <- standardize_yup(KITE, "up")
  expect_equal(unname(r$tip_input), c(0, 10), tolerance = 1e-8)
  expect_equal(unname(r$base_input), c(0, -2), tolerance = 1e-8)
  expect_length(r$warnings, 0)
})

test_that("an exact 45 degree tie resolves to the major axis", {
  r <- standardize_yup(rotate_deg(KITE, -45), "up")
  expect_equal(abs(r$axis), c(1, 1) / sqrt(2), tolerance = 1e-6)
  expect_equal(unname(r$tip_input), as.numeric(rotate_deg(matrix(c(0, 10), 1), -45)), tolerance = 1e-6)
})

test_that("a shape wider than long is oriented by the declared direction", {
  r <- standardize_yup(WIDE, "up")
  expect_equal(unname(r$tip_input), c(0, 3), tolerance = 1e-8)
  expect_lt(measure_outline(r$outline)$length_width_ratio, 1)
})

test_that("a near-circular outline triggers the axis warning", {
  th <- 2 * pi * (0:59) / 60
  r <- standardize_yup(cbind(cos(th), sin(th)), "up")
  expect_true(any(grepl("long axis could not be determined", r$warnings)))
})

test_that("a basal notch on the midline receives the base landmark", {
  r <- standardize_yup(NOTCHED, "up")
  expect_equal(unname(r$base_input), c(0, -1), tolerance = 1e-8)
  expect_false(any(grepl("more than two points", r$warnings)))
})

test_that("a midline that crosses a barb gives the lowest crossing and a warning", {
  r <- standardize_yup(BARBED, "up")
  expect_equal(unname(r$base_input[2]), -6, tolerance = 1e-8)
  expect_true(any(grepl("more than two points", r$warnings)))
})

test_that("a truncated tip triggers the blunt-tip warning and the demo outlines do not", {
  r <- standardize_yup(TRUNCATED, "up")
  expect_true(any(grepl("blunt or broken", r$warnings)))
  expect_equal(unname(r$tip_input), c(-1.5, 10), tolerance = 1e-8)   # tie: equal distance, smaller x
  for (f in DEMO_FILES) expect_length(process_image(load_working_image(f), direction = "left")$warnings, 0)
})

test_that("corresponding points are 200 distinct points with T at 1 and B at 101", {
  r <- standardize_yup(KITE_ASYM, "up")
  p <- r$points
  expect_equal(dim(p), c(200L, 2L))
  expect_equal(nrow(unique(round(p, 9))), 200L)
  expect_equal(unname(p[1, ]), unname(r$tip)); expect_equal(unname(p[101, ]), unname(r$base))
  same <- function(q) which(sqrt(rowSums(sweep(p, 2, q)^2)) < 1e-9)
  expect_identical(same(r$tip), 1L); expect_identical(same(r$base), 101L)
  seg <- sqrt(rowSums((p[c(2:200, 1), ] - p)^2))
  expect_true(all(seg > 0))
  # Spacing is constant in distance along the edge (the edges are polylines,
  # so straight-line spacing changes at corners). The position of each point
  # along its edge is found independently of resample_path().
  along <- function(path, q) {
    s0 <- c(0, cumsum(sqrt(rowSums(diff(path)^2))))
    best <- c(Inf, NA)
    for (i in seq_len(nrow(path) - 1L)) {
      v <- path[i + 1L, ] - path[i, ]
      t <- min(max(sum((q - path[i, ]) * v) / sum(v^2), 0), 1)
      dist <- sqrt(sum((path[i, ] + t * v - q)^2))
      if (dist < best[1] - 1e-12) best <- c(dist, s0[i] + t * sqrt(sum(v^2)))
    }
    best
  }
  edges <- list(r$outline[1:r$base_row, ], rbind(r$outline[r$base_row:nrow(r$outline), ], r$outline[1, ]))
  for (e in 1:2) {
    a <- t(apply(p[(e - 1) * 100 + 1:100, ], 1, function(q) along(edges[[e]], q)))
    expect_lt(max(a[, 1]), 1e-8)                               # every point lies on its edge
    total <- sum(sqrt(rowSums(diff(edges[[e]])^2)))
    expect_lt(max(abs(diff(c(a[, 2], total)) - total / 100)), 1e-8)
  }
  reg <- assign_regions(point_positions())
  expect_false(anyNA(reg)); expect_length(reg, 200L)
  expect_equal(point_positions()[c(1, 2, 100, 101, 102, 200)], c(0, 0.01, 0.99, 1, 0.99, 0.01))
  expect_equal(mirror_index()[c(1, 2, 100, 101, 102, 200)], c(1L, 200L, 102L, 101L, 100L, 2L))
  expect_equal(point_positions()[mirror_index()], point_positions())
})

test_that("arc-length spacing is constant within each edge on a smooth outline", {
  th <- 2 * pi * (0:1999) / 2000
  r <- standardize_yup(cbind(20 * cos(th), 50 * sin(th)), "up")
  seg <- sqrt(rowSums(diff(r$points[1:100, ])^2))
  expect_lt(diff(range(seg)) / mean(seg), 1e-3)
})

test_that("manual landmarks are snapped to vertices and rejected when they coincide", {
  r <- standardize_yup(KITE_ASYM, "up", tip = c(0.1, -9.8), base = c(-1.1, 2.1))   # y-up coords given to prepared outline
  expect_error(standardize_outline(prepare_outline(cbind(KITE[, 1], -KITE[, 2])), "up", tip = c(0, 10), base = c(0.01, 10)),
               "same outline point")
  m <- standardize_outline(prepare_outline(cbind(KITE[, 1], -KITE[, 2])), "left", tip = c(0.2, 9.9), base = c(1.2, -2))
  expect_equal(unname(m$tip_input), c(0, 10), tolerance = 1e-8)
  expect_equal(unname(m$base_input), c(1.2, -2), tolerance = 1e-8)
})

test_that("measurements match known shapes", {
  rect <- rbind(c(-20, 100), c(-20, 0), c(20, 0), c(20, 100))
  m <- measure_outline(rect)
  expect_equal(c(m$length_px, m$width_px, m$area_px2), c(100, 40, 4000), tolerance = 1e-8)
  th <- 2 * pi * (0:1999) / 2000
  e <- measure_outline(cbind(20 * cos(th), 50 * sin(th)))
  expect_equal(e$length_px, 100, tolerance = 0.01, scale = 1); expect_equal(e$width_px, 40, tolerance = 0.01, scale = 1)
  expect_equal(e$widest_position, 0.5, tolerance = 0.01, scale = 1)
  expect_lt(measure_outline(rbind(c(0, 10), c(-3, 0), c(3, 0)))$widest_position, 0.01)
  expect_equal(cross_section_width(U_SHAPE, -4), 6, tolerance = 1e-8)      # spans the notch
  mm <- measure_outline(rect, px_per_mm = 4)
  expect_equal(c(mm$length_mm, mm$width_mm, mm$area_mm2), c(25, 10, 250), tolerance = 1e-8)
  expect_true(is.na(measure_outline(rect)$length_mm))
})

test_that("a shared scale transfers between images resized by different factors", {
  ellipse_image <- function(w, h, a, b) {
    m <- outer(1:w, 1:h, function(x, y) ((x - w / 2) / a)^2 + ((y - h / 2) / b)^2 <= 1) * 0.8 + 0.1
    f <- tempfile(fileext = ".png"); imager::save.image(imager::as.cimg(m), f); f
  }
  # Same camera set-up: 20 original pixels per mm in both photographs.
  A <- load_working_image(ellipse_image(1200, 900, 300, 100))     # 30 mm long, not resized
  B <- load_working_image(ellipse_image(2400, 1800, 600, 200))    # 60 mm long, resized by 0.5
  expect_equal(A$factor, 1); expect_equal(B$factor, 0.5)
  len <- function(img) measure_outline(process_image(img, direction = "left")$std$outline)$length_px
  shared <- ppm_to_original(20, A$factor)                         # calibrated on A: 20 working px per mm
  expect_equal(len(A) / 20, 30, tolerance = 0.02)
  expect_equal(len(B) / ppm_to_working(shared, B$factor), 60, tolerance = 0.02)
  expect_equal(len(B) / 10, 60, tolerance = 0.02)                 # B's own two-click scale: 10 working px per mm
})

test_that("manual landmarks at the automatic positions give the same points and the same edge order", {
  kite <- rbind(KITE[1:4, ], c(0, -2), KITE[5:7, ])            # a vertex at the base position
  prep <- prepare_outline(cbind(kite[, 1], -kite[, 2]))
  auto <- standardize_outline(prep, "up")
  expect_equal(unname(auto$base_input), c(0, -2), tolerance = 1e-10)
  man <- standardize_outline(prep, "left", tip = auto$tip_input, base = auto$base_input)
  expect_equal(man$points, auto$points, tolerance = 1e-8)
  expect_identical(man$base_row, auto$base_row)
  # Edge 1 (points 2 to 100) is the left edge when the tip is up.
  for (r in list(auto, man)) { expect_lt(max(r$points[2:100, 1]), 0); expect_gt(min(r$points[102:200, 1]), 0) }
  asym <- standardize_yup(KITE_ASYM, "up")
  expect_lt(mean(asym$points[2:100, 1]), mean(asym$points[102:200, 1]))
})

test_that("a base point on either vertex of its edge is not duplicated", {
  no_zero <- function(r) all(sqrt(rowSums(diff(rbind(r$outline, r$outline[1, ]))^2)) > 1e-9)
  # Base on the second vertex of an edge with decreasing x (area centroid at 0, 0).
  second <- rbind(c(-sqrt(11), 1), c(-sqrt(11), -1), c(1, -1), c(1, -3), c(0, -4), c(2, -3), c(2, 3), c(0, 4), c(1, 3), c(1, 1))
  expect_equal(polygon_moments(second)$centroid, c(0, 0), tolerance = 1e-12)
  r <- standardize_yup(second, "up")
  expect_true(no_zero(r)); expect_equal(nrow(r$outline), nrow(second))
  expect_equal(unname(r$base), c(0, -4), tolerance = 1e-10); expect_equal(unname(r$outline[r$base_row, ]), c(0, -4), tolerance = 1e-10)
  # Its mirror image puts the base on the first vertex of an edge with increasing x.
  first <- cbind(-second[, 1], second[, 2])[nrow(second):1, ]
  m <- standardize_yup(first, "up")
  expect_true(no_zero(m)); expect_equal(nrow(m$outline), nrow(first)); expect_equal(unname(m$base), c(0, -4), tolerance = 1e-10)
  # The two are mirror images, so their points correspond under the mirror pairing.
  expect_equal(unname(cbind(-m$points[mirror_index(), 1], m$points[mirror_index(), 2])), unname(r$points), tolerance = 1e-8)
  for (k in c(2, 6)) expect_equal(standardize_yup(shift_start(second, k), "up")$points, r$points, tolerance = 1e-8)
  # Base inside an edge: one vertex is inserted.
  k <- standardize_yup(KITE, "up")
  expect_true(no_zero(k)); expect_equal(nrow(k$outline), nrow(KITE) + 1L)
})

test_that("the points where an outline meets a line are the same for a shape and its mirror image", {
  right <- rbind(c(0, 0), c(2, -1), c(2, 3), c(0, 2), c(1, 1))          # (0, 0) and (0, 2) touch x = 0 from the right
  left <- cbind(-right[, 1], right[, 2])[5:1, ]                          # the same vertices touch from the left
  for (p in list(right, left)) expect_equal(sort(vertical_crossings(p, 0)$y), c(0, 2))
  through <- rbind(c(-1, 0), c(0, -1), c(1, 0), c(0, 1))                 # the outline passes through x = 0 at two vertices
  cr <- vertical_crossings(through, 0)
  expect_equal(sort(cr$y), c(-1, 1)); expect_false(anyNA(cr$vertex))
  inside <- vertical_crossings(rbind(c(-1, -1), c(1, -1), c(1, 1), c(-1, 1)), 0.25)
  expect_equal(sort(inside$y), c(-1, 1)); expect_true(all(is.na(inside$vertex)))
})
