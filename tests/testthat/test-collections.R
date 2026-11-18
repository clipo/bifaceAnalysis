# Comparison of collections and locations.
# Specification: docs/specs/2026-11-17_collections_spec.md, section 10.

SEED <- 20261117

# Coefficient-like matrix with rows grouped into collections.
sim_groups <- function(sizes, p = 10, shift = 0) {
  g <- factor(rep(LETTERS[seq_along(sizes)], sizes))
  X <- matrix(stats::rnorm(sum(sizes) * p), ncol = p)
  X[, 1] <- X[, 1] + shift * as.integer(g)
  list(X = X, g = g)
}

# ---- labels -----------------------------------------------------------------

test_that("group labels from an outline file are trimmed and kept as text", {
  base <- process_image(load_working_image(DEMO_FILES[["artifact1"]]), direction = "left")
  labels <- c("A", " A ", "  ", "", "NA", "007")
  sp <- lapply(seq_along(labels), function(i)
    new_specimen(paste0("s", i), paste0("s", i), labels[i], base$std$outline * c(1 + i / 50, 1), base$std$base_row))
  f <- tempfile(fileext = ".csv")
  utils::write.csv(specimens_to_csv(sp), f, row.names = FALSE)
  loaded <- specimens_from_csv(f)
  expect_identical(unname(vapply(loaded, function(s) s$group, "")), c("A", "A", "", "", "NA", "007"))
  expect_identical(trim_label(c(" x ", NA, "NA")), c("x", "", "NA"))
})

# ---- partition --------------------------------------------------------------

test_that("the share between collections is SSB / SST and matches a linear model", {
  set.seed(SEED)
  d <- sim_groups(c(3, 5, 8), shift = 1)
  sst <- sum(sweep(d$X, 2, colMeans(d$X))^2)
  ssw <- sum(stats::resid(stats::lm(d$X ~ d$g))^2)
  expect_equal(share_between(d$X, d$g), (sst - ssw) / sst, tolerance = 1e-10)
  m <- collection_means(d$X, d$g)
  ssb <- sum(as.vector(table(d$g)) * rowSums(sweep(m, 2, colMeans(d$X))^2))
  expect_equal(ssb + ssw, sst, tolerance = 1e-10)

  same <- d$X[1:4, ]
  expect_equal(share_between(rbind(same, same), factor(rep(c("A", "B"), each = 4))), 0, tolerance = 1e-12)
  sep <- rbind(matrix(1, 3, 4), matrix(5, 3, 4))
  expect_equal(share_between(sep, factor(rep(c("A", "B"), each = 3))), 1, tolerance = 1e-12)
  expect_true(is.na(share_between(matrix(2, 6, 3), factor(rep(c("A", "B"), each = 3)))))
})

test_that("the share reference matches complete enumeration and restores the random state", {
  set.seed(SEED)
  d <- sim_groups(c(2, 2, 2), p = 4)
  perms <- as.matrix(expand.grid(rep(list(1:6), 6)))
  perms <- perms[apply(perms, 1, function(r) length(unique(r)) == 6L), ]
  labels <- unique(t(apply(perms, 1, function(p) as.character(d$g)[p])))
  expect_identical(nrow(labels), 90L)
  enum <- apply(labels, 1, function(l) share_between(d$X, factor(l)))
  expect_equal(mean(enum), 2 / 5, tolerance = 1e-10)

  set.seed(99); before <- .Random.seed
  ref <- share_reference(d$X, d$g)
  expect_identical(.Random.seed, before)
  expect_true(all(vapply(ref$reference$values, function(v) any(abs(enum - v) < 1e-12), logical(1))))
  expect_lt(abs(ref$reference$mean - mean(enum)), 0.02)
  expect_equal(ref$expected, 2 / 5)
  expect_identical(share_reference(d$X, d$g), ref)
})

# ---- distances --------------------------------------------------------------

test_that("raw and corrected distances match hand calculations", {
  A <- rbind(c(0, 0), c(2, 0))
  expect_equal(pair_distance(A, rbind(c(5, 0), c(7, 0))), c(raw = 5, square = 23, corrected = sqrt(23)))
  expect_equal(pair_distance(A, rbind(c(1, 0), c(1, 2))), c(raw = 1, square = -1, corrected = 0))
  set.seed(SEED)
  d <- sim_groups(c(4, 6), p = 5)
  pr <- pair_table(d$X, d$g, n_perm = 9)
  expect_equal(pr$raw_distance, sqrt(sum((colMeans(d$X[d$g == "A", ]) - colMeans(d$X[d$g == "B", ]))^2)), tolerance = 1e-12)
})

test_that("the corrected squared distance is unbiased and the raw and clipped values are not", {
  set.seed(SEED)
  run <- function(delta) {
    t(replicate(4000, {
      a <- matrix(stats::rnorm(5 * 12), 5); b <- matrix(stats::rnorm(8 * 12), 8); b[, 1] <- b[, 1] + delta
      pd <- pair_distance(a, b)
      c(q = pd[["square"]], raw2 = pd[["raw"]]^2)
    }))
  }
  se <- function(v) stats::sd(v) / sqrt(length(v))
  s2 <- run(2)
  expect_lt(abs(mean(s2[, "q"]) - 4), 3 * se(s2[, "q"]))
  expect_gt(mean(s2[, "raw2"]) - 4, 3 * se(s2[, "raw2"]))
  s0 <- run(0)
  expect_lt(abs(mean(s0[, "q"])), 3 * se(s0[, "q"]))
  clipped <- pmax(0, s0[, "q"])
  expect_gt(mean(clipped), 3 * se(clipped))
})

test_that("the pair reference pools only the two collections of the pair", {
  set.seed(SEED)
  d <- sim_groups(c(4, 5, 3), p = 3)
  d$X[d$g == "C", ] <- d$X[d$g == "C", ] + 1000
  pr <- pair_table(d$X, d$g)
  ab <- pr[pr$collection_a == "A" & pr$collection_b == "B", ]
  pool <- d$X[d$g %in% c("A", "B"), ]
  expect_lte(ab$reference_high, max(stats::dist(pool)))
})

test_that("leave-one-out variance ranges equal hand calculations", {
  X <- rbind(c(0, 0), c(2, 0), c(0, 4), c(3, 3))
  v <- vapply(1:4, function(i) { Y <- X[-i, ]; sum(sweep(Y, 2, colMeans(Y))^2) / 2 }, numeric(1))
  expect_equal(loo_variance_range(X), c(low = min(v), high = max(v)), tolerance = 1e-12)
  expect_equal(row_variance(X), sum(sweep(X, 2, colMeans(X))^2) / 3, tolerance = 1e-12)
  expect_equal(loo_variance_range(matrix(1, 3, 2)), c(low = 0, high = 0))
  wt <- within_table(rbind(matrix(1, 3, 2), X), factor(rep(c("A", "B"), c(3, 4))))
  expect_equal(wt$variance, c(0, row_variance(X)), tolerance = 1e-12)
})

# ---- summaries and shapes ---------------------------------------------------

test_that("collection summaries equal direct calculations and leave absent values empty", {
  g <- factor(c("A", "A", "A", "B", "B", "B"))
  scores <- cbind(PC1 = 1:6, PC2 = c(2, 2, 5, 0, 1, 2))
  meas <- data.frame(length_width_ratio = c(2, 2.5, 3, 1, 1.5, 2), widest_position = c(0.2, 0.3, 0.4, 0.5, 0.5, 0.2),
                     length_mm = c(10, 20, 30, 40, 50, 60))
  s <- collection_summary(g, scores, meas, c(A = 0.1, B = 0.2), located = "B")
  expect_equal(s$pc1_mean, c(2, 5)); expect_equal(s$pc2_mean, c(3, 1)); expect_true(all(is.na(s$pc3_mean)))
  expect_equal(s$length_width_mean, c(2.5, 1.5)); expect_equal(s$max_width_position_mean, c(0.3, 0.4))
  expect_equal(s$length_mm_mean, c(20, 50)); expect_identical(s$located, c(FALSE, TRUE))
  expect_equal(s$within_variance, c(0.1, 0.2))
  meas$length_mm[2] <- NA
  s1 <- collection_summary(g, scores[, 1, drop = FALSE], meas, c(A = 0.1, B = 0.2))
  expect_true(all(is.na(s1$length_mm_mean))); expect_true(all(is.na(s1$pc2_mean)))
  s0 <- collection_summary(g, matrix(0, 6, 0), meas, c(A = 0.1, B = 0.2))
  expect_true(all(is.na(s0$pc1_mean)))
})

test_that("a collection's mean shape equals the mean of its reconstructions", {
  set.seed(SEED)
  co <- matrix(stats::rnorm(5 * 40, sd = 0.05), 5); co[, 1] <- co[, 1] + 1; co[, 31] <- co[, 31] + 1
  m <- efa_inverse(colMeans(co), 200L)
  each <- Reduce(`+`, lapply(seq_len(5), function(i) efa_inverse(co[i, ], 200L))) / 5
  expect_equal(m, each, tolerance = 1e-10)
})

test_that("the variance text covers positive, mixed-zero, and all-zero variances", {
  cd <- list(sizes = table(factor(c("A", "A", "A", "B", "B", "B"))), eligible = c("A", "B"), unlabelled = 0L, small = integer(0))
  sh <- list(share = 0.3, expected = 0.2, reference = list(low = 0.1, high = 0.5))
  pr <- data.frame(collection_a = "A", collection_b = "B", corrected_distance = 0)
  txt <- function(v) paste(text_collections(cd, sh, pr, data.frame(collection = c("A", "B"), variance = v)), collapse = " ")
  expect_match(txt(c(1, 4)), "a ratio of 4.0 to 1")
  expect_match(txt(c(0, 4)), "The specimens of A do not differ")
  expect_match(txt(c(0, 0)), "No collection shows variation")
  expect_match(txt(c(1, 4)), "does not show that the collections have the same mean shape")
})

# ---- geography --------------------------------------------------------------

test_that("great-circle distances match known values", {
  expect_equal(haversine_km(0, 0, 0, 1), 111.195080, tolerance = 1e-6)
  expect_equal(haversine_km(0, 0, 90, 0), 10007.557221, tolerance = 1e-6)
  expect_equal(haversine_km(42.1, -75.9, 42.1, -75.9), 0)
  expect_equal(haversine_km(0, 179.5, 0, -179.5), 111.195080, tolerance = 1e-6)
})

test_that("the shape-geography correlation uses every assignment of positions", {
  set.seed(SEED)
  pos <- cbind(lat = c(40, 41, 42.5, 40.2, 43), lon = c(-75, -76, -74, -77.5, -75.5))
  geo <- geo_distance_matrix(pos[, 1], pos[, 2], LETTERS[1:5])
  sh <- as.matrix(stats::dist(matrix(stats::rnorm(10), 5))); dimnames(sh) <- dimnames(geo)
  res <- shape_geo_correlation(sh, geo, pos)
  expect_true(res$available)
  expect_equal(res$r, stats::cor(sh[lower.tri(sh)], geo[lower.tri(geo)]))
  perms <- as.matrix(expand.grid(rep(list(1:5), 5)))
  perms <- perms[apply(perms, 1, function(r) length(unique(r)) == 5L), ]
  enum <- apply(perms, 1, function(p) stats::cor(sh[lower.tri(sh)], geo[p, p][lower.tri(geo)]))
  expect_identical(res$n_assignments, 120L)
  expect_equal(sort(res$reference$values), sort(enum), tolerance = 1e-12)
  expect_identical(res$n_pairs, 10)

  same <- shape_geo_correlation(geo, geo, pos)
  expect_equal(same$r, 1); expect_equal(max(same$reference$values), 1)
  pos2 <- pos; pos2[2, ] <- pos2[1, ]
  geo2 <- geo_distance_matrix(pos2[, 1], pos2[, 2], LETTERS[1:5])
  expect_length(shape_geo_correlation(sh, geo2, pos2)$reference$values, 120L)

  flat <- matrix(rep(c(40, -75), each = 5), 5)
  expect_false(shape_geo_correlation(sh, geo_distance_matrix(flat[, 1], flat[, 2], LETTERS[1:5]), flat)$available)
  expect_false(shape_geo_correlation(sh[1:4, 1:4], geo[1:4, 1:4], pos[1:4, ])$available)
})

# ---- locations file ---------------------------------------------------------

write_lines_csv <- function(lines) { f <- tempfile(fileext = ".csv"); writeLines(lines, f); f }

test_that("a locations file is read with flexible headers and literal names", {
  loc <- read_locations(write_lines_csv(c(" Collection ,LATITUDE, longitude ,notes", " NA ,42.1,-75.9,x", "007, 41.5 ,-76.2,y")))
  expect_identical(loc$collection, c("NA", "007"))
  expect_equal(loc$latitude, c(42.1, 41.5)); expect_equal(loc$longitude, c(-75.9, -76.2))
})

test_that("each rejection condition of a locations file gives its message", {
  bad <- list(
    missing = c("collection,latitude", "A,1"),
    repeated = c("collection,latitude, Latitude ,longitude", "A,1,2,3"),
    empty = "collection,latitude,longitude",
    noname = c("collection,latitude,longitude", " ,1,2"),
    duplicate = c("collection,latitude,longitude", "A,1,2", " A ,3,4"),
    blank = c("collection,latitude,longitude", "A,,2"),
    text = c("collection,latitude,longitude", "A,north,2"),
    infinite = c("collection,latitude,longitude", "A,Inf,2"),
    lat = c("collection,latitude,longitude", "A,91,2"),
    lon = c("collection,latitude,longitude", "A,1,-181"),
    many = c("collection,latitude,longitude", sprintf("S%d,1,2", 1:201)))
  msgs <- c(missing = "Missing: longitude", repeated = "more than one column named latitude", empty = "no rows",
            noname = "must have a collection name", duplicate = "more than once: A", blank = "latitude in the locations file must be a number",
            text = "must be a number", infinite = "must be a number", lat = "between -90 and 90", lon = "between -180 and 180",
            many = "more than 200 rows")
  for (k in names(bad)) expect_error(read_locations(write_lines_csv(bad[[k]])), msgs[[k]], fixed = TRUE, info = k)
})

test_that("locations are matched by exact name and the template lists every label", {
  loc <- data.frame(collection = c("A", "C", "Z"), latitude = c(1, 2, 3), longitude = c(4, 5, 6))
  m <- match_locations(loc, c("A", "B", "C"))
  expect_identical(m$located$collection, c("A", "C")); expect_identical(m$unmatched_rows, "Z"); expect_identical(m$not_located, "B")
  expect_identical(nrow(match_locations(NULL, "A")$located), 0L)
  tp <- locations_template(c("B", " A", "", "A", "Q"), loc)
  expect_identical(names(tp), c("collection", "latitude", "longitude"))
  expect_identical(tp$collection, c("A", "B", "Q")); expect_equal(tp$latitude, c(1, NA, NA))
})

# ---- map ---------------------------------------------------------------------

test_that("the map window, background, and scale bar follow their rules", {
  w <- map_window(c(42, 42), c(-76, -76))                       # coincident positions: minimum spans
  expect_true(w$ok)
  expect_equal(diff(w$ylim), 0.5, tolerance = 1e-9)
  expect_equal(diff(w$xlim), 0.5 / cos(42 * pi / 180), tolerance = 1e-9)
  w2 <- map_window(c(40, 44), c(-80, -70))
  expect_equal(w2$ylim, c(42 - 2.8, 42 + 2.8), tolerance = 1e-9); expect_equal(w2$xlim, c(-75 - 7, -75 + 7), tolerance = 1e-9)
  expect_identical(map_layers(w2), c("world", "state"))
  expect_identical(map_layers(w), c("world", "state", "county"))
  expect_identical(map_layers(map_window(c(50, 52), c(5, 8))), "world")
  expect_identical(map_layers(map_window(c(48, 52), c(-130, -120))), c("world", "state"))   # straddles the box edge
  expect_true(map_window(c(79, 79.5), c(10, 20))$ok)
  expect_false(map_window(c(81, 70), c(10, 20))$ok)
  expect_false(map_window(c(10, 11), c(-170, 175))$ok)
  expect_false(map_window(10, 10)$ok)
  for (win in list(w, w2, map_window(c(79, 79.5), c(10, 20)))) {
    b <- scale_bar(win)
    width <- haversine_km(win$lat_mid, win$xlim[1], win$lat_mid, win$xlim[2])
    expect_identical(b$km, nice_length(width / 4))
    expect_true((b$km / 10^floor(log10(b$km))) %in% c(1, 2, 5))
    expect_equal(KM_PER_DEGREE * cos(b$y * pi / 180) * (b$x1 - b$x0), b$km, tolerance = 1e-9)   # length along its parallel
    expect_equal(b$y, win$ylim[1] + 0.08 * diff(win$ylim)); expect_equal(b$x0, win$xlim[1] + 0.05 * diff(win$xlim))
  }
  expect_identical(nice_length(30), 20); expect_identical(nice_length(40), 50); expect_identical(nice_length(80), 100)
})

test_that("the map builds with the window as its limits and the background layers it should have", {
  skip_if_not_installed("maps")
  loc <- data.frame(collection = c("A", "B", "C"), latitude = c(42, 42.3, 42.6), longitude = c(-76, -75.6, -75.9))
  w <- map_window(loc$latitude, loc$longitude)
  p <- plot_collection_map(loc, c(1, 2, 3), "Value", w)
  b <- expect_silent(ggplot2::ggplot_build(p))
  expect_equal(b$layout$panel_params[[1]]$x.range, w$xlim, tolerance = 1e-9)
  expect_equal(b$layout$panel_params[[1]]$y.range, w$ylim, tolerance = 1e-9)
  expect_length(p$layers, 3 + 4)                                  # world, state, county; points, labels, bar, bar label
  p2 <- plot_collection_map(transform(loc, longitude = longitude + 90), 1:3, "Value", map_window(loc$latitude, loc$longitude + 90))
  expect_length(p2$layers, 1 + 4)
})

# ---- figures ------------------------------------------------------------------

test_that("collection figures build for 2, 3, 6, and 9 collections with distinct colours", {
  for (G in c(2, 3, 6, 9)) {
    set.seed(SEED)
    d <- sim_groups(rep(3, G), p = 8, shift = 0.5)
    names_g <- levels(d$g)
    co <- d$X / 20; co[, 1] <- co[, 1] + 1; co[, 7] <- co[, 7] + 1
    shapes <- lapply(stats::setNames(names_g, names_g), function(g) efa_inverse(colMeans(co[d$g == g, , drop = FALSE]), 200L))
    pr <- pair_table(d$X, d$g, n_perm = 19); wt <- within_table(d$X, d$g)
    sh <- share_reference(d$X, d$g, n_perm = 19)
    plots <- list(plot_collection_means(shapes), plot_collection_share(sh), plot_pair_distances(pr),
                  plot_distance_matrix(pr, names_g), plot_within_variance(wt),
                  plot_shape_geo(transform(pr, distance_km = seq_len(nrow(pr)) * 10)))
    for (p in plots) expect_silent(ggplot2::ggplot_build(p))
    bm <- ggplot2::ggplot_build(plots[[1]])$data[[1]]
    expect_identical(length(unique(bm$group)), as.integer(G))
    expect_identical(length(unique(bm$colour)), as.integer(G))
    expect_identical(nrow(ggplot2::ggplot_build(plots[[3]])$data[[3]]), nrow(pr))
  }
  grid9 <- plot_outline_grid(stats::setNames(lapply(1:9, function(i) KITE * c(1 + i / 20, 1)), letters[1:9]), groups = letters[1:9])
  expect_identical(length(unique(ggplot2::ggplot_build(grid9)$data[[1]]$fill)), 9L)
  expect_length(unique(category_colours(letters[1:8])), 8L)
  expect_identical(unname(category_colours(letters[1:3])), OKABE_ITO[1:3])
})

test_that("collection membership excludes unlabelled specimens and small collections", {
  co <- matrix(seq_len(40), 10)
  cd <- collection_data(co, c("A", "A", "A", " B", "B ", "B", "", "C", "C", NA))
  expect_identical(cd$eligible, c("A", "B")); expect_identical(cd$rows, 1:6)
  expect_identical(cd$unlabelled, 2L); expect_identical(cd$small, c(C = 2L))
  expect_match(paste(text_collection_membership(cd), collapse = " "), "2 specimens without a group label")
})
