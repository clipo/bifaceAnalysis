test_that("Otsu separates two grey levels", {
  t <- otsu_threshold(c(rep(0.2, 500), rep(0.8, 500)))
  expect_gt(t, 0.2); expect_lt(t, 0.8)
  expect_equal(otsu_threshold(rep(0.4, 10)), 0.4)
})

test_that("opening removes specks narrower than the element and keeps larger regions", {
  m <- matrix(FALSE, 30, 30); m[5:20, 5:20] <- TRUE; m[25, 25] <- TRUE; m[27:28, 2:3] <- TRUE
  o <- open_square(m, 3L)
  expect_true(all(o[5:20, 5:20])); expect_false(o[25, 25]); expect_false(any(o[27:28, 2:3]))
  expect_identical(open_square(m, 0L), m)
  expect_true(all(!open_square(m, 3L) | m))            # opening never adds cells
})

test_that("the largest region not touching the border is selected and its holes are filled", {
  g <- matrix(0.1, 120, 100)
  disc <- function(cx, cy, r) outer(1:120, 1:100, function(x, y) (x - cx)^2 + (y - cy)^2 <= r^2)
  big <- disc(60, 55, 25)
  g[big] <- 0.9
  g[disc(60, 55, 6)] <- 0.1          # dark hole (ink)
  g[disc(20, 85, 7)] <- 0.9          # smaller light object
  g[, 1:8] <- 0.9                    # light band on one border (scale bar)
  s <- segment_image(g)
  expect_identical(s$polarity, "light")
  # Speck removal also removes the single cell at each of the disc's four
  # extremes, so the mask is the disc less at most those four cells.
  expect_false(any(s$mask & !big))
  expect_lte(sum(big & !s$mask), 4L)
  expect_true(all(s$mask[disc(60, 55, 6)]))             # hole filled
  expect_length(s$warnings, 0)
})

test_that("a region touching every side falls back to the largest region with a warning", {
  g <- matrix(0.1, 60, 60); g[, 20:40] <- 0.9
  s <- segment_image(g, polarity = "light")
  expect_false(is.null(s$mask))
  expect_true(any(grepl("touches the edge", s$warnings)))
  expect_null(segment_image(matrix(0.1, 40, 40), threshold = 0.5, polarity = "light")$mask)
})

test_that("the four demo images segment with automatic settings", {
  expect_length(DEMO_FILES, 4)
  for (n in names(DEMO_FILES)) {
    img <- load_working_image(DEMO_FILES[[n]])
    s <- segment_image(img$gray)
    expect_length(s$warnings, 0)
    expect_identical(s$polarity, "light", info = n)
    expect_false(any(border_values(s$mask)), info = n)
    expect_equal(nrow(component_table(!s$mask)$table), 1L, info = n)       # no holes
    # gray is indexed [x, y]; the scale bar is in the top rows (small y).
    expect_false(any(s$mask[, seq_len(floor(0.15 * ncol(s$mask)))]), info = n)
    ref <- REF$masks[[n]]
    expect_gte(sum(s$mask & ref) / sum(s$mask | ref), 0.99)
    o <- mask_outline(s$mask)
    expect_gte(nrow(o), 100)
    expect_lt(sqrt(sum((o[1, ] - o[nrow(o), ])^2)), 1.5)                   # closed: ends adjacent
  }
})

test_that("process_image reports unusable input with a message and no error", {
  img <- load_working_image(DEMO_FILES[[4]])
  expect_false(process_image(img, crop = c(0, 10, 0, 10))$ok)
  p <- process_image(img, crop = c(5, 130, 85, 165), direction = "left")
  expect_type(p$warnings, "character")
  ok <- process_image(img, direction = "left")
  expect_true(ok$ok); expect_false(ok$settings$cropped); expect_false(ok$settings$manual_landmarks)
})

test_that("a competing object of the other brightness class triggers the warning", {
  g <- matrix(0.1, 300, 200); g[151:300, ] <- 0.9      # dark half and light half, both on the border
  g[190:266, 60:136] <- 0.1                             # dark square on the light half: 5929 cells
  g[40:110, 60:130] <- 0.9                              # light square on the dark half: 5041 cells
  s <- segment_image(g, threshold = 0.5)
  expect_identical(s$polarity, "dark"); expect_identical(sum(s$mask), 5929L)
  expect_true(any(grepl("^1 other object", s$warnings)))
  g[40:110, 60:130] <- 0.1                              # remove the competitor
  expect_false(any(grepl("other object", segment_image(g, threshold = 0.5)$warnings)))
  # A hole inside the selected region is part of the artifact, not a competitor.
  h <- matrix(0.1, 200, 200); h[40:160, 40:160] <- 0.9; h[60:140, 60:140] <- 0.1
  sh <- segment_image(h, threshold = 0.5)
  expect_identical(sum(sh$mask), 121L * 121L); expect_length(sh$warnings, 0)
})

test_that("a crop outside the image gives a message", {
  img <- load_working_image(DEMO_FILES[["artifact4"]])
  expect_lt(img$width, 160)
  for (crop in list(c(160, 240, 70, 170), c(-50, -10, 0, 50), c(10, 100, 5000, 6000), c(NA, 50, 0, 50))) {
    p <- process_image(img, crop = crop, direction = "left")
    expect_false(p$ok); expect_match(p$warnings, "too small or lies outside")
  }
  p <- process_image(img, crop = c(-30, 1000, -30, 1000), direction = "left")     # larger than the image: clipped
  expect_true(p$ok); expect_identical(dim(p$mask), c(img$width, img$height)); expect_equal(p$mask_origin, c(0, 0))
})

test_that("an artifact cut by the crop edge is reported when a small interior object is selected instead", {
  img <- load_working_image(DEMO_FILES[["artifact1"]])
  whole <- process_image(img, direction = "left")
  bb <- apply(whole$outline_plot, 2, range)
  cut <- process_image(img, crop = c(mean(bb[, 1]), bb[2, 1] + 12, bb[1, 2] - 12, bb[2, 2] + 12), direction = "left")
  expect_true(any(grepl("may be cut off", cut$warnings)))
  if (!is.null(cut$mask)) expect_lt(sum(cut$mask), 0.25 * sum(whole$mask))    # what was outlined is not the artifact
  for (f in DEMO_FILES) expect_length(segment_image(load_working_image(f)$gray)$warnings, 0)
  g <- matrix(0.1, 120, 100); g[40:80, 30:70] <- 0.9; g[, 1:8] <- 0.9          # artifact larger than nothing on the border but the band
  expect_false(any(grepl("may be cut off", segment_image(g)$warnings)))
})
