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
