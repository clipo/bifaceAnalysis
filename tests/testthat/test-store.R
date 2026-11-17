demo_store <- function() {
  out <- lapply(names(DEMO_FILES), function(n) {
    p <- process_image(load_working_image(DEMO_FILES[[n]]), direction = "left")
    new_specimen(paste0(n, ".jpg"), n, if (n %in% c("artifact1", "artifact2")) "A" else "", p$std$outline, p$std$base_row,
                 if (n == "artifact1") 8.5 else NA_real_, p$settings)
  })
  stats::setNames(out, vapply(out, function(s) s$id, character(1)))
}

test_that("export then load reproduces the store", {
  st <- demo_store()
  f <- tempfile(fileext = ".csv")
  write.csv(specimens_to_csv(st), f, row.names = FALSE)
  expect_identical(names(read.csv(f)), CSV_COLUMNS)
  back <- merge_specimens(list(), specimens_from_csv(f))
  expect_identical(names(back), names(st))
  for (id in names(st)) {
    expect_identical(back[[id]]$name, st[[id]]$name); expect_identical(back[[id]]$group, st[[id]]$group)
    expect_equal(back[[id]]$px_per_mm, st[[id]]$px_per_mm)
    expect_equal(back[[id]]$outline, st[[id]]$outline, tolerance = 1e-9)
    expect_identical(back[[id]]$base_row, st[[id]]$base_row)
    expect_equal(back[[id]]$settings, st[[id]]$settings, tolerance = 1e-9)
  }
  expect_equal(specimen_measurements(back), specimen_measurements(st), tolerance = 1e-9)
  expect_equal(specimen_measurements(st)$length_mm[1], specimen_measurements(st)$length_px[1] / 8.5)
})

test_that("files that are not in the schema are rejected with a message", {
  df <- specimens_to_csv(demo_store()[1])
  wr <- function(d) { f <- tempfile(fileext = ".csv"); write.csv(d, f, row.names = FALSE); f }
  expect_error(specimens_from_csv(wr(df[, setdiff(names(df), "base_row")])), "not an outline file")
  v2 <- df; v2$schema_version <- 2L
  expect_error(specimens_from_csv(wr(v2)), "schema version")
  bad <- df; bad$x <- as.character(bad$x); bad$x[3] <- "abc"
  expect_error(specimens_from_csv(wr(bad)), "only numbers")
  short <- df[1:5, ]
  expect_error(specimens_from_csv(wr(short)), "incomplete outline")
  expect_error(specimens_from_csv(wr(data.frame(a = 1))), "not an outline file")
})

test_that("merging keeps unrelated specimens and replaces identical ones", {
  st <- demo_store()
  again <- merge_specimens(st, st)                       # identical contours: replaced in place
  expect_identical(names(again), names(st))
  other <- st[2]; other[[1]]$id <- names(st)[1]; other[[1]]$name <- st[[1]]$name   # same id and name, different contour
  merged <- merge_specimens(st, stats::setNames(other, names(st)[1]))
  expect_length(merged, 5L)
  nm <- vapply(merged, function(s) s$name, character(1))
  expect_identical(anyDuplicated(nm), 0L)
  expect_equal(merged[[1]]$outline, st[[1]]$outline)
  expect_identical(unique_specimen_name(" a ", list(list(name = "a"), list(name = "a (2)"))), "a (3)")
  expect_identical(unique_specimen_name("", list()), "specimen")
})

test_that("a loaded record with the same contour replaces name, group, scale, and settings", {
  st <- demo_store()
  f <- tempfile(fileext = ".csv"); write.csv(specimens_to_csv(st), f, row.names = FALSE)
  loaded <- specimens_from_csv(f)                        # contours differ from the originals only by CSV rounding
  loaded[[1]]$name <- "renamed"; loaded[[1]]$group <- "Z"; loaded[[1]]$px_per_mm <- 3; loaded[[1]]$settings$direction <- "up"
  m <- merge_specimens(st, loaded)
  expect_identical(names(m), names(st))                  # replaced in place, nothing duplicated
  expect_identical(m[[1]]$name, "renamed"); expect_identical(m[[1]]$group, "Z")
  expect_equal(m[[1]]$px_per_mm, 3); expect_identical(m[[1]]$settings$direction, "up")
  moved <- loaded[1]; moved[[1]]$outline[5, 1] <- moved[[1]]$outline[5, 1] + 1e-6      # a different contour
  expect_length(merge_specimens(st, moved), 5L)
})

test_that("two different specimens with the same name both survive a load", {
  st <- demo_store()
  df <- specimens_to_csv(st[1:2]); df$specimen <- "dup"
  f <- tempfile(fileext = ".csv"); write.csv(df, f, row.names = FALSE)
  m <- merge_specimens(list(), specimens_from_csv(f))
  expect_length(m, 2L)
  expect_identical(unname(vapply(m, function(s) s$name, character(1))), c("dup", "dup (2)"))
})

test_that("outlines that could not be analysed are rejected before they reach the store", {
  good <- specimens_to_csv(demo_store()[1])
  wr <- function(d) { f <- tempfile(fileext = ".csv"); write.csv(d, f, row.names = FALSE); f }
  tiny <- function(x, y, base_row) { d <- good[seq_along(x), ]; d$row <- seq_along(x); d$x <- x; d$y <- y; d$base_row <- base_row; d }
  bad <- list(
    infinite = { d <- good; d$x[5] <- Inf; d },
    huge = { d <- good; d$y[7] <- 1e200; d },
    repeated_vertex = { d <- good; d$x[9] <- d$x[8]; d$y[9] <- d$y[8]; d },
    row_gap = { d <- good; d$row[d$row > 20] <- d$row[d$row > 20] + 1; d },
    clockwise = { d <- good; n <- nrow(d); d$x <- d$x[c(1, n:2)]; d$y <- d$y[c(1, n:2)]; d$base_row <- n + 2L - d$base_row; d },
    collinear = tiny(rep(0, 10), c(9:5, 0:4) + 0.5, 6),
    tip_at_base = tiny(c(0, 1, 2, 3, 3, 0, 2, 1, 0, -1), c(5, 4, 3, 2, 1, 5, 0, 1, 2, 3), 6),
    tip_below_base = { d <- good; d$y <- -d$y; d$x <- -d$x; d },
    no_id = { d <- good; d$id <- ""; d },
    bad_scale = { d <- good; d$px_per_mm <- -2; d },
    mixed_metadata = { d <- good; d$base_row[3] <- d$base_row[3] + 1L; d })
  for (n in names(bad)) expect_error(specimens_from_csv(wr(bad[[n]])), "outline file", info = n)
  expect_null(outline_problem(demo_store()[[1]]$outline, demo_store()[[1]]$base_row))
  # Every accepted record supports the measurements and the analysis.
  ok <- specimens_from_csv(wr(specimens_to_csv(demo_store())))
  expect_true(all(is.finite(unlist(specimen_measurements(ok)[c("length_px", "width_px", "asymmetry", "widest_position")]))))
  expect_identical(shape_pca(efa_matrix(gpa(lapply(ok, function(s) s$points))$aligned, 12))$rank, 3L)
})
