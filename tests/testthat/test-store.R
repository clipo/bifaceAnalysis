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
