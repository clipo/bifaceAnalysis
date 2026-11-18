# Server logic of the app, run without a browser.
defaults <- list(image_source = "demo", auto_threshold = TRUE, threshold = 0.5, polarity = "auto", speck = 3,
                 direction = "left", mode = "crop", scale_mm = 10, scale_all = FALSE, colour_by = "group",
                 symmetrize = FALSE, auto_harmonics = TRUE, n_harmonics = 12, k = 2,
                 b_tip = 15, b_blade = 60, b_shoulder = 80, spec_name = "", spec_group = "")
run_app_test <- function(expr) {
  old <- setwd(APP_DIR); on.exit(setwd(old))
  shiny::testServer(shiny::shinyAppDir(APP_DIR), {{ expr }})
}
pdf_pages <- function(f) {
  length(grepRaw("/Type /Page /Parent", readBin(f, "raw", file.size(f)), fixed = TRUE, all = TRUE))
}
try_figure <- function(r) tryCatch({ p <- r(); ggplot2::ggplot_build(p); "plot" }, error = function(e) if (is_validation(e)) "message" else paste("ERROR:", conditionMessage(e)))
try_text <- function(r) tryCatch({ t <- r(); if (is.character(t) && length(t) > 0 && all(nzchar(t))) "text" else "EMPTY" }, error = function(e) if (is_validation(e)) "message" else paste("ERROR:", conditionMessage(e)))

test_that("the demonstration produces every figure, text, and download", {
  run_app_test({
    do.call(session$setInputs, defaults)
    session$setInputs(run_demo = 1)
    expect_identical(n_spec(), 4L)
    # The demonstration has no group labels, so the collection results give a message.
    for (f in names(figures)) expect_identical(try_figure(figures[[f]]), if (startsWith(f, "coll")) "message" else "plot", info = f)
    for (t in names(texts)) expect_identical(try_text(texts[[t]]), if (startsWith(t, "coll") && t != "coll_cautions") "message" else "text", info = t)

    pdf <- output$dl_report
    expect_identical(readBin(pdf, "raw", 4), charToRaw("%PDF"))
    expect_gte(pdf_pages(pdf), 10)

    meas <- read.csv(output$dl_meas)
    expect_identical(nrow(meas), 4L)
    expect_true(all(c("specimen", "group", "length_px", "width_px", "area_px2", "length_mm", "width_mm", "area_mm2",
                      "length_width_ratio", "widest_position", "px_per_mm", "asymmetry") %in% names(meas)))
    expect_true(all(is.na(meas$length_mm)))                 # demo has no scale: no millimetre values
    scores <- read.csv(output$dl_scores)
    expect_identical(nrow(scores), 4L); expect_identical(names(scores), c("specimen", "group", "cluster", "PC1", "PC2", "PC3"))
    coef <- read.csv(output$dl_coef)
    expect_identical(dim(coef), c(4L, 1L + 4L * n_harm()))
    outl <- read.csv(output$dl_outlines)
    expect_identical(names(outl), CSV_COLUMNS); expect_identical(length(unique(outl$id)), 4L)
    expect_identical(scores$specimen, meas$specimen); expect_identical(coef$specimen, meas$specimen)
  })
})

test_that("analysis settings do not alter saved specimens, and the symmetric option is reversible", {
  run_app_test({
    do.call(session$setInputs, defaults)
    session$setInputs(run_demo = 1)
    saved <- rv$specimens; meas <- measurements(); sc <- pca()$scores
    session$setInputs(auto_harmonics = FALSE, n_harmonics = 8, b_tip = 25, b_blade = 50, b_shoulder = 70, k = 3)
    expect_identical(n_harm(), 8L)
    session$setInputs(symmetrize = TRUE)
    expect_lt(asymmetry_score(aligned()$aligned[[1]], mirror_index()), 1e-8)
    expect_identical(rv$specimens, saved); expect_identical(measurements(), meas)
    session$setInputs(symmetrize = FALSE, auto_harmonics = TRUE, b_tip = 15, b_blade = 60, b_shoulder = 80, k = 2)
    expect_identical(pca()$scores, sc)
    session$setInputs(k = 99); expect_identical(k_used(), 3L)
    session$setInputs(k = NA); expect_identical(k_used(), 2L)
    session$setInputs(b_tip = 80, b_blade = 15, b_shoulder = 60)          # out of order
    expect_identical(boundaries(), c(0.15, 0.60, 0.80))
  })
})

test_that("saving replaces a specimen for the same image and advances", {
  run_app_test({
    do.call(session$setInputs, defaults)
    session$elapse(300)
    expect_true(processed()$ok)
    expect_identical(current_id(), "artifact1.jpg")
    session$setInputs(spec_name = "first", spec_group = " g1 ", save_spec = 1)
    expect_identical(n_spec(), 1L); expect_identical(rv$idx, 2L)
    expect_identical(unname(groups()), "g1")                  # trimmed when saved
    session$setInputs(prev_img = 1)
    session$elapse(300)
    session$setInputs(spec_name = "renamed", save_spec = 2)
    expect_identical(n_spec(), 1L)
    expect_identical(rv$specimens[["artifact1.jpg"]]$name, "renamed")
    expect_match(rv$last_saved, "^Replaced")
    session$setInputs(edit_which = "artifact1.jpg", edit_name = "again", edit_group = "  g2 ", edit_apply = 1)
    expect_identical(unname(groups()), "g2")                  # trimmed when edited
    session$setInputs(edit_remove = 1)
    expect_identical(n_spec(), 0L)
  })
})

test_that("with 0 to 3 specimens every output gives a message or a result, never an error", {
  run_app_test({
    do.call(session$setInputs, defaults)
    check <- function(n) {
      expect_identical(n_spec(), n)
      for (f in names(figures)) expect_false(grepl("^ERROR", try_figure(figures[[f]])), info = paste(n, f))
      for (t in names(texts)) expect_false(grepl("^ERROR|EMPTY", try_text(texts[[t]])), info = paste(n, t))
      expect_identical(readBin(output$dl_report, "raw", 4), charToRaw("%PDF"))
      for (d in c("dl_meas", "dl_scores", "dl_coef", "dl_outlines")) expect_gt(nrow(read.csv(output[[d]])), 0)
    }
    check(0L)
    session$setInputs(run_demo = 1)
    for (id in c("artifact4.jpg", "artifact3.jpg", "artifact2.jpg")) {
      session$setInputs(edit_which = id); session$setInputs(edit_remove = which(id == c("artifact4.jpg", "artifact3.jpg", "artifact2.jpg")))
      n <- n_spec(); check(n)
      if (n == 3L) {
        expect_identical(try_figure(figures$stab), "message"); expect_identical(try_figure(figures$dend), "message")
        expect_identical(try_figure(figures$morpho), "plot"); expect_identical(pca()$rank, 2L)
        expect_identical(ncol(read.csv(output$dl_scores)), 4L)           # specimen, group, PC1, PC2: no PC3
      }
      if (n == 2L) { expect_identical(try_figure(figures$morpho), "message"); expect_identical(try_figure(figures$stack), "plot") }
      if (n == 1L) { expect_identical(try_figure(figures$stack), "message"); expect_identical(try_figure(figures$meas), "plot") }
    }
  })
})

test_that("identical shapes give the 'identical' message", {
  run_app_test({
    do.call(session$setInputs, defaults)
    p <- process_image(load_working_image(file.path("www", "artifact1.jpg")), direction = "left")
    for (i in 1:4) rv$specimens[[paste0("copy", i)]] <- new_specimen(paste0("copy", i), paste0("copy", i), "", p$std$outline, p$std$base_row)
    expect_identical(n_spec(), 4L)
    msg <- tryCatch(figures$morpho(), error = function(e) conditionMessage(e))
    expect_match(msg, "identical in shape")
    expect_identical(try_figure(figures$rv), "message")
    expect_identical(readBin(output$dl_report, "raw", 4), charToRaw("%PDF"))
  })
})

test_that("an outline file loads into the running app", {
  run_app_test({
    do.call(session$setInputs, defaults)
    session$setInputs(run_demo = 1)
    f <- tempfile(fileext = ".csv"); file.copy(output$dl_outlines, f)
    for (id in names(rv$specimens)) rv$specimens[[id]] <- NULL
    expect_identical(n_spec(), 0L)
    session$setInputs(load_outlines = data.frame(name = "x.csv", datapath = f, stringsAsFactors = FALSE))
    expect_identical(n_spec(), 4L)
    expect_identical(try_figure(figures$morpho), "plot")
  })
})

test_that("crop, scale, and landmark interactions update the preview only", {
  run_app_test({
    do.call(session$setInputs, defaults)
    session$elapse(300)
    img <- current_image()
    auto <- processed()
    expect_true(auto$ok)
    expect_true(is.list(output$image_plot)); expect_true(is.list(output$std_plot)); expect_true(is.list(output$mask_plot))
    expect_identical(dim(auto$mask), c(img$width, img$height))

    # Crop: a box around the artifact keeps the outline; the crop is recorded.
    bb <- apply(auto$outline_plot, 2, range)
    session$setInputs(mode = "crop", img_brush = list(xmin = bb[1, 1] - 15, xmax = bb[2, 1] + 15, ymin = bb[1, 2] - 15, ymax = bb[2, 2] + 15))
    session$elapse(300)
    expect_true(processed()$ok); expect_true(processed()$settings$cropped)
    expect_equal(measure_outline(processed()$std$outline)$length_px, measure_outline(auto$std$outline)$length_px, tolerance = 0.01)
    session$setInputs(img_brush = list(xmin = 1, xmax = 5, ymin = 1, ymax = 5))      # too small: ignored
    expect_true(processed()$settings$cropped)
    session$setInputs(clear_crop = 1); session$elapse(300)
    expect_false(processed()$settings$cropped)

    # Scale: two clicks d working pixels apart for 25 mm give d / 25 px per mm.
    d <- round(img$width / 2)
    session$setInputs(mode = "scale", scale_mm = 25, img_click = list(x = 10, y = 30))
    expect_match(output$scale_text, "One point set")
    session$setInputs(img_click = list(x = 10 + d, y = 30))
    expect_equal(own_ppm(), d / 25); expect_equal(current_ppm(), d / 25)
    expect_equal(rv$shared_ppm, d / 25 / img$factor)
    session$setInputs(scale_mm = 50); expect_equal(own_ppm(), d / 50)
    session$setInputs(img_click = list(x = -20, y = 30))                              # outside the image: ignored
    expect_equal(own_ppm(), d / 50)

    # The shared scale applies to another image only when the box is ticked.
    session$setInputs(next_img = 1); session$elapse(300)
    expect_true(is.na(current_ppm()))
    session$setInputs(scale_all = TRUE)
    expect_equal(current_ppm(), rv$shared_ppm * current_image()$factor)
    session$setInputs(spec_name = "scaled", save_spec = 1)
    m <- measurements()
    expect_equal(m$length_mm, m$length_px / m$px_per_mm); expect_true(is.finite(m$length_mm))
    session$setInputs(prev_img = 1, clear_scale = 1); session$elapse(300)

    # Landmarks: tip then base, swapped relative to the automatic ones.
    session$setInputs(prev_img = 2); session$elapse(300)
    expect_identical(current_id(), "artifact1.jpg")
    auto <- processed()
    session$setInputs(mode = "landmarks", img_click = list(x = auto$base_plot[1], y = auto$base_plot[2]))
    expect_false(processed()$settings$manual_landmarks)                               # one landmark: still automatic
    session$setInputs(img_click = list(x = auto$tip_plot[1], y = auto$tip_plot[2]))
    session$elapse(300)
    man <- processed()
    expect_true(man$ok); expect_true(man$settings$manual_landmarks)
    expect_lt(sqrt(sum((man$tip_plot - auto$base_plot)^2)), 3)
    expect_lt(sqrt(sum((man$base_plot - auto$tip_plot)^2)), 3)
    session$setInputs(clear_landmarks = 1); session$elapse(300)
    expect_equal(processed()$tip_plot, auto$tip_plot)

    # Detection settings change the preview and are recorded.
    session$setInputs(auto_threshold = FALSE, threshold = 0.6, speck = 0, direction = "right"); session$elapse(300)
    p <- processed()
    expect_equal(p$settings$threshold, 0.6); expect_identical(as.integer(p$settings$speck), 0L)
    expect_lt(sqrt(sum((p$tip_plot - auto$base_plot)^2)), 0.25 * measure_outline(auto$std$outline)$length_px)
    session$setInputs(threshold = 0.98); session$elapse(300)
    expect_false(is.null(processed()$warnings))                                       # no error at an extreme threshold
    expect_identical(n_spec(), 1L)                                                    # none of this altered the saved specimen
  })
})

upload_input <- function(files) data.frame(name = names(files), size = file.size(files), type = "image/jpeg", datapath = unname(files), stringsAsFactors = FALSE)
figure_message <- function(r) tryCatch({ r(); "" }, error = function(e) conditionMessage(e))

test_that("Clear scale removes the effective scale of the image and leaves the shared scale for others", {
  run_app_test({
    do.call(session$setInputs, defaults)
    session$elapse(300)
    session$setInputs(mode = "scale", scale_mm = 10, scale_all = TRUE, img_click = list(x = 10, y = 30))
    session$setInputs(img_click = list(x = 110, y = 30))
    expect_equal(current_ppm(), 10)
    session$setInputs(clear_scale = 1)
    expect_true(is.na(own_ppm())); expect_true(is.na(current_ppm()))
    expect_match(output$scale_text, "Scale cleared for this image")
    session$setInputs(spec_name = "unscaled", save_spec = 1)
    expect_true(is.na(rv$specimens[["artifact1.jpg"]]$px_per_mm))
    expect_true(is.na(measurements()$length_mm))
    session$elapse(300)
    expect_identical(current_id(), "artifact2.jpg")
    expect_equal(current_ppm(), 10 * current_image()$factor)             # the shared scale still serves other images
    session$setInputs(prev_img = 1); session$elapse(300)
    expect_true(is.na(current_ppm()))                                    # and stays off for the cleared image
    session$setInputs(img_click = list(x = 10, y = 40)); session$setInputs(img_click = list(x = 60, y = 40))
    expect_equal(current_ppm(), 5)                                       # until a new scale is clicked on it
  })
})

test_that("an uploaded file never inherits the state or the specimen of an image with the same name", {
  up <- file.path(tempfile("up"), "artifact1.jpg"); dir.create(dirname(up))
  file.copy(DEMO_FILES[["artifact4"]], up)                               # different contents under a repeated name
  run_app_test({
    do.call(session$setInputs, defaults)
    session$elapse(300)
    session$setInputs(mode = "crop", img_brush = list(xmin = 0, xmax = 240, ymin = 70, ymax = 170))
    session$setInputs(spec_name = "demo one", save_spec = 1)
    demo_outline <- rv$specimens[["artifact1.jpg"]]$outline
    session$setInputs(images = upload_input(c(artifact1.jpg = up)), image_source = "upload")
    session$elapse(300)
    expect_identical(current_id(), "upload1/artifact1.jpg"); expect_identical(current_name(), "artifact1.jpg")
    expect_length(state(), 0L)                                           # no crop, scale, or landmarks carried over
    expect_true(processed()$ok); expect_false(processed()$settings$cropped)
    session$setInputs(img_brush = list(xmin = 2, xmax = 132, ymin = 60, ymax = 180))
    expect_false(is.null(state()$crop))
    session$setInputs(spec_name = "demo one", save_spec = 1)
    expect_identical(n_spec(), 2L)
    expect_identical(rv$specimens[["artifact1.jpg"]]$outline, demo_outline)   # the earlier specimen is untouched
    expect_identical(unname(vapply(rv$specimens, function(s) s$name, character(1))), c("demo one", "demo one (2)"))
    session$setInputs(images = upload_input(c(artifact1.jpg = up)))      # the same file uploaded again
    session$elapse(300)
    expect_identical(current_id(), "upload2/artifact1.jpg"); expect_length(state(), 0L)
    session$setInputs(image_source = "demo"); session$elapse(300)
    expect_identical(current_id(), "artifact1.jpg")
    expect_equal(state()$crop, c(0, 240, 70, 170))                       # the sample image keeps its own state
  })
})

test_that("an invalid pair of landmark clicks keeps the previous landmarks", {
  run_app_test({
    do.call(session$setInputs, defaults)
    session$elapse(300)
    auto <- processed()
    session$setInputs(mode = "landmarks", img_click = list(x = auto$base_plot[1], y = auto$base_plot[2]))
    expect_false(is.null(state()$pending_tip)); expect_false(processed()$settings$manual_landmarks)
    session$setInputs(img_click = list(x = auto$tip_plot[1], y = auto$tip_plot[2])); session$elapse(300)
    valid <- processed()
    expect_true(valid$ok); expect_true(valid$settings$manual_landmarks); expect_null(state()$pending_tip)
    v <- valid$outline_plot[40, ]
    session$setInputs(img_click = list(x = v[1] + 0.2, y = v[2])); session$setInputs(img_click = list(x = v[1], y = v[2] + 0.2))
    session$elapse(300)
    after <- processed()
    expect_true(after$ok); expect_null(state()$pending_tip)
    expect_equal(after$tip_plot, valid$tip_plot); expect_equal(after$std$points, valid$std$points)
    session$setInputs(spec_name = "kept", save_spec = 1)
    expect_identical(n_spec(), 1L); expect_true(rv$specimens[[1]]$settings$manual_landmarks)
  })
})

test_that("identical shapes give no components at 2, 3, and 5 specimens", {
  run_app_test({
    do.call(session$setInputs, defaults)
    p <- process_image(load_working_image(file.path("www", "artifact1.jpg")), direction = "left")
    for (i in 1:5) {
      rv$specimens[[paste0("copy", i)]] <- new_specimen(paste0("copy", i), paste0("copy", i), "", p$std$outline, p$std$base_row)
      if (i < 2) next
      for (f in names(figures)) expect_false(grepl("^ERROR", try_figure(figures[[f]])), info = paste(i, f))
      if (i >= 3) {
        expect_identical(pca()$rank, 0L, info = i)
        for (f in c("morpho", "scree", "pcs")) expect_match(figure_message(figures[[f]]), "identical in shape", info = paste(i, f))
        expect_identical(names(read.csv(output$dl_scores)), "message")
      }
      if (i >= 4) expect_match(texts$clust()[2], "no dashed line")
      expect_identical(readBin(output$dl_report, "raw", 4), charToRaw("%PDF"))
    }
  })
})

test_that("equal region boundaries withhold the region results with a message", {
  run_app_test({
    do.call(session$setInputs, defaults)
    session$setInputs(run_demo = 1)
    session$setInputs(b_tip = 60, b_blade = 60, b_shoulder = 80)
    expect_false(boundaries_ok())
    for (f in c("regions", "rvar", "rv")) expect_match(figure_message(figures[[f]]), "Two region boundaries are equal", info = f)
    expect_identical(try_text(texts$regions), "message")
    expect_identical(try_figure(figures$pvar), "plot"); expect_identical(try_figure(figures$morpho), "plot")
    expect_identical(readBin(output$dl_report, "raw", 4), charToRaw("%PDF"))
    session$setInputs(b_blade = 65)
    expect_true(boundaries_ok()); expect_identical(try_figure(figures$rv), "plot")
  })
})

test_that("figures that need more specimens say how many, and the outline file has one row per vertex", {
  run_app_test({
    do.call(session$setInputs, defaults)
    for (f in names(figures)) {
      if (startsWith(f, "coll")) expect_match(figure_message(figures[[f]]), "Collections are defined by group labels", info = f)
      else expect_match(figure_message(figures[[f]]), "needs at least [0-9] specimens\\. You have 0\\.", info = f)
    }
    session$setInputs(run_demo = 1)
    outl <- read.csv(output$dl_outlines)
    expect_identical(nrow(outl), sum(vapply(rv$specimens, function(s) nrow(s$outline), integer(1))))
    expect_identical(as.integer(table(outl$id)[names(rv$specimens)]), unname(vapply(rv$specimens, function(s) nrow(s$outline), integer(1))))
    # A rejected file leaves the store unchanged.
    before <- rv$specimens
    bad <- outl; bad$x[3] <- Inf
    f <- tempfile(fileext = ".csv"); write.csv(bad, f, row.names = FALSE)
    session$setInputs(load_outlines = data.frame(name = "bad.csv", datapath = f, stringsAsFactors = FALSE))
    expect_identical(rv$specimens, before)
    # A loaded record replaces the fields shown for the selected specimen.
    session$setInputs(edit_which = "artifact2.jpg")
    good <- outl; good$specimen[good$id == "artifact2.jpg"] <- "from file"; good$group[good$id == "artifact2.jpg"] <- "G"
    write.csv(good, f, row.names = FALSE)
    session$setInputs(load_outlines = data.frame(name = "good.csv", datapath = f, stringsAsFactors = FALSE))
    expect_identical(n_spec(), 4L)
    expect_identical(edited()$name, "from file"); expect_identical(edited()$group, "G")
  })
})

# ---- collections (docs/specs/2026-11-17_collections_spec.md) -----------------

# Outline file with 3 collections of 3 (A, B, C), one of 2 (D), and 2
# unlabelled specimens, built from the demonstration outlines by stretching
# their width.
collection_csv <- function() {
  base <- lapply(DEMO_FILES, function(f) process_image(load_working_image(f), direction = "left")$std)
  spec <- list(); k <- 0
  add <- function(src, stretch, group) {
    k <<- k + 1
    b <- base[[src]]
    spec[[k]] <<- new_specimen(sprintf("s%02d", k), sprintf("%s_%02d", src, k), group,
                               cbind(b$outline[, 1] * stretch, b$outline[, 2]), b$base_row, px_per_mm = 10)
  }
  for (st in c(0.9, 1, 1.12)) add("artifact1", st, "A")
  for (st in c(0.95, 1, 1.08)) add("artifact2", st, " B ")
  for (st in c(0.9, 1.05, 1.15)) add("artifact3", st, "C")
  for (st in c(1, 1.1)) add("artifact4", st, "D")
  add("artifact4", 0.92, ""); add("artifact1", 1.2, "")
  f <- tempfile(fileext = ".csv")
  utils::write.csv(specimens_to_csv(spec), f, row.names = FALSE)
  f
}
loc_file <- function(lines) { f <- tempfile(fileext = ".csv"); writeLines(lines, f); data.frame(name = "loc.csv", datapath = f, stringsAsFactors = FALSE) }
upload <- function(f) data.frame(name = "x.csv", datapath = f, stringsAsFactors = FALSE)

test_that("labelled collections produce every collection result, and labels do not change other results", {
  run_app_test({
    do.call(session$setInputs, defaults)
    session$setInputs(load_outlines = upload(collection_csv()))
    expect_identical(n_spec(), 13L)
    expect_identical(sort(unique(groups())), c("", "A", "B", "C", "D"))   # " B " was trimmed on loading
    expect_identical(coll()$eligible, c("A", "B", "C"))
    expect_identical(coll()$unlabelled, 2L); expect_identical(coll()$small, c(D = 2L))
    for (f in c("coll_means", "coll_share", "coll_pairs", "coll_matrix", "coll_within"))
      expect_identical(try_figure(figures[[f]]), "plot", info = f)
    expect_identical(try_figure(figures$coll_map), "message")
    expect_identical(try_text(texts$coll), "text")
    expect_true(coll_available())
    tab <- output$coll_pair_table
    expect_match(tab, "Corrected square", fixed = TRUE); expect_match(tab, "Corrected distance", fixed = TRUE)
    expect_match(paste(texts$coll(), collapse = " "), "2 specimens without a group label")

    before <- list(pca = pca(), tree = tree(), rv = rv_table())
    for (id in names(rv$specimens)) rv$specimens[[id]]$group <- ""
    expect_identical(list(pca = pca(), tree = tree(), rv = rv_table()), before)
    expect_identical(try_figure(figures$coll_means), "message")
    expect_false(coll_available())
  })
})

test_that("locations load, replace, rematch, and clear, and downloads have their stated columns", {
  run_app_test({
    do.call(session$setInputs, defaults)
    sm0 <- read.csv(output$dl_coll_summary, check.names = FALSE)
    expect_identical(nrow(sm0), 0L)
    expect_identical(names(sm0), c("collection", "n", "pc1_mean", "pc2_mean", "pc3_mean", "length_width_mean",
                                   "max_width_position_mean", "within_variance", "length_mm_mean", "located"))
    session$setInputs(load_outlines = upload(collection_csv()))
    session$setInputs(locations_file = loc_file(c("collection,latitude,longitude", "A,42.1,-75.9", "B,42.4,-76.3", "Z,40,-70")))
    expect_identical(coll_loc()$located$collection, c("A", "B"))
    expect_identical(coll_loc()$unmatched_rows, "Z"); expect_identical(coll_loc()$not_located, "C")
    expect_identical(try_figure(figures$coll_map), "plot")
    expect_identical(try_figure(figures$coll_geo), "message")          # 2 located collections, 5 needed
    expect_match(paste(texts$coll_geo(), collapse = " "), "at least 5 located collections")

    session$setInputs(locations_file = loc_file(c("collection,latitude", "A,1")))   # rejected: earlier locations kept
    expect_identical(coll_loc()$located$collection, c("A", "B"))
    session$setInputs(locations_file = loc_file(c("collection,latitude,longitude", "C,41,-77")))  # replaces
    expect_identical(coll_loc()$located$collection, "C")
    expect_identical(try_figure(figures$coll_map), "message")          # one located collection
    session$setInputs(locations_file = loc_file(c("collection,latitude,longitude", "A,42.1,-75.9", "B,42.4,-76.3", "C,41,-77")))

    sm <- read.csv(output$dl_coll_summary, check.names = FALSE, na.strings = "")
    expect_identical(sm$collection, c("A", "B", "C")); expect_identical(sm$n, c(3L, 3L, 3L))
    expect_true(all(sm$located)); expect_true(all(is.finite(sm$length_mm_mean)))
    expect_false(any(c("latitude", "longitude") %in% names(sm)))
    di <- read.csv(output$dl_coll_dist, check.names = FALSE, na.strings = "")
    expect_identical(names(di), c("collection_a", "collection_b", "n_a", "n_b", "raw_distance", "reference_mean",
                                  "reference_low", "reference_high", "corrected_square", "corrected_distance", "distance_km"))
    expect_identical(nrow(di), 3L); expect_true(all(is.finite(di$distance_km)))
    expect_equal(di$distance_km[1], haversine_km(42.1, -75.9, 42.4, -76.3))

    id_a <- names(rv$specimens)[vapply(rv$specimens, function(s) s$group == "A", logical(1))]
    for (id in id_a) rv$specimens[[id]]$group <- "A2"                   # renamed: loses its location
    expect_identical(coll_loc()$not_located, "A2"); expect_identical(coll_loc()$unmatched_rows, "A")
    tp <- read.csv(output$dl_loc_template, na.strings = "")
    expect_identical(tp$collection, c("A2", "B", "C", "D")); expect_true(is.na(tp$latitude[1]))
    di2 <- read.csv(output$dl_coll_dist, na.strings = "")
    expect_true(is.na(di2$distance_km[di2$collection_a == "A2" | di2$collection_b == "A2"][1]))

    session$setInputs(map_in_report = TRUE)             # the box is ticked by default in the page
    pages_with_map <- pdf_pages(output$dl_report)
    session$setInputs(map_in_report = FALSE)
    expect_lt(pdf_pages(output$dl_report), pages_with_map)
    session$setInputs(clear_locations = 1)
    expect_null(rv$locations)
    expect_identical(try_figure(figures$coll_map), "message")
  })
})

test_that("five located collections give the shape-geography comparison, and the map selector follows the PCA rank", {
  run_app_test({
    do.call(session$setInputs, defaults)
    session$setInputs(load_outlines = upload(collection_csv()))
    ids <- names(rv$specimens)
    labs <- c("A", "A", "A", "B", "B", "B", "C", "C", "C", "D", "D", "D", "")
    for (i in seq_along(ids)) rv$specimens[[ids[i]]]$group <- labs[i]
    session$setInputs(locations_file = loc_file(c("collection,latitude,longitude", "A,42.1,-75.9", "B,42.4,-76.3", "C,41,-77", "D,43,-74.5")))
    expect_identical(try_figure(figures$coll_geo), "message")          # 4 located
    expect_true(all(c("pc1_mean", "pc2_mean", "pc3_mean") %in% map_choices()))
    expect_identical(map_value(), "length_width_mean")
    session$setInputs(map_value = "pc2_mean")
    expect_identical(map_value(), "pc2_mean")
    expect_identical(try_figure(figures$coll_map), "plot")

    # A fifth collection of 3 specimens, located: the comparison is drawn.
    b <- rv$specimens[[ids[1]]]
    for (j in 1:3) rv$specimens[[paste0("e", j)]] <- new_specimen(paste0("e", j), paste0("extra_", j), "E",
                                                                  cbind(b$outline[, 1] * (0.8 + j / 20), b$outline[, 2]), b$base_row, 10)
    session$setInputs(locations_file = loc_file(c("collection,latitude,longitude", "A,42.1,-75.9", "B,42.4,-76.3",
                                                  "C,41,-77", "D,43,-74.5", "E,41.6,-75.1")))
    expect_identical(try_figure(figures$coll_geo), "plot")
    g <- coll_geo()
    expect_true(g$available); expect_identical(g$n_assignments, 120L); expect_identical(g$n_pairs, 10)
    expect_match(paste(texts$coll_geo(), collapse = " "), "5 located collections give 10 pairs")
    expect_true(all(is.finite(read.csv(output$dl_coll_dist, na.strings = "")$distance_km)))
  })
})

# Specimens made from one outline, stretched by the given factors.
stretched <- function(src, factors, group, prefix) {
  b <- process_image(load_working_image(DEMO_FILES[[src]]), direction = "left")$std
  out <- lapply(seq_along(factors), function(i)
    new_specimen(sprintf("%s%d", prefix, i), sprintf("%s_%d", prefix, i), group,
                 cbind(b$outline[, 1] * factors[i], b$outline[, 2]), b$base_row, px_per_mm = 10))
  stats::setNames(out, vapply(out, function(s) s$id, ""))
}

test_that("collections of identical specimens have zero variance, and the map selector follows the PCA rank", {
  run_app_test({
    do.call(session$setInputs, defaults)
    rv$specimens <- c(stretched("artifact1", c(0.9, 1, 1.1), "A", "a"), stretched("artifact2", c(0.9, 1, 1.1), "B", "b"),
                      stretched("artifact3", c(0.95, 1.1, 1.2), "C", "c"))
    expect_true(pca()$rank >= 3L)
    session$setInputs(map_value = "pc3_mean")
    expect_identical(map_value(), "pc3_mean")

    # Three copies of one outline and three of another: PCA rank 1.
    rv$specimens <- c(stretched("artifact1", c(1, 1, 1), "A", "a"), stretched("artifact2", c(1, 1, 1), "B", "b"))
    expect_identical(pca()$rank, 1L)
    expect_identical(coll_stats()$within$variance, c(0, 0))
    expect_identical(coll_stats()$within$loo_low, c(0, 0)); expect_identical(coll_stats()$within$loo_high, c(0, 0))
    txt <- paste(texts$coll(), collapse = " ")
    expect_match(txt, "No collection shows variation"); expect_false(grepl("a ratio of", txt))
    expect_false(any(c("pc2_mean", "pc3_mean") %in% map_choices()))
    expect_identical(map_value(), "length_width_mean")      # unavailable choice falls back to the default
    pr <- coll_stats()$pairs
    expect_equal(pr$corrected_square, pr$raw_distance^2)    # no sampling error to subtract

    # One constant and one varying collection.
    rv$specimens <- c(stretched("artifact1", c(1, 1, 1), "A", "a"), stretched("artifact2", c(0.9, 1, 1.1), "B", "b"))
    v <- coll_stats()$within$variance
    expect_identical(v[1], 0); expect_gt(v[2], 0)
    expect_match(paste(texts$coll(), collapse = " "), "The specimens of A do not differ")
  })
})

test_that("the report holds every collection section and caution, with and without the map", {
  run_app_test({
    do.call(session$setInputs, defaults)
    session$setInputs(load_outlines = upload(collection_csv()))
    headings <- function() vapply(report_content()$sections, function(x) x$heading, "")
    texts_of <- function() unlist(lapply(report_content()$sections, function(x) x$text))
    coll_heads <- c("Comparison of collections: mean shapes", "Comparison of collections: before drawing conclusions",
                    "Share of shape variation between collections", "Distances between collection means",
                    "Corrected distances between collection means", "Shape variance within collections")
    expect_true(all(coll_heads %in% headings()))
    expect_false(any(c("Map of collections", "Shape distance and geographic distance") %in% headings()))
    expect_true(all(COLLECTION_CAUTIONS %in% texts_of()))     # including the geographic caution, with no map
    # Checked by wording, so that removing a caution from the list is detected.
    for (w in c("Rounded or approximate coordinates change the geographic distances", "not confidence or credible intervals",
                "tick \"Analyse symmetric shapes\"", "not the directions that best separate the collections"))
      expect_true(any(grepl(w, texts_of(), fixed = TRUE)), info = w)

    session$setInputs(locations_file = loc_file(c("collection,latitude,longitude", "A,42.1,-75.9", "B,42.4,-76.3", "C,41,-77")))
    session$setInputs(map_in_report = TRUE)
    expect_true("Map of collections" %in% headings())
    expect_false("Shape distance and geographic distance" %in% headings())   # 3 located, 5 needed
    session$setInputs(map_in_report = FALSE)
    expect_false("Map of collections" %in% headings())
    expect_true(all(COLLECTION_CAUTIONS %in% texts_of()))

    lines <- readLines(output$dl_coll_dist)
    expect_identical(length(lines), 4L)                           # header and 3 pairs
    for (id in names(rv$specimens)) if (rv$specimens[[id]]$group == "C") rv$specimens[[id]]$group <- "C2"
    session$setInputs(locations_file = loc_file(c("collection,latitude,longitude", "A,42.1,-75.9", "B,42.4,-76.3")))
    dl <- readLines(output$dl_coll_dist)
    expect_true(all(grepl(",$", dl[grepl("C2", dl)])))             # distance_km written as an empty cell
    expect_false(any(grepl("NA", dl[-1])))
    for (id in names(rv$specimens)) rv$specimens[[id]]$group <- ""
    expect_identical(readLines(output$dl_coll_dist), paste(sprintf('"%s"', c("collection_a", "collection_b", "n_a", "n_b", "raw_distance",
                     "reference_mean", "reference_low", "reference_high", "corrected_square", "corrected_distance", "distance_km")), collapse = ","))
  })
})


test_that("whether the comparison is shown is decided from labels without the outline analysis", {
  run_app_test({
    do.call(session$setInputs, defaults)
    session$setInputs(load_outlines = upload(collection_csv()))
    expect_true(coll_available())
    # Unusable outlines make the analysis fail; availability must not depend on it.
    for (id in names(rv$specimens)) { rv$specimens[[id]]$outline <- rv$specimens[[id]]$outline[1:2, ]; rv$specimens[[id]]$points <- NULL }
    expect_error(coefs())
    expect_true(coll_available())
    for (id in names(rv$specimens)) rv$specimens[[id]]$group <- "one"
    expect_false(coll_available())
  })
})

test_that("only rendered outputs compute the collection statistics, never an observer", {
  run_app_test({
    do.call(session$setInputs, defaults)
    # Observers run even when their tab is hidden, so an observer that reads the
    # statistics recomputes the permutation references on every edit in step 3.
    # testServer renders every output regardless of visibility, so the check is
    # that every call to pair_table comes from inside an output's render
    # function. The app's R/ files are sourced above the reactives' environment.
    fenv <- environment(attr(coll_stats, "observable")$.origFunc)
    while (!exists("pair_table", envir = fenv, inherits = FALSE)) fenv <- parent.env(fenv)
    original <- get("pair_table", envir = fenv); calls <- 0L; outside <- 0L
    assign("pair_table", function(...) {
      calls <<- calls + 1L
      if (!any(vapply(sys.calls(), function(x) identical(x[[1]], as.name("renderFunc")), logical(1)))) outside <<- outside + 1L
      original(...)
    }, envir = fenv)
    on.exit(assign("pair_table", original, envir = fenv), add = TRUE)
    session$setInputs(load_outlines = upload(collection_csv()))
    session$setInputs(locations_file = loc_file(c("collection,latitude,longitude", "A,42.1,-75.9", "B,42.4,-76.3", "C,41,-77")))
    id <- names(rv$specimens)[1]
    session$setInputs(edit_which = id, edit_name = "renamed", edit_group = "B", edit_apply = 1)
    expect_identical(rv$specimens[[id]]$group, "B")
    expect_gt(calls, 0L)                  # the wrapper is the copy the app calls
    expect_identical(outside, 0L)
  })
})
