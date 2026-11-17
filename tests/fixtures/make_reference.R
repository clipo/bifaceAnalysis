# Generates tests/fixtures/reference.rds. Run from the repository root:
#   Rscript tests/fixtures/make_reference.R
# Needs Momocs 1.4.1 (CRAN archive) and geomorph. The tests read the committed
# file and do not need either package.
#
# Contents:
#   masks      artifact masks of the four demo images from this app's pipeline.
#              They were inspected against the photographs when generated; the
#              test that uses them guards against later changes in segmentation.
#   points     the 200 corresponding points of each demo specimen (pixels)
#   aligned    those points after this app's Procrustes alignment
#   momocs_coe Momocs::efourier(norm = FALSE, nb.h = 20) of `aligned`. Momocs
#              places t = 0 at the last vertex, so each polygon is passed with
#              its first vertex moved to the end.
#   momocs_ellipse  the same for a 3-by-1 ellipse at 2000 vertices, nb.h = 5
#   gpagen     geomorph::gpagen(Proj = FALSE) coordinates of `points`
suppressMessages({library(imager); library(Momocs); library(geomorph)})
for (f in list.files("bifaceAnalysis/R", full.names = TRUE)) source(f)

files <- list.files("bifaceAnalysis/www", pattern = "(?i)\\.jpe?g$", full.names = TRUE)
names(files) <- tools::file_path_sans_ext(basename(files))
masks <- lapply(files, function(f) segment_image(load_working_image(f)$gray)$mask)
points <- lapply(files, function(f) process_image(load_working_image(f), direction = "left")$std$points)
aligned <- gpa(points)$aligned

roll <- function(p) p[c(2:nrow(p), 1L), , drop = FALSE]
momocs <- function(shapes, h) {
  m <- efourier(Out(lapply(shapes, roll)), nb.h = h, norm = FALSE)$coe
  rownames(m) <- names(shapes)
  m
}
th <- 2 * pi * (0:1999) / 2000
ellipse <- cbind(3 * cos(th), sin(th))

reference <- list(
  masks = masks, points = points, aligned = aligned,
  momocs_coe = momocs(aligned, 20),
  ellipse = ellipse, momocs_ellipse = momocs(list(e = ellipse), 5)[1, ],
  gpagen = gpagen(simplify2array(points), Proj = FALSE, print.progress = FALSE)$coords,
  versions = c(Momocs = as.character(packageVersion("Momocs")), geomorph = as.character(packageVersion("geomorph")),
               imager = as.character(packageVersion("imager")), R = as.character(getRversion())),
  generated = as.character(Sys.Date())
)
saveRDS(reference, "tests/fixtures/reference.rds", version = 2)
cat("wrote tests/fixtures/reference.rds\n"); print(reference$versions)

# Overlay for visual inspection of the outlines and landmarks (committed).
png("tests/fixtures/mask_check.png", 1600, 420, res = 110); par(mfrow = c(1, 4), mar = c(0, 0, 1.5, 0))
for (n in names(files)) {
  img <- load_working_image(files[[n]]); p <- process_image(img, direction = "left")
  plot(NA, xlim = c(0, img$width), ylim = c(img$height, 0), asp = 1, axes = FALSE, xlab = "", ylab = "", main = paste(n, paste(p$warnings, collapse = "|")))
  rasterImage(img$raster, 0, img$height, img$width, 0); polygon(p$outline_plot, border = "red", lwd = 2)
  points(rbind(p$tip_plot, p$base_plot), pch = 4, col = "cyan", cex = 2, lwd = 3)
}
invisible(dev.off())
