# Load the app's functions and shared fixtures and shapes for the tests.
APP_DIR <- normalizePath(testthat::test_path("..", "..", "bifaceAnalysis"))
for (f in list.files(file.path(APP_DIR, "R"), full.names = TRUE)) source(f)
REF <- readRDS(testthat::test_path("..", "fixtures", "reference.rds"))
DEMO_FILES <- list.files(file.path(APP_DIR, "www"), pattern = "(?i)\\.jpe?g$", full.names = TRUE)
names(DEMO_FILES) <- tools::file_path_sans_ext(basename(DEMO_FILES))

rotate_deg <- function(xy, deg) {
  a <- deg * pi / 180
  xy %*% t(matrix(c(cos(a), sin(a), -sin(a), cos(a)), 2, 2))
}
shift_start <- function(xy, k) if (k == 0) xy else xy[c((k + 1):nrow(xy), 1:k), , drop = FALSE]
# Insert n - 1 evenly spaced points on every edge of a polygon.
densify <- function(xy, n = 10L) {
  nxt <- xy[c(2:nrow(xy), 1L), , drop = FALSE]
  do.call(rbind, lapply(seq_len(nrow(xy)), function(i)
    cbind(xy[i, 1] + (0:(n - 1)) / n * (nxt[i, 1] - xy[i, 1]), xy[i, 2] + (0:(n - 1)) / n * (nxt[i, 2] - xy[i, 2]))))
}
# Standardize a y-up polygon (prepare_outline expects image coordinates, y down).
standardize_yup <- function(xy, direction = "up", ...) standardize_outline(prepare_outline(cbind(xy[, 1], -xy[, 2])), direction, ...)

# Counterclockwise, tip up.
KITE <- rbind(c(0, 10), c(-3, 2), c(-1.2, 1.5), c(-1.2, -2), c(1.2, -2), c(1.2, 1.5), c(3, 2))
KITE_ASYM <- rbind(c(0, 10), c(-3, 2), c(-1.2, 1.5), c(-1.2, -2), c(1.2, -2), c(1.2, 1.5), c(3.6, 2.6))
WIDE <- rbind(c(0, 3), c(-5, 0), c(-4, -1), c(4, -1), c(5, 0))
NOTCHED <- rbind(c(0, 10), c(-3, 0), c(-2, -3), c(0, -1), c(2, -3), c(3, 0))
BARBED <- rbind(c(0, 10), c(-3, 0), c(-3, -6), c(2, -6), c(2, -4), c(-1, -4), c(-1, -2), c(3, -2), c(3, 0))
U_SHAPE <- rbind(c(0, 10), c(-3, 0), c(-3, -6), c(-1, -6), c(-1, -2), c(1, -2), c(1, -6), c(3, -6), c(3, 0))
TRUNCATED <- rbind(c(-1.5, 10), c(-3, 0), c(3, 0), c(1.5, 10))

is_validation <- function(e) inherits(e, "shiny.silent.error") || inherits(e, "validation")
