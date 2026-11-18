# Comparison of collections (specimens sharing a group label) and their
# locations. Statistics are descriptive. Every reference is a reassignment of
# the same fixed specimens or positions, not an uncertainty interval.
# Specification: docs/specs/2026-11-17_collections_spec.md.

MIN_COLLECTION_SIZE <- 3L
MIN_LOCATED_FOR_DISTANCE <- 5L
EARTH_RADIUS_KM <- 6371.0088
KM_PER_DEGREE <- pi * EARTH_RADIUS_KM / 180
MAX_MAP_LATITUDE <- 80
MAX_LOCATION_ROWS <- 200L
US_BOX <- list(lon = c(-125, -66), lat = c(24, 50))

# Group labels as used everywhere in the app: text, trimmed, never missing.
trim_label <- function(x) {
  x <- as.character(x)
  x[is.na(x)] <- ""
  trimws(x)
}

# Evaluate `expr` with the random number generator seeded, then restore the
# caller's random number state (or its absence).
with_seed <- function(seed, expr) {
  had <- exists(".Random.seed", envir = globalenv())
  old <- if (had) get(".Random.seed", envir = globalenv()) else NULL
  on.exit(if (had) assign(".Random.seed", old, envir = globalenv())
          else if (exists(".Random.seed", envir = globalenv())) rm(".Random.seed", envir = globalenv()))
  set.seed(seed)
  expr
}

summarize_reference <- function(values) {
  list(values = values, mean = mean(values),
       low = unname(stats::quantile(values, 0.025, names = FALSE)),
       high = unname(stats::quantile(values, 0.975, names = FALSE)))
}

# ---- collections ------------------------------------------------------------

# Rows of the full coefficient matrix that belong to eligible collections.
# `labels` are the group labels of the rows of `coefs`, in the same order.
# Returns list(X, group (factor), sizes, eligible, unlabelled, small, rows).
collection_data <- function(coefs, labels, min_size = MIN_COLLECTION_SIZE) {
  labels <- trim_label(labels)
  counts <- table(labels[nzchar(labels)])
  eligible <- sort(names(counts)[counts >= min_size])
  small <- counts[counts < min_size]
  rows <- which(labels %in% eligible)
  group <- factor(labels[rows], levels = eligible)
  list(X = coefs[rows, , drop = FALSE], group = group, sizes = table(group),
       eligible = eligible, unlabelled = sum(!nzchar(labels)),
       small = stats::setNames(as.integer(small), names(small)), rows = rows)
}

collection_means <- function(X, group) {
  m <- rowsum(X, as.character(group))[levels(group), , drop = FALSE]
  m / as.vector(table(group))
}

# Share of the total sum of squares that lies between collections, SSB / SST.
# NA when the rows do not differ.
share_between <- function(X, group) {
  if (no_variation(X)) return(NA_real_)
  centre <- colMeans(X)
  sst <- sum(sweep(X, 2, centre)^2)
  m <- collection_means(X, group)
  ssw <- sum((X - m[as.integer(group), , drop = FALSE])^2)
  (sst - ssw) / sst
}

# Observed share with the reference from reassigning specimens to collections
# at random, keeping the collection sizes.
share_reference <- function(X, group, n_perm = 999L, seed = 1L) {
  obs <- share_between(X, group)
  if (is.na(obs)) return(list(share = NA_real_, expected = NA_real_, reference = NULL))
  vals <- with_seed(seed, vapply(seq_len(n_perm), function(i) share_between(X, group[sample.int(length(group))]), numeric(1)))
  G <- nlevels(group); N <- nrow(X)
  list(share = obs, expected = (G - 1) / (N - 1), reference = summarize_reference(vals))
}

# Within-collection shape variance with divisor n - 1 (sum over variables).
row_variance <- function(X) sum(sweep(X, 2, colMeans(X))^2) / (nrow(X) - 1L)

within_variance <- function(X, group) {
  vapply(levels(group), function(g) row_variance(X[group == g, , drop = FALSE]), numeric(1))
}

# Smallest and largest within-collection variance with each specimen left out.
loo_variance_range <- function(X) {
  v <- vapply(seq_len(nrow(X)), function(i) row_variance(X[-i, , drop = FALSE]), numeric(1))
  c(low = min(v), high = max(v))
}

within_table <- function(X, group) {
  v <- within_variance(X, group)
  r <- t(vapply(levels(group), function(g) loo_variance_range(X[group == g, , drop = FALSE]), numeric(2)))
  data.frame(collection = levels(group), n = as.integer(table(group)), variance = unname(v),
             loo_low = unname(r[, 1]), loo_high = unname(r[, 2]), stringsAsFactors = FALSE)
}

# Raw and corrected distances between collection means for one pair of
# matrices, with Q the corrected squared distance.
pair_distance <- function(A, B) {
  d2 <- sum((colMeans(A) - colMeans(B))^2)
  q <- d2 - row_variance(A) / nrow(A) - row_variance(B) / nrow(B)
  c(raw = sqrt(d2), square = q, corrected = sqrt(max(0, q)))
}

# One row per pair of collections. The reference reassigns the specimens of
# the two collections between them at random, keeping the two sizes, and
# recomputes the raw distance.
pair_table <- function(X, group, n_perm = 999L, seed = 1L) {
  lv <- levels(group)
  pairs <- utils::combn(lv, 2, simplify = FALSE)
  with_seed(seed, do.call(rbind, lapply(pairs, function(p) {
    A <- X[group == p[1], , drop = FALSE]; B <- X[group == p[2], , drop = FALSE]
    d <- pair_distance(A, B)
    pool <- rbind(A, B); na <- nrow(A); n <- nrow(pool)
    ref <- vapply(seq_len(n_perm), function(i) {
      s <- sample.int(n)
      sqrt(sum((colMeans(pool[s[seq_len(na)], , drop = FALSE]) - colMeans(pool[s[-seq_len(na)], , drop = FALSE]))^2))
    }, numeric(1))
    r <- summarize_reference(ref)
    data.frame(collection_a = p[1], collection_b = p[2], n_a = na, n_b = nrow(B),
               raw_distance = d[["raw"]], reference_mean = r$mean, reference_low = r$low,
               reference_high = r$high, corrected_square = d[["square"]],
               corrected_distance = d[["corrected"]], stringsAsFactors = FALSE)
  })))
}

# Per-collection summary. `scores` are PC scores of the rows of X (may have
# 0 to 3 or more columns); `meas` are measurement rows in the same order.
collection_summary <- function(group, scores, meas, variance, located = NULL) {
  lv <- levels(group)
  pc_mean <- function(k) vapply(lv, function(g) if (ncol(scores) >= k) mean(scores[group == g, k]) else NA_real_, numeric(1))
  all_scaled <- all(is.finite(meas$length_mm))
  out <- data.frame(
    collection = lv, n = as.integer(table(group)),
    pc1_mean = unname(pc_mean(1)), pc2_mean = unname(pc_mean(2)), pc3_mean = unname(pc_mean(3)),
    length_width_mean = unname(vapply(lv, function(g) mean(meas$length_width_ratio[group == g]), numeric(1))),
    max_width_position_mean = unname(vapply(lv, function(g) mean(meas$widest_position[group == g]), numeric(1))),
    within_variance = unname(variance[lv]),
    length_mm_mean = if (all_scaled) unname(vapply(lv, function(g) mean(meas$length_mm[group == g]), numeric(1))) else NA_real_,
    located = if (is.null(located)) FALSE else lv %in% located,
    stringsAsFactors = FALSE)
  out
}

# ---- geography --------------------------------------------------------------

# Great-circle distance in km (haversine) between points given in degrees.
haversine_km <- function(lat1, lon1, lat2, lon2) {
  r <- pi / 180
  h <- sin((lat2 - lat1) * r / 2)^2 + cos(lat1 * r) * cos(lat2 * r) * sin((lon2 - lon1) * r / 2)^2
  2 * EARTH_RADIUS_KM * asin(pmin(1, sqrt(h)))
}

geo_distance_matrix <- function(lat, lon, names) {
  n <- length(lat)
  m <- outer(seq_len(n), seq_len(n), function(i, j) haversine_km(lat[i], lon[i], lat[j], lon[j]))
  dimnames(m) <- list(names, names)
  m
}

# Pearson correlation of the lower triangles of two distance matrices.
# NA when either set of distances does not vary.
matrix_correlation <- function(a, b) {
  x <- a[lower.tri(a)]; y <- b[lower.tri(b)]
  if (stats::sd(x) == 0 || stats::sd(y) == 0) return(NA_real_)
  stats::cor(x, y)
}

all_permutations <- function(n) {
  if (n == 1L) return(matrix(1L, 1, 1))
  p <- all_permutations(n - 1L)
  do.call(rbind, lapply(seq_len(n), function(i) cbind(i, matrix(c(seq_len(n))[-i][p], nrow = nrow(p)))))
}

# Matrix correlation between corrected shape distance and geographic distance
# with the reference from reassigning positions among the located collections:
# all assignments when there are at most 7, otherwise 999 random ones.
# Returns list(available, reason, r, reference, n_assignments, n_pairs).
shape_geo_correlation <- function(shape_d, geo_d, positions, n_perm = 999L, seed = 1L, min_located = MIN_LOCATED_FOR_DISTANCE) {
  G <- nrow(shape_d)
  distinct <- nrow(unique(round(positions, 12)))
  if (G < min_located)
    return(list(available = FALSE, reason = sprintf("The comparison of shape distance with geographic distance needs at least %d located collections. %d %s located.", min_located, G, if (G == 1) "is" else "are")))
  if (distinct < 2)
    return(list(available = FALSE, reason = "All located collections have the same position, so geographic distance does not vary."))
  r <- matrix_correlation(shape_d, geo_d)
  if (is.na(r))
    return(list(available = FALSE, reason = "Shape distance does not vary between the located collections, so no correlation can be computed."))
  perms <- if (G <= 7L) all_permutations(G) else with_seed(seed, t(vapply(seq_len(n_perm), function(i) sample.int(G), integer(G))))
  vals <- apply(perms, 1, function(p) matrix_correlation(shape_d, geo_d[p, p]))
  list(available = TRUE, r = r, reference = summarize_reference(vals), n_assignments = nrow(perms),
       n_pairs = G * (G - 1L) / 2L)
}

# ---- locations file ---------------------------------------------------------

# Read a locations file. Stops with a plain message when the file cannot be
# used. Returns data.frame(collection, latitude, longitude).
read_locations <- function(path) {
  df <- tryCatch(utils::read.csv(path, colClasses = "character", na.strings = character(0),
                                 check.names = FALSE, strip.white = FALSE, encoding = "UTF-8"),
                 error = function(e) NULL)
  if (is.null(df)) stop("This file could not be read as a CSV file.")
  norm <- tolower(trimws(names(df)))
  need <- c("collection", "latitude", "longitude")
  missing <- need[!need %in% norm]
  if (length(missing)) stop("The locations file must have the columns collection, latitude, and longitude. Missing: ", paste(missing, collapse = ", "), ".")
  repeated <- need[vapply(need, function(n) sum(norm == n) > 1L, logical(1))]
  if (length(repeated)) stop("The locations file has more than one column named ", paste(repeated, collapse = ", "), ".")
  if (nrow(df) == 0L) stop("The locations file has no rows.")
  if (nrow(df) > MAX_LOCATION_ROWS) stop("The locations file has more than ", MAX_LOCATION_ROWS, " rows.")
  col <- function(n) df[[which(norm == n)]]
  name <- trimws(col("collection"))
  if (any(!nzchar(name))) stop("Every row of the locations file must have a collection name.")
  if (anyDuplicated(name)) stop("The locations file names a collection more than once: ", paste(unique(name[duplicated(name)]), collapse = ", "), ".")
  num <- function(n) {
    v <- trimws(col(n))
    x <- suppressWarnings(as.numeric(v))
    if (any(!nzchar(v)) || anyNA(x) || any(!is.finite(x)))
      stop("Every ", n, " in the locations file must be a number.")
    x
  }
  lat <- num("latitude"); lon <- num("longitude")
  if (any(lat < -90 | lat > 90)) stop("Latitudes must lie between -90 and 90.")
  if (any(lon < -180 | lon > 180)) stop("Longitudes must lie between -180 and 180.")
  data.frame(collection = name, latitude = lat, longitude = lon, stringsAsFactors = FALSE)
}

# Match location rows to eligible collections by exact (trimmed) name.
match_locations <- function(locations, eligible) {
  if (is.null(locations))
    return(list(located = data.frame(collection = character(0), latitude = numeric(0), longitude = numeric(0)),
                unmatched_rows = character(0), not_located = eligible))
  hit <- locations$collection %in% eligible
  located <- locations[hit, , drop = FALSE]
  located <- located[order(match(located$collection, eligible)), , drop = FALSE]
  rownames(located) <- NULL
  list(located = located, unmatched_rows = locations$collection[!hit],
       not_located = setdiff(eligible, locations$collection))
}

# Template: one row per non-empty label, coordinates where a location is loaded.
locations_template <- function(labels, locations = NULL) {
  lv <- sort(unique(trim_label(labels)))
  lv <- lv[nzchar(lv)]
  i <- if (is.null(locations)) rep(NA_integer_, length(lv)) else match(lv, locations$collection)
  data.frame(collection = lv,
             latitude = if (is.null(locations)) rep(NA_real_, length(lv)) else locations$latitude[i],
             longitude = if (is.null(locations)) rep(NA_real_, length(lv)) else locations$longitude[i],
             stringsAsFactors = FALSE)
}

# ---- map geometry -----------------------------------------------------------

# Window for the map, or the reason it is withheld.
map_window <- function(lat, lon) {
  if (length(lat) < 2L) return(list(ok = FALSE, reason = "The map needs at least 2 located collections."))
  if (any(abs(lat) > MAX_MAP_LATITUDE))
    return(list(ok = FALSE, reason = sprintf("The map is not drawn because a location lies beyond %d degrees of latitude, where this map projection is strongly distorted. Distances and other results are not affected.", MAX_MAP_LATITUDE)))
  if (diff(range(lon)) > 180)
    return(list(ok = FALSE, reason = "The map is not drawn because the locations span more than 180 degrees of longitude, which this map cannot show across the 180 degree meridian. Distances and other results are not affected."))
  lat_mid <- mean(range(lat))
  lat_span <- max(1.4 * diff(range(lat)), 0.5)
  lon_span <- max(1.4 * diff(range(lon)), 0.5 / cos(lat_mid * pi / 180))
  lon_mid <- mean(range(lon))
  ylim <- c(max(-85, lat_mid - lat_span / 2), min(85, lat_mid + lat_span / 2))
  xlim <- c(max(-180, lon_mid - lon_span / 2), min(180, lon_mid + lon_span / 2))
  list(ok = TRUE, xlim = xlim, ylim = ylim, lat_mid = lat_mid)
}

# Background databases drawn for a window.
map_layers <- function(window) {
  overlaps_us <- window$xlim[1] <= US_BOX$lon[2] && window$xlim[2] >= US_BOX$lon[1] &&
    window$ylim[1] <= US_BOX$lat[2] && window$ylim[2] >= US_BOX$lat[1]
  c("world", if (overlaps_us) "state", if (overlaps_us && diff(window$ylim) < 3) "county")
}

nice_length <- function(target) {
  k <- floor(log10(target))
  cand <- c(1, 2, 5, 10) * 10^k
  cand[which.min(abs(cand - target))]
}

# Scale bar: length L in km and its endpoints. L is drawn along the parallel
# at 8 percent of the window height above the bottom, starting 5 percent of
# the width from the left.
scale_bar <- function(window) {
  width_km <- haversine_km(window$lat_mid, window$xlim[1], window$lat_mid, window$xlim[2])
  L <- nice_length(width_km / 4)
  y <- window$ylim[1] + 0.08 * diff(window$ylim)
  x0 <- window$xlim[1] + 0.05 * diff(window$xlim)
  list(km = L, x0 = x0, x1 = x0 + L / (KM_PER_DEGREE * cos(y * pi / 180)), y = y)
}

MAP_CACHE <- new.env()
map_polygons <- function(database) {
  if (is.null(MAP_CACHE[[database]])) MAP_CACHE[[database]] <- ggplot2::map_data(database)
  MAP_CACHE[[database]]
}
