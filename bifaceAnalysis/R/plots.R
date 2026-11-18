# One plotting function per figure. Each returns a ggplot object that is used
# for both the screen and the PDF report. Figures carry no titles.

OKABE_ITO <- c("#0072B2", "#E69F00", "#009E73", "#CC79A7", "#D55E00", "#56B4E9", "#F0E442", "#000000")
REGION_COLOURS <- c(Tip = "#D55E00", Blade = "#0072B2", Shoulder = "#009E73", Base = "#CC79A7")

theme_app <- function(base_size = 13) {
  ggplot2::theme_minimal(base_size = base_size, base_family = "sans") +
    ggplot2::theme(panel.grid.minor = ggplot2::element_blank(),
                   panel.grid.major = ggplot2::element_line(colour = "grey92"),
                   legend.position = "bottom")
}

theme_shape <- function(base_size = 13) {
  ggplot2::theme_void(base_size = base_size, base_family = "sans") +
    ggplot2::theme(legend.position = "bottom",
                   strip.text = ggplot2::element_text(size = base_size * 0.85))
}

# Colours for categories (group labels, collections, clusters): Okabe-Ito for
# up to 8 levels, discrete viridis for more, so that no colour repeats.
category_colours <- function(values) {
  lv <- levels(factor(values))
  cols <- if (length(lv) <= length(OKABE_ITO)) OKABE_ITO[seq_along(lv)] else scales::viridis_pal()(length(lv))
  stats::setNames(cols, lv)
}
category_scale <- function(values, aesthetic = "colour", name = NULL) {
  ggplot2::scale_discrete_manual(aesthetics = aesthetic, values = category_colours(values), name = name)
}

closed <- function(m) rbind(m, m[1, , drop = FALSE])

shapes_frame <- function(shapes) {
  do.call(rbind, lapply(names(shapes), function(n) {
    m <- closed(shapes[[n]])
    data.frame(specimen = n, x = m[, 1], y = m[, 2], stringsAsFactors = FALSE)
  }))
}

path_frame <- function(m) {
  m <- closed(m)
  data.frame(x = m[, 1], y = m[, 2])
}

has_labels <- function(groups) !is.null(groups) && any(nzchar(as.character(groups)))

# Each specimen's oriented outline, drawn to a common height.
plot_outline_grid <- function(outlines, groups = NULL) {
  norm <- lapply(outlines, function(o) {
    o <- sweep(o, 2, c(mean(range(o[, 1])), min(o[, 2])))
    o / max(o[, 2])
  })
  df <- shapes_frame(norm)
  df$specimen <- factor(df$specimen, levels = names(outlines))
  df$group <- if (has_labels(groups)) groups[match(df$specimen, names(outlines))] else "all"
  p <- ggplot2::ggplot(df, ggplot2::aes(x, y, group = specimen, fill = group)) +
    ggplot2::geom_polygon(colour = "grey20", linewidth = 0.4, alpha = 0.55) +
    ggplot2::facet_wrap(~specimen) + ggplot2::coord_equal() + theme_shape()
  if (has_labels(groups)) p + category_scale(df$group, "fill", "Group")
  else p + ggplot2::scale_fill_manual(values = c(all = "grey75"), guide = "none")
}

# All aligned shapes overlaid, with the mean shape in black.
plot_outline_stack <- function(aligned, mean_shape, groups = NULL) {
  df <- shapes_frame(aligned)
  df$group <- if (has_labels(groups)) groups[match(df$specimen, names(aligned))] else "all"
  p <- ggplot2::ggplot(df, ggplot2::aes(x, y)) +
    ggplot2::geom_path(ggplot2::aes(group = specimen, colour = group), linewidth = 0.5, alpha = 0.7) +
    ggplot2::geom_path(data = path_frame(mean_shape), colour = "black", linewidth = 1.1) +
    ggplot2::coord_equal() + theme_shape()
  if (has_labels(groups)) p + category_scale(df$group, "colour", "Group")
  else p + ggplot2::scale_colour_manual(values = c(all = "grey55"), guide = "none")
}

# Length against maximum width, in millimetres when every specimen has a scale.
plot_measurements <- function(meas, groups = NULL) {
  use_mm <- all(is.finite(meas$length_mm))
  df <- data.frame(specimen = meas$specimen,
                   length = if (use_mm) meas$length_mm else meas$length_px,
                   width = if (use_mm) meas$width_mm else meas$width_px)
  unit <- if (use_mm) "mm" else "pixels"
  df$group <- if (has_labels(groups)) groups else "all"
  p <- ggplot2::ggplot(df, ggplot2::aes(width, length, colour = group, label = specimen)) +
    ggplot2::geom_point(size = 3) +
    ggplot2::geom_text(vjust = -0.9, size = 3.5, family = "sans", show.legend = FALSE) +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = 0.12)) +
    ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = 0.12)) +
    ggplot2::labs(x = sprintf("Maximum width (%s)", unit), y = sprintf("Length (%s)", unit)) +
    theme_app()
  if (has_labels(groups)) p + category_scale(df$group, "colour", "Group")
  else p + ggplot2::scale_colour_manual(values = c(all = OKABE_ITO[1]), guide = "none")
}

# Share of the power beyond the first harmonic held by each harmonic (bars) and
# the cumulative share (line). The dashed line marks the number retained.
plot_harmonic_power <- function(power_table, n_used) {
  df <- power_table[-1, ]
  df$share <- diff(c(0, df$cumulative_excluding_first))
  ggplot2::ggplot(df, ggplot2::aes(harmonic)) +
    ggplot2::geom_col(ggplot2::aes(y = share), fill = OKABE_ITO[6], width = 0.8) +
    ggplot2::geom_line(ggplot2::aes(y = cumulative_excluding_first), colour = OKABE_ITO[5], linewidth = 0.8) +
    ggplot2::geom_point(ggplot2::aes(y = cumulative_excluding_first), colour = OKABE_ITO[5], size = 1.8) +
    ggplot2::geom_vline(xintercept = n_used + 0.5, linetype = "dashed", colour = "grey30") +
    ggplot2::scale_y_continuous(labels = function(v) paste0(round(100 * v), "%"), limits = c(0, 1.001)) +
    ggplot2::labs(x = "Harmonic", y = "Share of power beyond harmonic 1") +
    theme_app()
}

# One specimen rebuilt from increasing numbers of harmonics, over its outline.
plot_reconstruction <- function(shape, counts) {
  counts <- sort(unique(as.integer(counts)))
  full <- efa_forward(shape, max(counts))
  vec <- efa_vector(full)
  recs <- lapply(counts, function(h) efa_inverse(vec, 300L, h, full$x0, full$y0))
  names(recs) <- sprintf("%d harmonic%s", counts, ifelse(counts == 1, "", "s"))
  df <- shapes_frame(recs)
  df$specimen <- factor(df$specimen, levels = names(recs))
  ggplot2::ggplot(df, ggplot2::aes(x, y)) +
    ggplot2::geom_polygon(data = path_frame(shape), fill = "grey85", colour = NA) +
    ggplot2::geom_path(colour = OKABE_ITO[5], linewidth = 0.7) +
    ggplot2::facet_wrap(~specimen, nrow = 1) + ggplot2::coord_equal() + theme_shape()
}

plot_scree <- function(pca) {
  df <- data.frame(pc = seq_along(pca$explained), share = pca$explained)
  ggplot2::ggplot(df, ggplot2::aes(pc, share)) +
    ggplot2::geom_col(fill = OKABE_ITO[6], width = 0.7) +
    ggplot2::scale_x_continuous(breaks = df$pc) +
    ggplot2::scale_y_continuous(labels = function(v) paste0(round(100 * v), "%")) +
    ggplot2::labs(x = "Principal component", y = "Share of shape variance") + theme_app()
}

# Shape implied by a position in PC space (scores on the leading components),
# drawn centred (constant terms zero).
pc_shape <- function(pca, scores, n_points = 200L) {
  k <- seq_along(scores)
  efa_inverse(pca$center + as.numeric(pca$loadings[, k, drop = FALSE] %*% scores), n_points)
}

# PC scores with reconstructed shapes on a background grid.
plot_morphospace <- function(pca, colour_by = NULL, colour_name = NULL, axes = c(1L, 2L)) {
  sc <- pca$scores
  df <- data.frame(specimen = rownames(sc), a = sc[, axes[1]], b = sc[, axes[2]])
  df$colour <- if (has_labels(colour_by)) as.character(colour_by) else "all"
  pad <- function(r) { w <- diff(r); if (w == 0) w <- 1; r + c(-1, 1) * 0.18 * w }
  ra <- pad(range(df$a)); rb <- pad(range(df$b))
  ga <- seq(ra[1], ra[2], length.out = 5); gb <- seq(rb[1], rb[2], length.out = 4)
  cell <- min(diff(ga)[1], diff(gb)[1])
  bg <- do.call(rbind, lapply(ga, function(va) do.call(rbind, lapply(gb, function(vb) {
    s <- numeric(max(axes)); s[axes[1]] <- va; s[axes[2]] <- vb
    m <- closed(pc_shape(pca, s, 120L))
    m <- sweep(m, 2, colMeans(m))
    m <- m * 0.42 * cell / max(apply(m, 2, function(v) diff(range(v))))
    data.frame(id = paste(va, vb), x = va + m[, 1], y = vb + m[, 2])
  }))))
  lab <- function(i) sprintf("PC%d (%.0f%% of shape variance)", i, 100 * pca$explained[i])
  p <- ggplot2::ggplot() +
    ggplot2::geom_polygon(data = bg, ggplot2::aes(x, y, group = id), fill = "grey90", colour = "grey70", linewidth = 0.3) +
    ggplot2::geom_point(data = df, ggplot2::aes(a, b, colour = colour), size = 3) +
    ggplot2::geom_text(data = df, ggplot2::aes(a, b, label = specimen, colour = colour), vjust = -1, size = 3.5,
                       family = "sans", show.legend = FALSE) +
    ggplot2::coord_equal() + ggplot2::labs(x = lab(axes[1]), y = lab(axes[2])) + theme_app()
  if (has_labels(colour_by)) p + category_scale(df$colour, "colour", colour_name)
  else p + ggplot2::scale_colour_manual(values = c(all = OKABE_ITO[1]), guide = "none")
}

# Shapes at the mean and at -n_sd and +n_sd standard deviations on a component.
pc_extreme_shapes <- function(pca, pc, n_sd = 2) {
  s <- numeric(pc)
  lo <- s; lo[pc] <- -n_sd * pca$sdev[pc]
  hi <- s; hi[pc] <- n_sd * pca$sdev[pc]
  list(low = pc_shape(pca, lo), mean = pc_shape(pca, s), high = pc_shape(pca, hi))
}

# Mean, low-end, and high-end shapes for the first min(3, rank) components.
plot_pc_shapes <- function(pca) {
  pcs <- seq_len(min(3L, pca$rank))
  ends <- c("Low end (-2 SD)", "Mean", "High end (+2 SD)")
  df <- do.call(rbind, lapply(pcs, function(pc) {
    e <- pc_extreme_shapes(pca, pc)
    lab <- sprintf("PC%d (%.0f%%)", pc, 100 * pca$explained[pc])
    rbind(data.frame(pc = lab, end = ends[1], path_frame(e$low)),
          data.frame(pc = lab, end = ends[2], path_frame(e$mean)),
          data.frame(pc = lab, end = ends[3], path_frame(e$high)))
  }))
  df$end <- factor(df$end, levels = ends)
  ggplot2::ggplot(df, ggplot2::aes(x, y, colour = end, linetype = end, group = end)) +
    ggplot2::geom_path(linewidth = 0.9) +
    ggplot2::scale_colour_manual(values = stats::setNames(c(OKABE_ITO[1], "grey55", OKABE_ITO[5]), ends), name = NULL) +
    ggplot2::scale_linetype_manual(values = stats::setNames(c("solid", "dashed", "solid"), ends), name = NULL) +
    ggplot2::facet_wrap(~pc, nrow = 1) + ggplot2::coord_equal() + theme_shape()
}

# Leave-one-out agreement: the PC1-PC2 plane, and the direction of each
# component that is separable from its neighbours.
plot_stability <- function(stab) {
  df <- data.frame(left_out = names(stab$plane), quantity = "Plane of PC1 and PC2", value = stab$plane)
  for (a in names(stab$separable)[stab$separable])
    df <- rbind(df, data.frame(left_out = rownames(stab$axis), quantity = paste(a, "direction"), value = stab$axis[, a]))
  df$quantity <- factor(df$quantity, levels = unique(df$quantity))
  ggplot2::ggplot(df, ggplot2::aes(value, left_out)) +
    ggplot2::geom_segment(ggplot2::aes(x = 0, xend = value, yend = left_out), colour = "grey75") +
    ggplot2::geom_point(colour = OKABE_ITO[1], size = 2.5) +
    ggplot2::facet_wrap(~quantity, nrow = 1) + ggplot2::scale_x_continuous(limits = c(0, 1)) +
    ggplot2::labs(x = "Agreement with the full analysis (1 = unchanged)", y = "Specimen left out") + theme_app()
}

# Segment coordinates for drawing an hclust tree.
dendrogram_segments <- function(tree) {
  n <- length(tree$order)
  xpos <- numeric(n); xpos[tree$order] <- seq_len(n)
  node_x <- numeric(n - 1L); node_h <- tree$height
  segs <- vector("list", n - 1L)
  for (i in seq_len(n - 1L)) {
    ch <- tree$merge[i, ]
    cx <- ifelse(ch < 0, xpos[abs(ch)], node_x[pmax(ch, 1L)])
    chh <- ifelse(ch < 0, 0, node_h[pmax(ch, 1L)])
    node_x[i] <- mean(cx)
    segs[[i]] <- data.frame(x = c(cx[1], cx[2], cx[1]), xend = c(cx[1], cx[2], cx[2]),
                            y = c(chh[1], chh[2], node_h[i]), yend = rep(node_h[i], 3))
  }
  list(segments = do.call(rbind, segs), leaves = data.frame(x = seq_len(n), label = tree$labels[tree$order]))
}

# Height at which a horizontal cut of the tree gives k clusters: midway
# between the merges that leave k - 1 and k clusters. NA when those two merges
# are at the same height (within NUMERICAL_ZERO of the tallest merge), because
# no horizontal cut then gives exactly k clusters.
tree_cut_height <- function(tree, k) {
  h <- sort(tree$height, decreasing = TRUE)
  if (k < 2L || k > length(h)) return(NA_real_)
  if (h[k - 1L] - h[k] <= NUMERICAL_ZERO * max(h[1], .Machine$double.xmin)) return(NA_real_)
  mean(h[c(k - 1L, k)])
}

# Tree with leaves coloured by cluster and, when one exists, a dashed line at
# the height where a cut gives k clusters.
plot_dendrogram <- function(tree, k) {
  d <- dendrogram_segments(tree)
  cl <- stats::cutree(tree, k)
  d$leaves$cluster <- factor(cl[d$leaves$label])
  cut_h <- tree_cut_height(tree, k)
  p <- ggplot2::ggplot() +
    ggplot2::geom_segment(data = d$segments, ggplot2::aes(x = x, xend = xend, y = y, yend = yend), colour = "grey25") +
    ggplot2::geom_point(data = d$leaves, ggplot2::aes(x, 0, colour = cluster), size = 3) +
    ggplot2::scale_x_continuous(breaks = d$leaves$x, labels = d$leaves$label) +
    ggplot2::labs(x = NULL, y = "Distance at which groups merge") +
    category_scale(d$leaves$cluster, "colour", "Cluster") + theme_app() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
                   panel.grid.major.x = ggplot2::element_blank())
  if (is.finite(cut_h)) p <- p + ggplot2::geom_hline(yintercept = cut_h, linetype = "dashed", colour = "grey40")
  p
}

plot_silhouette <- function(sil, k) {
  sil$chosen <- sil$k == k
  ggplot2::ggplot(sil, ggplot2::aes(k, silhouette, fill = chosen)) +
    ggplot2::geom_col(width = 0.7) +
    ggplot2::scale_fill_manual(values = c(`TRUE` = OKABE_ITO[5], `FALSE` = OKABE_ITO[6]), guide = "none") +
    ggplot2::scale_x_continuous(breaks = sil$k) +
    ggplot2::labs(x = "Number of clusters", y = "Mean silhouette width") + theme_app()
}

# Mean shape with points coloured by region.
plot_regions <- function(mean_shape, regions) {
  df <- data.frame(x = mean_shape[, 1], y = mean_shape[, 2], region = regions)
  ggplot2::ggplot(df, ggplot2::aes(x, y)) +
    ggplot2::geom_path(data = path_frame(mean_shape), colour = "grey60") +
    ggplot2::geom_point(ggplot2::aes(colour = region), size = 1.8) +
    ggplot2::scale_colour_manual(values = REGION_COLOURS, name = "Region", drop = FALSE) +
    ggplot2::coord_equal() + theme_shape()
}

# Mean shape with points coloured by their variance across specimens.
plot_point_variance <- function(mean_shape, variance) {
  df <- data.frame(x = mean_shape[, 1], y = mean_shape[, 2], variance = variance)
  ggplot2::ggplot(df, ggplot2::aes(x, y)) +
    ggplot2::geom_path(data = path_frame(mean_shape), colour = "grey60") +
    ggplot2::geom_point(ggplot2::aes(colour = variance), size = 2.2) +
    ggplot2::scale_colour_viridis_c(name = "Variance at each point", option = "viridis") +
    ggplot2::coord_equal() + theme_shape() +
    ggplot2::theme(legend.key.width = grid::unit(1.4, "cm"),
                   legend.text = ggplot2::element_text(angle = 45, hjust = 1))
}

plot_region_variance <- function(region_variance) {
  df <- data.frame(region = factor(names(region_variance), levels = names(region_variance)),
                   variance = as.numeric(region_variance))
  df <- df[is.finite(df$variance), ]
  ggplot2::ggplot(df, ggplot2::aes(region, variance, fill = region)) +
    ggplot2::geom_col(width = 0.7) +
    ggplot2::scale_fill_manual(values = REGION_COLOURS, guide = "none") +
    ggplot2::labs(x = NULL, y = "Mean variance per point") + theme_app()
}

# Observed RV for each pair of regions beside its re-pairing reference.
plot_region_rv <- function(rv_table) {
  rv_table <- rv_table[is.finite(rv_table$rv), ]
  rv_table$pair <- paste(rv_table$region_a, "and", rv_table$region_b)
  rv_table$pair <- factor(rv_table$pair, levels = rev(rv_table$pair))
  ggplot2::ggplot(rv_table, ggplot2::aes(y = pair)) +
    ggplot2::geom_linerange(ggplot2::aes(xmin = reference_low, xmax = reference_high, colour = "ref"), linewidth = 4, alpha = 0.5) +
    ggplot2::geom_point(ggplot2::aes(x = reference_mean, colour = "ref"), shape = 124, size = 5) +
    ggplot2::geom_point(ggplot2::aes(x = rv, colour = "obs"), size = 3.5) +
    ggplot2::scale_colour_manual(values = c(obs = OKABE_ITO[5], ref = "grey55"), name = NULL,
                                 labels = c(obs = "Observed RV", ref = "RV with specimens re-paired at random (mean and central 95% of values)")) +
    ggplot2::scale_x_continuous(limits = c(0, 1)) +
    ggplot2::labs(x = "RV coefficient", y = NULL) + theme_app() +
    ggplot2::theme(legend.direction = "vertical")
}

# ---- collections ------------------------------------------------------------

# Mean shape of each collection, overlaid.
plot_collection_means <- function(shapes) {
  df <- shapes_frame(shapes)
  df$collection <- factor(df$specimen, levels = names(shapes))
  ggplot2::ggplot(df, ggplot2::aes(x, y, group = collection, colour = collection)) +
    ggplot2::geom_path(linewidth = 0.9) + ggplot2::coord_equal() +
    category_scale(df$collection, "colour", "Collection") + theme_shape()
}

# Observed value(s) beside a reassignment reference, one row per item.
reference_plot <- function(df, x_label, obs_label, ref_label, limits = NULL) {
  df$item <- factor(df$item, levels = rev(unique(df$item)))
  p <- ggplot2::ggplot(df, ggplot2::aes(y = item)) +
    ggplot2::geom_linerange(ggplot2::aes(xmin = low, xmax = high, colour = "ref"), linewidth = 4, alpha = 0.5) +
    ggplot2::geom_point(ggplot2::aes(x = mean, colour = "ref"), shape = 124, size = 5) +
    ggplot2::geom_point(ggplot2::aes(x = observed, colour = "obs"), size = 3.5) +
    ggplot2::scale_colour_manual(values = c(obs = OKABE_ITO[5], ref = "grey55"), name = NULL,
                                 labels = c(obs = obs_label, ref = ref_label)) +
    ggplot2::labs(x = x_label, y = NULL) + theme_app() + ggplot2::theme(legend.direction = "vertical")
  if (!is.null(limits)) p <- p + ggplot2::scale_x_continuous(limits = limits)
  p
}

plot_collection_share <- function(share) {
  r <- share$reference
  reference_plot(data.frame(item = "Between collections", observed = share$share, mean = r$mean, low = r$low, high = r$high),
                 "Share of shape variation between collections", "Observed share",
                 "Share with specimens reassigned to collections at random (mean and central 95% of values)", c(0, 1))
}

plot_pair_distances <- function(pairs) {
  reference_plot(data.frame(item = paste(pairs$collection_a, "and", pairs$collection_b), observed = pairs$raw_distance,
                            mean = pairs$reference_mean, low = pairs$reference_low, high = pairs$reference_high),
                 "Distance between collection means (raw)", "Observed distance",
                 "Distance with the two collections' specimens reassigned at random (mean and central 95% of values)",
                 c(0, max(c(pairs$raw_distance, pairs$reference_high))))
}

# Matrix of corrected distances between collection means.
plot_distance_matrix <- function(pairs, collections) {
  both <- rbind(data.frame(a = pairs$collection_a, b = pairs$collection_b, d = pairs$corrected_distance),
                data.frame(a = pairs$collection_b, b = pairs$collection_a, d = pairs$corrected_distance),
                data.frame(a = collections, b = collections, d = 0))
  both$a <- factor(both$a, levels = collections); both$b <- factor(both$b, levels = rev(collections))
  digits <- if (max(both$d) > 0) max(0, 2 - floor(log10(max(both$d)))) else 2
  ggplot2::ggplot(both, ggplot2::aes(a, b, fill = d)) + ggplot2::geom_tile(colour = "white") +
    ggplot2::geom_text(ggplot2::aes(label = formatC(d, format = "f", digits = digits)), size = 3.5, family = "sans",
                       colour = ifelse(both$d > 0.6 * max(both$d), "black", "white")) +
    ggplot2::scale_fill_viridis_c(name = "Corrected distance") + ggplot2::coord_equal() +
    ggplot2::labs(x = NULL, y = NULL) + theme_app() +
    ggplot2::theme(panel.grid.major = ggplot2::element_blank(), axis.text.x = ggplot2::element_text(angle = 30, hjust = 1))
}

# Within-collection variance with the range over leave-one-out values.
plot_within_variance <- function(within) {
  within$collection <- factor(within$collection, levels = rev(within$collection))
  ggplot2::ggplot(within, ggplot2::aes(y = collection, colour = collection)) +
    ggplot2::geom_linerange(ggplot2::aes(xmin = loo_low, xmax = loo_high), linewidth = 2, alpha = 0.5) +
    ggplot2::geom_point(ggplot2::aes(x = variance), size = 3.5) +
    category_scale(rev(levels(within$collection)), "colour", NULL) + ggplot2::guides(colour = "none") +
    ggplot2::expand_limits(x = 0) +
    ggplot2::labs(x = "Shape variance within the collection (bar: range with one specimen left out)", y = NULL) + theme_app()
}

# Corrected shape distance against geographic distance, one point per pair.
plot_shape_geo <- function(pairs_km) {
  pairs_km$label <- paste(pairs_km$collection_a, "-", pairs_km$collection_b)
  ggplot2::ggplot(pairs_km, ggplot2::aes(distance_km, corrected_distance)) +
    ggplot2::geom_point(colour = OKABE_ITO[1], size = 3) +
    ggplot2::geom_text(ggplot2::aes(label = label), vjust = -0.9, size = 3, family = "sans", colour = "grey30") +
    ggplot2::expand_limits(x = 0, y = 0) +
    ggplot2::labs(x = "Geographic distance (km)", y = "Corrected distance between collection means") + theme_app()
}

# Map of located collections coloured by a summary value.
plot_collection_map <- function(located, value, value_name, window) {
  located$value <- value
  bar <- scale_bar(window)
  layers <- map_layers(window)
  bg <- list(ggplot2::geom_polygon(data = map_polygons("world"), ggplot2::aes(long, lat, group = group),
                                   fill = "grey94", colour = "grey60", linewidth = 0.25))
  if ("state" %in% layers)
    bg <- c(bg, ggplot2::geom_polygon(data = map_polygons("state"), ggplot2::aes(long, lat, group = group),
                                      fill = NA, colour = "grey70", linewidth = 0.2))
  if ("county" %in% layers)
    bg <- c(bg, ggplot2::geom_polygon(data = map_polygons("county"), ggplot2::aes(long, lat, group = group),
                                      fill = NA, colour = "grey82", linewidth = 0.15))
  bar_df <- data.frame(x = bar$x0, xend = bar$x1, y = bar$y)
  ggplot2::ggplot() + bg +
    ggplot2::geom_point(data = located, ggplot2::aes(longitude, latitude, colour = value), size = 4) +
    ggplot2::geom_text(data = located, ggplot2::aes(longitude, latitude, label = collection),
                       vjust = -1.1, size = 3.5, family = "sans") +
    ggplot2::geom_segment(data = bar_df, ggplot2::aes(x = x, xend = xend, y = y, yend = y), linewidth = 1) +
    ggplot2::annotate("text", x = (bar$x0 + bar$x1) / 2, y = bar$y, label = sprintf("%s km (approximate)", format(bar$km, big.mark = ",", scientific = FALSE)),
                      vjust = -0.7, size = 3.2, family = "sans") +
    ggplot2::scale_colour_viridis_c(name = value_name) +
    ggplot2::coord_quickmap(xlim = window$xlim, ylim = window$ylim, expand = FALSE) +
    ggplot2::labs(x = "Longitude", y = "Latitude") + theme_app() +
    ggplot2::theme(panel.background = ggplot2::element_rect(fill = "#EAF2F8", colour = NA),
                   legend.key.width = grid::unit(1.5, "cm"))
}
