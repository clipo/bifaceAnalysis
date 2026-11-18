# Interpretation text. Each function returns a character vector of paragraphs
# built from the values it is given. Statements report what was computed and
# what the quantity measures. They make no archaeological claims and use no
# significance language.

fmt_pct <- function(x, digits = 0) paste0(formatC(100 * x, format = "f", digits = digits), "%")
fmt_num <- function(x, digits = 2) formatC(x, format = "f", digits = digits)
name_list <- function(x) {
  if (length(x) <= 1L) return(paste(x))
  paste(paste(x[-length(x)], collapse = ", "), "and", x[length(x)])
}

text_sample <- function(n) {
  out <- sprintf("This analysis uses %d specimen%s.", n, if (n == 1) "" else "s")
  if (n >= 3)
    out <- c(out, sprintf("With %d specimens, at most %d principal components can exist. Results from a small number of specimens describe those specimens and can change when specimens are added or removed.", n, n - 1L))
  out
}

text_outlines <- function(n, symmetrized) {
  c(text_sample(n)[1],
    "Each outline was turned so its tip points up, then divided at the tip and at the base into a left and a right edge. Each edge is represented by 100 points spaced evenly along its length, so the point 30 steps from the tip, for example, is 30 percent of the way along the edge on every specimen.",
    "The overlay shows all specimens after Procrustes alignment, which removes differences in position, size, and rotation. What remains is difference in shape. The black line is the mean shape. Where the thin lines spread apart, the specimens differ.",
    if (symmetrized) "The shapes shown were made symmetric by averaging each one with its mirror image. Measurements and saved outlines are not changed by this option." else
      "Left and right in these plots depend on which face of the artifact was photographed. Differences between the two edges are therefore not comparable from specimen to specimen unless every artifact was photographed with a consistent face up.")
}

text_measurements <- function(meas) {
  has_mm <- all(is.finite(meas$length_mm))
  out <- "Length is the extent of the outline along the tip-to-base axis. Maximum width is its extent at right angles to that axis. Both correspond to caliper measurements. The widest cross-section is located by measuring straight across the outline at 200 heights, and its height is given as a fraction of length above the base end."
  ratio <- sprintf("Length divided by width ranges from %s to %s.",
                   fmt_num(min(meas$length_width_ratio)), fmt_num(max(meas$length_width_ratio)))
  if (has_mm) {
    r <- range(meas$length_mm)
    out <- c(out, sprintf("Lengths range from %s to %s mm. %s", fmt_num(r[1], 1), fmt_num(r[2], 1), ratio))
  } else {
    out <- c(out, "At least one specimen has no scale, so lengths, widths, and areas are shown in image pixels. Pixel values depend on how each photograph was taken and cannot be compared between images. Ratios such as length divided by width do not depend on scale and can be compared.",
             ratio)
  }
  c(out, "The asymmetry score is the Procrustes distance between a shape and its own mirror image. It is 0 for a perfectly symmetric outline and increases as the two edges differ.",
    if (nrow(meas) > 1)
      sprintf("The most asymmetric specimen is %s (%s) and the least asymmetric is %s (%s).",
              meas$specimen[which.max(meas$asymmetry)], fmt_num(max(meas$asymmetry), 3),
              meas$specimen[which.min(meas$asymmetry)], fmt_num(min(meas$asymmetry), 3)))
}

text_harmonics <- function(power_table, n_used, automatic) {
  first <- power_table$cumulative_all[1]
  reached <- power_table$cumulative_excluding_first[n_used]
  c("Elliptical Fourier analysis describes a closed outline as a sum of harmonics. The first harmonic is an ellipse. Each later harmonic adds finer detail. The lower panel shows one specimen rebuilt from increasing numbers of harmonics.",
    sprintf("In these specimens the first harmonic holds %s of total power. The chart shows how the remaining power is divided among harmonics 2 and above.", fmt_pct(first, 1)),
    sprintf("The analysis uses %d harmonics (%s). Cumulative power excluding the first harmonic: these harmonics hold %s of the power beyond the first harmonic. Cumulative power including the first harmonic: they hold %s of total power. Harmonics that are left out carry the finest detail of the outline, which includes small real features and the pixel-level roughness of the traced edge.",
            n_used, if (automatic) "chosen automatically as the smallest number that reaches 99% of the power beyond the first harmonic, within the limits 6 to 30" else "set by hand",
            fmt_pct(reached, 1), fmt_pct(power_table$cumulative_all[n_used], 2)))
}

# Length/width ratio and widest position of a reconstructed (tip-up) shape.
shape_descriptors <- function(shape) {
  m <- measure_outline(shape)
  c(ratio = m$length_width_ratio, widest = m$widest_position)
}

text_pca <- function(pca, n) {
  ex <- pca$explained
  sc <- pca$scores
  ends <- function(i) sprintf("%s has the lowest score and %s the highest",
                              rownames(sc)[which.min(sc[, i])], rownames(sc)[which.max(sc[, i])])
  c("Principal component analysis finds the directions along which the specimens' shapes differ most. PC1 is the direction with the most variation, PC2 the direction with the most of what remains, and so on. Each point in the plot is one specimen. Specimens that plot close together have similar outlines.",
    sprintf("PC1 holds %s of the shape variance and PC2 holds %s, %s together. On PC1, %s. On PC2, %s.",
            fmt_pct(ex[1]), fmt_pct(ex[2]), fmt_pct(ex[1] + ex[2]), ends(1), ends(2)),
    "The grey shapes in the background show the outline that corresponds to each position in the plot. The distance between two specimens in this plot understates their full difference whenever later components hold a substantial share of the variance.",
    text_sample(n)[-1])
}

text_pc_shapes <- function(pca) {
  pcs <- seq_len(min(3L, pca$rank))
  desc <- vapply(pcs, function(pc) {
    e <- pc_extreme_shapes(pca, pc)
    lo <- shape_descriptors(e$low); hi <- shape_descriptors(e$high)
    sprintf("From the low end to the high end of PC%d, length divided by width goes from %s to %s, and the widest cross-section moves from %s to %s of the length above the base.",
            pc, fmt_num(lo["ratio"]), fmt_num(hi["ratio"]), fmt_pct(lo["widest"]), fmt_pct(hi["widest"]))
  }, character(1))
  c("Each panel shows the mean shape (dashed) and the shapes that correspond to scores 2 standard deviations below and above the mean on one component. Compare the two solid outlines to see which parts of the outline move along that component.",
    desc,
    "These two numbers summarize only part of each shape change. Look at the outlines for changes in the shoulders, barbs, and base. The plotted shapes are positions along a component, not real specimens, and with few specimens they can look unlike any artifact in the sample.")
}

text_stability <- function(stab) {
  worst <- which.min(stab$plane)
  out <- c("To show how much the components depend on individual specimens, the analysis was repeated with each specimen left out in turn, keeping the same alignment and number of harmonics. The plot compares the directions of the components in each repeat with those of the full analysis. A value near 1 means the direction, or the plane, was nearly unchanged when the specimen named on that row was left out. A value well below 1 means it changed substantially. These values do not measure change in the scores or in the shares of variance.",
           sprintf("For the plane containing PC1 and PC2, which is what the scatter plot shows, the lowest agreement is %s (when %s is left out) and the median is %s.",
                   fmt_num(min(stab$plane, na.rm = TRUE)), names(stab$plane)[worst], fmt_num(stats::median(stab$plane, na.rm = TRUE))))
  for (a in names(stab$separable)) {
    if (stab$separable[[a]]) {
      v <- stab$axis[, a]
      out <- c(out, sprintf("For the direction of %s alone, the lowest agreement is %s (when %s is left out).",
                            a, fmt_num(min(v, na.rm = TRUE)), rownames(stab$axis)[which.min(v)]))
    } else {
      out <- c(out, sprintf("%s holds a share of variance too close to a neighbouring component's for its direction to be distinguished reliably, so no separate value is shown for it. Components with similar shares can exchange order or rotate into each other when a specimen is removed. Read such components together, as a plane, and not one at a time.", a))
    }
  }
  out
}

text_clusters <- function(tree, k, sil, groups = NULL) {
  sizes <- table(stats::cutree(tree, k))
  out <- c("Hierarchical clustering (Ward's method) joins the most similar outlines first and then joins groups step by step. In the tree, the height at which two branches join is the dissimilarity between them.",
           if (is.finite(tree_cut_height(tree, k))) "Cutting the tree at the dashed line gives the number of clusters you chose." else
             "Several branches join at the same height here, so no single horizontal cut gives the number of clusters you chose and no dashed line is drawn. The division among those tied branches is arbitrary.",
           sprintf("With %d clusters, the cluster sizes are %s.", k, name_list(as.integer(sizes))))
  if (nrow(sil) > 0) {
    best <- sil$k[which.max(sil$silhouette)]
    out <- c(out, sprintf("Mean silhouette width measures how much closer specimens are to their own cluster than to the nearest other cluster, on a scale from -1 to 1. Values near 1 indicate compact, well separated clusters. Values near 0 or below indicate clusters that overlap. Here it is highest at %d clusters (%s) and is %s at your choice of %d.",
                          best, fmt_num(max(sil$silhouette)), fmt_num(sil$silhouette[sil$k == k]), k))
  }
  out <- c(out, "Clustering always produces clusters, including from data with no real groups. Treat clusters as a description of similarity in outline shape, to be compared with other information about the artifacts.")
  if (has_labels(groups))
    out <- c(out, "The table counts the specimens in each cluster by the group label you entered. It shows whether specimens that share a label were placed in the same cluster at this number of clusters.")
  out
}

text_regions <- function(region_variance, rv_table, n, boundaries) {
  rvv <- region_variance[is.finite(region_variance)]
  b <- sort(boundaries)
  out <- c(sprintf("Each edge is divided by position along its length from the tip: tip from 0 to %s, blade to %s, shoulder to %s, and base from there to the base landmark. The same point numbers are used on every specimen. These divisions are proportional. They do not locate the actual shoulder or haft of each artifact, so check the region plot against your specimens and adjust the boundaries if needed.",
                   fmt_pct(b[1]), fmt_pct(b[2]), fmt_pct(b[3])))
  if (length(rvv) >= 2 && min(rvv) > 0)
    out <- c(out, sprintf("Variance per point is highest in the %s region and lowest in the %s region, a ratio of %s to 1. This measures how much specimens differ from one another in that part of the outline after alignment.",
                          tolower(names(rvv)[which.max(rvv)]), tolower(names(rvv)[which.min(rvv)]), fmt_num(max(rvv) / min(rvv), 1)))
  out <- c(out, "Procrustes alignment fits whole outlines to each other, so a large difference in one part of the outline is partly spread to the other parts. The variance at a point reflects that overall fit and not only the local difference at that point.")
  if (!is.null(rv_table)) rv_table <- rv_table[is.finite(rv_table$rv), , drop = FALSE]
  if (!is.null(rv_table) && nrow(rv_table) > 0) {
    top <- rv_table[which.max(rv_table$rv), ]
    above <- rv_table$rv > rv_table$reference_high
    out <- c(out,
             "The RV coefficient measures how much the coordinates of two regions vary together across specimens in a linear way. It ranges from 0 (no linear association between the two sets of coordinates) to 1 (the specimens differ from one another in the same pattern in both regions). A value of 0 does not rule out other kinds of relation between the regions.",
             sprintf("The highest observed value is %s, between the %s and the %s.", fmt_num(top$rv), tolower(top$region_a), tolower(top$region_b)),
             sprintf("RV tends to be large when specimens are few, including when the pairing of specimens between two regions is arbitrary. The grey bars show the values obtained when each specimen's coordinates for one region were paired at random with another specimen's coordinates for the other region. With %d specimens, %s. An observed value inside the grey bar is no larger than values that arise from arbitrary pairing.",
                     n, if (any(above)) sprintf("the observed value lies above the grey bar for %d of %d pairs of regions", sum(above), nrow(rv_table)) else "no observed value lies above its grey bar"),
             if (n <= 5) sprintf("With %d specimens there are only %d distinct ways to pair them, so the grey bars are coarse.", n, as.integer(factorial(n))),
             "The grey bar depends on the number of specimens and on the pattern of variation within each of the two regions, all of which are held fixed. It is not a margin of error for the observed value. The random pairing reuses the same aligned coordinates, so it does not remove shared variation that the alignment of whole outlines itself produces. Treat these values as descriptions, not as evidence of how the artifacts were made or used.")
  }
  out
}
