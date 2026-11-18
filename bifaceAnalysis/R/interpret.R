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

# ---- collections ------------------------------------------------------------

COLLECTION_CAUTIONS <- c(
  "A difference between collections is a difference between these samples. It can arise from the age of the artifacts, raw material, resharpening and breakage, how each collection was gathered, or local practice. The analysis does not distinguish among these.",
  "The grey reference ranges show what reassigning these same specimens produces. They are not confidence or credible intervals, and a value outside its range is not evidence of a real difference between populations or places.",
  "Left and right in an outline depend on which face was photographed. Unless every collection was photographed with a consistent face up, tick \"Analyse symmetric shapes\" before comparing collections.",
  "Collections photographed or outlined by different people, at different image resolutions, or with different detection settings can differ for that reason alone. Edge roughness from image resolution affects the higher harmonics.",
  "The PC axes describe variation among all specimens. They are not the directions that best separate the collections.",
  "Every saved specimen, including unlabelled ones and those in collections too small to compare, is part of the alignment, the automatic choice of the number of harmonics, and the PCA. Adding, removing, or editing any specimen can therefore change these results."
)

COLLECTIONS_HELP <- "Collections are defined by group labels. Give each specimen the name of its collection, for example its site, in step 2 or step 3. Each collection needs at least 3 specimens, and at least 2 collections are needed for a comparison. A file of collection locations is optional and adds a map."

text_collection_membership <- function(cd) {
  sizes <- sprintf("%s (%d)", names(cd$sizes), as.integer(cd$sizes))
  out <- sprintf("The comparison uses %d collections with %d specimens in total: %s.", length(cd$eligible), sum(cd$sizes), name_list(sizes))
  left <- character(0)
  if (cd$unlabelled > 0) left <- c(left, sprintf("%d specimen%s without a group label", cd$unlabelled, if (cd$unlabelled == 1) "" else "s"))
  if (length(cd$small)) left <- c(left, sprintf("%d specimen%s in collections with fewer than %d specimens (%s)", sum(cd$small),
                                                if (sum(cd$small) == 1) "" else "s", MIN_COLLECTION_SIZE, name_list(names(cd$small))))
  if (length(left)) out <- c(out, sprintf("Left out of the collection results: %s.", name_list(left)))
  out
}

text_collections <- function(cd, share, pairs, within) {
  out <- c(text_collection_membership(cd),
           "The plot shows the mean outline of each collection, rebuilt from the mean of its specimens' Fourier coefficients.")
  if (is.finite(share$share)) {
    r <- share$reference
    out <- c(out, sprintf("Of the total variation in outline shape among these specimens, %s lies between the collection means and the rest within collections. If specimens were assigned to collections at random with the same sizes, the expected share would be %s (the number of collections minus 1, divided by the number of specimens minus 1). Random reassignment gave shares from %s to %s (central 95%% of values).",
                          fmt_pct(share$share, 1), fmt_pct(share$expected, 1), fmt_pct(r$low, 1), fmt_pct(r$high, 1)))
  } else out <- c(out, "The specimens in these collections have identical outlines, so no share of variation can be computed.")
  np <- nrow(pairs)
  far <- pairs[which.max(pairs$corrected_distance), ]; near <- pairs[which.min(pairs$corrected_distance), ]
  zero <- pairs[pairs$corrected_distance == 0, , drop = FALSE]
  out <- c(out,
           "The raw distance between two collection means is larger on average than the distance between the means of the populations they were drawn from, because each mean carries sampling error. The corrected distance subtracts the expected size of that error. It assumes that specimens are independent of one another and does not account for the alignment and the number of harmonics having been fitted to these same specimens.",
           if (np > 1) sprintf("The largest corrected distance is between %s and %s, and the smallest between %s and %s.", far$collection_a, far$collection_b, near$collection_a, near$collection_b),
           if (nrow(zero)) sprintf("The corrected distance is 0 for %s: their raw distance is no larger than sampling error of the expected size produces. This does not show that the collections have the same mean shape.",
                                   name_list(paste(zero$collection_a, "and", zero$collection_b))),
           sprintf("There %s %d pair%s of collections. Pairs that share a collection are not independent of one another, and with several pairs an observed distance can fall outside its grey range through sampling alone.",
                   if (np == 1) "is" else "are", np, if (np == 1) "" else "s"))
  v <- stats::setNames(within$variance, within$collection)
  out <- c(out, if (all(v > 0)) sprintf("Shape variance within collections is largest in %s and smallest in %s, a ratio of %s to 1. Variance from 3 to 5 specimens is imprecise. The bars show how much it changes when one specimen is left out; they are not intervals.",
                                        names(v)[which.max(v)], names(v)[which.min(v)], fmt_num(max(v) / min(v), 1))
         else if (any(v > 0)) sprintf("The specimens of %s do not differ in outline, so no ratio of variances is given.", name_list(names(v)[v == 0]))
         else "No collection shows variation in outline shape among its specimens.")
  c(out, "Colour the PCA plot by group label to see the specimens of each collection in the space of the principal components.")
}

text_collection_geo <- function(geo) {
  if (!isTRUE(geo$available)) return(geo$reason)
  r <- geo$reference
  G <- (1 + sqrt(1 + 8 * geo$n_pairs)) / 2
  c(sprintf("Each point is a pair of collections: %d located collections give %d pairs. Pairs that share a collection are not independent, so the points are not independent observations.", as.integer(G), as.integer(geo$n_pairs)),
    sprintf("The correlation between corrected shape distance and geographic distance is %s. Reassigning the positions among the collections (%d assignments%s) gave values from %s to %s (central 95%%), with a mean of %s.",
            fmt_num(geo$r), as.integer(geo$n_assignments), if (geo$n_assignments <= 5040) ", all possible" else ", chosen at random",
            fmt_num(r$low), fmt_num(r$high), fmt_num(r$mean)),
    "With a few locations, any pattern here can arise from the particular collections sampled. A few locations cannot establish a spatial pattern. Rounded or approximate coordinates change the geographic distances.")
}

text_collection_map <- function(window, not_located, unmatched) {
  out <- if (isTRUE(window$ok)) c(
    "Each point is a located collection, coloured by the selected summary. The colours show means of samples, and a mean from a small collection can differ from the mean of the population it came from by more than the differences between collections.",
    "The coastlines and boundaries are low resolution and are for orientation only. A small window may show no boundaries, and outside the contiguous United States the map shows country borders only. The scale bar is correct along its own line of latitude and approximate elsewhere. Labels of nearby collections can overlap.")
  else window$reason
  if (length(not_located)) out <- c(out, sprintf("Collections without a location: %s.", name_list(not_located)))
  if (length(unmatched)) out <- c(out, sprintf("Rows of the locations file that match no collection in the comparison: %s. Names must match the group labels exactly.", name_list(unmatched)))
  out
}
