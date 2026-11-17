# Design specification: rewrite of the biface morphometrics Shiny app

Date: 2026-11-17
Branch: `app-rewrite`
Replaces: `bifaceAnalysis/app.R` at commit f8eeed2 (the version deployed at
https://lipolab.shinyapps.io/bifaceAnalysis/)
Status: design, not yet built. To be audited before implementation.

## 1. Purpose and audience

The app teaches outline-based shape analysis of bifaces and projectile points to
student analysts. A student photographs or uploads artifact images, extracts
outlines, and examines shape variation in the assemblage. The rewrite has three
goals set by the author: make the app more robust, make it simpler to use, and
give students more information about what the results mean.

## 2. Defects in the current app that the design must remove

Each item was confirmed by reading `bifaceAnalysis/app.R` at f8eeed2. Line
numbers refer to that file.

Results that are wrong:

- D1. Harmonic power (lines 1562 to 1575) groups coefficient columns as
  `(4*(h-1)+1):(4*h)`. Momocs orders the columns A1..An, B1..Bn, C1..Cn, D1..Dn,
  so each bar sums coefficients from unrelated harmonics.
- D2. Region variance, variance heatmap, region correlations, and PC1
  contributions (lines 1668 to 1947) compare outline point i across specimens.
  Points are produced by `coo_interpolate(n = 100)` from a contour whose start
  point is arbitrary (line 788), so point i is not the same location on
  different specimens. `custom_mshape` averages non-corresponding points.
- D3. The "modularity index" (lines 1903 to 1925) is described as the covariance
  ratio of Adams (2016) [CITE-CHECK: Adams 2016]. It is a sum of correlations
  divided by a sum of mean distances and is not that statistic.
- D4. Geometric region detection applies curvature from the last processed
  specimen to the mean shape (lines 1239 to 1251). The left blade index range is
  computed with the same expression as the right blade (lines 279 to 292).
  `find_inflection_points` (lines 90 to 112) recurses without a stopping rule.
- D5. Manual mode collects blade-junction points that are never used, and applies
  point indices chosen on one specimen to the mean shape (lines 1125 to 1233).
- D6. The "Thin-Plate Splines" tab computes no thin-plate spline (lines 1632 to
  1665).

Robustness defects:

- D7. Each press of "Process Image" appends an outline (lines 803 to 808).
  Reprocessing an image duplicates the specimen. No outline can be removed.
- D8. Defaults ("darker than background", tip "up") do not match the demo
  images, which show light artifacts on a dark background with the tip to the
  left.
- D9. The largest contour is chosen by point count, not area (lines 766 to 777).
  Ink catalog numbers on the demo artifacts produce holes in the mask.
- D10. Click markers use image coordinates as CSS pixels (lines 590 to 606).
- D11. Analysis fails with one specimen, with fewer than three specimens, and
  when the cluster count exceeds the specimen count.
- D12. `efourier(norm = TRUE)` after `coo_align()` (lines 1077 to 1085) does not
  fix which end is the tip, so near-symmetric outlines can be rotated 180
  degrees relative to each other.
- D13. Momocs was archived from CRAN on 2026-03-18. A redeploy cannot install it
  from CRAN.
- D14. The PDF handler duplicates about 300 lines of plotting code. The README
  describes exports and a script that the app does not contain.

## 3. Decisions already made by the author (2026-11-17)

1. Elliptical Fourier analysis is implemented in base R inside the app and
   validated numerically against Momocs 1.4.1. Momocs is not a runtime
   dependency.
2. Region analysis is rebuilt on corresponding points. The invented modularity
   index is dropped. The covariance ratio is not computed.
3. Work is done on branch `app-rewrite` of `clipo/bifaceAnalysis`. Nothing is
   pushed or deployed without the author's instruction.
4. Additions: interpretation text beside each result, scaled linear
   measurements, group labels, and data exports.

## 4. Scope

In scope: the Shiny app in `bifaceAnalysis/`, its tests, and its README.

Out of scope: the command-line script `biface_analysis.R` and the PDF and CSV
outputs at the repository root. That script shares defects D1, D2, and D12. The
README will state this. The script is not modified in this work.

Not included, with reasons: per-specimen shoulder landmarks (adds two clicks per
specimen and does not apply to lanceolate forms), the covariance ratio (author
decision 2), and thin-plate spline deformation grids (PC shape reconstructions
in section 7.4 serve the teaching purpose more directly).

## 5. Dependencies and layout

Runtime packages: shiny, bslib, imager, ggplot2. imager 1.0.8 is current on CRAN
(published 2025-12-23). dplyr and Momocs are removed.

```
bifaceAnalysis/            deployed app directory
  app.R                    UI and server wiring only
  R/image.R                thresholding, mask cleaning, contour extraction
  R/outline.R              orientation, landmarks, resampling, measurements
  R/efa.R                  elliptical Fourier forward, inverse, power
  R/shape.R                Procrustes alignment, PCA, clustering, regions
  R/plots.R                one plotting function per figure (screen and PDF)
  R/interpret.R            interpretation text built from computed values
  www/                     demo images
tests/testthat/            unit tests (outside the deployed directory)
tests/fixtures/            Momocs 1.4.1 reference values and generating script
docs/specs/, docs/audits/
```

Tests and the Momocs fixture script sit outside `bifaceAnalysis/` so that
rsconnect's dependency scan does not record testthat or Momocs for deployment.

## 6. Processing pipeline

All functions in this section are pure functions of their inputs, with no Shiny
dependency, so each can be unit tested.

### 6.1 Image to binary mask

1. Load JPEG or PNG. If the longer side exceeds 1200 pixels, downscale to 1200.
   All later coordinates, including scale calibration, are in the downscaled
   image.
2. Optional crop from a drag rectangle (Shiny brush). Default is the whole image.
3. Convert to grayscale. Compute a global threshold by Otsu's method on a
   256-bin histogram (implemented in base R). The student can override the
   threshold with a slider.
4. Polarity. The artifact is assumed not to touch the crop border. The class
   (above or below threshold) holding the majority of border pixels is
   background. The student can override.
5. Morphological opening with a student-set radius (default 2 pixels) removes
   small specks.
6. Fill holes: every background component that does not touch the border
   becomes foreground. This removes ink labels and specular highlights.
7. Keep the connected component with the largest pixel area.
8. Extract the single closed boundary of that component.

Quality checks, shown as plain warnings and never as silent failures: the mask
touches the crop border; the mask covers less than 1 percent or more than 90
percent of the crop; the boundary has fewer than 100 points; no foreground
found.

### 6.2 Orientation and landmarks

Coordinates are converted to a y-up system. Traversal direction is set to
counterclockwise using the sign of the shoelace area.

1. The student declares the tip direction for the image (up, down, left, right).
   The setting carries over to the next image.
2. Compute the two principal axes of the outline polygon from its area second
   moments. Choose the axis whose direction is within 45 degrees of the declared
   tip direction. This is the long axis for most points and the short axis for
   points that are wider than long. Orient the unit vector u toward the declared
   tip direction.
3. If the ratio of the principal second moments is below 1.1, the axes are
   poorly determined. Use the declared direction itself as u and show a warning.
4. Tip landmark T: the outline point with the largest projection on u.
5. Midline: the line through the area centroid parallel to u. Base landmark B:
   the intersection of the midline with the outline that is most distant from T
   in the direction of -u.
6. Manual override: the student can click T and B on the outline preview. u then
   runs from B to T.
7. Rotate so that u points up.

### 6.3 Corresponding points

T and B divide the outline into two edges. Each edge is resampled by linear
interpolation to 100 points equally spaced in arc length, starting at T for the
first edge and at B for the second, giving 200 points with T at index 1 and B at
index 101. Point k therefore lies at the same relative arc-length position
between tip and base on every specimen.

Stated limitation (shown in the app): this correspondence is geometric. Equal
arc-length fractions are not guaranteed to fall on the same morphological
feature. A point with long barbs has more edge length in the shoulder zone than
a lanceolate point. Which face was photographed is arbitrary, so left and right
are not comparable across specimens, and all region results pool both edges.

### 6.4 Measurements

Computed from the oriented full-resolution outline, in pixels and in millimetres
when a scale is set: length (extent along u), maximum width (extent
perpendicular to u), length/width, position of maximum width as a fraction of
length from the base, area (shoelace), perimeter, and an asymmetry score. The
asymmetry score is the Procrustes distance between the 200-point configuration
and its mirror image about the tip-base axis with the two edges relabelled.

Scale calibration: the student clicks two points on the scale bar and enters the
distance in millimetres. Calibration is stored per image. A checkbox applies the
current calibration to all images, for batches photographed at a fixed distance.
Without calibration, only pixel values and ratios are reported and the app says
so.

### 6.5 Specimen store

Specimens are keyed by image name. Saving an outline for an image that already
has one replaces it. Each specimen holds: name (editable), group label
(editable, optional), the 200 corresponding points in pixel units, the
full-resolution oriented outline, measurements, calibration, and the processing
settings used. Specimens can be removed. The store can be downloaded as a CSV of
outline coordinates and reloaded later, because shinyapps.io sessions time out.
CSV is used instead of RDS because `readRDS` on an uploaded file is unsafe.

## 7. Analysis

Minimum sample sizes are enforced with a plain message, not an error: 2
specimens for the outline overlays, 3 for PCA, 4 for clustering and region
covariation.

### 7.1 Alignment

Generalized Procrustes alignment of the 200-point configurations: translate to
the centroid, scale to unit centroid size, rotate (no reflection) to the
iteratively updated mean by least squares until the change in summed squared
distance is below 1e-10. The mean shape is rotated so its tip-base axis is
vertical. Implemented in base R. Points are treated as fixed, not sliding,
semilandmarks.

An optional setting analyses symmetrized shapes (the average of each
configuration and its relabelled mirror image). Default is off. The app explains
that with the photographed face arbitrary, asymmetric variation has an arbitrary
sign.

### 7.2 Elliptical Fourier analysis

Coefficients a_n, b_n, c_n, d_n for n = 1..H are computed from the aligned closed
200-point polygons with the Kuhl and Giardina (1982) formulas
[CITE-CHECK: Kuhl and Giardina 1982; the PDF is not in the repository]. With
segment increments dx_p, dy_p, dt_p = sqrt(dx_p^2 + dy_p^2), cumulative length
t_p, and perimeter T:

    a_n = T / (2 n^2 pi^2) * sum_p (dx_p/dt_p) [cos(2 n pi t_p / T) - cos(2 n pi t_{p-1} / T)]
    b_n = T / (2 n^2 pi^2) * sum_p (dx_p/dt_p) [sin(2 n pi t_p / T) - sin(2 n pi t_{p-1} / T)]

c_n and d_n use dy_p in place of dx_p.

No first-ellipse normalization is applied. Size, position, rotation, and start
point are already standardized by sections 6.2, 6.3, and 7.1. This removes
defect D12.

Harmonic power for harmonic n is (a_n^2 + b_n^2 + c_n^2 + d_n^2) / 2, averaged
over specimens. Cumulative power is reported both including and excluding the
first harmonic, labelled as such. The default H is the smallest number of
harmonics reaching 99 percent of cumulative power excluding the first harmonic,
bounded to 6..30. The student can override H. A reconstruction panel shows one
selected specimen rebuilt with 1, 2, 4, 8, 16, and H harmonics over its original
outline.

Validation: (i) forward coefficients match `Momocs::efourier(norm = FALSE)` on
the same coordinates to 1e-8, from a committed fixture; (ii) inverse of forward
reproduces the input polygon within a stated tolerance at H = 100; (iii) an
ellipse with semi-axes 3 and 1 sampled at 2000 points gives first-harmonic
amplitudes 3 and 1 and higher harmonics below 1e-4; (iv) power computed by the
app equals the same formula applied to the Momocs coefficient matrix with its
documented column order.

### 7.3 PCA and stability

`prcomp` on the specimen-by-4H coefficient matrix, centered and unscaled.
Reported: variance explained per component and a scree plot. Stability is
described without inference: each specimen is left out in turn, the PCA is
refitted, and the absolute cosine between the refitted and full-data loading
vector is recorded for PC1 and PC2. The minimum and median are shown, with text
explaining that values near 1 mean the axis does not depend on any single
specimen.

### 7.4 Shape change along components

For PC1 to PC3, the mean coefficient vector plus and minus 2 standard deviations
of the scores along the component is passed through the inverse transform and
drawn over the mean shape. The morphospace plot shows PC scores with
reconstructed shapes on a background grid, coloured by group label or cluster.
This replaces the tab named "Thin-Plate Splines".

### 7.5 Clustering

Ward (D2) hierarchical clustering on Euclidean distances between coefficient
vectors. The dendrogram is shown. The student sets k between 2 and n - 1. Mean
silhouette width is shown for each k in that range as a descriptive guide. When
group labels exist, a cluster-by-group table is shown. No test statistic or
p-value is reported.

### 7.6 Regions

Regions are intervals of relative arc-length position from tip (0) to base (1),
applied to both edges: tip, blade, shoulder, base, with default boundaries 0.15,
0.60, 0.80, adjustable by the student. Membership is by point index and is
identical on every specimen. Boundaries are drawn on the mean shape and on each
specimen so the student can see where they fall.

Outputs:

- Per-point variance: the sum of the x and y variances of the Procrustes
  coordinates at each point, drawn as colour on the mean shape.
- Per-region variance: the mean per-point variance within the region.
- Between-region covariation: the RV coefficient (Escoufier 1973)
  [CITE-CHECK: Escoufier 1973] between the coordinate blocks of each pair of
  regions. For each pair, the specimens of one block are permuted 999 times at a
  fixed seed. The mean and the 2.5 to 97.5 percent range of the permuted RV are
  drawn next to the observed value as a reference for what RV looks like at this
  sample size with no association. No p-value is reported.

Reason a frequentist reference is used here: a Bayesian model of covariance
between coordinate blocks of 30 to 180 variables with fewer than 30 specimens
is not tractable in an interactive classroom app, and the quantity is
descriptive.

Stated limitations (shown in the app): RV is biased upward when specimens are
few relative to variables, which the permutation reference makes visible. A
single Procrustes alignment of the whole outline spreads variation across
regions and can induce covariation between them
[CITE-CHECK: source for superimposition-induced covariance]. Region boundaries
are proportional, not anatomical.

## 8. Interface

Four steps in a navigation bar, plus a Help page.

1. Images. Choose demo images or upload JPEG or PNG files. A "Run demo" button
   processes the four demo images with automatic settings and opens Results.
2. Outline. One image at a time, with previous and next controls. The image,
   mask, and oriented outline with T, B, and region boundaries update when any
   setting changes; there is no separate "process" button. Controls: drag to
   crop, threshold (automatic by default), polarity (automatic by default),
   speck removal, tip direction, optional scale calibration, specimen name,
   group label. Warnings from section 6.1 appear here. "Save outline" adds or
   replaces the specimen.
3. Specimens. A table of saved specimens with thumbnails of the outlines,
   measurements, and group labels. Edit name or group, remove a specimen,
   download or reload the outline CSV.
4. Results. Sub-tabs: Outlines, Measurements, Harmonics, PCA, Shape change,
   Clusters, Regions. Results recompute when specimens or settings change; there
   is no "Run analysis" button. Each sub-tab has a "How to read this" panel
   whose text is generated from the current values (section 9).

Downloads: PDF report (figures, interpretation text, settings, and package
versions), and CSV files for measurements, PC scores, Fourier coefficients, and
outline coordinates.

Figures use ggplot2, the Okabe-Ito palette for categories and viridis for
continuous values, a sans-serif font, and no in-figure titles. The same plotting
function produces the screen figure and the PDF page.

## 9. Interpretation text

`R/interpret.R` returns text built from computed values. Each statement reports
a value and what it measures. Text does not make archaeological claims (function,
chronology, cultural affiliation) and does not use significance language.
Examples of content: the share of variance on PC1 and PC2 and the specimens at
each extreme; what shape difference the PC1 reconstructions show, described by
measured differences in length/width and position of maximum width between the
plus and minus shapes; which region has the highest per-point variance and by
what factor; the observed RV next to its permutation range; the leave-one-out
stability values; and a statement of sample size with the number of components
that can exist (n - 1).

## 10. Testing and acceptance

Unit tests (testthat), run from the repository root:

- Otsu threshold on a synthetic bimodal image; polarity detection; hole filling;
  largest-component selection with a second smaller object present.
- Segmentation of the four demo images yields one closed outline each with no
  warnings under automatic settings and tip direction "left".
- Orientation: a synthetic triangle-plus-stem shape rotated by arbitrary angles
  with the matching declared direction returns the same T, B, and resampled
  points (within tolerance) after orientation. A wider-than-long shape is
  oriented by the declared direction. A near-circular shape triggers the
  warning.
- Resampling: 200 points, T at 1, B at 101, equal spacing within each edge.
- Measurements on a rectangle and an ellipse of known dimensions; calibration
  conversion.
- EFA validation items (i) to (iv) in section 7.2.
- Procrustes alignment matches `Momocs::fgProcrustes` or `geomorph::gpagen`
  (fixed semilandmarks) from a committed fixture, up to a global rotation.
- PCA scores match `prcomp` on Momocs coefficients for the same input, up to
  sign.
- RV equals 1 for identical blocks and matches a direct trace formula
  computation on random data.
- Edge cases: 1, 2, and 3 specimens produce messages, not errors; k above n - 1
  is clamped; re-saving a specimen replaces it.

App-level test (shinytest2): load the demo, run it, confirm each Results sub-tab
renders and the PDF and CSV downloads are non-empty.

Acceptance for the build milestone: all tests pass; the app runs locally on the
demo images from a clean session; `rsconnect::appDependencies("bifaceAnalysis")`
lists no package that is absent from CRAN.

## 11. Open points for the auditor and author

- Whether 100 points per edge and fixed (non-sliding) semilandmarks are adequate
  for the region analysis, given that sliding would need a much larger
  implementation.
- Whether Procrustes-aligned coordinates or axis-oriented, centroid-size-scaled
  coordinates are the better input to the Fourier analysis.
- Whether the base landmark rule in 6.2 step 5 behaves acceptably on
  basally notched and concave-base points.
