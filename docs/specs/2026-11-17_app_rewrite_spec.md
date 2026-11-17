# Design specification: rewrite of the biface morphometrics Shiny app

Date: 2026-11-17
Branch: `app-rewrite`
Replaces: `bifaceAnalysis/app.R` at commit f8eeed2 (the version deployed at
https://lipolab.shinyapps.io/bifaceAnalysis/)
Status: design, revision 2 (2026-11-17), not yet built. Revision 1 was audited
by Codex (stem 2026-11-17_spec). Revision 2 resolves findings F1 to F12; see
`docs/audits/codex_audit_2026-11-17_spec_ledger.md` and section 12.

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
- D11. The analysis has no guards on specimen count or cluster count. With one
  specimen `hclust` stops. With two specimens it stops at the default cluster
  count of 3 and runs at 2. With one specimen the region plot's `1:(n-1)` loop
  indexes out of range. Components beyond the data's rank are reported as if
  they carried information.
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
(published 2025-12-23). dplyr and Momocs are removed. Silhouette width,
Procrustes alignment, Fourier analysis, and the RV coefficient are implemented
in base R as specified below.

Test-only packages: testthat; cluster (a recommended package shipped with R,
used to check the silhouette implementation). Fixture generation only: Momocs
1.4.1 (CRAN archive) and geomorph, used by `tests/fixtures/make_reference.R` to
write committed reference values. Tests read the committed values and do not
need Momocs. Deployment: rsconnect.

```
bifaceAnalysis/            deployed app directory
  app.R                    UI and server wiring only
  R/image.R                thresholding, mask selection, contour extraction
  R/outline.R              orientation, landmarks, resampling, measurements
  R/efa.R                  elliptical Fourier forward, inverse, power
  R/shape.R                Procrustes alignment, PCA, clustering, regions
  R/store.R                specimen records, CSV schema, image loading
  R/plots.R                one plotting function per figure (screen and PDF)
  R/interpret.R            interpretation text built from computed values
  R/report.R               PDF report
  www/                     demo images
tests/testthat/            unit and server tests (outside the deployed directory)
tests/fixtures/            reference values, reference masks, generating script
docs/specs/, docs/audits/
```

Tests and the fixture script sit outside `bifaceAnalysis/` so that rsconnect's
dependency scan does not record testthat, cluster, or Momocs for deployment.

## 6. Processing pipeline

All functions in this section are pure functions of their inputs, with no Shiny
dependency, so each can be unit tested.

### 6.1 Image to binary mask

Coordinates. "Original pixels" are pixels of the file. "Working pixels" are
pixels of the image after the resize in step 1. The resize factor f (working
divided by original, at most 1) is kept for each image. Every outline and
measurement is in working pixels of its own image.

1. Load JPEG or PNG. Drop any alpha channel. If the longer side exceeds 1200
   pixels, resize by linear interpolation so the longer side is 1200.
2. Optional crop from a drag rectangle. Default is the whole image.
3. Convert to grayscale. Compute a global threshold by Otsu's method on a
   256-bin histogram over [0, 1] (base R). The student can override it.
4. Form two candidate classes: pixels above the threshold ("light") and at or
   below it ("dark"). The student can restrict the search to one class.
5. Within each candidate class, apply a morphological opening with a square
   structuring element of side s pixels (student-set, default 3, 0 disables),
   then label 4-connected components.
6. Select the artifact: the component with the largest pixel area among those
   that do not touch the border of the crop, over all candidate classes. This
   replaces the border-majority polarity rule of revision 1, which failed on
   three of the four demo images because the scale bar and light margins occupy
   the border (audit finding F1). If every component touches the border, select
   the largest component overall and warn.
7. Fill holes: every 4-connected region outside the selected component that
   does not reach the border becomes part of it.
8. Pad the mask with one background pixel on every side and extract the
   iso-contour at level 0.5. Keep the longest contour. It is closed because the
   padded mask cannot touch the border.

Warnings shown to the student, never silent: no component found; every
component touched the border (step 6 fallback); the mask covers less than 1
percent or more than 90 percent of the crop; another candidate component has at
least one quarter of the selected area; the contour has fewer than 100 points
(the outline is then rejected and cannot be saved).

### 6.2 Orientation and landmarks

The contour is converted to a y-up system, consecutive duplicate vertices are
removed, and traversal is set to counterclockwise by the sign of the shoelace
area.

1. The student declares the tip direction d for the image (up, down, left,
   right). The setting carries over to the next image.
2. Compute the principal axes of the polygon from its area second moments. If
   the ratio of the larger to the smaller principal moment is below 1.1, set
   u = d and warn that the axis is poorly determined. Otherwise choose the
   principal axis with the larger absolute cosine to d. On an exact tie (both
   axes at 45 degrees) choose the major axis. Orient u so that u . d > 0.
3. Rotate so that u points up.
4. Tip T: the vertex with the largest y. If several vertices lie within 1e-9 of
   the maximum, choose the one nearest the midline of step 5; if still tied,
   the one with the smaller x.
5. Midline: the vertical line through the area centroid. An edge from vertex i
   to i+1 crosses the midline at x0 when x_i <= x0 < x_{i+1} or
   x_{i+1} <= x0 < x_i (half-open, so a vertex on the line is counted once).
   Base B: the crossing with the smallest y. B is inserted into the polygon as
   a new vertex by linear interpolation along its edge. If there are more than
   two crossings the app warns that the base may be misplaced.
6. Blunt-tip warning: if the cross-section width (section 6.4) at 5 percent of
   length below T exceeds 35 percent of the maximum width, the app warns that
   the tip end is blunt or broken and that T is a geometric choice. The outline
   can still be saved (author decision, 2026-11-17).
7. Manual override: the student clicks the tip, then the base. Each click is
   snapped to the nearest outline vertex. u runs from the base vertex to the
   tip vertex. If both clicks snap to the same vertex the app reports an error
   and keeps the previous landmarks. Manual landmarks are cleared when the crop
   changes.

T, B, and the resampled points do not depend on which vertex the contour
tracer started from, because every rule above refers to geometry only.

### 6.3 Corresponding points

T and B divide the polygon into two directed paths in counterclockwise order:
edge 1 from T to B (the left edge when the tip is up) and edge 2 from B to T.
Each path is sampled by linear interpolation at arc-length fractions j/100 of
its own length, j = 0..99, which includes the path's start and excludes its
end. This gives 200 distinct points: T at index 1, B at index 101.

Position from the tip, used for regions: index i in 1..100 has position
(i - 1)/100; index 101 has position 1; index i in 102..200 has position
(201 - i)/100. Mirror pairing: index i in 2..100 pairs with 202 - i; indices 1
and 101 pair with themselves.

Stated limitation (shown in the app): this correspondence is geometric. Equal
arc-length fractions are not guaranteed to fall on the same morphological
feature. A point with long barbs has more edge length in the shoulder zone than
a lanceolate point. Which face was photographed is arbitrary, so left and right
are not comparable across specimens. Region results pool both edges, which
gives each region the same definition on every specimen but does not remove
variation caused by the arbitrary face. Whether 100 points per edge and fixed
(non-sliding) points are adequate has not been tested and is stated as such.

### 6.4 Measurements

Computed from the oriented working-resolution contour (all contour vertices,
not the 200 resampled points):

- Length: extent along the axis (largest minus smallest y).
- Maximum width: extent at right angles to the axis (largest minus smallest x).
  Length and maximum width correspond to caliper measurements.
- Cross-section width at height y: the distance between the two outermost
  crossings of the horizontal line at y with the contour. It spans any notch.
- Position of the widest cross-section: cross-section width is evaluated at 200
  levels, the midpoints of 200 equal slices of the length. The level with the
  largest width is reported as a fraction of length above the lowest contour
  point. Ties go to the level nearest the base.
- Length divided by maximum width.
- Area by the shoelace formula.
- Asymmetry score: the partial Procrustes distance between the 200-point
  configuration and its mirror image with points relabelled by the mirror
  pairing of section 6.3.

Perimeter is not reported because a pixel contour's length depends on image
resolution.

Scale. The student clicks two points on the scale bar and enters the distance
in millimetres, which gives pixels per millimetre in working pixels of that
image. A specimen stores the working-pixel value that applied when it was
saved. The shared scale (checkbox "use this scale for images that have none of
their own") is held in original pixels per millimetre (working value divided
by f) and converted for each image by that image's f, so it is correct for
images resized by different factors. Cropping is a translation and does not
affect scale. Without a scale, lengths and areas are reported in pixels and the
app says that pixel values cannot be compared between images.

### 6.5 Specimen store

Three kinds of state are kept separate. Preview state (crop, threshold, tip
direction, clicks) belongs to the image being viewed and changes nothing that
is saved. A saved specimen is fixed when "Save outline" is pressed and changes
only if its image is saved again, it is renamed or relabelled, or it is
removed. Analysis settings (harmonics, symmetric shapes, clusters, region
boundaries) change results only and never alter a saved specimen.

Each specimen has a stable id (the image file name, or the id read from a
loaded file), a display name (unique; a suffix is added on collision), a group
label, the oriented contour with the row of the base landmark, pixels per
millimetre, and the processing settings used. The 200 points and all
measurements are derived from the contour and are not stored separately.

Persistence. The store is downloaded as a CSV, schema version 1, one row per
contour vertex, with columns: schema_version, id, specimen, group, px_per_mm,
base_row, threshold, polarity, speck, direction, manual_landmarks, cropped,
row, x, y. Loading a file restores every field above. A loaded specimen whose
id matches an existing specimen replaces it when the contours are identical
and otherwise is added under a new id. A file with a missing column, a
different schema version, or non-numeric coordinates is rejected with a
message. CSV is used instead of RDS because `readRDS` on an uploaded file is
unsafe.

## 7. Analysis

Minimum sample sizes are enforced with a plain message, not an error: 2
specimens for the outline overlays, harmonic power, and region variance; 3 for
PCA; 4 for clustering, the stability check, and region covariation.

### 7.1 Alignment

Generalized Procrustes alignment of the 200-point configurations: translate to
the centroid, scale to unit centroid size, rotate (no reflection) by least
squares to the iteratively updated unit-size mean, starting from the first
specimen, until the summed squared distance to the mean changes by less than
1e-10. The aligned set is then rotated as a whole so that the mean shape's
base-to-tip direction is vertical. Points are fixed, not sliding.

Symmetric shapes (optional, default off, analysis only). For each specimen the
centred unit-size configuration X and its relabelled mirror image M, rotated
onto X, are averaged, and the average is rescaled to unit centroid size. The
alignment above is then run on the averages. Saved contours and all
measurements, including the asymmetry score, are unaffected. Switching the
option off restores the previous results exactly.

### 7.2 Elliptical Fourier analysis

Coefficients a_n, b_n, c_n, d_n for n = 1..H are computed from the aligned closed
200-point polygons with the Kuhl and Giardina (1982) formulas
[CITE-CHECK: Kuhl and Giardina 1982; the PDF is not in the repository]. The
polygon is closed by joining point 200 to point 1. The parameter t is
cumulative chord length with t = 0 at point 1 (the tip). With segment p running
from point p to point p+1, increments dx_p, dy_p, dt_p = sqrt(dx_p^2 + dy_p^2),
t_p the parameter at the end of segment p, and perimeter T:

    a_n = T / (2 n^2 pi^2) * sum_p (dx_p/dt_p) [cos(2 n pi t_p / T) - cos(2 n pi t_{p-1} / T)]
    b_n = T / (2 n^2 pi^2) * sum_p (dx_p/dt_p) [sin(2 n pi t_p / T) - sin(2 n pi t_{p-1} / T)]

c_n and d_n use dy_p in place of dx_p. The constant terms are the arc-length
means x0 = (1/T) sum_p dt_p (x_p + x_{p+1})/2 and likewise y0. They are not
zero for a vertex-centred polygon and are kept for reconstruction. The inverse
is

    x(t) = x0 + sum_n [a_n cos(2 n pi t / T) + b_n sin(2 n pi t / T)]
    y(t) = y0 + sum_n [c_n cos(2 n pi t / T) + d_n sin(2 n pi t / T)]

Reconstructions of a specimen use its own x0, y0 and are drawn over its
aligned outline. Shapes reconstructed from a position in PC space (section 7.4)
use x0 = y0 = 0 and are drawn centred.

No first-ellipse normalization is applied. Size, position, rotation, traversal
direction, and start point are already standardized by sections 6.2, 6.3, and
7.1. This removes defect D12.

Harmonic power for harmonic n is (a_n^2 + b_n^2 + c_n^2 + d_n^2) / 2, averaged
over specimens, always computed from a reference spectrum of 30 harmonics.
Cumulative power is reported both including and excluding the first harmonic,
labelled as such. The automatic H is the smallest number of harmonics reaching
99 percent of cumulative power excluding the first harmonic, bounded to 6..30.
If the power beyond the first harmonic is zero (every outline a circle), H is
6. The student can set H by hand between 2 and 30. A reconstruction panel shows
one selected specimen rebuilt with the sorted distinct values of 1, 2, 4, 8,
and H harmonics.

Validation: (i) forward coefficients match `Momocs::efourier(norm = FALSE)` to
1e-8 on a demo outline. Momocs places t = 0 at the last vertex, so the Momocs
input is the same polygon with its first vertex moved to the end. Reference
values are committed. (ii) A circle of radius 2 sampled at 2000 vertices gives
a_1 = d_1 = 2 within 1e-4 and all other coefficients below 1e-8 in absolute
value. (iii) A 3-by-1 ellipse sampled at 2000 vertices equally spaced in angle
matches committed Momocs values, including a_1 = 2.6245 and a_3 = 0.2418 to
four decimals; the first harmonic is not 3 and 1 because the parameter is arc
length (audit finding F2). (iv) The inverse at H = 100, evaluated at the
parameter values of the original vertices, reproduces a 200-point demo
configuration with a maximum vertex error below 0.5 percent of its centroid
size. (v) Power equals the formula applied to the Momocs coefficient matrix in
its A1..AH, B1..BH, C1..CH, D1..DH order.

### 7.3 PCA and stability

`prcomp` on the specimen-by-4H coefficient matrix, centered and unscaled. The
number of usable components r is the number whose variance exceeds 1e-12 of the
total, at most n - 1. Only those are reported, plotted, or exported. If the
total variance is zero (identical shapes) the app says so and shows no PCA.
Reported: share of variance per component and a scree plot.

Stability is described without inference, for n >= 4 and r >= 2. Each specimen
is left out in turn and the PCA refitted on the remaining rows of the same
coefficient matrix. H is held at its full-data value and the alignment is not
refitted. Two quantities are reported:

- Plane agreement (primary): the product of the cosines of the two principal
  angles between the plane of full-data PC1 and PC2 and the plane of the
  refitted PC1 and PC2. It is unchanged if the two components exchange order or
  rotate within the plane.
- Direction agreement for PC1 and for PC2: the absolute cosine between the
  full-data and refitted loading vectors of the same rank. It is reported for a
  component only when that component is separable, defined as its variance
  being at least 1.5 times that of the next component and, for PC2, at most
  two-thirds that of PC1. Otherwise the app states that the component is not
  separable from its neighbour and reports the plane only.

### 7.4 Shape change along components

For each of the first min(3, r) components, the mean coefficient vector plus
and minus 2 standard deviations of the scores along the component is passed
through the inverse transform and drawn over the mean shape. The morphospace
plot (r >= 2) shows PC scores with reconstructed shapes on a background grid,
coloured by group label or cluster. This replaces the tab named "Thin-Plate
Splines".

### 7.5 Clustering

Ward (D2) hierarchical clustering on Euclidean distances between coefficient
vectors, for n >= 4. The dendrogram is shown. The student sets k; values
outside 2..n-1 are clamped. Mean silhouette width is shown for each k in that
range as a descriptive guide. For specimen i with mean distance a to the other
members of its cluster and smallest mean distance b to another cluster, the
width is (b - a) / max(a, b). It is 0 for a specimen alone in its cluster and
0 when a = b = 0. The implementation is checked against `cluster::silhouette`.
When group labels exist, a cluster-by-group table is shown. No test statistic
or p-value is reported.

### 7.6 Regions

Regions are intervals of position from the tip (section 6.3) with boundaries
b1 < b2 < b3 (defaults 0.15, 0.60, 0.80, adjustable): tip [0, b1), blade
[b1, b2), shoulder [b2, b3), base [b3, 1]. Each point belongs to exactly one
region, and membership is the same on every specimen. Boundaries entered out of
order are sorted. Boundaries are drawn on the mean shape and on the outline
preview.

Outputs:

- Per-point variance: the sum of the x and y variances of the aligned
  coordinates at each point, drawn as colour on the mean shape.
- Per-region variance: the mean per-point variance within the region.
- Between-region covariation (n >= 4): the RV coefficient (Escoufier 1973)
  [CITE-CHECK: Escoufier 1973]. For column-centred blocks X (n by p) and Y
  (n by q) of aligned coordinates, RV = tr(XX'YY') / sqrt(tr(XX'XX') tr(YY'YY')).
  If either block has zero variance, RV is not available. Reference: the rows
  of Y are re-paired with the rows of X at random 999 times (seed 1), holding
  the aligned coordinates fixed, and the mean and the 2.5 to 97.5 percent range
  of the resulting RV values are drawn beside the observed value. No p-value is
  reported.

The reference is a re-pairing of fixed coordinates (author decision,
2026-11-17). It shows how large RV is at this sample size when the pairing of
the two regions is arbitrary. It is not an uncertainty interval, and it does
not remove covariation produced by aligning all regions together, because that
alignment is held fixed. The app states both points. The statistic is
descriptive, which the shared conventions permit (exception c); no claim is
made about the feasibility of a Bayesian alternative.

Stated limitations (shown in the app): RV is large when specimens are few
relative to variables. A single alignment of the whole outline spreads
variation across regions. Region boundaries are proportional, not anatomical.

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

Unit tests (testthat), run from the repository root with
`Rscript tests/testthat.R`. Tolerances are absolute unless stated.

Segmentation
- Otsu on a synthetic image with half its pixels at 0.2 and half at 0.8 returns
  a threshold strictly between 0.2 and 0.8.
- A synthetic image with a light disc containing a dark hole, a smaller light
  disc, and a light band along one border: the mask equals the filled large
  disc exactly; the smaller disc and the band are excluded.
- Each of the four demo images, unchanged, with automatic settings: no
  warnings; artifact class "light"; the mask does not touch the border; the
  complement of the mask is one 4-connected component (no holes); no mask pixel
  lies in the top 15 percent of rows of the image as displayed (the scale bar);
  intersection over union with the committed reference mask is at least 0.99.
  Reference masks are produced once by this pipeline, inspected against the
  photographs, and committed.

Orientation and points
- A kite-with-stem polygon rotated by 0, 37, 90, 180, and 265 degrees, with the
  declared direction nearest the rotated tip, and with the contour start vertex
  shifted by 0, 7, and 23: all cases return the same 200 points within 1e-8.
- An exact 45 degree case selects the major axis. A shape wider than long is
  oriented by the declared direction. A regular 60-gon triggers the
  poorly-determined-axis warning. A shape with a basal notch on the midline
  places B at the notch. A shape whose midline crosses a barb gives more than
  two crossings, places B at the lowest, and warns. A truncated triangle
  triggers the blunt-tip warning; the four demo outlines do not.
- Resampling of an asymmetric polygon: 200 distinct points, T only at index 1,
  B only at index 101, every segment length positive, spacing constant within
  each edge to 1e-8, each point in exactly one region.
- Manual landmarks on the same vertex give an error.

Measurements and scale
- A 40 by 100 rectangle: length 100, maximum width 40, area 4000, each within
  1e-8. An ellipse with semi-axes 50 and 20 at 2000 vertices: length 100 and
  width 40 within 0.01; widest position within 0.01 of 0.5. A triangle with its
  base down: widest position below 0.01. A notched shape: cross-section width
  spans the notch.
- Scale: the same synthetic image at resize factors 1 and 0.5 gives lengths in
  millimetres equal within 2 percent, both with its own two-click scale and
  with the shared scale set on the other image.

Fourier, alignment, PCA
- Items (i) to (v) of section 7.2.
- Alignment matches `geomorph::gpagen(Proj = FALSE)` reference coordinates for
  the demo configurations within 1e-5 after one global rotation.
- PCA scores equal `prcomp` scores of the committed Momocs coefficient matrix
  for the same aligned input within 1e-8 up to sign.
- Symmetrized shapes have asymmetry below 1e-10.
- Stability: for the rows (1.01, 0), (-1.01, 0), (0, 1), (0, -1), plane
  agreement is 1 within 1e-8 for every left-out row and neither component is
  reported as separable. With three specimens no stability result and no PC3
  are produced. Identical shapes produce the "identical" message.

Clustering and regions
- Silhouette widths equal `cluster::silhouette` within 1e-10 on random data
  for k = 2..5, on a partition with a singleton, and on data with duplicates.
- RV of a non-constant block with itself is 1 within 1e-12; RV equals the
  covariance-matrix trace form within 1e-12 on random data; RV is NA for a
  constant block; the reference is identical on repeat calls.

Store and downloads
- Saving a specimen under an existing id replaces it. Export then load into an
  empty store reproduces ids, names, groups, scale, settings, contours (1e-9),
  and measurements (1e-9). Loading a file with a missing column, a wrong schema
  version, or text in x is rejected with a message. Two different specimens
  with the same name both survive a load.
- Changing harmonics, region boundaries, clusters, or the symmetric option
  leaves saved contours and measurements identical; switching the symmetric
  option on and off restores identical PC scores.

Server test (`shiny::testServer` on the app directory): run the demo; every
figure builds; every interpretation text is non-empty; the PDF starts with
`%PDF` and has at least 10 pages; each CSV parses with its documented columns
and one row per specimen (or per vertex for outlines); with 0, 1, 2, and 3
specimens every output returns a message and no error; k = 99 is clamped.

shinytest2 is not used: its dependency websocket does not compile against the
OpenSSL 4.1 installed on the development machine (checked 2026-11-17). Browser
behaviour (drag to crop, clicks for scale and landmarks, tab rendering) is
checked by hand against a written checklist in `tests/manual_checklist.md`, and
that check is reported as manual.

Acceptance for the build milestone: all tests pass; the manual checklist is
completed on the demo images; `rsconnect::appDependencies("bifaceAnalysis")`
lists no package that is absent from CRAN.

## 11. Open points

- Whether 100 points per edge and fixed (non-sliding) points are adequate for
  the region analysis has not been tested. A sensitivity check at 50 and 200
  points per edge is a possible later addition.
- The citations marked CITE-CHECK must be verified against the sources before
  any of this text is reused in a publication. The app's own text names
  methods without attributing specific claims to those sources.

## 12. Revision history

Revision 2, 2026-11-17, after Codex audit 2026-11-17_spec. F1: artifact
selection rule replaced (6.1). F2: ellipse test corrected, circle test added
(7.2). F3: inverse, constant terms, start point, Momocs phase, reference
spectrum, zero-power case specified (7.2). F4: tie rules, base insertion,
click snapping, blunt-tip warning (6.2). F5: sampling fractions, positions,
mirror pairing, region intervals (6.3, 7.6). F6: plane agreement, separability
rule, rank gating (7.3, 7.4). F7: width definitions and scale transfer through
original pixels (6.1, 6.4). F8: state separation, CSV schema, symmetrization
(6.5, 7.1). F9: RV formula, reference relabelled, justification changed (7.6).
F10: silhouette specified, dependency sets listed (5, 7.5). F11: tolerances and
content checks (10). F12: D11 reworded (2).
