# Design specification: comparison of collections by location

Date: 2026-11-17
Branch: `app-rewrite`
Extends: `docs/specs/2026-11-17_app_rewrite_spec.md` (revision 2 with build
notes a to k), as built at commit 8d387bf. Section numbers written as "main
6.3" refer to that document.
Status: design, revision 1, not yet audited and not yet built.

## 1. Purpose

The author asked for an optional feature that compares collections of bifaces
from different locations, and that shows the results on a map when the
coordinates of the locations are supplied. The feature must not be required:
the app must work as it does now for a user who has one group of bifaces and no
locations file.

The feature answers three descriptive questions for a student:

1. How different are the mean outline shapes of the collections?
2. How much of the shape variation lies between collections, and is that more
   than random grouping of the same specimens produces?
3. Do collections that are closer together on the ground have more similar
   mean shapes?

## 2. Decisions made by the author (2026-11-17)

- Scale of use: a class exercise, about 3 to 6 locations with 5 to 20 bifaces
  each. Thresholds are set for that case.
- Map: a static map with coastline and boundary outlines from the `maps`
  package, identical in the app and the PDF report. No interactive map and no
  requests to outside servers.
- Inference: descriptive in the app. The app reports effect sizes with
  random-reassignment references and exports scores by collection. No model is
  fitted in the app. This uses exception (c) of the shared conventions
  (descriptive methods). A multilevel model fitted in Stan to the exported
  scores is the route to any published inference and is outside this
  specification.
- The live app is not changed by this work. The feature is built on the branch
  and deployed to the test address only.

## 3. Scope

In scope: a "Collections" sub-tab in step 4, a locations file, six figures,
interpretation text, two CSV downloads, report pages, tests, and checklist
items.

Out of scope: any change to image processing, outlines, alignment, Fourier
analysis, or the existing tabs. Ordinations that maximize separation between
groups (canonical variates, between-group PCA) are excluded, because with more
variables than specimens they show separation between groups that are random
samples of one population [CITE-CHECK: source for this property of
between-group PCA; no PDF is in the repository]. A demonstration dataset is out
of scope because the repository holds four sample images and no set of
provenanced collections (section 12).

## 4. Definitions and inputs

### 4.1 Collections

A collection is the set of saved specimens that share a group label (main 6.5),
compared exactly after removing leading and trailing spaces. Specimens with an
empty label belong to no collection.

An eligible collection has at least 3 specimens. The comparison needs at least
2 eligible collections. Specimens that are unlabelled or in a collection with
fewer than 3 specimens are left out of every collection statistic, and the tab
states how many were left out and names the collections concerned.

All specimens, including those left out here, remain in the alignment, the
Fourier analysis, and the PCA of the other tabs. The collection statistics
therefore use the same aligned shapes, the same number of harmonics H, the same
symmetric-shape option, and the same PC axes as the rest of step 4. A
consequence, stated in the tab, is that adding or removing an unlabelled
specimen can change the collection results slightly through the alignment.

Shape data. X is the N by 4H matrix of Fourier coefficients (main 7.2) of the N
specimens in eligible collections. Row i is x_i. Collection g has n_g specimens
and mean vector m_g. The grand mean of the N rows is c.

### 4.2 Locations file (optional)

A CSV file with a header row and the columns `collection`, `latitude`,
`longitude`. Column names are matched without regard to case or surrounding
spaces. Other columns are ignored. Latitude and longitude are decimal degrees
on WGS84, north and east positive.

The file is rejected with a message, and any previously loaded locations are
kept, when: a required column is missing; the file has no rows or more than 200
rows; a collection name is empty or repeated; a coordinate is missing, not a
number, or not finite; a latitude is outside [-90, 90]; or a longitude is
outside [-180, 180]. The file size limit of the outline file (main 6.5) applies.

Rows are matched to collections by exact name after trimming. The tab lists
rows that match no collection and eligible collections that have no row. A
located collection is an eligible collection with a row. Two collections may
have the same coordinates.

Locations are held for the session only and are not written to the outline
file. A "Clear locations" button removes them. A "Download a locations
template" button returns a CSV with the three columns and one row per current
collection name, with empty coordinates.

The tab states that the file is sent to the server for the session, and that
coordinates of protected sites should be rounded or replaced with an
approximate position before upload.

## 5. Statistics

All statistics are computed in `R/collections.R` from X, the collection labels,
and the locations. Random reassignments use 999 draws and seed 1, and restore
the caller's random number state, as `region_rv` does (main 7.6).

### 5.1 Mean shapes

For each eligible collection the mean coefficient vector m_g is passed through
the inverse transform (main 7.2) with x0 = y0 = 0. The mean shapes are drawn
overlaid, one colour per collection.

### 5.2 Share of shape variation between collections

    SST = sum_i |x_i - c|^2
    SSW = sum_g sum_{i in g} |x_i - m_g|^2
    SSB = sum_g n_g |m_g - c|^2 = SST - SSW

The reported value is SSB / SST. It is not available when X has no variation
(main build note i).

Reference: the N specimens are reassigned to collections at random, keeping the
collection sizes, and SSB / SST is recomputed. The mean and the 2.5 to 97.5
percent range of the reassigned values are shown beside the observed value.
With random grouping the expected share is (G - 1) / (N - 1) for G collections,
which is far from zero in small samples: 2 / 14 = 0.14 for 3 collections of 5. The tab
states this number for the current data. No p-value is reported.

### 5.3 Distance between collection means

For collections g and h:

    raw:        d_gh  = |m_g - m_h|
    v_g         = sum_{i in g} |x_i - m_g|^2 / (n_g - 1)
    corrected:  dc_gh = sqrt(max(0, d_gh^2 - v_g / n_g - v_h / n_h))

The raw distance is larger on average than the distance between the population
means, because each sample mean carries sampling error of expected squared size
v_g / n_g. The quantity under the square root is an unbiased estimate of the
squared distance between population means when the specimens of each collection
are independent draws from that collection's population. dc_gh is 0 when the
raw distance is no larger than sampling error of that size produces. The square
root of an unbiased estimate is not itself unbiased, and the tab calls dc_gh an
approximate correction.

Reference for each pair: the n_g + n_h specimens of the two collections are
reassigned between them at random, keeping the two sizes, and the raw distance
is recomputed. The mean and the 2.5 to 97.5 percent range are shown beside the
observed raw distance. The tab states the number of pairs and that, with a 95
percent range and several pairs, an observed value can fall outside its range
through sampling alone.

Outputs: a pair table (collections, sizes, raw distance, reference mean and
range, corrected distance), a figure of observed raw distances against their
reference ranges in the form of the region RV figure, and a matrix figure of
corrected distances. Distances are in units of the coefficients of
unit-centroid-size shapes and are comparable only within one analysis.

No tree of collections is drawn. With 3 to 6 collections a tree adds nothing to
the matrix and suggests a branching order that the data cannot support.

### 5.4 Variation within each collection

v_g of section 5.3 is the within-collection shape variance. Its expected value
does not depend on n_g, so collections of different sizes can be compared. Its
precision does depend on n_g. To show this, v_g is recomputed with each
specimen of the collection left out in turn, and the smallest and largest of
those values are drawn as a range around v_g. The range describes sensitivity
to single specimens. It is not an uncertainty interval, and the tab says so.

### 5.5 Collection summaries

For each eligible collection: n_g; the mean score on each PC of the existing
PCA (main 7.3) for up to the first 3 components; the mean of length divided by
width; the mean position of maximum width (main 6.4); v_g; and the mean length
in mm. Mean length is given only when every specimen in every eligible
collection has a scale.

### 5.6 Shape distance against geographic distance

Requires at least 5 located collections with at least two distinct positions.
With fewer, the tab states the requirement and the number available. Three
locations give 3 pairs and four give 6, which cannot show a relationship.

Geographic distance is the great-circle distance by the haversine formula on a
sphere of radius 6371.0088 km:

    h = sin^2((lat2 - lat1) / 2) + cos(lat1) cos(lat2) sin^2((lon2 - lon1) / 2)
    distance = 2 R asin(min(1, sqrt(h)))

The figure plots dc_gh against distance in km, one labelled point per pair. The
reported value r is the Pearson correlation between the two sets of pairwise
distances, the matrix correlation of the Mantel test [CITE-CHECK: Mantel 1967;
no PDF is in the repository]. It is not available when either set has no
variation.

Reference: the positions are reassigned among the located collections and r is
recomputed. All G! assignments are used when G is at most 7 (5040 assignments),
and 999 random assignments otherwise. The mean and the 2.5 to 97.5 percent
range are reported with r. The tab states the number of pairs and the number of
assignments used, and that with 5 locations the reference has at most 120
values. No p-value is reported.

## 6. Map

Shown when at least 2 collections are located. Collections that are eligible
and not located are named below the map.

The map is a ggplot2 figure. Each located collection is a point at its
longitude and latitude, labelled with its name above the point, and coloured on
the viridis scale by a quantity the student selects from section 5.5: mean PC1
score (default), mean PC2 score, mean length divided by width, mean position of
maximum width, within-collection variance, number of specimens, and mean length
in mm when available.

Window. The bounding box of the located collections is extended by 20 percent
of its span on each side, with a minimum span of 0.5 degrees of latitude and
the equivalent ground distance in longitude at the middle latitude. The aspect
ratio is set by `coord_quickmap`.

Background. Country polygons from the `world` database of the `maps` package
are drawn as a light fill with thin grey boundaries. State boundaries from the
`state` database are added when the window overlaps the contiguous United
States (longitude -125 to -66, latitude 24 to 50), and county boundaries from
the `county` database when, in addition, the window spans less than 3 degrees
of latitude. A scale bar in km is drawn, with a length of 1, 2, or 5 times a
power of ten nearest to one quarter of the window's width at the middle
latitude.

Limits, stated below the map. The outlines are low resolution and are for
orientation only. Outside the contiguous United States a window smaller than a
country shows no boundary. Collections on both sides of the 180 degree meridian
are not drawn correctly, although their distances in section 5.6 are correct.
Labels of collections that are close together can overlap.

The author chose the static map knowing it has no terrain, rivers, or roads.
Section 12 records the consequence for small study areas.

## 7. Interface

A sub-tab "Collections" is added to step 4 after "Regions". It is always
present. With fewer than 2 eligible collections it shows only a short
explanation: group labels entered in step 2 or step 3 define collections, each
needs at least 3 specimens, and a locations file is optional. The other tabs
and the sidebar are unchanged, and no existing control gains a meaning it did
not have.

With 2 or more eligible collections the tab shows, in this order:

1. A sentence giving the collections in the analysis with their sizes, and any
   specimens left out (section 4.1).
2. Mean shapes (5.1) beside the "How to read this" text for the tab.
3. Share between collections with its reference (5.2).
4. Pair distances with references, the matrix of corrected distances, and the
   pair table (5.3).
5. Within-collection variance with leave-one-out ranges (5.4).
6. A "Locations" card: the file input, the template download, the clear button,
   the privacy statement (4.2), the list of unmatched names, the selector for
   the mapped quantity, the map (6), and the distance figure (5.6) or the
   sentence saying why it is withheld.

The PCA tab's "Colour specimens by: group label" option already colours the
morphospace by collection. The tab's text points to it.

Downloads added to the sidebar: "Collection summary (CSV)" with the columns of
section 5.5, and "Collection distances (CSV)" with the pair table of 5.3 and
the distance in km where both collections are located. Neither file contains
coordinates. Both are empty files with a header when the comparison is not
available.

Report. When the comparison is available, the PDF gains one page for each of
the figures above with its text. The map page is included only when the box
"Include the map in the PDF report" in the Locations card is ticked (default
ticked). The card states that the map shows the positions of the collections.

Categorical colours use the Okabe-Ito palette for up to 8 collections. With
more than 8 the discrete viridis scale is used, because the present
`category_scale` repeats colours beyond 8 levels.

## 8. Interpretation text

Built in `R/interpret.R` from computed values, under the rules of main
section 9. Content:

- The collections compared, their sizes, and the specimens left out.
- The observed share between collections, the share expected from random
  grouping at these sizes, and the reference range.
- The pair with the largest and the pair with the smallest corrected distance,
  and the pairs whose corrected distance is 0.
- The collections with the largest and smallest within-collection variance and
  their ratio, with the caution about small collections.
- For the distance figure: r, its reference range, the number of pairs, and the
  statement that a few locations cannot establish a spatial pattern.

Standing cautions, shown in full on the tab and in the report:

- A difference between collections is a difference between these samples. It
  can arise from the age of the artifacts, raw material, resharpening and
  breakage, how each collection was gathered, or local practice. The analysis
  does not distinguish among these.
- Left and right in an outline depend on which face was photographed. Unless
  every collection was photographed face-consistently, tick "Analyse symmetric
  shapes" before comparing collections.
- Collections photographed or outlined by different people, at different image
  resolutions, or with different detection settings can differ for that reason
  alone. Edge roughness from image resolution affects the higher harmonics.
- The PC axes describe variation among all specimens. They are not the
  directions that best separate the collections.

## 9. Dependencies and layout

One runtime package is added: `maps`. It is not installed on the development
machine and must be installed before the build. It is used through
`ggplot2::map_data`. New file `bifaceAnalysis/R/collections.R` holds the
statistics, the locations file reader, and the distance function. Figures go in
`R/plots.R` and text in `R/interpret.R`. New test file
`tests/testthat/test-collections.R`.

## 10. Testing and acceptance

Unit tests, tolerances absolute unless stated.

Partition
- For simulated coefficient matrices with 3 collections of unequal size,
  SSB + SSW equals SST to 1e-10 relative, and SSB / SST equals the value from
  `summary(lm(X ~ factor))` sums of squares pooled over columns.
- Collections that are identical copies of one set of rows give a share of 0.
  Collections with no within-collection variation and different means give 1.
- The mean of the reassigned shares over 999 draws is within 0.02 of
  (G - 1) / (N - 1) for a simulated matrix with no collection effect.
- The reference is identical on repeated calls, and the caller's random number
  state is unchanged by the call.

Distances
- Raw distance equals the Euclidean distance between row means.
- Corrected distance: over 4000 simulated pairs of samples (sizes 5 and 8, 12
  variables, known mean difference of length 2 and of length 0), the mean of
  the quantity under the square root is within 3 Monte Carlo standard errors of
  4 and of 0. The mean of the raw squared distance exceeds 4 by more than 3
  standard errors in the first case.
- A negative value under the root is reported as 0.
- Leave-one-out ranges contain v_g for a collection of identical rows (all
  zero) and bracket the values computed by hand for a 4-specimen example.

Geography
- Haversine: (0, 0) to (0, 1) is 111.195 km within 0.001. The equator to the
  pole is 10007.543 km within 0.001. A point to itself is 0. Two points across
  the 180 degree meridian at the equator, longitudes 179.5 and -179.5, are
  111.195 km apart.
- Matrix correlation: equals `cor` of the lower triangles. With 5 located
  collections all 120 assignments are used and include the identity. When shape
  distance is set equal to geographic distance, r is 1 and is the largest
  reference value.
- Not available (no error) when all positions coincide.

Locations file
- A valid file with mixed-case headers, extra columns, and surrounding spaces
  is read. Each rejection condition of section 4.2 returns its message and
  leaves earlier locations in place. Unmatched names on each side are listed.

Thresholds
- 1 eligible collection, 2 eligible collections, a collection of 2 specimens,
  unlabelled specimens, 4 located collections, and 5 located collections each
  give the outputs and messages of sections 4 to 6, with no error.

Figures and report
- Each new plotting function returns a ggplot object that builds without error
  or warning for 2, 3, 6, and 9 collections, and for a window inside and
  outside the contiguous United States. The map test is skipped with a message
  when `maps` is not installed, and the full suite is run with it installed
  before acceptance.
- The report has the expected number of pages with and without locations and
  with the map box unticked.

Server tests (`shiny::testServer`)
- The existing outputs are unchanged when group labels are added: values of the
  PCA, cluster, and region reactives are identical before and after.
- Loading, replacing, and clearing a locations file; the two downloads have the
  stated columns and no coordinate columns.

Existing tests continue to pass unchanged.

Manual checklist items are added for: the explanation with no collections; the
tab with 3 labelled collections built from the sample images by repeated
upload; loading a locations file and a rejected file; the map inside and
outside the United States; the selector for the mapped quantity; the template
download; the map box in the report; and the narrow layout of the tab.

Acceptance: all automated tests pass, the checklist items are run on the test
deployment, and the Codex build audit has no blocking finding open.

## 11. Performance

With 6 collections of 20 specimens and H = 30, X is 120 by 120. The partition
reference is 999 computations of group means. The pair references are 15 pairs
of 999 each. The matrix correlation uses at most 5040 assignments of a 7 by 7
matrix. Each is expected to finish in well under one second, and the build
must confirm a total under 2 seconds on the development machine for this case.
The statistics are computed in one reactive that depends on the coefficient
matrix, the labels, and the locations, so that changing the mapped quantity
does not repeat them.

## 12. Open points

1. No example data. The four sample images cannot form collections. The tests
   use simulated coefficient matrices, and the checklist uses the sample images
   uploaded repeatedly under different labels, which exercises the interface
   and not the archaeology. A demonstration needs provenanced images from the
   author.
2. Map context for small areas. A class project with sites a few kilometres
   apart gets county lines in the United States and an empty background
   elsewhere. The scale bar and labels remain. An interactive map with tiles
   would give context and was not chosen.
3. Two citations are flagged and unverified: Mantel 1967, and a source for the
   behaviour of group-separating ordinations with more variables than
   specimens.
4. The corrected distance assumes that the specimens of a collection are
   independent draws. Specimens from one cache or one knapper are not, and the
   correction is then too small. The tab states the assumption.

## 13. Revision history

Revision 1, 2026-11-17. First draft, after the author's four decisions in
section 2.
