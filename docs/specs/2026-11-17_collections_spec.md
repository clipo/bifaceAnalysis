# Design specification: comparison of collections by location

Date: 2026-11-17
Branch: `app-rewrite`
Extends: `docs/specs/2026-11-17_app_rewrite_spec.md` (revision 2 with build
notes a to k), as built at commit 8d387bf. Section numbers written as "main
6.3" refer to that document.
Status: design, revision 2 (2026-11-17), not yet built. Revision 1 was audited
by Codex (stem 2026-11-17_spec_2). Revision 2 resolves findings F1 to F9; see
`docs/audits/codex_audit_2026-11-17_spec_2_ledger.md` and section 13.

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
- After the spec audit: the map is withheld beyond 80 degrees of latitude and
  across the 180 degree meridian (section 6); group labels are trimmed of
  surrounding spaces everywhere in the app (section 4.1); a new locations file
  replaces the old one, and locations are matched to collections by current
  name (section 4.2).

## 3. Scope

In scope: a "Collections" sub-tab in step 4, a locations file, seven figures,
interpretation text, two CSV downloads, report pages, tests, and checklist
items. Two changes to existing behaviour are permitted, and no others: group
labels are trimmed wherever they enter the app (section 4.1), and figures that
colour by group label use the palette rule of section 7 when there are more
than 8 labels.

Out of scope: any other change to image processing, outlines, alignment,
Fourier analysis, or the existing tabs. Ordinations that maximize separation
between groups (canonical variates, between-group PCA) are excluded, because
with more variables than specimens they show separation between groups that
are random samples of one population [CITE-CHECK: source for this property of
between-group PCA; no PDF is in the repository]. A demonstration dataset is out
of scope because the repository holds four sample images and no set of
provenanced collections (section 12).

## 4. Definitions and inputs

### 4.1 Collections and labels

Group labels are trimmed of leading and trailing white space wherever they
enter the app: when a specimen is saved in step 2, when a label is edited in
step 3, and when an outline file is loaded. A label that is empty after
trimming is empty. Labels are text and are never converted to numbers or to
missing values: a label "NA" or "007" is kept as written. Before this change
the outline file reader passed labels through untrimmed, so the existing PC
score export and group colouring could show a label with surrounding spaces
as a separate group. After it, every tab and export uses the same labels.

A collection is the set of saved specimens that share a non-empty label,
compared exactly after trimming.

An eligible collection has at least 3 specimens. The comparison needs at least
2 eligible collections. Specimens that are unlabelled or in a collection with
fewer than 3 specimens are left out of every collection statistic, and the tab
states how many were left out and names the collections concerned.

All specimens, including those left out here, remain in the alignment, the
Fourier analysis, and the PCA of the other tabs. The collection statistics
use the aligned shapes, the number of harmonics H, the symmetric-shape option,
the coefficient matrix, and the PC axes computed from all saved specimens.
Consequently, adding, removing, or editing any specimen, including one in no
eligible collection, can change the collection results through three routes:
the alignment, the automatic choice of H (main 7.2), and the PCA. The size of
the change is not bounded. The tab states these dependencies without a
qualifier of size.

Shape data. X is the N by 4H matrix of Fourier coefficients (main 7.2) of the N
specimens in eligible collections, taken from the full coefficient matrix.
Row i is x_i. Collection g has n_g specimens and mean vector m_g. The grand
mean of the N rows is c. G is the number of eligible collections.

### 4.2 Locations file (optional)

A CSV file, UTF-8, with a header row and the columns `collection`, `latitude`,
`longitude`. The file is read with every column as text and with no value
treated as missing, so that collection names are kept as written. Header names
are compared after trimming and conversion to lower case. The file is rejected
when a required normalized name occurs more than once (for example both
`latitude` and ` Latitude `). Other columns are ignored. Latitude and longitude
are converted from text to numbers after trimming. They are decimal degrees on
WGS84, north and east positive.

The file is rejected with a message, and the locations already loaded are kept
unchanged, when: a required column is missing or repeated; the file has no
data rows or more than 200; a collection name is empty after trimming or is
repeated after trimming; a coordinate is empty, does not convert to a number,
or is not finite; a latitude is outside [-90, 90]; or a longitude is outside
[-180, 180]. Files are limited to 40 MB, the limit the app already sets for
all uploads (`shiny.maxRequestSize` in `app.R`).

A file that is accepted replaces the whole locations table. Rows are not
merged with an earlier file. Locations are matched to collections by exact
trimmed name each time the collections change, so a collection renamed in
step 3 loses its location until the file or the name is corrected. The tab
lists file rows that match no eligible collection and eligible collections
that have no row. A located collection is an eligible collection with a
matching row. Two collections may have the same coordinates.

Locations are held for the session only and are not written to the outline
file. A "Clear locations" button removes them. A "Download a locations
template" button returns a CSV with the three columns and one row for every
non-empty label among the saved specimens, eligible or not, in alphabetical
order. Coordinates are filled for labels that have a loaded location and are
empty otherwise.

The tab states that the file is sent to the server for the session, and that
coordinates of protected sites should be rounded or replaced with an
approximate position before upload. It also states that approximate positions
change the geographic distances of section 5.6.

## 5. Statistics

All statistics are computed in `R/collections.R` from X, the collection labels,
and the locations. Random reassignments use 999 draws and seed 1, and restore
the caller's random number state, as `region_rv` does (main 7.6). Every
reference describes the variation produced by reassigning these fixed
specimens. It is not an uncertainty interval for a population quantity, and
the app never calls it one.

### 5.1 Mean shapes

For each eligible collection the mean coefficient vector m_g is passed through
the inverse transform (main 7.2) with x0 = y0 = 0. Because the transform is
linear, this equals the mean of the specimens' reconstructions at the same H.
The mean shapes are drawn overlaid, one colour per collection.

### 5.2 Share of shape variation between collections

    SST = sum_i |x_i - c|^2
    SSW = sum_g sum_{i in g} |x_i - m_g|^2
    SSB = sum_g n_g |m_g - c|^2 = SST - SSW

The reported value is SSB / SST. It is not available when X has no variation
(main build note i).

Reference: the N specimens are reassigned to collections at random, keeping the
collection sizes, and SSB / SST is recomputed. The mean and the 2.5 to 97.5
percent range of the reassigned values are shown beside the observed value.
With random grouping of fixed rows into fixed sizes, the expected share is
(G - 1) / (N - 1), which is far from zero in small samples: 2 / 14 = 0.14 for
3 collections of 5. The tab states this number for the current data. No
p-value is reported.

### 5.3 Distance between collection means

For collections g and h:

    raw distance:      d_gh = |m_g - m_h|
    within variance:   v_g  = sum_{i in g} |x_i - m_g|^2 / (n_g - 1)
    corrected square:  Q_gh = d_gh^2 - v_g / n_g - v_h / n_h
    corrected:         dc_gh = sqrt(max(0, Q_gh))

Assumptions. The coefficient vectors of each collection are independent draws
from that collection's population, the two collections are independent of each
other, the representation (alignment, H, and the coefficient space) is fixed,
and second moments are finite. Under these assumptions the raw squared
distance d_gh^2 is larger on average than the squared distance between the
population means by v_g / n_g + v_h / n_h in expectation, and Q_gh is an
unbiased estimate of that squared distance. Equal covariance matrices are not
required.

Two operations make dc_gh differ from the unbiased Q_gh. Setting a negative
Q_gh to 0 raises the average, and taking the square root changes it. When the
population means are equal, Q_gh is negative in about half of samples, and the
average of max(0, Q_gh) is positive. A corrected distance of 0 means that the
raw distance is no larger than sampling error of the expected size produces.
It does not show that the populations have the same mean shape.

In this app the alignment, the automatic H, and the coefficients are fitted to
the same specimens, and specimens from one cache or one knapper are not
independent. The correction does not account for either, and the tab states
this.

Reference for each pair: the n_g + n_h specimens of the two collections are
reassigned between them at random, keeping the two sizes, and the raw distance
is recomputed. The mean and the 2.5 to 97.5 percent range are shown beside the
observed raw distance. The tab states the number of pairs, that pairs that
share a collection are not independent of one another, and that with a 95
percent range and several pairs an observed value can fall outside its range
through sampling alone.

Outputs: a pair table (collections, sizes, raw distance, reference mean and
range, corrected square, corrected distance), a figure of observed raw
distances against their reference ranges in the form of the region RV figure,
and a matrix figure of corrected distances. Distances are in units of the
coefficients of unit-centroid-size shapes and are comparable only within one
analysis.

No tree of collections is drawn. With 3 to 6 collections a tree adds nothing to
the matrix and suggests a branching order that the data cannot support.

### 5.4 Variation within each collection

v_g of section 5.3 is the within-collection shape variance. Under the
assumptions of section 5.3 its expected value does not depend on n_g, so
collections of different sizes can be compared. Its precision does depend on
n_g. To show this, v_g is recomputed with each specimen of the collection left
out in turn, and the smallest and largest of those values are drawn as a range
around v_g. The range describes sensitivity to single specimens. It is not an
uncertainty interval, and the tab says so. A collection of 3 specimens gives
leave-one-out values computed from 2 specimens each.

Ratio in the text. When every v_g is positive the text gives the ratio of the
largest to the smallest. When some but not all are 0 the text names the
collections whose specimens do not differ and gives no ratio. When all are 0
the text says that no collection shows variation.

### 5.5 Collection summaries

For each eligible collection: n_g; the mean score on each existing PC of the
PCA (main 7.3), for up to the first 3 components; the mean of length divided by
width; the mean position of maximum width (main 6.4); v_g; and the mean length
in mm. A PC column is empty when the PCA has fewer components than its number
(rank r of main 7.3, which can be 0, 1, or 2 with few or similar specimens).
Mean length is given only when every specimen in every eligible collection has
a scale, and is empty otherwise.

### 5.6 Shape distance against geographic distance

Requires at least 5 located collections with at least two distinct positions.
This threshold is an operational choice for a classroom tool and not a derived
limit. With fewer, the tab states the requirement and the number available.
The figure has one point per pair of collections, not per collection: 5
collections give 10 pairs, and pairs that share a collection are not
independent.

Geographic distance is the great-circle distance by the haversine formula on a
sphere of radius R = 6371.0088 km, with latitudes and longitudes converted from
degrees to radians:

    h = sin^2((lat2 - lat1) / 2) + cos(lat1) cos(lat2) sin^2((lon2 - lon1) / 2)
    distance = 2 R asin(min(1, sqrt(h)))

The figure plots dc_gh against distance in km, one labelled point per pair. The
reported value r is the Pearson correlation between the two sets of pairwise
distances, the matrix correlation of the Mantel test [CITE-CHECK: Mantel 1967;
no PDF is in the repository]. It is not available when either set has no
variation.

Reference: the positions are reassigned among the located collections and r is
recomputed. All G_L! assignments of positions to the G_L located collections
are used when G_L is at most 7 (5040 assignments), including the observed
assignment, and 999 random assignments otherwise. Assignments are counted as
permutations of collections, so when positions coincide some assignments give
the same arrangement and are counted with their multiplicity. The mean and the
2.5 to 97.5 percent range are reported with r. The tab states the number of
pairs and the number of assignments used, and that with 5 locations the
reference has 120 values. No p-value is reported.

## 6. Map

Shown when at least 2 collections are located and none of the conditions
below withholds it. Collections that are eligible and not located are named
below the map.

Withholding. The map is withheld, with a sentence giving the reason, when any
located collection has a latitude beyond 80 degrees north or south, or when
the longitudes of the located collections span more than 180 degrees (a set
that the window rule below would draw across the wrong side of the globe).
Distances in section 5.6 and all other results are unaffected.

The map is a ggplot2 figure. Each located collection is a point at its
longitude and latitude, labelled with its name above the point, and coloured on
the viridis scale by a quantity the student selects from section 5.5. The
choices are the mean score of each PC that exists (up to 3), mean length
divided by width (the default), mean position of maximum width,
within-collection variance, number of specimens, and mean length in mm when
available. When the selected quantity becomes unavailable (for example a PC
that no longer exists after a specimen is removed) the selection returns to
the default.

Window, computed in this order:

1. Take the minimum and maximum latitude and longitude of the located
   collections. Let lat_mid be the middle latitude.
2. Latitude span: the larger of 1.4 times the latitude range and 0.5 degrees.
   Longitude span: the larger of 1.4 times the longitude range and
   0.5 / cos(lat_mid) degrees. These are centred on the middle of the ranges.
   Coincident positions give ranges of 0 and therefore the minimum spans.
3. Clip the latitude limits to [-85, 85] and the longitude limits to
   [-180, 180].

The window is applied with `coord_quickmap(xlim, ylim, expand = FALSE)`. The
polygon data are not filtered or clipped before drawing, so polygons that
cross the edge of the window keep their shape.

Background. Country polygons from the `world` database of the `maps` package,
read with `ggplot2::map_data`, are drawn as a light fill with thin grey
boundaries. State boundaries from the `state` database are added when the
window overlaps the box of the contiguous United States (longitude -125 to -66,
latitude 24 to 50). County boundaries from the `county` database are added
when the window overlaps that box and its latitude span is below 3 degrees.

Scale bar. The window width in km is the haversine distance between the left
and right edges at lat_mid. The bar length L is the value of the form 1, 2, or
5 times a power of ten nearest to one quarter of that width. The bar is drawn
along the parallel at 8 percent of the window height above the bottom edge,
starting 5 percent of the window width from the left edge, and its length in
degrees of longitude is L / (111.19508 cos(lat_bar)), where lat_bar is the
latitude of the bar. It is labelled "L km (approximate)".

Limits, stated below the map. The outlines are low resolution and are for
orientation only. A small window may show no boundaries, and outside the
contiguous United States it shows country borders only. The scale bar is
correct along its own parallel and approximate elsewhere in the window.
Labels of collections that are close together can overlap. The colours show
means of samples, and a mean from a small collection can differ from the
population mean by much more than the differences between collections.

The author chose the static map knowing it has no terrain, rivers, or roads.
Section 12 records the consequence for small study areas.

## 7. Interface

A sub-tab "Collections" is added to step 4 after "Regions". It is always
present. With fewer than 2 eligible collections it shows only a short
explanation: group labels entered in step 2 or step 3 define collections, each
needs at least 3 specimens, and a locations file is optional. The other tabs
and the sidebar are unchanged except as permitted in section 3.

With 2 or more eligible collections the tab shows, in this order:

1. A sentence giving the collections in the analysis with their sizes, and any
   specimens left out (section 4.1).
2. Mean shapes (5.1) beside the "How to read this" text for the tab.
3. Share between collections with its reference, drawn as a figure in the form
   of the region RV figure with one row (5.2).
4. Pair distances with references, the matrix of corrected distances, and the
   pair table (5.3).
5. Within-collection variance with leave-one-out ranges (5.4).
6. A "Locations" card: the file input, the template download, the clear button,
   the privacy statement (4.2), the list of unmatched names, the selector for
   the mapped quantity, the map (6) or the sentence saying why it is withheld,
   and the distance figure (5.6) or the sentence saying why it is withheld.

The seven figures are: mean shapes, share, pair distances, distance matrix,
within-collection variance, map, and shape distance against geographic
distance.

The PCA tab's "Colour specimens by: group label" option colours the
morphospace by the same trimmed labels. The tab's text points to it.

Downloads added to the sidebar. Missing values are written as empty cells.
When the comparison is not available each file has its header row only.
Neither file contains coordinates.

- "Collection summary (CSV)", one row per eligible collection, columns in this
  order: `collection`, `n`, `pc1_mean`, `pc2_mean`, `pc3_mean`,
  `length_width_mean`, `max_width_position_mean`, `within_variance`,
  `length_mm_mean`, `located`. The PC columns are always present and empty
  where the component does not exist. `located` is TRUE or FALSE.
- "Collection distances (CSV)", one row per pair, columns in this order:
  `collection_a`, `collection_b`, `n_a`, `n_b`, `raw_distance`,
  `reference_mean`, `reference_low`, `reference_high`, `corrected_square`,
  `corrected_distance`, `distance_km`. `distance_km` is empty unless both
  collections are located.

Report. When the comparison is available, the PDF gains one section for each
of the figures that is shown, with its text. A section may continue onto text
pages, as existing sections do (`report_section`). The map section is included
only when the box "Include the map in the PDF report" in the Locations card is
ticked (default ticked). The card states that the map shows the positions of
the collections.

Colours. A single palette function replaces `category_scale` for every figure
that colours by group label or collection, for both the colour and the fill
aesthetics: the Okabe-Ito palette for up to 8 levels, and the discrete viridis
scale for more than 8. The present `category_scale` repeats colours beyond 8
levels (9 levels give 8 distinct colours).

## 8. Interpretation text

Built in `R/interpret.R` from computed values, under the rules of main
section 9. Content:

- The collections compared, their sizes, and the specimens left out.
- The observed share between collections, the share expected from random
  grouping at these sizes, and the reference range.
- The pair with the largest and the pair with the smallest corrected distance,
  and the pairs whose corrected distance is 0, with the statement that 0 does
  not show equal mean shapes.
- The within-collection variance statement of section 5.4, with the caution
  about small collections.
- For the distance figure: r, its reference range, the number of collections
  and of pairs, that pairs sharing a collection are not independent, and that
  a few locations cannot establish a spatial pattern.

Standing cautions, shown in full on the tab and in the report:

- A difference between collections is a difference between these samples. It
  can arise from the age of the artifacts, raw material, resharpening and
  breakage, how each collection was gathered, or local practice. The analysis
  does not distinguish among these.
- The reference ranges show what reassigning these same specimens produces.
  They are not confidence or credible intervals, and a value outside its range
  is not evidence of a real difference between populations or places.
- Left and right in an outline depend on which face was photographed. Unless
  every collection was photographed face-consistently, tick "Analyse symmetric
  shapes" before comparing collections.
- Collections photographed or outlined by different people, at different image
  resolutions, or with different detection settings can differ for that reason
  alone. Edge roughness from image resolution affects the higher harmonics.
- The PC axes describe variation among all specimens. They are not the
  directions that best separate the collections.
- Map colours are sample means. Rounded or approximate coordinates change the
  geographic distances.

## 9. Dependencies and layout

One runtime package is added: `maps`. It is not installed on the development
machine and must be installed before the build. It is used through
`ggplot2::map_data`. New file `bifaceAnalysis/R/collections.R` holds the
statistics, the locations file reader, and the distance function. Figures go in
`R/plots.R` and text in `R/interpret.R`. The trimming of labels is added to
`specimens_from_csv` in `R/store.R`. New test file
`tests/testthat/test-collections.R`.

## 10. Testing and acceptance

Unit tests, tolerances absolute unless stated. Every simulation uses
`set.seed(20261117)` at its start and independent standard normal values
unless stated.

Labels
- An outline file with labels "A", " A ", "  ", "", "NA", and "007" loads with
  labels "A", "A", "", "", "NA", and "007". The first two form one collection.
  Labels saved in step 2 and edited in step 3 are trimmed (server test).
- PCA, cluster, and region values are identical before and after group labels
  are added to the same specimens (server test).

Partition
- For simulated coefficient matrices with 3 collections of sizes 3, 5, and 8
  and 10 variables, SSB + SSW equals SST to 1e-10 relative, and SSB / SST
  equals the ratio of the between and total sums of squares from
  `lm(X ~ factor)` pooled over columns.
- Collections that are identical copies of one set of rows give a share of 0.
  Collections with no within-collection variation and different means give 1.
- For collections of sizes 2, 2, and 2 the 90 distinct assignments are
  enumerated independently in the test. Every reassigned share returned is one
  of the enumerated values, and the mean of 999 reassigned shares is within
  0.02 of the enumerated mean, which equals (G - 1) / (N - 1) = 0.4 to 1e-10.
- The reference is identical on repeated calls, and the caller's random number
  state is unchanged by the call.

Distances
- Raw distance equals the Euclidean distance between row means to 1e-12.
- Hand examples: collection A rows (0, 0), (2, 0) and collection B rows
  (5, 0), (7, 0) give d^2 = 25, v_A = v_B = 2, Q = 23. Collection B rows
  (1, 0), (1, 2) give d^2 = 1, v_B = 2, Q = -1, and corrected distance 0.
- Simulation, 4000 replicates, sizes 5 and 8, 12 variables, population mean
  difference of length 2 (on the first variable) and of length 0: the mean of
  Q is within 3 Monte Carlo standard errors of 4 and of 0. The mean of the raw
  squared distance exceeds 4 by more than 3 standard errors in the first case.
  The mean of max(0, Q) exceeds 0 by more than 3 standard errors in the
  second case. An implementation without the correction, or with divisor n_g
  in v_g, fails at least one of these.
- The pair reference pools only the two collections of the pair: for a fixture
  in which a third collection is far from both, every reassigned raw distance
  is no larger than the largest distance obtainable from the pair's own rows.
- Leave-one-out: for a 4-specimen collection the four leave-one-out values are
  computed by hand in the test, and the returned smallest and largest equal
  them to 1e-12. A collection of identical rows gives 0, 0, and 0.

Summaries and shapes
- Each column of the summary equals a value computed directly in the test from
  the PC scores and measurements for a fixture, to 1e-12. With PCA rank 1 the
  PC2 and PC3 columns are empty. Mean length in mm is empty when one specimen
  lacks a scale.
- The mean shape equals the mean of the specimens' reconstructions at the same
  H to 1e-10.
- Variance ratio text for the three cases of section 5.4.

Geography
- Haversine: (0, 0) to (0, 1) is 111.195080 km within 1e-6. The equator to the
  pole, (0, 0) to (90, 0), is 10007.557221 km within 1e-6. A point to itself is
  0. Longitudes 179.5 and -179.5 at the equator are 111.195080 km apart within
  1e-6. An implementation that omits the conversion to radians fails.
- Matrix correlation equals `cor` of the lower triangles. With 5 located
  collections the 120 assignments are enumerated independently in the test and
  the returned reference values equal them as a sorted vector. When shape
  distance is set equal to geographic distance, r is 1 and is the largest
  reference value. With two positions coincident the reference still has 120
  values.
- Not available (no error) when all positions coincide.

Locations file
- A valid file with mixed-case headers, surrounding spaces, and extra columns
  is read, with names "NA" and "007" kept as text. Each rejection condition of
  section 4.2, including a repeated normalized header, returns its message and
  leaves earlier locations unchanged. An accepted second file replaces the
  first. A renamed collection becomes unmatched. Unmatched names on each side
  are listed. The template has one row per non-empty label with the stated
  columns.

Map
- With `maps` installed: for fixtures of coincident positions, a single
  positive span, positions at 79 degrees, positions inside, outside, and
  straddling the contiguous United States box, and a window under and over 3
  degrees of latitude, the built plot's limits equal the window rule to 1e-9
  and the layers present equal the background rule. The scale bar's endpoints
  are L km apart by the haversine distance within 0.1 percent, and L is the
  value given by the rule.
- Positions at 81 degrees, and longitudes spanning more than 180 degrees,
  withhold the map with its message.
- The map test is skipped with a message when `maps` is not installed. The full
  suite is run with it installed before acceptance.

Figures and report
- Each new plotting function returns a ggplot object that builds without error
  or warning for 2, 3, 6, and 9 collections. The built data have one mean shape
  per collection, one point per pair, and, for 9 collections, 9 distinct
  colours in both the colour and the fill scales.
- Unavailable states give their messages and no error: PCA rank 0 and 1, all
  within-collection variances 0, fewer than 5 located collections, and a
  withheld map.
- The report's sections and page count match the figures shown, with and
  without locations and with the map box unticked, allowing text continuation
  pages.

Server tests (`shiny::testServer`)
- Loading, replacing, and clearing a locations file. The two downloads have the
  stated columns in order, empty cells for missing values, a header only when
  the comparison is not available, and no coordinate columns.

Existing tests continue to pass unchanged, apart from tests updated for label
trimming.

Manual checklist items are added for: the explanation with no collections; the
tab with 3 labelled collections built from the sample images by repeated
upload; loading a locations file, a replacement file, and a rejected file;
renaming a located collection; the map inside and outside the United States
and a withheld map; the selector for the mapped quantity; the template
download; the map box in the report; and the narrow layout of the tab.

Acceptance: all automated tests pass, the checklist items are run on the test
deployment, and the Codex build audit has no blocking finding open.

## 11. Performance

With 6 collections of 20 specimens and H = 30, X is 120 by 120. The partition
reference is 999 computations of group means. The pair references are 15 pairs
of 999 each. The matrix correlation uses at most 5040 assignments of a 7 by 7
matrix. The build must confirm a total under 2 seconds on the development
machine for this case. The statistics are computed in one reactive that
depends on the coefficient matrix, the labels, and the locations, so that
changing the mapped quantity does not repeat them.

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
4. The corrected distance assumes independent specimens in a fixed
   representation (section 5.3). Neither holds exactly in the app, and the
   size of the resulting error is not established.

## 13. Revision history

Revision 1, 2026-11-17. First draft, after the author's four decisions in
section 2.

Revision 2, 2026-11-17. Resolves the Codex spec audit (stem 2026-11-17_spec_2).
F1: equator-to-pole test value corrected to 10007.557221 km and radian
conversion stated. F2: window order, clipping, withholding beyond 80 degrees
and across the meridian, unfiltered polygons, and scale-bar rule defined;
boundary limitation reworded. F3: PC selectors depend on rank, with a default
fallback and empty values; variance-ratio text for zero variances. F4: labels
trimmed everywhere and kept as text; one palette function for colour and fill.
F5: dependencies on alignment, automatic H, and PCA named; size qualifier
removed. F6: Q defined separately, assumptions stated, clipping and square-root
effects described, and the meaning of a corrected 0 stated. F7: header,
text-reading, replacement, rematching, template, file-size, and CSV schema
rules defined. F8: seeded and specified simulations, exact hand and enumerated
oracles, and tests for summaries, shapes, plot data, map geometry, and
unavailable states. F9: the threshold of 5 stated as operational; cautions on
fixed-data references, dependent pairs, sample means, and approximate
coordinates added. The figure count is corrected from six to seven.

## Build notes

(a) Section 10, scale bar test (2026-11-17). Revision 2 asked for the bar's
endpoints to be L km apart by the haversine distance within 0.1 percent. The
bar is drawn along a parallel with the length rule of section 6, and the
haversine distance between two points on one parallel is the great-circle
distance, which is shorter than the distance along the parallel. For a bar of
1000 km at 60 degrees the difference is about 0.3 percent, so the test as
written would fail the rule it checks. The test instead checks that the length
along the parallel, 111.19508 cos(lat_bar) times the longitude extent, equals L
to 1e-9. The rule and the label "approximate" are unchanged.
(b) Section 7, report (2026-11-17). The standing cautions of section 8 form
their own report section, "Comparison of collections: before drawing
conclusions", placed after the mean shapes.
(c) Section 11, performance (2026-11-17). Six collections of 20 specimens with
H = 30 took 0.46 s for the share, pair, and within-collection statistics, and
7 located collections (5040 assignments) took 0.07 s, on the development
machine.

