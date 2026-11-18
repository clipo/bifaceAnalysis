# Biface outline analysis

A Shiny app for teaching outline-based morphometrics with photographs of
bifaces and projectile points. It finds the outline of one artifact per image,
turns it tip up, takes measurements, and compares the outline shapes of the
saved specimens. Every result is shown with a note on how to read it.

The design, the methods, and the defects of the earlier version that this
rewrite corrects are documented in
`docs/specs/2026-11-17_app_rewrite_spec.md`.

## Running the app

```r
install.packages(c("shiny", "bslib", "ggplot2", "imager", "maps"))
shiny::runApp("bifaceAnalysis")     # from the repository root
```

The app does not use Momocs, which was archived from CRAN on 2026-03-18. The
elliptical Fourier and Procrustes code is in base R and is tested against
reference values from Momocs 1.4.1 and geomorph (see Tests).

## Using the app

1. **Images.** Use the four sample images or upload JPEG or PNG files. The
   demonstration button processes the sample images and opens the results.
2. **Outline.** State which way the tip points. Check the outline, tip, and
   base that the app finds. Drag a box to crop if other objects interfere. Set
   a scale by clicking the two ends of a known distance. Set the tip and base
   by hand if the automatic markers are misplaced. Name the specimen and save
   it.
3. **Specimens.** Review, rename, label, or remove saved specimens. Download
   the outlines as a CSV file and load that file later to continue.
4. **Results.** Outlines, measurements, harmonics, principal components, shape
   change along components, clusters, regions, and collections. A PDF report
   and CSV files of measurements, PC scores, Fourier coefficients, and
   collection summaries and distances can be downloaded.

## Comparing collections (optional)

The Collections tab compares groups of specimens, for example collections from
different sites. A collection is the set of specimens that share a group label.
Each collection needs at least 3 specimens, and at least 2 collections are
needed. The tab shows the mean shape of each collection, the share of shape
variation that lies between collections, distances between collection means
(raw, and corrected for sampling error), and the variation within each
collection. Each is shown beside the values produced by reassigning the same
specimens at random. These are descriptions of the samples, not tests.

A map is added when a locations file is loaded: a CSV file with the columns
`collection`, `latitude`, and `longitude` (decimal degrees, WGS84, north and
east positive), one row per collection, with names matching the group labels.
The tab offers a template with the current labels. With 5 or more located
collections, shape distance is also plotted against geographic distance. The
map uses the outlines of the `maps` package and draws no terrain. Locations are
kept for the session only and are not written to the outline file. The design
is in `docs/specs/2026-11-17_collections_spec.md`.

Analyses need a minimum number of specimens: 2 for the overlay, harmonics, and
region variance; 3 for principal components; 4 for clustering, the
leave-one-out check, and covariation between regions.

## Files

```
bifaceAnalysis/
  app.R            interface and reactive wiring
  R/image.R        thresholding, artifact selection, contour
  R/outline.R      orientation, tip and base, corresponding points, measurements
  R/efa.R          elliptical Fourier analysis
  R/shape.R        Procrustes alignment, symmetry, PCA, clusters, regions
  R/collections.R  comparison of collections, locations file, map geometry
  R/store.R        image loading, specimen records, CSV schema
  R/plots.R        one function per figure
  R/interpret.R    the "how to read this" text
  R/report.R       PDF report
  www/             sample images
tests/
  testthat.R       runs all tests
  testthat/        unit tests and server tests
  fixtures/        reference values and the script that generates them
  manual_checklist.md   browser checks that are not automated
```

## Tests

```r
Rscript tests/testthat.R            # from the repository root
```

The tests read `tests/fixtures/reference.rds` and need only the app's packages
plus testthat and cluster. Regenerating the reference values
(`Rscript tests/fixtures/make_reference.R`) needs Momocs 1.4.1 from the CRAN
archive and geomorph.

Interactions that require a browser (dragging a crop box, clicking on the
image, layout) are listed in `tests/manual_checklist.md`.

## Limits

The outline records plan shape only. Breaks and resharpening alter outlines
and are not distinguished from original form. Left and right depend on which
face was photographed. Results from a few specimens can change when one
specimen is added or removed, and the app reports how much. The Help page of
the app describes each method and its limits.

## Other files in this repository

The command-line script `biface_analysis.R` and the PDF and CSV files in the
repository root were produced by the earlier Momocs-based pipeline. They were
not changed by the rewrite and are not consistent with this app. The script
shares defects D1, D2, and D12 described in the specification and should not
be used for analysis without correction.
