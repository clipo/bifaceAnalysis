---
date: 2026-11-17
mode: spec
title:
scope:
  - docs/specs/2026-11-17_collections_spec.md
  - docs/specs/2026-11-17_app_rewrite_spec.md
  - bifaceAnalysis/R/shape.R
  - bifaceAnalysis/R/efa.R
  - bifaceAnalysis/R/store.R
  - bifaceAnalysis/R/plots.R
  - bifaceAnalysis/R/interpret.R
  - bifaceAnalysis/R/report.R
  - bifaceAnalysis/app.R
note: Extension spec for an optional collections-by-location feature; the app it extends is built at this commit
prior:
review:
commit: 5481f0b2f346604195c2d4c6ebf070c6a158eaae
model: gpt-6-astra
effort: high
codex_version: codex-cli 0.169.0
started: 2026-11-17T19:58:21-05:00
message_sha256: 40fd85b9c0a1329e193861206105dd04389dcf6294b461a53d9852048d9d0dce
session_id: 01a1ec29-72d7-75f3-a441-e01da97c124f
sandbox: workspace-write [workdir, $TMPDIR]
finished: 2026-11-17T20:06:09-05:00
exit_status: ok
changed_outside_stem:
---

Mode: spec
Stem: 2026-11-17_spec_2
Audited commit: 5481f0b2f346604195c2d4c6ebf070c6a158eaae
Date: 2026-11-17
Repository root: the parent of your working directory (..); the scope paths in the brief are relative to it, so read src/x.md as ../src/x.md. Your own files go under the working directory.

## Project

A Shiny app in R that teaches outline-based shape analysis of bifaces to students: image segmentation, 200 corresponding outline points, Procrustes alignment, elliptical Fourier analysis, PCA, clustering, and region statistics. The app was rewritten against `docs/specs/2026-11-17_app_rewrite_spec.md` and is built at the audited commit. The document under audit, `docs/specs/2026-11-17_collections_spec.md`, specifies an optional extension that compares collections of bifaces (specimens sharing a group label) and, when a file of coordinates is supplied, maps the results. Nothing in it is built yet.

## Read in this order

1. `docs/specs/2026-11-17_collections_spec.md`, in full. This is the document under audit.
2. `docs/specs/2026-11-17_app_rewrite_spec.md`: sections 6.3 to 6.5, 7, 8, 9, and build notes (a) to (k) at the end. The new spec refers to these as "main 6.3" and so on.
3. `bifaceAnalysis/R/shape.R`: `gpa`, `no_variation`, `shape_pca`, `rv_coefficient`, `region_rv` (the pattern the new references are to follow).
4. `bifaceAnalysis/R/efa.R`: the coefficient layout and the inverse transform used for mean shapes.
5. `bifaceAnalysis/R/store.R`: `new_specimen`, `specimen_measurements`, `specimens_from_csv` (group labels, file limits).
6. `bifaceAnalysis/app.R`: the reactives that produce groups, the coefficient matrix, and the PCA; the Results tab layout; the download handlers; the report assembly.
7. `bifaceAnalysis/R/plots.R` (`category_scale`, `plot_region_rv`), `R/interpret.R`, and `R/report.R` as needed for the claims the new spec makes about them.

## Check

1. Section 5.2. Is the statement that the expected share SSB / SST under random reassignment with fixed sizes is (G - 1) / (N - 1) correct, and is 2 / 14 = 0.14 the right value for 3 collections of 5? Check by derivation and by simulation.
2. Section 5.3. Is `d_gh^2 - v_g / n_g - v_h / n_h` an unbiased estimate of the squared distance between population means under the stated assumption, with v_g defined with the n_g - 1 divisor? Is the description of the truncated square root accurate? Is the pair reference (reassignment within the pair, applied to the raw distance) a valid reference for the quantity shown beside it?
3. Section 5.4. Is the claim that the expected value of v_g does not depend on n_g correct, and is the leave-one-out range described accurately as a sensitivity range and not an interval?
4. Section 5.6. Are the haversine formula, the radius, and the four test values in section 10 (111.195 km, 10007.543 km, 0, and the meridian case) correct? Is the reference by reassigning positions among collections well defined, including the count of 120 assignments for 5 locations and the handling of coincident positions? Is the threshold of 5 located collections justified by the text?
5. Section 4.1. The collection statistics use the coefficient matrix, alignment, H, and PCA computed from all specimens, including those in no eligible collection. Does the built code make that possible without changing existing results, and does anything in the existing reactives contradict the spec's description (for example how groups, H, or the symmetric option reach the coefficient matrix)?
6. Section 6. Can the map be produced as specified with ggplot2 and the `maps` package: the `world`, `state`, and `county` databases through `ggplot2::map_data`, `coord_quickmap` with a zoomed window, the minimum window, and the scale bar rule? Identify any case in which the specified window, background, or scale bar is undefined or wrong (single position, positions at high latitude, a window crossing the edge of the contiguous United States box). `maps` may not be installed in your environment. If it is not, say so under not_examined.
7. Sections 4.2 and 7. Is the locations file specification complete enough to implement and test without further decisions: header matching, rejection conditions, matching of names, behaviour when a file is replaced, and what the downloads contain? Does any existing behaviour of group labels (trimming, empty labels, renaming in step 3, labels in a loaded outline file) conflict with "compared exactly after removing leading and trailing spaces"?
8. Section 7. Does the claim that `category_scale` repeats colours beyond 8 levels hold, and is the stated replacement sufficient for every figure that colours by collection?
9. Section 10. For each test, is the expected value correct and is the test able to fail if the implementation is wrong? Name any statistic or rule in sections 4 to 7 that has no test. Is the simulation test of the corrected distance specified tightly enough (number of draws, tolerance) to be reproducible and to detect an implementation that omits the correction?
10. Sections 2, 5, and 8. Does any wording describe a random-reassignment range as an uncertainty interval, a test, or evidence of a real difference, contrary to the stated descriptive intent? Is any required caution missing that a student would need in order not to over-read a map of 3 to 6 collections?
11. Anything in the specification that contradicts itself, contradicts the main specification, or cannot be built from the text without returning to the author.

## Note from the author

Extension spec for an optional collections-by-location feature; the app it extends is built at this commit

## Milestone focus

This is a specification or design written before any implementation. Audit whether it is complete enough to build from without going back to the author, whether its stated assumptions hold for the data, methods, and tools it names, whether any decision is unjustified or contradicts another, and whether anything it promises cannot be delivered with the stated tools.

## Rules

- Read every file in scope in full from the working tree. Do not rely on summaries the authors wrote about their own work when the underlying file is in scope.
- Treat every claim as unverified until you have checked it against the files or the tools.
- Classify each finding as wrong, unsupported, or plausible but not demonstrated (not_demonstrated). Give a concrete repair and an acceptance check that would show the repair holds.
- List what you checked and found correct under held. List what you did not examine, with the reason, under not_examined. List the questions you would have asked the authors under questions.
- Put scripts and logs under evidence/2026-11-17_spec_2/, anything derived from project data (extracts, copies, page images, large tables) under evidence/2026-11-17_spec_2/data/, and your prose report at audits/2026-11-17_spec_2_report.md, all relative to your working directory. Files under data/ are not committed, so every value your report or a finding cites must appear in the report itself or in a tracked log. List every file you wrote under evidence_files.
- Record `git rev-parse HEAD` in the report.
- If the sandbox blocks a command, record it under not_examined rather than working around it.
- Do not report style preferences. Do not pad the list. If you find nothing, say so and return an empty findings list.
- Respond only with JSON matching the provided schema.
