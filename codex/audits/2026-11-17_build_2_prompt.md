---
date: 2026-11-17
mode: build
title:
scope:
  - docs/specs/2026-11-17_collections_spec.md
  - bifaceAnalysis/R/collections.R
  - bifaceAnalysis/R/plots.R
  - bifaceAnalysis/R/interpret.R
  - bifaceAnalysis/R/store.R
  - bifaceAnalysis/R/shape.R
  - bifaceAnalysis/R/efa.R
  - bifaceAnalysis/R/report.R
  - bifaceAnalysis/app.R
  - tests/testthat.R
  - tests/testthat/helper-app.R
  - tests/testthat/test-collections.R
  - tests/testthat/test-server.R
  - tests/manual_checklist.md
  - bifaceAnalysis/README.md
note: Build of the optional collections-by-location feature against spec revision 2; not yet deployed or run on real data
prior:
review:
commit: beedd0415ebbfad8c2e2c7c2f39c99c652c40290
model: gpt-6-astra
effort: high
codex_version: codex-cli 0.169.0
started: 2026-11-17T21:36:59-05:00
message_sha256: b766e95fac232da19d9da61645d311c6ecaf13d9b4209f4b85ac4919fad5b305
session_id: 01a1ec83-c12b-7840-8bc5-2213ef4b5db9
sandbox: workspace-write [workdir, $TMPDIR]
finished: 2026-11-17T21:47:02-05:00
exit_status: ok
changed_outside_stem:
---

Mode: build
Stem: 2026-11-17_build_2
Audited commit: beedd0415ebbfad8c2e2c7c2f39c99c652c40290
Date: 2026-11-17
Repository root: the parent of your working directory (..); the scope paths in the brief are relative to it, so read src/x.md as ../src/x.md. Your own files go under the working directory.

## Project

A Shiny app in R that teaches outline-based shape analysis of bifaces to students. The commit under audit adds an optional comparison of collections (specimens sharing a group label), with an optional locations file, a static map from the `maps` package (3.4.4, installed on the development machine), and a comparison of shape distance with geographic distance. It was built against `docs/specs/2026-11-17_collections_spec.md` revision 2, which was audited by you earlier (stem 2026-11-17_spec_2; ledger `docs/audits/codex_audit_2026-11-17_spec_2_ledger.md`). The spec's build notes (a) to (c) record deviations made during the build. Nothing has been deployed or run on real collections. The author's mutation run of 17 deliberate errors is not in the repository; do not rely on it.

## Read in this order

1. `docs/specs/2026-11-17_collections_spec.md`, in full, including build notes (a) to (c).
2. `bifaceAnalysis/R/collections.R`, in full.
3. `bifaceAnalysis/app.R`: the section headed "collections", the Collections tab in the UI, the figures and texts lists, the three new download handlers, and the report assembly in `dl_report`.
4. `bifaceAnalysis/R/plots.R`: `category_colours`, `category_scale`, and every function after the "collections" heading.
5. `bifaceAnalysis/R/interpret.R`: everything after the "collections" heading.
6. `bifaceAnalysis/R/store.R`: `specimens_from_csv` (label trimming).
7. `tests/testthat/test-collections.R` and the collections tests at the end of `tests/testthat/test-server.R`, with `tests/testthat/helper-app.R`.
8. `tests/manual_checklist.md` (section "Collections") and `bifaceAnalysis/README.md`.
9. `R/shape.R`, `R/efa.R`, and `R/report.R` as needed to check calls into them.

## Check

1. Does each statistic in `collections.R` implement the formula of spec section 5 exactly: SSB/SST, the share reference, d, v_g, Q, dc, the pair reference (pooling only the pair), leave-one-out ranges, the summary columns, haversine, the matrix correlation, and the complete or random assignments? Reproduce at least the share, Q, and the correlation on an independent fixture.
2. Do the references restore the caller's random number state in every path, including the shape-geography reference when it uses complete enumeration and when it uses 999 random assignments?
3. Section 4.1: are the collection statistics computed from the full-data coefficient matrix, H, and PCA, with only eligible rows used, and are rows, labels, scores, and measurements aligned in the same order in `coll_summary` and the downloads?
4. Section 4.2: does `read_locations` implement every rejection condition and the text-reading rule, and does the app replace the table only on an accepted file, rematch by current name, and build the template from all non-empty labels?
5. Section 6: does `map_window` implement the withholding conditions and the three ordered window steps, `map_layers` the background rule, and `scale_bar` the length and placement rule? Is build note (a), which replaces the spec's haversine check of the scale bar with a check along the parallel, correct?
6. Section 7: do the two CSV downloads have exactly the stated columns in order, empty cells for missing values, a header only when the comparison is unavailable, and no coordinates? Does the map selector offer only available quantities and fall back to the default? Is the map page in the report controlled by the tick box?
7. Section 8: is every standing caution present on the tab and in the report, and does any generated sentence describe a reference range as an interval, a test, or evidence of a real difference, or call the clipped distance unbiased?
8. Are there inputs that make any new reactive, plot, text, or download raise an error instead of giving a message: for example rank 0 or 1 PCA, all within-collection variances 0, identical collection means, coincident locations, a collection name with special characters, more than 8 collections, or a locations file whose names match no collection?
9. Do the tests in `test-collections.R` and the new server tests check what section 10 of the spec requires, with correct expected values, and could any of them pass with a wrong implementation? Name any statistic, rule, or output that has no test.
10. Did the two permitted changes to existing behaviour (label trimming on load, palette beyond 8 levels) change anything else in the existing tabs, downloads, or tests?
11. Run `Rscript tests/testthat.R` from the repository root and report the result. If the sandbox prevents it, record that under not_examined.

## Note from the author

Build of the optional collections-by-location feature against spec revision 2; not yet deployed or run on real data

## Milestone focus

This is an implementation reported as complete and not yet run on real data. Audit whether the code does what its specification and docstrings say, whether the tests would fail if the code were wrong, whether any parameter, seed, path, or unit is inconsistent between files, and whether running it on real data would produce a result the author would misread.

## Rules

- Read every file in scope in full from the working tree. Do not rely on summaries the authors wrote about their own work when the underlying file is in scope.
- Treat every claim as unverified until you have checked it against the files or the tools.
- Classify each finding as wrong, unsupported, or plausible but not demonstrated (not_demonstrated). Give a concrete repair and an acceptance check that would show the repair holds.
- List what you checked and found correct under held. List what you did not examine, with the reason, under not_examined. List the questions you would have asked the authors under questions.
- Put scripts and logs under evidence/2026-11-17_build_2/, anything derived from project data (extracts, copies, page images, large tables) under evidence/2026-11-17_build_2/data/, and your prose report at audits/2026-11-17_build_2_report.md, all relative to your working directory. Files under data/ are not committed, so every value your report or a finding cites must appear in the report itself or in a tracked log. List every file you wrote under evidence_files.
- Record `git rev-parse HEAD` in the report.
- If the sandbox blocks a command, record it under not_examined rather than working around it.
- Do not report style preferences. Do not pad the list. If you find nothing, say so and return an empty findings list.
- Respond only with JSON matching the provided schema.
