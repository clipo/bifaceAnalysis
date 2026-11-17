---
date: 2026-11-17
mode: build
title:
scope:
  - docs/specs/2026-11-17_app_rewrite_spec.md
  - bifaceAnalysis/app.R
  - bifaceAnalysis/R
  - tests/testthat.R
  - tests/testthat
  - tests/fixtures/make_reference.R
  - tests/fixtures/reference.rds
  - tests/manual_checklist.md
  - bifaceAnalysis/README.md
  - bifaceAnalysis/www
note: Rewrite of the Shiny app against spec revision 2; build notes at the end of the spec list three known differences
prior:
review:
commit: 5b0ec11c1b7af780c9d37e79b6f2037332c96d4b
model: gpt-6-astra
effort: high
codex_version: codex-cli 0.169.0
started: 2026-11-17T16:00:53-05:00
message_sha256: a7faae5a5d4bb0fadcfec1aaca64a43fbef354b3cd890277c153b5e5cf4b5b7f
session_id: 01a1eb50-0bac-7ed0-92b7-7a69c0d6f649
sandbox: workspace-write [workdir, $TMPDIR]
finished: 2026-11-17T16:18:20-05:00
exit_status: ok
changed_outside_stem:
---

Mode: build
Stem: 2026-11-17_build
Audited commit: 5b0ec11c1b7af780c9d37e79b6f2037332c96d4b
Date: 2026-11-17
Repository root: the parent of your working directory (..); the scope paths in the brief are relative to it, so read src/x.md as ../src/x.md. Your own files go under the working directory.

## Project

An R Shiny app for teaching outline-based morphometrics of bifaces and projectile points to students. It segments one artifact per photograph, standardizes the outline (tip up, 200 corresponding points anchored at tip and base), takes measurements, and compares shapes by Procrustes alignment, base-R elliptical Fourier analysis, PCA, Ward clustering, and region statistics, with generated interpretation text and a PDF report. This build replaces an earlier single-file app and was written against revision 2 of the specification, which was itself revised after a prior audit (2026-11-17_spec).

## Read in this order

1. `docs/specs/2026-11-17_app_rewrite_spec.md`: sections 5 to 10 define the required behaviour and tests. The "Build notes" at the end of section 12 list three differences the author already knows about.
2. `bifaceAnalysis/R/image.R` (`segment_image`, `open_square`, `component_table`, `mask_outline`), then `bifaceAnalysis/R/outline.R` (`axis_direction`, `standardize_outline`, `points_from_outline`, `measure_outline`).
3. `bifaceAnalysis/R/efa.R` (`efa_forward`, `efa_inverse`, `harmonic_power`, `choose_harmonics`) and `bifaceAnalysis/R/shape.R` (`gpa`, `symmetrize_shape`, `shape_pca`, `pca_stability`, `mean_silhouette`, `rv_coefficient`, `region_rv`).
4. `bifaceAnalysis/R/store.R` (`load_working_image`, `process_image`, `specimens_to_csv`, `specimens_from_csv`, `merge_specimens`), `bifaceAnalysis/R/plots.R`, `bifaceAnalysis/R/interpret.R`, `bifaceAnalysis/R/report.R`.
5. `bifaceAnalysis/app.R`: the server function, in particular the scale handling (`own_ppm`, `current_ppm`, `rv$shared_ppm`), `processed`, the save and edit observers, `figures`, `texts`, and the download handlers.
6. `tests/testthat/` (helper and six test files), `tests/fixtures/make_reference.R`, `tests/manual_checklist.md`, `bifaceAnalysis/README.md`.

The tests run from the repository root with `Rscript tests/testthat.R`. `tests/fixtures/reference.rds` holds values from Momocs 1.4.1 and geomorph. Regenerating it needs those packages. The tests do not.

## Check

1. Does each function in `bifaceAnalysis/R/` do what the specification section it implements requires? Name every place where code and specification differ, beyond the three build notes.
2. Coordinate conventions: `load_working_image` returns a grey matrix indexed [x, y] with y down, `prepare_outline` flips y, and `process_image` maps the tip, base, and outline back to image coordinates for plotting, with and without a crop. Is any transposition, flip, or crop offset wrong? Verify on a demo image in `bifaceAnalysis/www` that `tip_plot` and `base_plot` fall on the outline at the ends the declared direction implies.
3. `standardize_outline`: are the tie rules, the half-open midline crossing rule, the insertion of the base point, and the manual-landmark branch correct for every vertex arrangement, including a base point that falls on a vertex, an edge that wraps from the last vertex to the first, and a tip chosen among several vertices at the maximum? Does the manual branch produce the same left and right edge assignment as the automatic branch?
4. `efa_forward` and `efa_inverse`: are the coefficients and constant terms those of Kuhl and Giardina's formulas with an arc-length parameter and t = 0 at the first point, and is the inverse consistent with the forward transform? Is the comparison with Momocs in `tests/fixtures/make_reference.R` (first vertex moved to the end) a valid comparison or one that would hide a phase error?
5. `gpa`: does it converge to the same alignment as a standard generalized Procrustes fit, and is the final rotation of the whole set applied consistently to the specimens and the mean? Could the reflection guard in `rotate_onto` change a shape?
6. `pca_stability`, `mean_silhouette`, `rv_coefficient`, `region_rv`: are the formulas correct, including the separability rule, singleton clusters, constant blocks, and restoration of the caller's random number state?
7. Scale: is the conversion between working and original pixels correct in `app.R` (`ppm_to_original`, `ppm_to_working`), and is the scale stored with a saved specimen in the units that `measure_outline` divides by? Can a scale set on one image be applied to another image without the user asking for it, or persist after "Clear scale"?
8. State in `app.R`: can any sequence of user actions (changing image source, uploading new files with the same names, replacing a specimen, editing names, loading a CSV, changing analysis settings) alter a saved specimen unintentionally, attach a specimen to the wrong image, leave stale crop or landmark state, or raise an unhandled error? Check behaviour at 0, 1, 2, 3, and 4 specimens and when all shapes are identical.
9. `specimens_from_csv` and `merge_specimens`: can a file that passes validation produce a specimen that later causes an error, and does the merge rule match specification section 6.5?
10. Interpretation text in `interpret.R` and the Help text in `app.R`: is every statement true of the quantity the code computes? Identify any sentence a student would read as supporting a conclusion the analysis does not support, and any statement that uses the language of statistical significance.
11. Tests: for each test in `tests/testthat/`, would it fail if the code it covers were wrong? Identify tests that compare a function with itself, tests whose tolerance is too loose to detect an error, and requirements in specification section 10 that have no test. Mutate the code where practical and report which mutations the tests miss.
12. Are the three build notes accurate, and is the weakened disc test in `test-image.R` still sufficient to detect a wrong artifact selection?
13. Run `Rscript tests/testthat.R` if the sandbox allows it and report the result.

## Note from the author

Rewrite of the Shiny app against spec revision 2; build notes at the end of the spec list three known differences

## Milestone focus

This is an implementation reported as complete and not yet run on real data. Audit whether the code does what its specification and docstrings say, whether the tests would fail if the code were wrong, whether any parameter, seed, path, or unit is inconsistent between files, and whether running it on real data would produce a result the author would misread.

## Rules

- Read every file in scope in full from the working tree. Do not rely on summaries the authors wrote about their own work when the underlying file is in scope.
- Treat every claim as unverified until you have checked it against the files or the tools.
- Classify each finding as wrong, unsupported, or plausible but not demonstrated (not_demonstrated). Give a concrete repair and an acceptance check that would show the repair holds.
- List what you checked and found correct under held. List what you did not examine, with the reason, under not_examined. List the questions you would have asked the authors under questions.
- Put scripts and logs under evidence/2026-11-17_build/, anything derived from project data (extracts, copies, page images, large tables) under evidence/2026-11-17_build/data/, and your prose report at audits/2026-11-17_build_report.md, all relative to your working directory. Files under data/ are not committed, so every value your report or a finding cites must appear in the report itself or in a tracked log. List every file you wrote under evidence_files.
- Record `git rev-parse HEAD` in the report.
- If the sandbox blocks a command, record it under not_examined rather than working around it.
- Do not report style preferences. Do not pad the list. If you find nothing, say so and return an empty findings list.
- Respond only with JSON matching the provided schema.
