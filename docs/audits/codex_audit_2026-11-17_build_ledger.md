# Codex Audit Ledger: bifaceAnalysis, build (2026-11-17_build)

**Date:** 2026-11-17
**Prompt:** codex/audits/2026-11-17_build_prompt.md
**Findings file:** codex/audits/2026-11-17_build_findings.json
**Report:** codex/audits/2026-11-17_build_report.md
**Scope:** docs/specs/2026-11-17_app_rewrite_spec.md; bifaceAnalysis/app.R; bifaceAnalysis/R; tests/testthat.R; tests/testthat; tests/fixtures/make_reference.R; tests/fixtures/reference.rds; tests/manual_checklist.md; bifaceAnalysis/README.md; bifaceAnalysis/www
**Audited commit:** 5b0ec11c1b7af780c9d37e79b6f2037332c96d4b
**Model and effort:** gpt-6-astra, high. Codex codex-cli 0.169.0. Session 01a1eb50-0bac-7ed0-92b7-7a69c0d6f649. Effective sandbox: workspace-write [workdir, $TMPDIR].
**Run:** 2026-11-17T16:00:53-05:00 to 2026-11-17T16:18:20-05:00, exit status ok. No file under `codex/` outside the run's stem changed.

## Summary

| Verdict | Count |
|---|---|
| Confirmed, fixed and verified in tree | 12 |
| Confirmed, deferred with owner | 1 (F12, Carl Lipo) |
| Confirmed, declined by the author | 0 |
| Refuted | 0 |
| Undetermined | 0 |

F11 is counted as fixed. One part of it (global detection controls) was resolved by a clarifying spec note, not a code change.

Blocking findings open: 1 (F12, the manual browser checklist). Milestone proceeds: no, until the checklist is completed or the author accepts the build without it.

How the checks were made. Each finding was reproduced by a script written for this ledger (not Codex's scripts) before any repair, repaired, and then re-run. Tests were added for each repair. Seventeen deliberately wrong versions of the code (the five that Codex reported as passing, and twelve more aimed at the repairs) were run against the test files. Sixteen are detected. The seventeenth (removing the finite-coordinate check in `specimens_from_csv`) is not detected because `outline_problem` rejects the same input at a second place with a message. Final state: `Rscript tests/testthat.R` completes with no failure in seven test files.

## Findings

### F1. CSV validation accepts specimens that fail downstream

| Field | Value |
|---|---|
| Severity (Codex) | major |
| Classification (Codex) | wrong |
| Location | bifaceAnalysis/R/store.R, `specimens_from_csv` |
| Codex claim | Outlines with infinite coordinates, coincident tip and base, or no area pass validation and cause errors later. |
| Check performed | Wrote CSVs with an infinite x, a ten-vertex contour whose base row repeats the tip position, and a collinear contour. All three were accepted before the repair. Measuring the first then failed with "infinite or missing values in 'x'". |
| Acceptance check performed | Run. Eleven malformed files (infinite, 1e200, repeated vertex, row gap, clockwise, collinear, tip at base, tip below base, empty id, negative scale, inconsistent metadata) are each rejected with a message. A rejected file leaves the store unchanged (server test). Accepted demo records give finite measurements and a rank 3 PCA. |
| Verdict | Confirmed |
| Resolution | Fixed. Added `outline_problem()` and checks of coordinates, row numbers, ids, scale, and per-specimen metadata. The whole file is validated before the store is changed. |
| Tree evidence | `tests/testthat/test-store.R`, test "outlines that could not be analysed are rejected before they reach the store"; `test-server.R`, "A rejected file leaves the store unchanged". Mutant "no outline validation in CSV": 5 failures. |
| Resolved on | 2026-11-17 |

### F2. Filename-only state attaches old processing to new images

| Field | Value |
|---|---|
| Severity (Codex) | major |
| Classification (Codex) | wrong |
| Location | bifaceAnalysis/app.R (image identity and state); bifaceAnalysis/R/store.R, `process_image` |
| Codex claim | State and saved specimens are keyed by file name, so an upload with a repeated name inherits crop and landmarks, can replace an unrelated specimen, and a stale crop raises "subscript out of bounds". |
| Check performed | `process_image(artifact4, crop = c(160, 240, 70, 170))` on the 134-pixel-wide image raised "subscript out of bounds" before the repair. Read `current_id` and confirmed it returned the file name. |
| Acceptance check performed | Run as a server test: after cropping and saving sample image artifact1, an upload of different contents named artifact1.jpg has id `upload1/artifact1.jpg`, empty state, and saves as a second specimen without changing the first. A second upload of the same file gets `upload2/...` and empty state. Returning to the sample images restores the sample's own crop. Crops outside the image return a message. Loading a CSV refreshes the fields of the selected record. |
| Verdict | Confirmed |
| Resolution | Fixed, following the author's decision of 2026-11-17 ("new identity"): each upload receives new ids. Crop bounds are checked before indexing. The name, group, and edit fields follow the saved record. The part of the finding about global threshold, polarity, and speck controls is resolved by spec build note (e), which states that these controls apply to the image shown. They do not alter saved specimens (held by the auditor). |
| Tree evidence | `test-server.R`, "an uploaded file never inherits the state or the specimen of an image with the same name"; `test-image.R`, "a crop outside the image gives a message". Mutants "upload id = file name" (9 failures) and "crop bounds unchecked" (1 error). |
| Resolved on | 2026-11-17 |

### F3. Clear scale leaves an effective calibration

| Field | Value |
|---|---|
| Severity (Codex) | major |
| Classification (Codex) | wrong |
| Location | bifaceAnalysis/app.R, scale handling |
| Codex claim | After "Clear scale", the shared scale still applies to the image when the shared-scale box is ticked. |
| Check performed | Read `current_ppm`: it fell back to `rv$shared_ppm` whenever the image had no scale of its own and the box was ticked, and clearing changed only the image's clicks. The existing server test pressed Clear scale without asserting the result. |
| Acceptance check performed | Run as a server test: after clearing, `current_ppm()` is NA, the scale text says so, and a specimen saved then has no millimetre values. |
| Verdict | Confirmed |
| Resolution | Fixed, following the author's decision of 2026-11-17 ("opt this image out"): clearing removes the image's scale and stops the shared scale from applying to it until a new scale is clicked on it. Other images still receive the shared scale. The button and checkbox labels state this. |
| Tree evidence | `test-server.R`, "Clear scale removes the effective scale of the image and leaves the shared scale for others". Mutant "clear scale keeps shared": 5 failures. |
| Resolved on | 2026-11-17 |

### F4. PDF pages discard text as the store grows

| Field | Value |
|---|---|
| Severity (Codex) | major |
| Classification (Codex) | wrong |
| Location | bifaceAnalysis/R/report.R |
| Codex claim | The cover page does not paginate. With 40 specimens only part of the processing records appears. |
| Check performed | Wrote a cover with 40 records and extracted the text with pdftotext: 1 page, 21 of 40 records before the repair. |
| Acceptance check performed | Run at the supported maximum: 500 records plus 30 interpretation paragraphs on a figure page. Every record, both end markers, and every paragraph are found by text extraction. Wrapped lines are no wider than the text area by measured string width, including a 400-character word. |
| Verdict | Confirmed |
| Resolution | Fixed. Text is wrapped to the measured page width and continued on "(continued)" pages. A figure page keeps at least 4.5 inches for the figure. |
| Tree evidence | `tests/testthat/test-report.R` (two tests). Mutant "report not paginated": 6 failures. The text-extraction part of the test is skipped where pdftotext is not installed. It ran here. |
| Resolved on | 2026-11-17 |

### F5. Tests miss incorrect implementations

| Field | Value |
|---|---|
| Severity (Codex) | major |
| Classification (Codex) | not demonstrated |
| Location | tests/testthat/test-shape.R, test-store.R, test-outline.R, test-server.R |
| Codex claim | Five wrong implementations pass the tests: plane agreement always 1, PC2 never separable, observed RV substituted for its reference, loaded metadata ignored for an identical contour, and reversed traversal in the manual-landmark branch. |
| Check performed | Read the cited tests and confirmed each gap (for example the RV reference test compared the function with itself). |
| Acceptance check performed | Run. Tests were added with expected values computed independently of the code: plane agreement from the normals of the planes in three dimensions; separability for three stated variance sets; the reference against the exact distribution over all 24 pairings (mean 1/3); changed metadata on an identical contour; manual landmarks at the automatic positions giving the same 200 points with the left edge first. All five mutants now fail: 8, 1, 3, 4, and 3 failures. |
| Verdict | Confirmed |
| Resolution | Fixed. Further gaps from the report's per-test table were also closed: two specimens with the same name through a CSV load, the row count of the outline file, visible messages at zero specimens, and the server tests named under F2, F3, and F9. One difference from spec section 10 remains by design and is recorded as build note (j): the scale test uses two images at the same original resolution, not one image at two resolutions. |
| Tree evidence | New tests in `test-shape.R`, `test-store.R`, `test-outline.R`, `test-server.R`; mutation results above. |
| Resolved on | 2026-11-17 |

### F6. Interpretation text states more than the computations establish

| Field | Value |
|---|---|
| Severity (Codex) | major |
| Classification (Codex) | wrong |
| Location | bifaceAnalysis/R/interpret.R; bifaceAnalysis/app.R Help text |
| Codex claim | RV of 0 is described as "unrelated", the reference as "for sample size only", point 30 as 30 percent, every outline as an elongated oval, agreement near 1 as "unchanged", the dashed cut as always giving k clusters, and labels as predicting cluster purity. |
| Check performed | Computed RV(x, x squared) = 0 for x = -2..2. Computed the exact mean reference at n = 4: 1/3 for a rank-one block and 1 for an isotropic block, so the reference depends on the blocks and not only on n. `point_positions()[30]` is 0.29. Read each sentence at the cited lines. |
| Acceptance check performed | Read back the revised sentences. `grep` for "unrelated", "sample size only", "elongated", "mostly one label", "was unchanged" in `interpret.R` and `app.R` returns nothing. The tied-height case is tested at 4 and 5 identical specimens (the text says no dashed line is drawn). |
| Verdict | Confirmed (all seven sentences) |
| Resolution | Fixed. RV is described as linear covariation of coordinates and 0 as not ruling out other relations. The reference is described as depending on the number of specimens and the pattern of variation within each region. The point example reads "30 steps from the tip". The oval explanation is removed. Agreement is described as applying to directions and planes, not scores or variance shares. The dashed line is drawn and described only when a single cut exists. The label sentence describes the table only. |
| Tree evidence | `bifaceAnalysis/R/interpret.R` (`text_outlines`, `text_harmonics`, `text_stability`, `text_clusters`, `text_regions`); `plots.R` `tree_cut_height`; `test-server.R` "identical shapes give no components at 2, 3, and 5 specimens". |
| Resolved on | 2026-11-17 |

### F7. Competing objects of the other brightness class give no warning

| Field | Value |
|---|---|
| Severity (Codex) | minor |
| Classification (Codex) | wrong |
| Location | bifaceAnalysis/R/image.R, `segment_image` |
| Codex claim | The count of other large objects uses only the selected class. |
| Check performed | Built an image with a dark square of 5,929 cells and a light square of 5,041 cells, each surrounded by the other class. The dark square was selected with no warning before the repair. |
| Acceptance check performed | Run: the case warns ("1 other object"); removing the competitor removes the warning; the demo images and the synthetic disc still give no warning. |
| Verdict | Confirmed |
| Resolution | Fixed. Competitors are counted across every searched class under the same border rule. A region inside the filled mask is a hole of the artifact and is not counted, which is tested. |
| Tree evidence | `test-image.R`, "a competing object of the other brightness class triggers the warning". Mutant: 1 failure. |
| Resolved on | 2026-11-17 |

### F8. Base insertion mishandles a crossing at an edge's second vertex

| Field | Value |
|---|---|
| Severity (Codex) | minor |
| Classification (Codex) | wrong |
| Location | bifaceAnalysis/R/outline.R, `vertical_crossings`, `standardize_outline` |
| Codex claim | The base can coincide with the second vertex of its edge, which duplicated that vertex. A vertex that touches the midline can be counted twice. |
| Check performed | Ran Codex's ten-vertex polygon: one zero-length segment in the returned contour before the repair. Checking the mirror image of the same polygon then showed a further defect that Codex did not report: the half-open rule counted the touching vertex twice from the right and not at all from the left, so the mirror image received a different base (y = -1 instead of -4). |
| Acceptance check performed | Run: second-vertex and first-vertex cases have one base vertex and no zero-length segment; the two mirror images give points that correspond under the mirror pairing within 1e-8; shifting the start vertex changes nothing; touching vertices are returned once from either side. The closing-edge case cannot occur in the automatic branch, because the contour starts at the tip and the base is the lowest point. |
| Verdict | Confirmed, and extended |
| Resolution | Fixed. The crossing rule now returns every vertex on the midline and every point where an edge passes through it. A base on a vertex reuses the vertex. This replaces the rule in spec section 6.2 and is recorded as build note (d). The author should note this change to a rule that the spec audit had reviewed. |
| Tree evidence | `test-outline.R`, "a base point on either vertex of its edge is not duplicated" and "the points where an outline meets a line are the same for a shape and its mirror image". Mutant "base vertex not reused": 5 failures. Demo outlines are unchanged (reference fixture tests pass at 1e-10). |
| Resolved on | 2026-11-17 |

### F9. Invalid manual landmarks discard the previous valid landmarks

| Field | Value |
|---|---|
| Severity (Codex) | minor |
| Classification (Codex) | wrong |
| Location | bifaceAnalysis/app.R, landmark clicks |
| Codex claim | Two replacement clicks that snap to the same vertex leave the preview unusable, where spec section 6.2 requires the previous landmarks to be kept. |
| Check performed | Read the click handler: the first replacement click overwrote the stored tip and cleared the base before any validation. |
| Acceptance check performed | Run as a server test: after a valid manual pair, two clicks near one vertex leave the preview valid, with the same tip and the same 200 points, and the specimen can be saved. |
| Verdict | Confirmed |
| Resolution | Fixed. The first click is held as a pending tip. The pair replaces the landmarks only when the two clicks snap to different outline points. Otherwise a notice is shown. |
| Tree evidence | `test-server.R`, "an invalid pair of landmark clicks keeps the previous landmarks". Mutant: 4 failures, 1 error. |
| Resolved on | 2026-11-17 |

### F10. RV is unavailable for some blocks that vary

| Field | Value |
|---|---|
| Severity (Codex) | minor |
| Classification (Codex) | wrong |
| Location | bifaceAnalysis/R/shape.R, `rv_coefficient` |
| Codex claim | The zero-variance guard is an absolute cutoff, so RV depends on units. |
| Check performed | `rv_coefficient(x, x)` for `matrix(1:20, 4)` was 1, and NA after multiplying by 1e-8, before the repair. |
| Acceptance check performed | Run: RV is unchanged within 1e-10 under factors 1e-8, 1e-3, and 1e6; constant blocks, including ones whose value is not exactly representable, give NA. |
| Verdict | Confirmed |
| Resolution | Fixed. Each centred block is divided by its root sum of squares. A block has no variation when that quantity is at most 1e-10 of the root sum of squares of the uncentred values (build note (i)). The aligned coordinates the app passes were never near the old cutoff, so no reported value changes. |
| Tree evidence | `test-shape.R`, "RV does not depend on the units of the blocks". Mutant: 3 failures. |
| Resolved on | 2026-11-17 |

### F11. Required interface information is missing

| Field | Value |
|---|---|
| Severity (Codex) | minor |
| Classification (Codex) | not demonstrated |
| Location | bifaceAnalysis/app.R (Outline page, region controls, Help); `plot_harmonic_power`; `interpret.R` |
| Codex claim | No mask preview; cumulative power including the first harmonic is not reported; two stated limitations are absent from Help; equal region boundaries are allowed. |
| Check performed | Compared spec sections 6.3, 7.2, 7.6, and 8 with the interface. All four gaps are present. |
| Acceptance check performed | Run where testable: `output$mask_plot` renders and the page contains the mask panel; equal boundaries give a message on the three region figures and leave other figures unaffected (server test). The two text additions were read back. |
| Verdict | Confirmed |
| Resolution | Fixed. Added a "Detected region (mask)" panel. The harmonics text reports cumulative power both excluding and including the first harmonic, labelled. Help states that pooling both edges does not remove variation from which face was photographed, and that the adequacy of 100 fixed points per edge has not been tested. Boundaries must differ (author decision, 2026-11-17, build note (h)). |
| Tree evidence | `app.R` (`mask_plot`, `boundaries_ok`, Help text); `interpret.R` `text_harmonics`; `test-server.R`, "equal region boundaries withhold the region results with a message". Mutant: 5 failures. |
| Resolved on | 2026-11-17 |

### F12. The manual acceptance condition is unmet

| Field | Value |
|---|---|
| Severity (Codex) | major |
| Classification (Codex) | not demonstrated |
| Location | tests/manual_checklist.md; spec section 10 |
| Codex claim | The browser checklist is not completed, the checklist overstated what the server tests cover, and CRAN availability of the dependencies was not established. |
| Check performed | The checklist status is "NOT YET COMPLETED". Three attempts to load the local app in Chrome during the build failed. The dependency check was re-run with network access: 55 dependencies, all present in the CRAN index of 25,774 packages retrieved on 2026-11-17. |
| Acceptance check performed | Partly. Dependencies: verified. Browser checklist: not run. |
| Verdict | Confirmed |
| Resolution | Deferred to Carl Lipo: run `shiny::runApp("bifaceAnalysis")` and complete `tests/manual_checklist.md` with date and browser. The checklist's statement of coverage was corrected, and items were added for the mask panel, scale clearing, landmark rollback, same-name upload, and equal boundaries. |
| Tree evidence | `tests/manual_checklist.md`, status paragraph and items 3, 8a, 9, 9a, 20a, 20b. |
| Resolved on | Open |

### F13. Three identical specimens produce a spurious component

| Field | Value |
|---|---|
| Severity (Codex) | major |
| Classification (Codex) | wrong |
| Location | bifaceAnalysis/R/shape.R, `shape_pca` |
| Codex claim | Centring identical non-integer rows leaves a residue that is counted as variance, so three identical shapes give rank 1. |
| Check performed | Identical copies of one demo specimen: rank 0 at n = 2 and 4, and rank 1 at n = 3 (variance 1.5e-36) and n = 5 (3.2e-37) before the repair. |
| Acceptance check performed | Run: rank 0 at n = 2 to 6 for identical rows; in the app at 3 and 5 identical specimens the PCA figures give the "identical in shape" message and the scores file holds no scores; rows that differ by one part in a million still give rank 2. |
| Verdict | Confirmed |
| Resolution | Fixed. `shape_pca` reports no components when the rows do not differ by the rule in build note (i). |
| Tree evidence | `test-shape.R`, "identical non-integer rows give no components at any count"; `test-server.R`, "identical shapes give no components at 2, 3, and 5 specimens". Mutant: 1 error. |
| Resolved on | 2026-11-17 |

## Held by the auditor

- HEAD matched the requested commit, all files in scope were read, and the test suite passed.
- Demo overlays were inspected. Tip and base markers lie on the contours within 5.69e-14 pixels and match the declared direction. Crop coordinates agree exactly at a fixed threshold.
- Independent Fourier quadrature agrees within 3.15e-11. Momocs agrees within 5.86e-17 with the rolled input and differs by 0.0026 without it, so the phase adjustment is valid.
- Alignment agrees with the geomorph fixture within 1.28e-6 after one rotation. The final rotation is consistent, and the reflection guard preserves shape.
- Resampling, measurements, symmetrization, silhouette, region partition, the RV trace formula, and random-state restoration pass, subject to the edge cases reported.
- Scale conversion and saved scale units are correct. With the box unticked there is no cross-image fallback. Analysis settings do not alter saved specimens.
- Normal CSV round trips, merge branching, sample-size guards, cluster-count clamping, and demo downloads work.
- Build notes (a) and (b) are accurate. The disc test still fails when the smaller disc is selected.
- No significance or p-value language was found in student-facing text.

The report also notes that tie-breaking rounds to nine decimals. A comment in `outline.R` (`axis_direction`) now states this policy.

## Not examined by the auditor

- Browser interaction, client-side field updates, layout, and deployment: open, see F12.
- CRAN availability: covered on 2026-11-17 with network access, all 55 dependencies present.
- One mutation edit blocked by the auditor's shell: no action needed. The mutants were re-run for this ledger.
- The Kuhl and Giardina and Escoufier papers: open. The PDFs are not in the repository. The CITE-CHECK flags in the spec remain.
- Real collections, exhaustive malformed-file fuzzing, performance at maximum size, and all PDF layouts: open. The report was tested at 500 records for text preservation only.
- The legacy command-line script and root outputs: out of scope by spec section 4.
- Project CLAUDE.md: absent from the repository. Closed.

## Questions and answers

1. Should Clear scale clear the shared calibration globally or opt the current image out? Author, 2026-11-17: opt this image out. Implemented under F3.
2. Should uploads with repeated file names replace existing identities or receive new ids? Author, 2026-11-17: new identity. Implemented under F2.
3. What tolerance defines identical contours for merging? Set to a maximum coordinate difference of 1e-9 pixels with the same base row and vertex count, which is the round-trip tolerance of spec section 10, and recorded as build note (f). Author, 2026-11-17: accepted.
4. Are equal region boundaries allowed? Author, 2026-11-17: no, require distinct values. Implemented under F11.
5. Which browser and platform will be used for manual acceptance? Open. The author's own browser. See F12.

No follow-up round was run.

## Sign-off

2026-11-17. The author (Carl Lipo) signed off on this ledger with F12 open and deferred to him. He accepted build notes (d) the midline rule, (e) global detection controls with per-image crop, scale, and landmarks, and (f) the merge tolerance. The repairs are in commit daa5cf1. The app is not to be used on real images or deployed until the manual browser checklist in `tests/manual_checklist.md` is completed. No repair-check audit of commit daa5cf1 has been run.
