# Build audit, 2026-11-17

Audited working tree HEAD: `5b0ec11c1b7af780c9d37e79b6f2037332c96d4b`, matching the brief. The specification, all eight R modules, the complete app, helper and six test files, fixture generator, manual checklist, and README were read. The shared CLAUDE.md was read. The project CLAUDE.md does not exist at the stated location or in HEAD. The prior specification resolution ledger was also read. Its decisions about geometric landmarks, analysis-only symmetrization, full CSV persistence, and the fixed-alignment RV reference are retained.

The build is not ready for acceptance. `Rscript tests/testthat.R` completed successfully, but targeted checks found accepted inputs that break analysis, state that transfers to unrelated images, a calibration that survives Clear scale, PDF content lost at realistic specimen counts, and a spurious PCA component for three identical specimens. Five deliberately incorrect implementations pass their relevant test files. The required browser checklist remains incomplete.

## Findings

### F1. CSV validation accepts specimens that fail downstream

**Major, wrong.** `bifaceAnalysis/R/store.R:139-158` checks numeric type and missing coordinates, then row count and base-row range. It does not require finite coordinates or valid derived configurations. The numeric test is `!is.numeric(df$x) || !is.numeric(df$y) || anyNA(df$x) || anyNA(df$y)`.

`checks.R` imports an outline with an infinite coordinate successfully. Measurements and analysis then fail with `infinite or missing values in 'x'`. Finite coordinates of order 1e200 also pass and overflow downstream. A finite ten-vertex contour with coincident tip and base positions at different rows passes, then GPA/Fourier/PCA fails with `missing value where TRUE/FALSE needed`. A collinear outline passes and reports `widest_position = NaN`. These are accepted store records, not rejection messages from the upload handler.

Repair: validate finite and numerically usable coordinates, nonzero area and centroid size, positive edge lengths after duplicate handling, distinct geometric tip/base, row numbering, nonempty ids, and consistent per-specimen metadata. Validate all derived points before committing a loaded batch. Acceptance: these files are rejected with a useful message and leave the store unchanged; accepted records build every applicable measurement, figure, and export.

### F2. Filename-only state attaches old processing to new images

**Major, wrong.** `bifaceAnalysis/app.R:261-315` resets only `rv$idx` when the file list changes, returns `f$name[...]` as identity, and looks up state and saved specimens using that name. It does not distinguish demo files, uploads, or replacement uploads with the same name. `store.R:45-52` then constructs potentially descending index sequences from stale crop bounds.

The server check changed from demo artifact1 to an upload named `artifact1.jpg` whose actual contents were artifact4. Its width changed to 134, but the old crop `0,240,70,170` remained and the full state compared identical. Another replacement upload retained the same state. A stale crop `160,240,70,170` on the narrower image raises `subscript out of bounds`. Saving also addresses the existing specimen under that filename, so unrelated content can replace it. Threshold, polarity, and speck controls are global, although section 6.5 assigns preview processing state to the image. The name/group refresh and edit-field refresh depend on id changes, not replacement or renaming under an unchanged id.

Repair: distinguish image content/source from display filename, specify deliberate replacement behavior, reset or restore preview state on actual identity changes, and refresh editable fields when the selected record changes. Validate intersected crop bounds before indexing. Acceptance: switching sources, uploading changed contents with repeated names, replacing a selected record by CSV, and returning after renaming cannot inherit another image's crop, clicks, calibration, or stale record fields. Out-of-image crops return a message.

### F3. Clear scale leaves an effective calibration

**Major, wrong.** `bifaceAnalysis/app.R:350-373` clears only `state()$scale`. Every finite `own_ppm()` writes `rv$shared_ppm`, even while sharing is unchecked. `current_ppm()` falls back to that retained value when sharing is checked.

Server reproduction: before Clear scale, own/current/shared were `10,10,10`. Afterwards they were `NA,10,10`. A specimen saved after clearing therefore still receives millimetre measurements. Sharing unchecked does correctly prevent ordinary cross-image fallback; the conversion arithmetic itself is correct.

Repair: define and implement clearing local and shared calibration so the effective scale is removed as the control promises. If local-only clearing is desired, expose that distinction and a separate shared-clear control. Acceptance: after the applicable clear action `current_ppm()` is unavailable, scale text says so, and newly saved records have no millimetre scale until explicitly calibrated or opted into an existing shared scale.

### F4. PDF pages silently discard text as the store grows

**Major, wrong.** `bifaceAnalysis/R/report.R:3-24` assigns text height from line count within a fixed ten-inch layout. `app.R:705-718` puts all specimen processing records on one cover page. No overflow handling or pagination exists.

A cover with forty ordinary processing records and start/end markers produced 124 wrapped lines requiring 20.66 inches. Text extraction from the resulting PDF contains only 21 specimen records, numbers 009 through 029. Both markers are absent. The file is a valid PDF, so the current header/page-count tests do not detect this loss.

Repair: paginate cover records and interpretation paragraphs within measured page capacity. Acceptance: at forty and at the supported maximum specimen count, all settings, first/last records, and complete interpretation text are visible and extractable, without overlap or clipping.

### F5. Tests miss materially incorrect implementations

**Major, not demonstrated.** `tests/testthat/test-shape.R:44-53,89-100`, `test-store.R:41-53`, `test-outline.R:101-108`, and `test-server.R:135-197` do not establish their broad correctness claims.

Mutations applied only in an auditor-owned R environment produced these results:

| Mutation | Test file | Failures | Errors |
|---|---|---:|---:|
| Return plane agreement 1 for every leave-one-out fit | test-shape.R | 0 | 0 |
| Never report PC2 as separable | test-shape.R | 0 | 0 |
| Replace every RV reference summary with observed RV | test-shape.R | 0 | 0 |
| Ignore loaded metadata when the contour is identical | test-store.R | 0 | 0 |
| Reverse traversal only in the manual-landmark branch | test-outline.R | 0 | 0 |
| Select the smaller disc in the synthetic selection test | test-image.R | 3 | 0 |
| Make symmetrization a no-op | test-shape.R | 4 | 0 |

The RV repeat-call check compares the function with itself. The merge check `merge_specimens(st, st)` supplies no changed metadata, so retaining the old record passes. The manual test checks snapped endpoints, not edge assignment. The server test invokes Clear scale without asserting its effect. See the per-test assessment below for additional gaps and the limited meaning of reference/self comparisons.

Repair: add independent expected values and adversarial state transitions, including non-unit plane agreement, positive PC2 separability, a known permutation distribution, changed metadata on identical contours, and manual/automatic correspondence equality. Acceptance: all five surviving mutants fail for their intended defect, with the normal implementation passing.

### F6. Interpretation makes conclusions the computations do not establish

**Major, wrong.** `bifaceAnalysis/R/interpret.R:134` calls RV zero evidence that region shapes are “unrelated”. At line 139 the reference is “for sample size only”. Neither statement follows from the calculation. For `x = -2:2` and `y = x^2`, RV is zero despite a deterministic relationship. At the same n=4, the exact mean over all 24 permutations is 1/3 for a rank-one block paired with itself, and 1 for an isotropic centered block. The reference depends on the fixed coordinate structure as well as sample size.

Other literal overstatements: line 22 and `app.R:215` place point 30 at 30 percent, whereas `point_positions()[30]` is 0.29. Line 51 asserts every outline is an elongated oval without testing that claim. Line 90 describes agreement near 1 as an unchanged result, though agreement measures loading directions/planes, not unchanged scores or variances. Line 107 says a horizontal dashed cut gives the chosen k even when tied merge heights make that impossible, including identical shapes. Line 116 promises predominantly one label per cluster whenever labels correspond to shape, independent of the requested k or number of labels.

Repair: describe RV as linear coordinate covariation and its reference as conditional on fixed blocks and n. Restrict stability language to directions/planes, correct indexing, remove unmeasured geometric explanations and guaranteed label purity, and qualify or omit the horizontal cut with tied heights. Acceptance: the quadratic and fixed-n reference examples cannot be read as independence or a sample-size-only effect; tied-height trees and point indices are described correctly. No explicit statistical-significance or p-value claim was found in student-facing text.

### F7. Competing objects in the opposite polarity are not warned about

**Minor, wrong.** `bifaceAnalysis/R/image.R:99-104` computes `n_other` from the winning class's filtered table, rather than the candidate set searched across classes. Section 6.1 requires a warning when another candidate is at least one quarter of the selected area.

The independent synthetic case has a selected interior dark component of 5,929 pixels and an interior light component of 5,041 pixels. It returns no warning. This is independent of the ambiguity about whether border-rejected components count as candidates.

Repair: count competitors across all eligible classes after selecting the winner. Acceptance: this opposite-polarity case warns, while the supplied demo and synthetic-disc behavior remains correct.

### F8. Base insertion mishandles a crossing at the second endpoint

**Minor, wrong.** `bifaceAnalysis/R/outline.R:179-184` says the half-open rule permits the base to coincide with the first vertex “but never with its second”. For an edge with decreasing x it can coincide with the second. The code tests only the first endpoint and duplicates the second endpoint in the saved contour.

`edge_checks.R` supplies a simple counterclockwise polygon with vertices `(-sqrt(11),1), (-sqrt(11),-1), (1,-1), (1,-3), (0,-4), (2,-3), (2,3), (0,4), (1,3), (1,1)`. Its area centroid is (0,0). Standardization inserts the base (0,-4) twice consecutively. Ordinary basal-notch/manual comparisons pass, and the current resampler tolerates the duplicate, so this reproduction does not establish wrong 200-point coordinates. It establishes an invalid claim and a zero-length edge introduced into the returned contour. A tangent vertex can also count twice under the specified half-open inequalities; “a vertex ... counted once” is not universally true.

Repair: detect both endpoints, reuse an existing base vertex, and document whether tangent contacts are counted or coalesced. Acceptance: first-endpoint, second-endpoint, closing-edge, and tangent cases contain one base vertex and no zero-length contour segments, with consistent left/right paths and cyclic-start behavior.

### F9. Invalid manual landmarks discard the previous valid landmarks

**Minor, wrong.** `bifaceAnalysis/app.R:339-341,387-395` replaces state with the new clicks before validating them. Section 6.2 requires a same-vertex error while keeping previous landmarks.

After a valid manual pair, two new clicks snapping to the same vertex change `processed()$ok` from TRUE to FALSE with the same-point error. The previous valid pair is not retained. The unit test checks that the pure function throws, which does not test server rollback.

Repair: retain the last valid pair and treat replacement clicks as pending until validation succeeds. Acceptance: the error is displayed, but the prior valid preview and ability to save remain intact.

### F10. RV is unavailable for some nonconstant blocks

**Minor, wrong.** `bifaceAnalysis/R/shape.R:175` rejects `den < 1e-24`. This is an absolute magnitude threshold, not the specified zero-variance rule, and breaks scale invariance.

For `x = matrix(1:20,4)`, RV(x,x)=1; multiplying both blocks by 1e-8 returns NA. Repair: normalize centered blocks before forming Gram matrices and detect zero variation robustly without an arbitrary unit-dependent cutoff. Acceptance: self-RV and cross-RV are unchanged across safe positive rescalings; truly constant blocks remain NA.

### F11. Required interface information is missing

**Minor, not demonstrated.** Spec sections 6.3, 7.2, 7.6, and 8 are not fully implemented by the Outline, region-control, and Help UI in `app.R:122-131,171-173,215-233`, `plots.R` function `plot_harmonic_power`, and `interpret.R`.

There is no mask preview, only the image/contour and standardized outline. Only cumulative power excluding harmonic 1 is shown as a series; including-first cumulative power is computed but not reported as required. The Help text does not state that pooling edges fails to remove arbitrary-face variation or that adequacy of fixed 100-point edges has not been tested. Region sliders permit equal boundaries, although the specification requires b1 < b2 < b3; sorting alone produces empty regions.

Repair: supply the required views and limitations and either enforce distinct boundaries or explicitly revise the empty-region contract. Acceptance: all required quantities and cautions are visible in the app, and equal-boundary inputs have a defined, tested result.

### F12. The build's manual acceptance condition is unmet

**Major, not demonstrated.** `tests/manual_checklist.md:7-11` says “Status: NOT YET COMPLETED”; spec section 10 requires completion on the demo images. Build note (c) accurately discloses the omission. The stronger checklist claim that every interaction is covered server-side is not accurate: uploads/source replacement, real client input updates, invalid-landmark rollback, and scale clearing are not established by its server tests.

Repair: complete and date the browser checklist after fixes, including actual upload, pointer/crop behavior, field synchronization, PDF inspection, and narrow layouts. Acceptance: browser/version and observed results are recorded for all items, with no unresolved errors. The dependency scan returned 55 CRAN-labelled packages without Momocs, but its CRAN index request failed, so current availability was not independently established.

### F13. Three identical specimens produce a spurious PCA component

**Major, wrong.** `bifaceAnalysis/R/shape.R:77-85` uses `total > 0` after `prcomp` centering, then only a threshold relative to that total. Floating-point centering of three identical noninteger coefficient rows leaves a tiny common residual that is treated as variance. `app.R:546-550,600-606,624,675-679` then permits a scree and shape-change plot and an exported PC1 instead of the identical-shape message.

Three copies of the same processed demo specimen have maximum coefficient-row difference exactly 0 but return rank 1 and retained variance 1.508674e-36. Morphospace says `These specimens have 1` component. The scree and shape-change figures build and the shape-change interpretation appears. Four identical copies return rank 0, which is the only identical-shape sample count covered by the existing server test.

Repair: recognize exactly identical rows before PCA and use stable centering, with a documented numerical-zero policy for nearly identical inputs. Acceptance: three, four, and other counts of identical noninteger configurations report no variation, produce no PCA components/shape-change figures, and export no spurious scores; genuinely small nonzero variation remains analyzable according to the specified rank rule.

## Checks that held

The test command exited 0. All six test files completed without reported test failures or skips (`tests.log`). The process emitted the environment warning `OMP: Warning #179: Function Can't set size of /tmp file failed`; it did not prevent the tests.

The displayed demo images were visually inspected in `data/overlay.png`. Their masks follow the artifacts and exclude the top scale bars. Coordinates use [x,y] with y down until preparation, and the inverse landmark transform correctly restores image coordinates. For artifact1, tip=(5,97.5), base=(226,113.9229); artifact2, tip=(11,105.5), base=(264,123.8687); artifact3, tip=(13,113.5), base=(184.3328,124.8328); artifact4, tip=(5,123.5), base=(119,116.4917). All tips are at the declared left end. Marker distances to contour segments were at most 5.684342e-14 pixels. With a fixed threshold, artifact1's crop and full-image markers agree exactly. Recomputing Otsu after cropping changes the mask slightly, explaining marker changes of 0.04812449, 0.01621129, 0.1352301, and 1 pixel in the four automatic-threshold comparisons. Those changes are not offset errors.

Opening, 4-connectivity, area-based selection, hole filling, padded contour extraction, ordinary axis selection, tip tie-breaking examples, landmark snapping, asymmetric resampling, width profiles and shoelace measurements pass their specified examples. The weakened disc test still rejects selecting the smaller object: its mutation causes three failures. Note (a) is accurate: the independent run omits exactly four disc cells, includes zero outside cells, and leaves zero hole cells unfilled. Note (b) correctly distinguishes distance along a polyline from straight chord distance; its independent projection check is meaningful.

The Fourier implementation follows the supplied arc-length equations, including the closing edge, t=0 at the first vertex, constants, coefficient ordering, and power. Independent midpoint quadrature over 1,000 subdivisions per segment differs from the coefficients by at most 3.145939e-11. Against installed Momocs 1.4.1, rolling the input as specified gives maximum error 5.854692e-17; leaving it unrolled gives 0.00258335. The roll aligns different phase conventions; it does not hide an implementation phase error. Circle, ellipse, inverse-at-vertices, power and constant-term tests pass.

GPA agrees with the geomorph fixture after one global rotation to maximum coordinate error 1.275032e-06, below 1e-5. The final rotation is applied to both specimens and mean. The reflection guard uses an orthogonal determinant-one matrix; on a reflected target it changes the optimal orientation, not shape. Interpoint distances changed by at most 5.551115e-17. Symmetrization, PCA centering/unscaled fitting, rank gating, variance shares and the implemented full-data separability rule agree with the spec. The maximum-iteration cap is a finite implementation safeguard, not evidence that convergence always occurs on arbitrary input.

Silhouette formulas agree with cluster for ordinary partitions, singletons, duplicates and zero distances. RV's trace numerator/denominator are correct apart from F10. Region assignments, variance formulas and the fixed-coordinate permutation procedure agree with the spec. RNG restoration works both with an existing seed (project test) and when no seed existed (independent check). Seed 1 and 999 permutations are the defaults.

Scale algebra is correct: original ppm=working ppm/f, working ppm=original ppm*f. Saved records use working ppm and measurements divide by those units. No ordinary shared fallback occurs while the sharing checkbox is false. Analysis settings leave saved contours, measurements and scale unchanged; symmetrization off restores the prior PCA scores. The supplied server tests cover 0,1,2,3,4 specimens, four identical shapes, normal downloads, replacement, removal and k clamping. Four-copy identical-shape PCA and RV return the expected availability messages. Independent checks of every figure and interpretation reactive at 1..4 identical specimens raised no unhandled errors, but found incorrect PCA availability at n=3 (F13).

CSV round trips preserve the normal fixture records and settings. The merge branching implements the intended same-id/same-contour replacement and different-contour new-id rule on ordinary data, including unique display names. Its `all.equal(..., tolerance=1e-9)` is approximate and normally relative, not literal contour identity. A precise identity tolerance remains an author specification detail. Cached `points` are kept in the in-memory record despite section 6.5 saying points are not stored separately, but are derived when constructing/importing a specimen, never exported, and no inconsistency from that cache was demonstrated.

Other literal precision differences: axis cosines and tip distances are rounded to nine decimals before tie-breaking (`outline.R:57,167`), and widths before choosing the widest level (`outline.R:225`). Thus near-ties can be treated as exact ties. These deterministic numeric policies should be stated if retained. They are not grounds to reopen the accepted geometric-landmark decision.

## Per-test assessment

These assessments cover every `test_that` in scope. Passing means the stated examples passed, not that every implementation error is detected. Except the explicit mutation results, sensitivity judgments below are from reading assertions.

| File and test start line | What the test establishes and what it can miss |
|---|---|
| image 1 | Two-level and constant Otsu cases. Any threshold in the gap passes; a generally wrong threshold function could satisfy these cases. |
| image 7 | Opening removes small blocks and preserves a larger one, zero disables. No even-sized element or shifted-origin oracle. |
| image 15 | Strong large-disc selection and hole-fill check despite four-cell allowance. Smaller-disc mutant fails. |
| image 33 | Border fallback warning and empty forced-light mask. Does not test cross-class competitor warnings. |
| image 41 | All demo masks satisfy border, hole, scale-bar and 0.99 reference-IoU conditions. Masks are frozen outputs from the same pipeline, so this is regression evidence; visual inspection is the independent selection check. |
| image 60 | Tiny crop rejection and normal process success. The nontrivial crop only asserts character warnings, not correct selection/coordinates. |
| outline 1 | Rotation/start-shift and densification invariance on one asymmetric polygon. Baseline is the same standardizer, so a consistently wrong invariant result could pass. |
| outline 12 | Independent expected tip/base for kite. |
| outline 19 | Correct major-axis choice in the exact 45-degree fixture. Tolerance 1e-6 is looser than the 1e-8 rotation requirement, but a wrong perpendicular-axis choice is much larger and would fail. |
| outline 25 | Declared-direction behavior for a wide polygon. |
| outline 31 | Near-circle warning exists, not precise fallback axis/points. |
| outline 37 | Base at the notch and no crossing warning for this ordinary notch. Does not exercise decreasing-x second endpoint. |
| outline 43 | Lowest-y base and warning for the barb fixture. |
| outline 49 | Blunt-tip warning, tied tip choice and no demo warnings. |
| outline 56 | Distinct samples, endpoints, positive segments, independent along-edge spacing, positions and mirror pairing. Endpoint equality also compares fields produced by one function; independent path projection supplies the stronger check. |
| outline 94 | Straight spacing approximately constant for a smooth ellipse, tolerance 1e-3 relative. This is supplementary; it is not the exact polyline-spacing oracle. |
| outline 101 | Snap two points and reject coincident snaps. Reversed manual traversal survives. The first assigned manual result is unused. No server rollback assertion. |
| outline 110 | Independent rectangle, ellipse, triangle, notch and scale measurements. Strong basic oracles. |
| outline 125 | Correct unit transfer at resize factors 1 and 0.5 for objects expected to be 30 and 60 mm. It does not execute two-click UI calibration or the exact same-image/same-physical-length experiment specified in section 10. Its tolerance .02 uses testthat's relative comparison, appropriate for the stated two-percent target. |
| efa 1 | Independent committed Momocs coefficients/order, absolute tolerance1e-8. |
| efa 8 | Analytic sampled-circle first and higher coefficients. |
| efa 15 | Independent Momocs ellipse and rounded a1/a3 values. |
| efa 21 | Forward/inverse consistency at original vertex parameters, not an independent inverse oracle. Constant and forward tests constrain compensating errors. |
| efa 31 | Independent coefficient grouping and power calculation against Momocs matrix. |
| efa 41 | H bounds, minimal qualifying H, exact zero residual vector and lower bound. Does not exercise floating near-zero residual spectra broadly. |
| efa 55 | Independent nonzero arc-length constant for triangle and repeated-point rejection. |
| shape 1 | Independent geomorph fit after one common rotation, size and mean orientation. Equality to REF$aligned is additionally same-implementation regression evidence. |
| shape 15 | Independent Momocs matrix and base prcomp scores. Taking absolute values cellwise is weaker than checking one sign per component, but detects ordinary numerical differences. |
| shape 23 | Rank at n=3, exact identical and line data; null stability when unavailable. |
| shape 34 | Symmetry and unit size; original examples asymmetric. No-op mutant fails four assertions. |
| shape 44 | Unit plane agreement under order exchange, inseparability, PC1 separability. Always-one planes and always-false PC2 separability survive. |
| shape 55 | Independent cluster silhouette oracle, singletons, duplicates, zero distances and k range. Does not independently test Ward D2 tree construction. |
| shape 71 | Exact partition counts and boundary membership, including unsorted boundaries. No duplicate-boundary contract. |
| shape 78 | Independent RV trace identity, self-RV, constant block and scale factor2. Does not catch the absolute small-scale cutoff. |
| shape 89 | RNG restoration, repeatability, bounds, pair count and one point variance. Reference distribution has no independent oracle; observed-as-reference mutant survives. |
| store 10 | Complete normal round trip and measurement consistency. Comparing measurements computed by the same function before/after establishes persistence, not measurement truth. |
| store 28 | Required column/version/text-x/short rejection. No finite/degenerate-geometry checks. |
| store 41 | Different-contour collision and unique names. Same-store remerge cannot distinguish replacement from ignoring loaded metadata. |
| server 16 | Every demo figure builds, text is nonempty, PDF header and at least10 pages, parsed CSV dimensions/basic columns. No text truth, PDF overflow, or exact outline-row-count assertion. |
| server 43 | Store immutability, reversible symmetrization, H/k controls and sorted boundaries. |
| server 62 | Save/replace/advance, explicit rename/group/remove. Inputs are supplied programmatically; client update synchronization is not tested. |
| server 83 | Figures/text/downloads at 0..3, key count gates, no PC3 at n=3. `is_validation` also accepts silent `req()` errors as a “message”, without asserting a visible explanatory sentence. |
| server 109 | Four identical demo copies give PCA/RV availability messages and a PDF. Other identical-shape outputs are not individually checked here. |
| server 122 | Normal CSV loads into empty running store and PCA builds. No duplicate-name CSV or replacement-metadata assertion. |
| server 135 | Happy-path crop, scale, shared fallback, saved units, landmarks and detection controls. Clear scale has no effect assertion; small crop ignored. No same-filename upload/source transition, failed-landmark rollback, or browser pointer check. |

No numerical tolerance inspected is so large that it would plausibly conceal a gross axis, phase or unit error in the independent numerical comparisons. Most weaknesses are missing expected behavior, self-comparison, or narrow fixtures. Default `expect_equal` tolerances are relative unless `scale=1` is supplied, despite the spec's absolute-tolerance default. Most numerical fixture tests explicitly use absolute scale where essential. The 0.99 mask IoU and inverse 0.005 centroid-size tolerance are intentional specification tolerances, not auditor-selected relaxations.

Section 10 is broadly represented, but the exact same-image resize/calibration experiment, two different same-name specimens through CSV loading, exact per-vertex CSV row-count check, genuinely visible small-n messages, and the manual acceptance run are absent or weaker than specified. Outside section 10's enumerated examples, the first/second/wrapping base endpoints, nontrivial stability planes, PC2 separability, independent permutation reference, upload identity, clear-scale behavior, malformed numeric geometry and large reports need regression coverage.

## Limits, questions, and records

No browser session or actual pointer interaction was run. The incomplete manual checklist remains open. No real archaeological upload collection, deployment, command-line legacy pipeline, root CSV/PDF outputs, exhaustive malformed-file fuzzing, or performance/maximum-size run was examined. The PDF stress test examined text preservation, not all screen/page layouts. Kuhl and Giardina and Escoufier originals were not supplied or opened; equations were checked against the supplied mathematical definitions, numerical integration, and installed/reference implementations, not attributed to independently inspected papers. Point-count adequacy and anatomical correspondence remain the explicitly deferred research questions.

`rsconnect::appDependencies("../bifaceAnalysis")` returned 55 packages labelled CRAN, with no Momocs or geomorph, but reported `cannot open URL 'https://cloud.r-project.org/src/contrib/PACKAGES'`. Current remote availability was not verified and the failed network request was not worked around. A shell command attempting a Python heredoc edit of the auditor's mutation harness was blocked with `can't create temp file for here document: operation not permitted`; that shell command was not retried. Normal auditor-owned file edits used the provided patch tool. No git command changed the working tree, index, or history.

Questions for the authors: should Clear scale clear the shared calibration globally or explicitly opt the current image out? Should replacement uploads with the same filename deliberately replace a saved identity, or receive new ids? What exact numerical tolerance defines “identical” contours for merge? Are equal region boundaries meant to be allowed? What browser/platform will be used to complete manual acceptance? These decisions do not excuse the demonstrated failure states.

Versions: R4.6.2 (2026-10-29), shiny1.14.0, bslib0.12.0, imager1.0.8, ggplot2 4.0.3, testthat3.3.2, cluster2.1.8.3, rsconnect1.11.2, Momocs1.4.1, geomorph4.1.2; pdftotext/Poppler26.11.0. Reference fixture metadata reports Momocs1.4.1, geomorph4.1.2, imager1.0.8 and R4.6.2. Momocs was used for an independent audit comparison, not by the app or its tests. The inherited shared conventions and prior ledger were contextual records, not proof of implementation correctness.

Evidence: `evidence/2026-11-17_build/{checks,edge_checks,state_checks,mutations,report_checks}.R` and matching `.log` files; `tests.log`, `versions.log`, `dependencies.log`. Derived files are `data/overlay.png`, `data/degenerate.csv`, `data/large_report.pdf`, and `data/large_report.txt`. Every numeric value relied on above is printed in this report or a tracked log. The temporary CSV is overwritten across adversarial cases; `checks.R` and its log define all cases reproducibly. No historical audit report was edited.
