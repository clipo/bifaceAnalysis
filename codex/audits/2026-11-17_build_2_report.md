# Collections build audit

Date: 2026-11-17. Mode: build. Stem: `2026-11-17_build_2`.

Audited working-tree HEAD: `beedd0415ebbfad8c2e2c7c2f39c99c652c40290`, matching the requested commit. The scoped implementation files had no working-tree differences from HEAD. No project file, git index, or history was changed.

The main statistical formulas and reactive alignment are correct. The full test suite passes. Acceptance is blocked by F1: exactly constant collections can produce an enormous reported variance ratio from floating-point residue. The other findings concern an incorrect description of random assignments, colliding pair labels, incomplete acceptance tests, and omitted interface content.

## Findings

### F1. Constant collections produce a spurious variance ratio

Severity: blocking. Classification: wrong.

Location: `bifaceAnalysis/R/collections.R:82`, `bifaceAnalysis/R/interpret.R:188-192`.

The implementation computes `sum(sweep(X, 2, colMeans(X))^2) / (nrow(X) - 1L)` and chooses the ratio sentence with `all(v > 0)`. It does not use the existing numerical-zero handling. This violates the zero-variance cases in spec 5.4.

Reproduction through the actual server: two collections containing three exact copies each of two synthetic elliptical outlines yield PCA rank 1. Within each collection the maximum absolute difference from its first coefficient row is exactly 0, and `no_variation(X)` is TRUE. Nevertheless, A's variance is `4.5623235653455283e-66` and B's is `2.8888949165808538e-34`. The generated sentence reports “a ratio of 63320693396766192788713188098048.0 to 1.” With six copies of one shape, rank 0, both variances are `2.82118644197349e-37` and the text reports a ratio of 1 instead of no variation. These values and sentences are in `evidence/2026-11-17_build_2/degenerate.log`.

Repair: make constant/numerically constant inputs return zero variance under the app's documented numerical-zero policy, including leave-one-out subsets, and use that result consistently in the text, summaries, and pair correction.

Acceptance: test non-integer coefficient rows and actual duplicate outlines through the server. Both constant collections must have zero variance and zero leave-one-out ranges, with the all-zero sentence and no ratio. A mixed constant/varying fixture must use the mixed-zero sentence. Retain the ordinary-data formula and Q checks.

### F2. Random geography assignments are called complete enumeration

Severity: major. Classification: wrong.

Location: `bifaceAnalysis/R/interpret.R:201-203`; compare `bifaceAnalysis/R/collections.R:193`.

The text chooses “all possible” when `geo$n_assignments <= 5040`. The random branch produces 999 assignments, which also satisfies that condition. Both an independent 8-collection fixture and an actual 9-collection server fixture report “999 assignments, all possible.” There are 40320 possible collection permutations for eight collections. The underlying random reference values are correct.

Repair: return an explicit enumeration/random method field with the computed reference and use it in the tab and report text.

Acceptance: 5 and 7 located collections say all possible with 120 and 5040 assignments; 8 and 9 say chosen at random with 999. Check both rendered text and PDF content. Evidence: `numeric.log` and `server.log`.

### F3. Distinct pairs can share one plotted row

Severity: major. Classification: wrong.

Location: `bifaceAnalysis/R/plots.R:325,344-346`.

Pair identities are constructed with `paste(collection_a, "and", collection_b)`, then converted to a factor using `unique(df$item)`. Valid collection names A, B and C, A and B, and C give six pairs but only five plotted y positions. Pairs (A, B and C) and (A and B, C) both become “A and B and C”; their observations and reference ranges overlap on one row. This is not merely overlapping text: the data share a categorical position. The CSV retains the distinct pairs.

Repair: use a unique pair identifier as the plotted category, separate from the display label. Quote or otherwise disambiguate the names in the displayed label.

Acceptance: this four-name fixture has six distinct y positions and identifiable labels, with each observation and reference matched to its own pair. Evidence: `numeric.log`.

### F4. The acceptance tests do not establish several required results

Severity: major. Classification: not_demonstrated.

Location: `tests/testthat/test-collections.R:97-105,165-192,232-266`; `tests/testthat/test-server.R:428-464`; spec section 10.

The only numeric assertion on a pair reference is `expect_lte(ab$reference_high, max(stats::dist(pool)))`. This is a loose bound on one quantile, not the specified check of every reassigned distance against the pair's attainable distances. An in-memory replacement of all pair-reference means, lower limits, and upper limits with zero passed every test in `test-collections.R`. The mutation was run by this audit; the author's unavailable mutation run was not used. Evidence: `check_numeric.R` and the BEGIN/END MUTATION block of `numeric.log`.

Other untested or incompletely tested contracts are listed below. Most notably, no repository geography test exercises the 999-random branch, and the test named “the map selector follows the PCA rank” never removes a component or checks fallback. The report test checks only that unticking the box reduces page count, not which sections/text disappeared or the expected count.

Repair: add independent numeric reference oracles and the missing branch, state, and output checks. Do not use the implementation's own helper as the sole expected result.

Acceptance: the zero-reference mutation must fail; random-reference, RNG, rank-transition, exact-zero, plot-key, report-content, and missing-value regressions must fail their corresponding tests. The corrected suite must then pass with maps installed.

### F5. Required cautions disappear in supported states

Severity: minor. Classification: wrong.

Location: `bifaceAnalysis/R/interpret.R:148-155,189-204`; `bifaceAnalysis/app.R:973-981`; spec sections 5.4 and 8.

`COLLECTION_CAUTIONS` omits the standing statement that map colours are sample means and rounded/approximate coordinates change distances. Those statements occur only in conditional map/geography text. In the two-collection report fixture, turning the map off produces a 23-page report with neither statement; no geography figure is available to carry the coordinate statement. With no locations the same omissions remain. Section 8 requires the standing cautions in full in the report and tab. Separately, the small-sample and leave-one-out sensitivity cautions are inside `all(v > 0)` and disappear in the mixed-zero case, although a varying collection and its bar remain.

Repair: put the standing statements in the unconditional caution vector, and keep the variance precision/sensitivity explanation outside the ratio's three-way branch.

Acceptance: verify all required cautions in tab text and extracted PDF text with no locations, a withheld map, the map box off, and positive/mixed/all-zero variances. Evidence: `server.log`; PDFs and text extracts under this stem's `data/` directory.

### F6. The displayed pair table omits Q

Severity: minor. Classification: wrong.

Location: `bifaceAnalysis/app.R:855-863`; spec section 5.3, Outputs.

The specified pair table includes the corrected square. The rendered table includes raw distance, the reference, and corrected distance, but no `corrected_square`. The CSV does include it. A student cannot inspect the signed quantity that is clipped to obtain a zero distance in the tab.

Repair: add a clearly labelled corrected-square column to the pair table.

Acceptance: a pair with negative Q shows the negative value and corrected distance 0 in both the table and CSV.

### F7. The unavailable tab is not restricted to its explanation

Severity: minor. Classification: wrong.

Location: `bifaceAnalysis/app.R:213-236,846-852`; spec section 7.

Only `coll_intro` switches to `COLLECTIONS_HELP`. All plot containers, the pair-table card, Locations controls, and the full caution card remain unconditionally in the UI. The figure reactives each emit the same explanation. Section 7 says that fewer than two eligible collections show only the short explanation; the build has no note authorizing this deviation.

Repair: condition the comparison content on eligibility, leaving one explanation when unavailable.

Acceptance: browser or DOM checks with zero, one eligible, and two eligible collections confirm a single explanation in the first two states and the complete comparison in the third.

## Held

- Independent synthetic fixture: SST `160.2222222222222`, SSB `77.55555555555554`, SSW `82.66666666666666`, share `0.4840499306518723`. All agree with the implementation. The 999-draw share reference was independently regenerated at seed 1: mean `0.2467877030983826`, lower `0.0429958391123436`, upper `0.5880721220527045`; analytical expected share `0.25`.
- On the same fixture, independent Q values for AB, AC, BC are `32.33333333333333`, `29.33333333333334`, and `-11.66666666666667`. Corrected distances are `5.686240703077327`, `5.416025603090641`, and 0. Divisors are n minus 1 for variance and n for each mean's correction. Raw distances, within-pair pooling, reference summaries, and clipping agree. The first pair's reference was regenerated independently and checked against all 20 three-of-six assignments; largest attainable raw distance `6.009252125773315`.
- Ordinary-data leave-one-out ranges, column means, summary columns, all-eligible scale gating, and inverse-transform linearity match the specified calculations. The defect in F1 concerns numerical constant inputs, not a changed variance divisor.
- Great-circle distances agree with an independent 3D unit-vector/atan2 calculation. The lower-triangle correlation for five collections is `0.347822001268706`; all 120 reference values agree with independent enumeration. Seven collections return 5040 assignments. Eight return 999 random assignments whose values agree with an independently seeded calculation; observed r `0.7079475323054896`. Coincident positions retain permutation multiplicity in the repository test.
- Share, pair, complete-geography (5 and 7), and random-geography (8) references restore a present RNG state, preserve its absence, and repeat identically. Early returns for zero share, too few locations, coincident positions, and zero shape distances also preserve it. The seed helper restores state after an intentional error.
- The app computes alignment, symmetry, H, coefficients, and PCA on all specimens. Eligible row indices select coefficients, full-data scores, and measurements in the same order. An interleaved fixture with nine eligible collections, three unlabelled rows, and one singleton verified every exported summary value directly; H was 7 and PCA rank 7. Labels NA, 007, commas, quotes, angle brackets, ampersands, and percent signs survived this check. The pair-label ambiguity is separately identified in F3.
- Locations are read as text without missing-value coercion, normalized headers are checked for missing/repeated required names, and all specified rejection categories are implemented. Acceptance replaces the table only after parsing succeeds. Rematching, unmatched names, all-nonempty-label templates, clearing, and rejection retention agree with the contract. The upload limit remains 40 MB. No coordinates enter either new result CSV.
- Both result CSVs have the specified column order, empty missing cells, and header-only unavailable output. Rank 0 and rank 1 summaries leave absent PC columns empty. One missing eligible scale withholds length-mm means and that map choice. The selector falls back from unavailable PC3 to length/width in the additional server checks.
- Map windows agree with independent padding/minimum-span/clipping calculations for coincident, latitude-only, longitude-only, high-latitude, clipped, US, non-US, and straddling fixtures. The unfiltered polygon/background code and built limits agree with the rule. Scale-bar candidate selection, parallel length, and placement agree. Build note (a) is correct: a 1000 km parallel segment at 60 degrees has haversine endpoint distance `996.9156797055513` km, a `0.3084320294448717` percent difference. A spherical constant with more precision than the spec's printed 111.19508 is used consistently.
- All seven new figures build for the nine-collection server fixture; repository tests also cover six non-map figures at 2, 3, 6, and 9 collections. Complete lack of matching locations, coincident locations, rank 0/1, identical collection means, and a withheld map produced messages or figures rather than new reactive/download errors in the additional checks. Numerical zero interpretation still fails as F1 describes.
- The map checkbox actually removes the map section: 24 pages with the map, 23 without, and 23 without locations for the report fixture. Text continues on additional pages. Build note (b)'s separate cautions section is present. No examined sentence calls the clipped distance unbiased or calls a reassignment range a population interval/test/evidence of a real difference.
- The source diff shows one label-trimming line in the outline reader and the shared palette replacement in the existing plotting code. Up to eight levels retain their original colours; nine distinct colour/fill levels pass. The same helper also changes cluster colours above eight. No numerical changes to existing analysis paths were found, and the existing-label invariance server test passes. Existing server assertions were adjusted for the newly added collection outputs, not weakened for the existing numerical results.
- The 120 by 120 performance fixture took `0.437` seconds for share, pair, and within statistics on this audit machine, below two seconds. Map quantity changes do not invalidate these statistics. The exact author timings in build note (c) were not independently certified.

## Test coverage and execution

`Rscript tests/testthat.R` ran from the repository root, with TMPDIR placed under this audit's writable tree. Exit status 0. Collections, efa, image, outline, report, server, shape, and store tests completed with no reported failures or skips. Maps 3.4.4 was installed, so map tests ran. Full output: `evidence/2026-11-17_build_2/test-suite.log`.

The supplied hand expectations, haversine constants, enumeration counts, and simulation targets are correct. The tests exercise every main statistic at least once, but do not meet every section-10 acceptance requirement. In addition to F4:

- No pair-reference mean/lower/upper oracle, pair-reference RNG check, geography RNG check, absent-seed check, or geography random-branch test is present.
- The five-location enumeration is tested, but the seven-location branch boundary and complete enumeration at seven are not independently checked by repository tests.
- Rank 0/1 are tested in the summary helper. The existing identical-outline server test stops at five unlabelled specimens, so it never exercises an eligible zero-variance comparison. No eligible rank transition, zero-variance collection text through the analysis pipeline, or missing-scale map-selector transition is tested.
- No direct server oracle verifies all collection CSV numeric values against shuffled/interleaved full-data rows. Unit tests check summary values, while server tests mostly check names, sizes, presence, and one geographic distance.
- No raw CSV assertion distinguishes empty missing cells from other missing-value encodings across both files. The unavailable distances download is not checked. Server rejection retention is exercised only for a missing-column file, although unit tests exercise each rejection category.
- Save/edit server tests use g1 and g2, without surrounding spaces. They do not establish trimming at those entry points. The source does trim them.
- Map tests lack complete built-plot coverage of single-axis spans, clipping, high latitude, and straddling. Several expected results call the implementation's own window/length helpers. Map plots are not included in the 2/3/6/9-collection figure loop. Pair y-position identity, plotted numeric values, full report headings/cautions, exact page counts, and every unavailable-state message are not tested.

The audit's additional scripts fill several of these checks for this commit, but they do not replace maintained regression tests. The all-zero pair-reference mutation passing the collections unit file is direct evidence of a weak oracle, not an inference from test names.

## Prior ledger reconciliation

Read `docs/audits/codex_audit_2026-11-17_spec_2_ledger.md`. Its nine resolutions were specification repairs, not proof that the implementation meets them. Latitude/meridian withholding, trimming, replacement/rematching, shared full-data analysis, fixed-data references, Q versus clipped distance, and the operational five-location threshold were treated as settled decisions. F1 is new implementation evidence against the zero-variance contract associated with prior F3. F4 identifies incomplete implementation of the stronger testing contract associated with prior F8. No prior report or ledger was edited and no author decision is reopened.

## Limits and not examined

- The project-root `CLAUDE.md` is absent (file read and root git tree lookup); the global `~/.claude/CLAUDE.md` was read.
- Browser layout and manual checklist C1 to C11 were not run. The checklist says “Status: NOT YET RUN.” No deployment or real-collection analysis was performed. Browser acceptance therefore remains unestablished independently of F1.
- No use was made of the author's non-repository mutation run. Only the explicitly logged in-memory mutation was run by this audit.
- No external literature, Momocs archival-date claim, real-data validity, or fitted-representation unbiasedness claim was independently established. The last remains an acknowledged design limitation, not a new finding.
- PDFs were generated and their text/page counts checked; pages were not visually inspected for overlap. An auxiliary full-report run emitted two removed-row drawing warnings, recorded in `server.log`. They were not traced to an existing plot. Rendering the new non-map collection figures for the three degenerate fixtures separately emitted no such warnings (`degenerate.log`). No claim of warning-free full-report rendering is made.
- A final Python here-document command intended to check the report for em dashes and verify that every listed evidence file exists was blocked before Python ran: `zsh:1: can't create temp file for here document: operation not permitted`. That optional report check was not rerun or worked around. See `environment-check.log`. R also printed `OMP: Warning #179: Function Can't set size of /tmp file failed:` at startup; the suite and independent scripts still completed successfully. That diagnostic was left as recorded.
- Files outside the brief were not substantively audited. Runtime dependencies/fixtures were executed by the requested test suite and the app. Specified app, plot, interpret, and checklist sections were read in full; collections, tests, README, and the needed shape/EFA/report implementations were read in full.

## Environment

R 4.6.2 (2026-10-29); testthat 3.3.2; shiny 1.14.0; bslib 0.12.0; ggplot2 4.0.3; maps 3.4.4; scales 1.4.0; imager 1.0.8; cluster 2.1.8.3. Base statistics, utilities, graphics, and grid are from this R installation. PDF text extraction used pdftotext/Poppler 26.11.0. Package versions are recorded in `numeric.log`.

## Question for the authors

Who will run and record Collections checklist C1 to C11 on the test deployment after the repairs?

## Evidence files written

All paths are relative to `codex/`:

- `audits/2026-11-17_build_2_report.md`
- `evidence/2026-11-17_build_2/check_numeric.R`
- `evidence/2026-11-17_build_2/numeric.log`
- `evidence/2026-11-17_build_2/check_server.R`
- `evidence/2026-11-17_build_2/server.log`
- `evidence/2026-11-17_build_2/check_degenerate.R`
- `evidence/2026-11-17_build_2/degenerate.log`
- `evidence/2026-11-17_build_2/test-suite.log`
- `evidence/2026-11-17_build_2/environment-check.log`
- `evidence/2026-11-17_build_2/data/commit-diff.txt`
- `evidence/2026-11-17_build_2/data/with-map.pdf`
- `evidence/2026-11-17_build_2/data/with-map.txt`
- `evidence/2026-11-17_build_2/data/without-map.pdf`
- `evidence/2026-11-17_build_2/data/without-map.txt`
- `evidence/2026-11-17_build_2/data/no-locations.pdf`
- `evidence/2026-11-17_build_2/data/no-locations.txt`

Disposable R/test temporary files were confined to `evidence/2026-11-17_build_2/tmp/`; durable cited outputs are listed above. Every value used in a finding appears in this report or a tracked log, not solely in an ignored data file.
