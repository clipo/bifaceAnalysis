# Codex Audit Ledger: bifaceAnalysis, build (2026-11-17_build_2)

**Date:** 2026-11-17
**Prompt:** codex/audits/2026-11-17_build_2_prompt.md
**Findings file:** codex/audits/2026-11-17_build_2_findings.json
**Report:** codex/audits/2026-11-17_build_2_report.md
**Scope:** docs/specs/2026-11-17_collections_spec.md; bifaceAnalysis/R/collections.R; bifaceAnalysis/R/plots.R; bifaceAnalysis/R/interpret.R; bifaceAnalysis/R/store.R; bifaceAnalysis/R/shape.R; bifaceAnalysis/R/efa.R; bifaceAnalysis/R/report.R; bifaceAnalysis/app.R; tests/testthat.R; tests/testthat/helper-app.R; tests/testthat/test-collections.R; tests/testthat/test-server.R; tests/manual_checklist.md; bifaceAnalysis/README.md
**Audited commit:** beedd0415ebbfad8c2e2c7c2f39c99c652c40290
**Model and effort:** gpt-6-astra, high. Codex codex-cli 0.169.0. Session 01a1ec83-c12b-7840-8bc5-2213ef4b5db9. Effective sandbox: workspace-write [workdir, $TMPDIR].
**Run:** 2026-11-17T21:36:59-05:00 to 2026-11-17T21:47:02-05:00, exit status ok. No file under `codex/` outside the run's stem changed. The runner warned that the report links to files under `evidence/2026-11-17_build_2/data/`, which are not committed. The values on which the findings rest (the variances and ratio of F1, the pair-label collision of F3, the page counts) are stated in the report text, so the findings can be checked without those files.

## Summary

| Verdict | Count |
|---|---|
| Confirmed, fixed and verified in tree | 7 |
| Confirmed, deferred with owner | 0 |
| Confirmed, declined by the author | 0 |
| Refuted | 0 |
| Undetermined | 0 |

Blocking findings open: 0. Milestone proceeds: deployment to the test address and the collections browser checklist (C1 to C11). Replacement of the live app follows the author's instruction of 2026-11-17 (automatic once the checklist passes), after the author's sign-off on this ledger.

How the checks were made. The author confirmed all seven findings for checking on 2026-11-17. Each was reproduced by reading the cited lines and, where numeric, by a script written for this ledger. Repairs are in commit 66d2573. Every repair has a test. Twelve deliberately wrong versions of the repaired code were then run against the tests: seven against `tests/testthat/test-collections.R` and five against `tests/testthat/test-server.R`, with the unmutated files passing first. All twelve are detected, including the mutation Codex reported as passing (every pair-reference value replaced by zero). One of them, removal of the geographic caution, was at first not detected, because the test compared the report text with the caution list itself. The test now checks the wording of four cautions, and the mutant is detected. Final state: `Rscript tests/testthat.R` completes with no failure in eight test files, with `maps` 3.4.4 installed.

Two process errors occurred during these checks and are recorded here. A first mutation run of the original build reported 11 of 17 mutants undetected, because the harness saved backups as `.bak` files inside `bifaceAnalysis/R/` and the test helper sourced every file in that folder. The helper now sources only `.R` files, and the rerun detected all 17. During the repair mutations, a newly written server test failed on unmutated code. The run was stopped, which left `interpret.R` mutated. The file was restored from the saved copy, and each repaired line was confirmed present by reading the file before the server mutants were rerun.

## Findings

### F1. Identical specimens receive a rounding-residue variance and an absurd ratio

| Field | Value |
|---|---|
| Severity (Codex) | blocking |
| Classification (Codex) | wrong |
| Location | bifaceAnalysis/R/collections.R `row_variance`; bifaceAnalysis/R/interpret.R `text_collections` |
| Codex claim | Collections of identical specimens get variances of about 1e-34 and 1e-66 instead of 0, and the text reports a ratio of about 6.3e31 to 1. |
| Check performed | Three identical non-integer rows (a random row multiplied by pi) gave `row_variance` = 7.4e-32 while `no_variation` returned TRUE. The function did not apply the app's numerical-zero rule (main spec build note i). The text chose the ratio sentence with `all(v > 0)`, so residues produced a ratio. |
| Acceptance check performed | Run. Duplicate non-integer rows now give variance 0 and leave-one-out range 0 to 0 (unit test). Through the app, two collections of three copies each give variances exactly 0, the all-zero sentence, and no ratio; one constant and one varying collection give the mixed sentence (server test). The ordinary formula and the Q hand examples still pass. Mutant removing the rule: detected. |
| Verdict | Confirmed |
| Resolution | Fixed. `row_variance` returns 0 when `no_variation` holds. Q, the summaries, and the text use the same value. |
| Tree evidence | `bifaceAnalysis/R/collections.R`, `row_variance`; tests "duplicate non-integer rows have zero variance..." and "collections of identical specimens have zero variance...". |
| Resolved on | 2026-11-17 |

### F2. Random geography assignments were described as all possible

| Field | Value |
|---|---|
| Severity (Codex) | major |
| Classification (Codex) | wrong |
| Location | bifaceAnalysis/R/interpret.R `text_collection_geo` |
| Codex claim | The text says "all possible" whenever there are at most 5040 assignments, which includes the 999 random ones used for 8 or more collections. |
| Check performed | Read: the condition was `geo$n_assignments <= 5040`, and 999 is less than 5040. |
| Acceptance check performed | Run. `shape_geo_correlation` returns `method` ("complete" or "random") and `n_located`. With 7 collections: 5040 assignments, values equal to an independent enumeration by Heap's algorithm, text "all possible". With 8: 999 assignments equal to an independent seeded draw, random state restored, text "chosen at random" and "8 located collections give 28 pairs". Mutants labelling the random branch complete, deciding by count, and not restoring the random state: all detected. |
| Verdict | Confirmed |
| Resolution | Fixed. |
| Tree evidence | `bifaceAnalysis/R/collections.R`, `shape_geo_correlation`; `bifaceAnalysis/R/interpret.R`, `text_collection_geo`; test "seven locations use all 5040 assignments and eight use 999 random ones". |
| Resolved on | 2026-11-17 |

### F3. Distinct pairs could share one row of the pair figure

| Field | Value |
|---|---|
| Severity (Codex) | major |
| Classification (Codex) | wrong |
| Location | bifaceAnalysis/R/plots.R `reference_plot` |
| Codex claim | Rows were keyed by label, so the pairs (A, B and C) and (A and B, C) both became "A and B and C". |
| Check performed | Read: `factor(df$item, levels = rev(unique(df$item)))` merges equal labels. |
| Acceptance check performed | Run. Collections named "A", "B and C", "A and B", and "C" give 6 pairs, 6 distinct rows, 6 distinct labels (coinciding labels are numbered), and each row's observed value equals its own raw distance. Mutant keying by label: detected. |
| Verdict | Confirmed |
| Resolution | Fixed. Rows are keyed by position, and labels are applied separately. |
| Tree evidence | `bifaceAnalysis/R/plots.R`, `reference_plot`; test "pairs whose labels coincide keep separate rows with their own values". |
| Resolved on | 2026-11-17 |

### F4. The tests did not establish several section-10 requirements

| Field | Value |
|---|---|
| Severity (Codex) | major |
| Classification (Codex) | not demonstrated |
| Location | tests/testthat/test-collections.R; tests/testthat/test-server.R |
| Codex claim | Replacing every pair-reference value with zero passed the unit tests; the random geography branch, rank transitions, report contents, raw CSV cells, and trimming at save and edit were not tested. |
| Check performed | Read: the only pair-reference check was an upper bound on one quantile; the selector test never removed a component; the report test compared page counts only; the save and edit tests used labels without spaces. |
| Acceptance check performed | Run. Added: an independent calculation of every pair-reference mean and limit; random-state restoration for pair and geography references, with and without an existing seed; the 7- and 8-location branches; the label-collision plot; a rank transition through the app with fallback of the map selector; zero-variance collections through the app; report sections with and without locations and the tick box, and caution wording; raw distance CSV lines with empty cells and a header-only file; trimming at save and edit; map plots for 2, 3, 6, and 9 collections; a clipped window and a one-axis span drawn as computed. The zero-reference mutant is detected. |
| Verdict | Confirmed |
| Resolution | Fixed. |
| Tree evidence | Tests after the heading "repairs after the build audit" in `tests/testthat/test-collections.R`, and the last three tests in `tests/testthat/test-server.R`; `report_content` in `bifaceAnalysis/app.R`, separated from the download handler so the sections can be tested. |
| Resolved on | 2026-11-17 |

### F5. Required cautions disappeared in some states

| Field | Value |
|---|---|
| Severity (Codex) | minor |
| Classification (Codex) | wrong |
| Location | bifaceAnalysis/R/interpret.R |
| Codex claim | The geographic cautions were missing from the report without a map section, and the variance caution appeared only when every variance was positive. |
| Check performed | Read: the geographic caution was only in the map text, and the variance caution was inside `if (all(v > 0))`. |
| Acceptance check performed | Run. The geographic caution is in `COLLECTION_CAUTIONS`, which the report includes whenever the comparison is available. The report holds it with no locations and with the map box unticked. The variance caution appears in all three variance cases. Mutants removing either caution: detected. |
| Verdict | Confirmed |
| Resolution | Fixed. |
| Tree evidence | `bifaceAnalysis/R/interpret.R`, `COLLECTION_CAUTIONS` and `text_collections`. |
| Resolved on | 2026-11-17 |

### F6. The pair table omitted the corrected square

| Field | Value |
|---|---|
| Severity (Codex) | minor |
| Classification (Codex) | wrong |
| Location | bifaceAnalysis/app.R `coll_pair_table` |
| Codex claim | The on-screen table lacked Q, which spec section 5.3 lists. |
| Check performed | Read: no `corrected_square` column. The CSV already had it. |
| Acceptance check performed | Run. The rendered table contains "Corrected square" (server test). With identical collections, the corrected square equals the squared raw distance. Mutant removing the column: detected. |
| Verdict | Confirmed |
| Resolution | Fixed. |
| Tree evidence | `bifaceAnalysis/app.R`, `output$coll_pair_table`. |
| Resolved on | 2026-11-17 |

### F7. The unavailable tab showed more than the explanation

| Field | Value |
|---|---|
| Severity (Codex) | minor |
| Classification (Codex) | wrong |
| Location | bifaceAnalysis/app.R, Collections tab |
| Codex claim | With fewer than 2 eligible collections the tab showed empty cards and the Locations card as well as the explanation. |
| Check performed | Read: the cards were built unconditionally. |
| Acceptance check performed | Partly run. The comparison cards are inside a `conditionalPanel` on `coll_available`, which is TRUE with 3 eligible collections and FALSE after labels are removed (server test). Mutant forcing TRUE: detected. What the browser shows is checklist item C1, not yet run. |
| Verdict | Confirmed |
| Resolution | Fixed in code. The browser check is part of C1. |
| Tree evidence | `bifaceAnalysis/app.R`, `conditionalPanel("output.coll_available", ...)` and `coll_available`. |
| Resolved on | 2026-11-17 |

## Held by the auditor

- The full test suite passed with `maps` 3.4.4.
- Independent partition values (SST 160.222, SSB 77.556, SSW 82.667, share 0.4840) and the seeded share reference matched.
- Independent Q values (32.333, 29.333, -11.667), raw distances, divisors, clipping, pair-only pooling, and reference summaries matched.
- Haversine distances matched an independent vector calculation. The five-collection correlation (0.3478) and its 120 reference values, and an eight-collection random reference, matched.
- References preserved existing and absent random states on every path checked.
- Full-data alignment, H, coefficients, and PCA feed the eligible rows correctly, including interleaved rows and special-character labels.
- The locations rules, text parsing, replacement, rematching, template, clearing, CSV schemas, empty cells, header-only files, and absence of coordinates matched the contract.
- Map geometry, backgrounds, and the scale-bar rule matched. Build note (a) is correct: 1000 km along the parallel at 60 degrees is 996.916 km by haversine.
- Degenerate states produced figures or messages without errors (F1 was an error in text only).
- The map tick box removed the map section (24 pages with, 23 without). Nine-level palettes were distinct. The performance fixture took 0.437 s.

## Not examined by the auditor

- The browser checklist C1 to C11, deployment, and real collections. Open: C1 to C11 are run by Claude on the test address next.
- The author's mutation run, which is not in the repository. Closed: this ledger records the mutation runs made for the repairs.
- External literature, the README's archival date for Momocs, and unbiasedness after fitted alignment and harmonic selection. Open: the two CITE-CHECK flags remain, and the unbiasedness limitation is stated in the spec and the app.
- Visual inspection of PDF pages, and two removed-row warnings in an auxiliary report run that Codex did not trace. Open until the report download is checked in the browser checklist.
- The project-root `CLAUDE.md` is absent. Closed.
- An optional report-validation step blocked by the sandbox. Closed: the report contents are now covered by the server test of `report_content`.

## Questions and answers

1. Who will run and record checklist C1 to C11 on the test deployment after the repairs? Answered by the author, 2026-11-17: Claude, in Chrome, recording the date, browser, results, and limits.

No answer led to a follow-up round.

## Sign-off

Pending: the author's sign-off on this ledger.
