# Codex Audit Ledger: bifaceAnalysis, spec (2026-11-17_spec_2)

**Date:** 2026-11-17
**Prompt:** codex/audits/2026-11-17_spec_2_prompt.md
**Findings file:** codex/audits/2026-11-17_spec_2_findings.json
**Report:** codex/audits/2026-11-17_spec_2_report.md
**Scope:** docs/specs/2026-11-17_collections_spec.md; docs/specs/2026-11-17_app_rewrite_spec.md; bifaceAnalysis/R/shape.R; bifaceAnalysis/R/efa.R; bifaceAnalysis/R/store.R; bifaceAnalysis/R/plots.R; bifaceAnalysis/R/interpret.R; bifaceAnalysis/R/report.R; bifaceAnalysis/app.R
**Audited commit:** 5481f0b2f346604195c2d4c6ebf070c6a158eaae (specification revision 1)
**Model and effort:** gpt-6-astra, high. Codex codex-cli 0.169.0. Session 01a1ec29-72d7-75f3-a441-e01da97c124f. Effective sandbox: workspace-write [workdir, $TMPDIR].
**Run:** 2026-11-17T19:58:21-05:00 to 2026-11-17T20:06:09-05:00, exit status ok. No file under `codex/` outside the run's stem changed.

## Summary

| Verdict | Count |
|---|---|
| Confirmed, fixed and verified in tree | 9 |
| Confirmed, deferred with owner | 0 |
| Confirmed, declined by the author | 0 |
| Refuted | 0 |
| Undetermined | 0 |

Blocking findings open: 0. Milestone proceeds: yes, to the build, once the author signs off on specification revision 2.

How the checks were made. The author confirmed all nine findings for checking on 2026-11-17 and answered three questions (map limits, label trimming, locations replacement). Each finding was checked by reading the cited code and specification lines and, where numeric, by an R script written for this ledger and not taken from Codex's evidence. The script ran the distance-correction simulation (seed 11, 4000 replicates, sizes 5 and 8, 12 variables), computed the pole distance, read a CSV with repeated headers and the names "NA" and "007", counted palette colours for 9 levels, and computed the PCA rank of three specimens of which two are identical. All repairs are to the specification. Nothing is built. The revised lines were read back from the tree with `grep` in a separate command after writing.

## Findings

### F1. The equator-to-pole test value is wrong

| Field | Value |
|---|---|
| Severity (Codex) | blocking |
| Classification (Codex) | wrong |
| Location | collections spec revision 1, section 10, Geography |
| Codex claim | With R = 6371.0088 km the equator-to-pole distance is 10007.557 km, not 10007.543 km, so the test would reject a correct implementation. |
| Check performed | Computed pi R / 2 = 10007.5572 km and pi R / 180 = 111.1951 km in R. The stated 10007.543 differs by 0.014 km, outside the stated tolerance of 0.001. |
| Acceptance check performed | Run for the two values above. The one-degree, identical-point, and meridian values agree with the formula. |
| Verdict | Confirmed |
| Resolution | Fixed. Section 10 gives 10007.557221 km and 111.195080 km, both within 1e-6, and requires that an implementation omitting the radian conversion fails. Section 5.6 states the conversion. |
| Tree evidence | Spec revision 2, lines 503 to 509 (Geography tests) and section 5.6 ("converted from degrees to radians"). |
| Resolved on | 2026-11-17 |

### F2. The map window, background, and scale bar are not defined for all accepted locations

| Field | Value |
|---|---|
| Severity (Codex) | major |
| Classification (Codex) | unsupported |
| Location | collections spec revision 1, section 6 |
| Codex claim | The minimum longitude span, aspect, and scale bar fail near the poles; the order of padding and minimum span is unstated; scale-bar placement is undefined; the boundary limitation is misstated. |
| Check performed | Computed the minimum longitude span 0.5 / cos(lat): 0.5 degrees at the equator, 2.88 at 80 degrees, 28.65 at 89 degrees, 286.48 at 89.9 degrees. Above about 89 degrees the rule produces a window wider than any useful map, and at 90 degrees it is undefined. Read revision 1: the order of steps, clipping, and bar position were not stated. "Outside the contiguous United States a window smaller than a country shows no boundary" omits the case of a small window inside the United States with no county line, so the wording was inexact. |
| Acceptance check performed | Not runnable: `maps` is not installed and nothing is built. The revised tests in section 10 (Map) are the acceptance check for the build. |
| Verdict | Confirmed |
| Resolution | Fixed, with the author's decision of 2026-11-17 (withhold beyond 80 degrees). Section 6 now withholds the map beyond 80 degrees and when longitudes span more than 180 degrees, gives the window in three ordered steps with clipping, requires unfiltered polygons with `coord_quickmap(..., expand = FALSE)`, defines the scale-bar length, position, and conversion, and rewords the limitation. |
| Tree evidence | Spec revision 2, section 6 ("Withholding", "Window, computed in this order", "Scale bar", "Limits"). |
| Resolved on | 2026-11-17 |

### F3. Map selectors and text are undefined when components or variances are missing

| Field | Value |
|---|---|
| Severity (Codex) | major |
| Classification (Codex) | unsupported |
| Location | collections spec revision 1, sections 5.5, 6, 8 |
| Codex claim | PC1 and PC2 may not exist, and a largest-to-smallest variance ratio is undefined when a variance is 0. |
| Check performed | Ran `shape_pca` on three specimens of which two are identical: rank 1, so PC2 does not exist. Read `shape_pca` (shape.R lines 90 to 113): components beyond the rank are dropped. Revision 1 made "mean PC1 score" the default and required a ratio with no zero case. |
| Acceptance check performed | Run for the rank case. The zero-variance case follows from the definition of v_g. |
| Verdict | Confirmed |
| Resolution | Fixed. Selectors list only existing components, the default is mean length divided by width, an unavailable selection returns to the default, PC summary columns are empty where absent, and section 5.4 gives text for positive, mixed-zero, and all-zero variances. |
| Tree evidence | Spec revision 2, sections 5.4 ("Ratio in the text"), 5.5, 6 (selector paragraph), 7 (summary CSV columns). |
| Resolved on | 2026-11-17 |

### F4. Labels are not trimmed everywhere, and the palette repeats beyond 8 levels

| Field | Value |
|---|---|
| Severity (Codex) | major |
| Classification (Codex) | wrong |
| Location | store.R line 221; app.R lines 586, 656, 750; plots.R `category_scale` |
| Codex claim | The existing app does not use the trimmed collection identity everywhere, so colouring and exports can disagree with the collections tab, and the palette rule does not reach every scale. |
| Check performed | Read app.R lines 509 and 567: labels are trimmed when saved in step 2 and edited in step 3. Read store.R line 221: `specimens_from_csv` passes `d$group[1]` to `new_specimen` untrimmed. `groups()` (app.R line 586) returns labels unchanged, and they feed colouring (line 656) and the PC score export (line 750). A loaded label " A " is therefore a separate group from "A". Counted `rep(OKABE_ITO, length.out = 9)`: 8 distinct colours. Revision 1 said the PCA tab "already colours the morphospace by collection", which is wrong for labels loaded with spaces. |
| Acceptance check performed | Run for the code reading and the colour count. The load test of "A", " A ", "  ", "", "NA", "007" is specified for the build. |
| Verdict | Confirmed |
| Resolution | Fixed, with the author's decision of 2026-11-17 (trim everywhere). Section 4.1 trims labels at every entry point, including the outline file reader, and keeps them as text. Section 3 permits this change and the palette change to existing tabs. Section 7 replaces `category_scale` with one function for colour and fill. |
| Tree evidence | Spec revision 2, sections 3, 4.1, 7 ("Colours"), 9 (store.R change), 10 (Labels). |
| Resolved on | 2026-11-17 |

### F5. Unlabelled specimens affect collection results by more routes than stated

| Field | Value |
|---|---|
| Severity (Codex) | major |
| Classification (Codex) | wrong |
| Location | collections spec revision 1, section 4.1; app.R lines 624 to 638 |
| Codex claim | An unlabelled specimen changes results through automatic H and the PCA as well as alignment, and the change need not be slight. |
| Check performed | Read app.R lines 628 to 637: `aligned`, `power`, `n_harm` (through `choose_harmonics(power())`), `coefs`, and `pca` are computed from all specimens. A change in automatic H changes the dimension of the coefficient space. Codex's evidence records H changing from 7 to 17 after one shape was added. That figure was not reproduced, but the dependency it illustrates is established by the code. |
| Acceptance check performed | Code path traced. The numeric example is Codex's and was not rerun. |
| Verdict | Confirmed |
| Resolution | Fixed. Section 4.1 names alignment, automatic H, and PCA as dependencies and removes "slightly". |
| Tree evidence | Spec revision 2, section 4.1, paragraph beginning "All specimens". |
| Resolved on | 2026-11-17 |

### F6. The corrected distance was described as unbiased when only the quantity under the root is

| Field | Value |
|---|---|
| Severity (Codex) | major |
| Classification (Codex) | wrong |
| Location | collections spec revision 1, section 5.3 |
| Codex claim | Revision 1 conflates the unbiased squared quantity with the clipped square root, understates assumptions, and does not say that a corrected 0 is not evidence of equal populations. |
| Check performed | Simulation (seed 11, 4000 replicates, sizes 5 and 8, 12 standard normal variables). Equal means: mean Q = -0.024 (Monte Carlo SE 0.027), mean max(0, Q) = 0.657, mean raw d^2 = 3.879 (expected 12/5 + 12/8 = 3.9), Q negative in 55 percent of samples. Mean difference of length 2: mean Q = 3.933 (SE 0.045), within 3 SE of 4; mean raw d^2 = 7.842. Revision 1 called dc "approximate" but called the clipped root's input unbiased without separating clipping. Its sentence "the quantity under the square root" was correct for Q and ambiguous for max(0, Q). |
| Acceptance check performed | Run for Q and max(0, Q) as above. |
| Verdict | Confirmed |
| Resolution | Fixed. Section 5.3 defines Q separately, states independence within and between collections, a fixed representation, and finite second moments, describes the upward effect of clipping, states that a corrected 0 does not show equal mean shapes, and states that fitted alignment and H are not accounted for. Section 10 tests Q, raw d^2, and max(0, Q) under a seeded design. |
| Tree evidence | Spec revision 2, section 5.3 ("Assumptions", "Two operations"), section 10 (Distances). |
| Resolved on | 2026-11-17 |

### F7. The locations file and download contracts are incomplete

| Field | Value |
|---|---|
| Severity (Codex) | major |
| Classification (Codex) | unsupported |
| Location | collections spec revision 1, sections 4.2, 7; app.R line 10 |
| Codex claim | Repeated normalized headers, label parsing, replacement, rematching, template contents, the size limit, and the CSV schemas are undefined. |
| Check performed | Read a CSV with header `collection, Latitude ,latitude,longitude` and rows named "NA" and "007" using `read.csv`: both latitude columns are kept, and the collection column becomes integer with values NA and 7. Read app.R line 10: the upload limit is 40 MB, which revision 1 referred to only as "the file size limit of the outline file". Revision 1 did not define replacement, rematching after a rename, template rows, or column lists. |
| Acceptance check performed | Run for header and parsing behaviour. The remaining contracts are specification gaps confirmed by reading. |
| Verdict | Confirmed |
| Resolution | Fixed, with the author's decision of 2026-11-17 (replace and rematch). Section 4.2 reads all columns as text with no missing-value strings, rejects repeated normalized headers, states the 40 MB limit, defines atomic replacement, rematching by current name, and template rows. Section 7 lists both CSV schemas in order with empty cells for missing values. |
| Tree evidence | Spec revision 2, sections 4.2 and 7 ("Downloads"). |
| Resolved on | 2026-11-17 |

### F8. The acceptance tests could pass with wrong implementations

| Field | Value |
|---|---|
| Severity (Codex) | major |
| Classification (Codex) | not demonstrated |
| Location | collections spec revision 1, section 10 |
| Codex claim | Simulations lack a generating distribution and seed, leave-one-out checks only bracket values, figures only need to build, and summaries, mean shapes, plot data, map geometry, and unavailable states are untested. |
| Check performed | Read revision 1 section 10 against each statistic of sections 4 to 7. Confirmed: no seed or distribution; "bracket" instead of equality; no test of section 5.5 values, mean reconstruction, colour counts, map limits, or rank-0 and rank-1 states. Checked the new oracles: for sizes 2, 2, 2 there are 6! / (2! 2! 2!) = 90 assignments and (G - 1)/(N - 1) = 2/5. Hand examples: A = (0,0),(2,0), B = (5,0),(7,0) gives d^2 = 25, v = 2, Q = 23; B = (1,0),(1,2) gives d^2 = 1, v_B = 2, Q = -1. Both recomputed by hand. |
| Acceptance check performed | Partly. Oracles verified by hand. Whether each test detects its target mutation is checked at build time by mutation runs. |
| Verdict | Confirmed |
| Resolution | Fixed. Section 10 fixes seed 20261117 and the generating distribution, uses exact hand and enumerated oracles, requires leave-one-out equality to 1e-12, and adds tests for labels, summaries, mean shapes, plot data, 9-level colours, map geometry, unavailable states, and report sections. |
| Tree evidence | Spec revision 2, section 10. |
| Resolved on | 2026-11-17 |

### F9. The location threshold and the cautions overstate what the design supports

| Field | Value |
|---|---|
| Severity (Codex) | minor |
| Classification (Codex) | unsupported |
| Location | collections spec revision 1, sections 5.6, 8 |
| Codex claim | Five locations is operational rather than derived, and cautions on fixed-data references, dependent pairs, sample means on the map, and approximate coordinates are missing. |
| Check performed | Read revision 1: section 5.6 justified 5 by saying 3 or 4 locations "cannot show a relationship", which is a judgement and not a derivation. Section 8 had no caution distinguishing references from intervals, none on dependence among pairs that share a collection, none on map colours as sample means, and none on rounded coordinates, although section 4.2 recommends rounding. |
| Acceptance check performed | Read. Text checks are specified for the build. |
| Verdict | Confirmed |
| Resolution | Fixed. Section 5.6 states the threshold as operational and counts pairs, sections 5 and 8 add the reference, dependence, sample-mean, and coordinate cautions, and section 4.2 notes the effect of approximate positions. |
| Tree evidence | Spec revision 2, sections 4.2, 5 (opening paragraph), 5.3, 5.6, 6 ("Limits"), 8. |
| Resolved on | 2026-11-17 |

## Held by the auditor

- E(SSB/SST) = (G - 1)/(N - 1) for fixed rows and fixed collection sizes, including unequal sizes, by derivation, simulation, and enumeration. 2/14 rounds to 0.14.
- The untruncated squared-distance correction is unbiased for independent coefficient vectors in a fixed representation with finite second moments, without equal covariance matrices.
- The n - 1 divisor gives a within-collection variance whose expectation does not depend on sample size. Leave-one-out ranges are correctly described as sensitivity ranges.
- Within-pair reassignment of raw distances matches the quantity displayed beside its reference.
- The haversine formula is correct in radians. The one-degree, identical-point, and meridian values hold.
- Five locations give 120 assignments including the identity, and coincident positions keep their multiplicities.
- Existing reactives can supply the shared alignment, symmetry option, H, coefficients, and PCA without group labels changing numerical results.
- The coefficient layout and zero-offset inverse transform support mean-shape reconstruction.
- `category_scale` repeats beyond 8 levels. A shared fallback for colour and fill can cover the group-coloured figures.
- The report architecture supports shared plotting functions and text continuation pages. A zoomed `coord_quickmap` window built on a synthetic polygon.

## Not examined by the auditor

- `maps` is not installed, so real map geometry was not tested. Open: it must be installed before the build, and the map tests of section 10 cover it.
- The project-root `CLAUDE.md` is absent. Recorded as absent and closed.
- The extension is unbuilt, so no server, browser, PDF, deployment, or performance check was run. Open until the build, the build audit, and the checklist.
- Prior ledgers, image data, and unlisted files were excluded by the brief. Closed: the exclusion was intended.
- Flagged literature attributions were not verified. Open: the two CITE-CHECK flags in section 12 remain.
- Independence and unbiasedness after fitted alignment and harmonic selection were not established. Open: stated in sections 5.3 and 12 as an unquantified limitation.
- No command was blocked by the sandbox. Closed.

## Questions and answers

1. What latitude and extent limits should the map support, and where is the scale bar valid? Answered by the author, 2026-11-17: withhold the map beyond 80 degrees of latitude and across the 180 degree meridian, with the scale bar labelled approximate. Applied in section 6.
2. Should existing plots and score exports adopt canonical labels and consistent colours? Answered by the author, 2026-11-17: trim labels everywhere. Applied in sections 3, 4.1, and 7.
3. What replaces an unavailable PC selection, and what are the CSV schemas and missing-value conventions? Proposed by Claude and not objected to by the author, 2026-11-17: the selector lists existing components and falls back to mean length divided by width; fixed column lists with empty cells. Applied in sections 6 and 7.
4. Should uploads replace the locations table, with later edits rematched by current name? Answered by the author, 2026-11-17: replace and rematch. Applied in section 4.2.
5. Will population claims be qualified as applying to independent vectors in a fixed representation? Proposed by Claude and not objected to by the author, 2026-11-17: yes. Applied in section 5.3.
6. Is the share a numeric output, fixing six figures, and should report acceptance allow continuation pages? Proposed by Claude and not objected to by the author, 2026-11-17: the share is drawn as a figure (seven figures in total), and continuation pages are allowed. Applied in sections 7 and 10.

No answer led to a follow-up round.

## Sign-off

2026-11-17. The author (Carl Lipo) signed off on this ledger and on specification revision 2, and instructed that the build proceed.
