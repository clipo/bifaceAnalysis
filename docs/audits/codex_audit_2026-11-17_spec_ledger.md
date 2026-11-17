# Codex Audit Ledger: bifaceAnalysis, spec (2026-11-17_spec)

**Date:** 2026-11-17
**Prompt:** codex/audits/2026-11-17_spec_prompt.md
**Findings file:** codex/audits/2026-11-17_spec_findings.json
**Report:** codex/audits/2026-11-17_spec_report.md
**Scope:** docs/specs/2026-11-17_app_rewrite_spec.md; bifaceAnalysis/app.R; bifaceAnalysis/README.md; bifaceAnalysis/www
**Audited commit:** a85cb88d0e3f966e553b23e0be5958612b7a38e7
**Model and effort:** gpt-6-astra, high. Codex codex-cli 0.169.0. Session 01a1eb09-bcf0-7ee0-8a0f-1d5bacd77627. Effective sandbox: workspace-write [workdir, $TMPDIR].
**Run:** 2026-11-17T14:44:06-05:00 to 2026-11-17T14:54:57-05:00, exit status ok.

Checks by Claude are in `docs/audits/2026-11-17_spec_checks.R` with output in
`docs/audits/2026-11-17_spec_checks.log`. Prototype results described as
"scratchpad" were produced outside the repository and are not tracked; they are
superseded by the tests required in spec section 10.

## Summary

| Verdict | Count |
|---|---|
| Confirmed, fixed and verified in tree | 12 |
| Confirmed, deferred with owner | 0 |
| Confirmed, declined by the author | 0 |
| Refuted | 0 |
| Undetermined | 0 |

Blocking findings open: 0. Milestone proceeds: yes. Author sign-off: 2026-11-17.

All fixes are to the specification (revision 2). Most of Codex's acceptance
checks describe tests of code that does not exist yet. They were written into
spec section 10 and will be run at the build milestone.

## Findings

### F1. Automatic segmentation fails on the demo images

| Field | Value |
|---|---|
| Severity (Codex) | blocking |
| Classification (Codex) | wrong |
| Location | spec 6.1, 8, 10 (rev 1 lines 126-134, 317-318, 361-362) |
| Codex claim | The border-majority polarity rule does not yield one closed artifact outline on demo images 1, 2, and 4. |
| Check performed | Ran the rev 1 rule on the four JPEGs (`docs/audits/2026-11-17_spec_checks.R`, section F1). Border fractions above threshold 0.525, 0.521, 0.311, 0.540 (Codex: 0.525, 0.521, 0.311, 0.540). Images 1, 2, 4 take the artifact as dark, and in all four the foreground touches the border. |
| Acceptance check performed | Not runnable in the tree at spec stage. A scratchpad prototype of the replacement rule gave one closed contour per image, class light, no warnings, ink holes filled (visually inspected). The reference-mask comparison is now a required test in spec section 10. |
| Verdict | Confirmed |
| Resolution | Fixed in the spec: section 6.1 steps 4 to 8 replace the rule with largest component not touching the crop border, over both classes; connectivity, structuring element, and border handling are specified. |
| Tree evidence | Spec line 149: "Select the artifact: the component with the largest pixel area among those that do not touch the border of the crop". Read back with grep after the edit. |
| Resolved on | 2026-11-17 |

### F2. Ellipse test expectation is wrong

| Field | Value |
|---|---|
| Severity (Codex) | blocking |
| Classification (Codex) | wrong |
| Location | spec 7.2 item (iii) (rev 1 lines 247-253) |
| Codex claim | A 3-by-1 ellipse does not give first-harmonic amplitudes 3 and 1 under arc-length parameterization. |
| Check performed | Momocs 1.4.1 `efourier(norm = FALSE)` on a 2000-vertex 3-by-1 ellipse: A1 2.624521, D1 1.099272, A3 0.241791 (checks log, section F2). Circle of radius 2: A1 = D1 = 1.999988, higher harmonics 0. |
| Acceptance check performed | Run in part: independent base-R prototype agreed with Momocs to 4.97e-14 on a demo outline (scratchpad). Formal fixture tests are in spec section 10. |
| Verdict | Confirmed |
| Resolution | Fixed in the spec: item (iii) now expects the Momocs values (a_1 = 2.6245, a_3 = 0.2418); a circle test (ii) covers the single-harmonic case. Equations unchanged. |
| Tree evidence | Spec line 350: "matches committed Momocs values, including a_1 = 2.6245 and a_3 = 0.2418". Read back with grep after the edit. |
| Resolved on | 2026-11-17 |

### F3. Inverse transform and conventions unspecified

| Field | Value |
|---|---|
| Severity (Codex) | major |
| Classification (Codex) | unsupported |
| Location | spec 7.2, 7.4 |
| Codex claim | Inverse, constant terms, start point, Momocs phase, reconstruction tolerance, and the zero-residual-power case are not specified. |
| Check performed | Read rev 1 section 7.2: only forward equations given; "within a stated tolerance" with no value. Arc-length mean of a vertex-centred 3-4-5 triangle computed as (0.1667, 0), not zero (checks log, F3). Momocs phase confirmed: prototype matches Momocs only when the first vertex is moved to the end. |
| Acceptance check performed | Not runnable at spec stage; the listed cases (circle, H below 16, zero residual power, Momocs phase) are tests in section 10 and 7.2. |
| Verdict | Confirmed |
| Resolution | Fixed in the spec: section 7.2 gives the inverse, x0 and y0, t = 0 at the tip, closure, the Momocs comparison convention, a fixed 30-harmonic reference spectrum, the zero-power rule, and a 0.5 percent reconstruction tolerance. |
| Tree evidence | Spec line 323: "x(t) = x0 + sum_n [a_n cos(2 n pi t / T) + b_n sin(2 n pi t / T)]"; line 354 tolerance. Read back with grep after the edit. |
| Resolved on | 2026-11-17 |

### F4. Orientation and landmark ties unspecified

| Field | Value |
|---|---|
| Severity (Codex) | major |
| Classification (Codex) | unsupported |
| Location | spec 6.2, 10 |
| Codex claim | Axis ties, tip ties, click projection, base insertion, and coincident landmarks are undefined. |
| Check performed | Read rev 1 section 6.2 steps 2 to 6: none of these cases is stated. |
| Acceptance check performed | Not runnable at spec stage; the fixtures Codex lists are tests in section 10 ("Orientation and points"). |
| Verdict | Confirmed |
| Resolution | Fixed in the spec: section 6.2 now states the 45 degree rule, tip tie rule, half-open crossing rule, base insertion, click snapping, same-vertex error, and a blunt-tip warning. Author decision 2026-11-17: blunt or broken tips are placed geometrically with a warning. |
| Tree evidence | Spec line 178: "On an exact tie (both axes at 45 degrees) choose the major axis". Read back with grep after the edit. |
| Resolved on | 2026-11-17 |

### F5. Edge sampling ambiguous

| Field | Value |
|---|---|
| Severity (Codex) | major |
| Classification (Codex) | unsupported |
| Location | spec 6.3, 7.6 |
| Codex claim | Endpoint handling is unspecified; sampling both ends of both edges gives 198 distinct points. |
| Check performed | Sampled both edges of a square with both endpoints included: 200 rows, 198 distinct (checks log, F5). |
| Acceptance check performed | Not runnable at spec stage; the 200-distinct-points test is in section 10. |
| Verdict | Confirmed |
| Resolution | Fixed in the spec: section 6.3 specifies fractions j/100, j = 0..99 per directed edge, positions from the tip, mirror pairing; section 7.6 gives half-open region intervals. |
| Tree evidence | Spec line 208: "its own length, j = 0..99, which includes the path's start and excludes its end". Read back with grep after the edit. |
| Resolved on | 2026-11-17 |

### F6. Same-rank leave-one-out cosines misreport stability

| Field | Value |
|---|---|
| Severity (Codex) | major |
| Classification (Codex) | wrong |
| Location | spec 7.1, 7.3, 7.4 |
| Codex claim | When PC1 and PC2 exchange order the same-rank cosines are 0 while the plane is unchanged; PC3 is promised at n = 3. |
| Check performed | Reproduced Codex's fixture: same-rank cosines 0 and 0, principal-angle cosines 1 and 1. Three specimens give 2 components (checks log, F6). |
| Acceptance check performed | Not runnable at spec stage; the order-swap, three-specimen, and identical-shape cases are tests in section 10. |
| Verdict | Confirmed |
| Resolution | Fixed in the spec: section 7.3 makes plane agreement primary, reports direction agreement only for separable components (variance ratio 1.5), freezes H, states the alignment is not refitted, and gates all outputs by rank; 7.4 uses min(3, r) components. |
| Tree evidence | Spec line 371: "Plane agreement (primary): the product of the cosines of the two principal angles". Read back with grep after the edit. |
| Resolved on | 2026-11-17 |

### F7. Width position and scale transfer undefined

| Field | Value |
|---|---|
| Severity (Codex) | major |
| Classification (Codex) | unsupported |
| Location | spec 6.1, 6.4 |
| Codex claim | Bounding width has no longitudinal position; a scale shared across images resized by different factors is wrong; "full-resolution" contradicts the downscaling. |
| Check performed | Read rev 1 sections 6.1 and 6.4: all three inconsistencies are present as described. |
| Acceptance check performed | Not runnable at spec stage; rectangle, ellipse, notched, and resize-factor tests are in section 10. |
| Verdict | Confirmed |
| Resolution | Fixed in the spec: section 6.4 separates maximum width (transverse extent) from cross-section width and defines the widest position, ties, and notches (my choice, see question 2); the shared scale is held in original pixels and converted by each image's resize factor; "working resolution" replaces "full resolution". |
| Tree evidence | Spec line 252: "is held in original pixels per millimetre (working value divided by f)". Read back with grep after the edit. |
| Resolved on | 2026-11-17 |

### F8. Persistence and state separation unspecified

| Field | Value |
|---|---|
| Severity (Codex) | major |
| Classification (Codex) | unsupported |
| Location | spec 6.5, 7.1, 7.6, 8 |
| Codex claim | CSV content, ids, symmetrization, and what recomputes on a settings change are not defined. |
| Check performed | Read rev 1 section 6.5: the store holds more than the "CSV of outline coordinates" could restore; symmetrization is one sentence. |
| Acceptance check performed | Not runnable at spec stage; export, reload, duplicate-name, and settings-invariance tests are in section 10. |
| Verdict | Confirmed |
| Resolution | Fixed in the spec: section 6.5 separates preview, saved, and analysis state, defines schema version 1 with full restore (author decision 2026-11-17), stable ids, and collision behaviour; section 7.1 specifies symmetrization as analysis only (author decision 2026-11-17). |
| Tree evidence | Spec line 273: "The store is downloaded as a CSV, schema version 1, one row per contour vertex". Read back with grep after the edit. |
| Resolved on | 2026-11-17 |

### F9. RV reference mislabelled; Bayesian justification unsupported

| Field | Value |
|---|---|
| Severity (Codex) | major |
| Classification (Codex) | unsupported |
| Location | spec 7.6 |
| Codex claim | The reference is conditional on a fixed joint alignment; the intractability claim is unsupported. |
| Check performed | Read rev 1 section 7.6: the text calls the reference "no association" without the condition, and asserts intractability with no evidence. `~/.claude/CLAUDE.md` exception (c) covers descriptive statistics. |
| Acceptance check performed | Not runnable at spec stage; RV formula, self-RV, constant-block, and fixed-seed tests are in section 10. |
| Verdict | Confirmed |
| Resolution | Fixed in the spec: section 7.6 gives the trace formula, the NA rule, labels the reference as re-pairing of fixed coordinates (author decision 2026-11-17), says it is not an uncertainty interval or a correction for alignment, and cites exception (c). |
| Tree evidence | Spec line 426: "The reference is a re-pairing of fixed coordinates (author decision, 2026-11-17)". Read back with grep after the edit. |
| Resolved on | 2026-11-17 |

### F10. Dependencies incomplete

| Field | Value |
|---|---|
| Severity (Codex) | minor |
| Classification (Codex) | unsupported |
| Location | spec 5, 7.5, 10 |
| Codex claim | Silhouette has no stated source; test and deployment dependencies are not listed. |
| Check performed | Read rev 1 section 5: only four runtime packages named. `cluster::silhouette` conventions checked: singleton 0 (checks log, F10). |
| Acceptance check performed | Not runnable at spec stage; the comparison with `cluster::silhouette` is a test in section 10. |
| Verdict | Confirmed |
| Resolution | Fixed in the spec: section 7.5 specifies a base-R silhouette with singleton and zero-distance rules; section 5 lists runtime, test-only, fixture, and deployment packages. |
| Tree evidence | Spec line 399: "The implementation is checked against `cluster::silhouette`". Read back with grep after the edit. |
| Resolved on | 2026-11-17 |

### F11. Acceptance tests not objective

| Field | Value |
|---|---|
| Severity (Codex) | major |
| Classification (Codex) | not demonstrated |
| Location | spec 7.2, 10 |
| Codex claim | Several tests lack tolerances, fixtures, or content checks. |
| Check performed | Read rev 1 section 10: "within tolerance" appears without values; downloads are checked only for being non-empty. |
| Acceptance check performed | Partly run: section 10 was rewritten; whether every assertion is deterministic is settled only when the tests exist and run (build milestone). |
| Verdict | Confirmed |
| Resolution | Fixed in the spec: section 10 rewritten with fixtures, numeric tolerances, parsed CSV and PDF checks, and explicit small-sample cases. shinytest2 replaced by `shiny::testServer` plus a written manual checklist, with the reason stated. |
| Tree evidence | Spec line 537 and section 10 as a whole. Read back with grep after the edit. |
| Resolved on | 2026-11-17 |

### F12. Defect D11 overstated

| Field | Value |
|---|---|
| Severity (Codex) | minor |
| Classification (Codex) | unsupported |
| Location | spec 2 (D11); bifaceAnalysis/app.R:1088-1104 |
| Codex claim | The current pipeline does not fail unconditionally below three specimens. |
| Check performed | Ran the current sequence with Momocs 1.4.1 on its `bot` data: n=1 stops in hclust; n=2, k=3 stops in cutree; n=2, k=2 runs (checks log, F12). |
| Acceptance check performed | Run: matches Codex's small_n.log for n=1 and n=2 at k=2 and k=3. The interactive observer itself was not exercised. |
| Verdict | Confirmed |
| Resolution | Fixed in the spec: D11 reworded to state the missing guards and the cases separately. |
| Tree evidence | Spec line 57: "With two specimens it stops at the default cluster count of 3 and runs at 2". Read back with grep after the edit. |
| Resolved on | 2026-11-17 |

## Held by the auditor

- All scoped files and images read in full; HEAD matched the requested commit; no project CLAUDE.md.
- D1 confirmed from Momocs 1.4.1 source: coefficient columns are A1..AH, B1..BH, C1..CH, D1..DH. (Claude also confirmed this from `colnames(efourier(...)$coe)`.)
- D2 to D6 supported by the code, with qualifications: D4 has a successful return but no bounded termination on failure; D5 blade clicks are displayed and checked for completeness but construct no regions.
- D7 to D10 supported. D12 lacks anatomical tip anchoring. D14 duplicates plotting in the PDF handler and the README promises exports absent from app.R.
- D13 archive date and the imager version confirmed from CRAN.
- The Fourier coefficient equations and power formula are correct and agree with Momocs within 3.47e-18. Aligned coordinates are valid inputs; first-ellipse normalization is unnecessary when scale, rotation, traversal, and phase are fixed.
- The farthest-midline-intersection rule resolves ordinary multiple crossings. No demo image exercises the 45 degree or wider-than-long cases.
- The geometric-correspondence limitation is accurate. It does not establish homology or the adequacy of 100 points per edge, and pooling edges does not remove face-related variation. (Spec section 6.3 now says this.)
- RV(X, X) = 1 reproduced; small-sample inflation illustrated; alignment-induced covariance is a valid limitation.
- Descriptive permutation references are permitted by the shared conventions. All runtime features are implementable with base R and the stated packages.

## Not examined by the auditor

- The rewrite, browser behaviour, deployment, downloads, and the rsconnect scan: open until the build milestone. shinytest2 is unavailable on this machine (websocket fails to compile against OpenSSL 4.1), so browser behaviour will be checked by hand and reported as manual.
- Kuhl and Giardina (1982), Escoufier (1973), and Adams (2016) were not obtained. Open. The spec keeps CITE-CHECK flags on all three; the equations were verified numerically against Momocs, not against the papers.
- Momocs source URL and the PMC copy of Klingenberg (2009) were unreachable for the auditor. Closed for Momocs (installed 1.4.1 was used). The Klingenberg reference is not cited in the spec.
- The command-line script `biface_analysis.R`, root outputs, and the deployed bundle: outside the brief. Claude confirmed before the audit that the deployed bundle's app.R is byte-identical to `bifaceAnalysis/app.R` at f8eeed2. The script remains out of scope (spec section 4).
- Anatomical-correspondence validation and point-count sensitivity: open; listed in spec section 11.
- The auditor's final self-check script was blocked by the sandbox. Closed: the runner validated the findings JSON against its schema and reported status ok.

## Questions and answers

1. Should ambiguous or broken tips require manual landmarks, or remain geometric extremes with a warning? Author, 2026-11-17: geometric placement with a warning. Spec 6.2 steps 4 and 6.
2. Which cross-sectional width and location should be reported for tied widths or notched outlines? Not put to the author as a separate question. Claude chose: maximum width is the transverse extent; the widest position uses the outermost crossings at each level (spanning notches), with ties going to the level nearest the base. Spec 6.4. Author, 2026-11-17: accepted as defined, including dropping perimeter.
3. Should symmetrization affect analysis only? Author, 2026-11-17: yes, analysis only. Spec 7.1.
4. Must CSV reload restore the complete store? Author, 2026-11-17: yes, including scale and processing settings. Spec 6.5.
5. Should the RV reference remain a re-pairing of fixed coordinates? Author, 2026-11-17: yes, labelled as such. Spec 7.6.

No follow-up round was requested.

## Other events

- Codex wrote only under `codex/`. `changed_outside_stem` in the prompt frontmatter is empty.
- The report links to four images under `codex/evidence/2026-11-17_spec/data/`, which is not committed. The values they illustrate are in the tracked `checks.log`.
- Changes Claude made to the spec that no finding asked for: perimeter was dropped from the measurements (resolution dependent), and the speck-removal default changed from radius 2 to a square of side 3.
