# Specification audit, 2026-11-17

Audited working-tree HEAD: `a85cb88d0e3f966e553b23e0be5958612b7a38e7`, matching the brief. Read-only git diff found no changes to the scoped spec, app, README or images. The brief identifies `f8eeed2` as the deployed version, which was not separately checked. No project files were edited and no git mutation commands were used.

The spec is not ready to build. Two demonstrated contradictions block its acceptance gate: the automatic segmentation procedure fails on three demo images, and its ellipse test rejects correct arc-length Fourier coefficients. The remaining findings specify missing or incorrect contracts that affect results and reproducibility. The author's decisions to replace Momocs, retain descriptive region analysis and omit the covariance ratio remain accepted.

## Execution and limits

Read `~/.claude/CLAUDE.md` first, then the complete spec, complete 2237-line app, four JPEGs and complete README in order. The project CLAUDE.md is absent as the brief stated. Evidence scripts ran in R 4.6.2 (2026-10-29) with Momocs 1.4.1, imager 1.0.8, shiny 1.14.0, bslib 0.12.0, ggplot2 4.0.3, cluster 2.1.8.3 and testthat 3.3.2. Full loaded-package versions are in checks.log. Python 3.9.6 and its standard library generated this report and JSON. NumPy, SciPy, Pillow and scikit-image were unavailable; no installation was attempted. shinytest2 was unavailable. The runtime emitted an OpenMP temporary-file warning but computations completed. A final here-document evidence-integrity assertion was sandbox-blocked and was not worked around; see Not examined. Earlier JSON syntax validation passed.

Reproduce from `codex/`: `Rscript evidence/2026-11-17_spec/checks.R > evidence/2026-11-17_spec/checks.log 2>&1`, then the analogous command for `small_n.R`. `make_report.py` regenerates this report and the structured finding ledger. Derived mask previews are under `data/`; all reported numbers are also in this report or tracked logs.

## Defect-by-defect verdict

| Item | Verdict and checked location in app.R |
| --- | --- |
| D1 | Confirmed, 1562-1575. Installed efourier.Out names and fills coefficient blocks A, B, C, D. |
| D2 | Confirmed, 19-44, 788, 1668-1947. Indexwise means and region comparisons lack anatomical or tip/base reindexing. Alignment changes coordinates, not the indexing. |
| D3 | Confirmed, 1903-1925. Numerator is a sum of correlations; denominator a sum of mean distances. Rescaling distances changes this ratio, which is not a covariance ratio. |
| D4 | Confirmed, 90-112, 279-292, 1239-1251. There is a successful-return condition but no iteration/depth/failure bound. The left blade uses the same forward traversal pattern as the right, with a different shoulder endpoint. |
| D5 | Confirmed for region construction, 1125-1233. Blade clicks are displayed and required by the completeness guard but do not define a region. Indices are carried from the current specimen to the mean. |
| D6 | Confirmed, 1632-1665. The tab overlays actual specimens at PCA extrema, with no TPS fit or deformation grid. |
| D7 | Confirmed for successful presses, 803-808. Append-only lists; reset clears previews, not stored outlines; no removal control. |
| D8 | Confirmed, 377-403 and all four photographs. Defaults darker/up disagree with light/left demos. |
| D9 | Confirmed, 742-777. Thresholding and opening do not fill enclosed ink holes. Largest contour is chosen by nrow. |
| D10 | Confirmed, 590-606. Coordinates from the plot click are inserted directly into CSS left/top pixels. |
| D11 | Qualified, 1067-1104. n=1 and k>n fail. n=2,k=2 passes the isolated PCA/clustering path but is rank deficient; default k=3 fails. See F12. |
| D12 | Confirmed structurally, 1077-1085. coo_align uses SVD and normalized EFA has no anatomical tip input. No endpoint correspondence is enforced. |
| D13 | Confirmed for normal CRAN installation by the CRAN package page. Archive sources remain available. |
| D14 | PDF duplication and missing CSV/RDS handlers confirmed, 1950-2237. README promises these exports and analysis_script.R. That separate script's existence was outside the permitted file reads. |

## Numerical and image results

The coefficient equations are correct, as follows by integrating the piecewise-linear contour by parts twice, or by direct comparison with installed Momocs 1.4.1. The sign and factor T/(2*pi^2*n^2) match. Power is the sum of squared sine/cosine coefficients divided by two, by orthogonality. The independent 2000-vertex 3-by-1 ellipse gives first-harmonic amplitudes 2.62452424251 and 1.09927297350. Its third-harmonic x cosine coefficient is 0.241790947734; first, third and fifth powers are 4.04826428490, 0.0375992257941 and 0.00343731058546. Maximum Momocs discrepancy is 3.46944695195e-18. The equations and the ellipse acceptance expectation cannot both hold.

The prototype used native-size grayscale JPEGs, the first maximizing Otsu bin split, a disk of radius 2, 4-connected components, explicit opening border handling, hole filling and area selection. It used no crop and no downscaling, because every image is already below 1200 pixels.

| Image | Dimensions | Otsu | Light border fraction | Automatic result |
| --- | --- | --- | --- | --- |
| artifact1.jpg | 244 x 176 | 0.34765625 | 439/836 = 0.525119617225 | Background selected, area 37485, fraction 0.872880961252, border contact, 2 open contours of 410 and 4 vertices |
| artifact2.jpg | 278 x 182 | 0.41796875 | 477/916 = 0.520742358079 | Background selected, area 45140, fraction 0.892165388568, border contact, 2 open contours of 446 and 4 vertices |
| artifact3.jpg | 194 x 186 | 0.375 | 235/756 = 0.310846560847 | Artifact selected, area 10374, fraction 0.287495843033, no border contact, 1 closed contour of 725 vertices |
| artifact4.jpg | 134 x 170 | 0.44140625 | 326/604 = 0.539735099338 | Background selected, area 19340, fraction 0.848990342406, border contact, 2 open contours of 289 and 4 vertices |

The top ruler and light margins make border-majority background estimation fail. A diagnostic run forcing light foreground, otherwise identical, produced one contour without border contact for each image. Artifact/ruler component areas were 14836/5252, 15095/5219, 10374/4076 and 5920/3271. Thus area selection separates the ruler once polarity is correct, and hole filling removes enclosed marks. This diagnostic is not evidence that the automatic requirement passed. Silhouette accuracy against hand-traced outlines was not assessed.

On these forced-light artifact contours, principal-moment ratios were 5.45762584711, 8.15069631739, 2.78218000218 and 2.74029876743. Axial length/width ratios were 2.19711803625, 2.61580446263, 1.53308801939 and 1.48652344348. Each had a unique maximal projection and one base-side midline intersection. The demos have lateral constrictions and irregular/blunt-looking tips, but none exercises the exact-45-degree, multiple-base-crossing or wider-than-long cases. The photographs alone do not establish whether an apex is broken. For a basal notch, selecting the most negative intersection is a determined geometric rule; it may place B inside the notch rather than on the basal-most extent. It must not be presented as inferred anatomical homology.

Inclusive sampling of both endpoints gives 200 rows and only 198 distinct points. The half-open convention in F5 removes this problem. Fixed proportional correspondence is usable for descriptive variation along the chosen sampling map, but does not justify anatomical region claims. More landmarks or sliding are not automatically required. Adequacy of 100 points per edge remains a sensitivity question, especially for long barbs and small notches. Pooling both edges avoids separate left/right reports but does not make results invariant to independent reflections of specimens; symmetrization is the specified way to remove signed asymmetric variation.

GPA-aligned outlines are legitimate EFA inputs. Fixed T-based phase and traversal avoid needing first-ellipse normalization. Individual tip/base axes need not remain vertical after GPA, which is acceptable for least-squares shape alignment. Rotating the mean to vertical must apply the same final rotation to every aligned specimen. Measurements on saved pixel outlines remain distinct from dimensionless aligned reconstructions. Symmetrization order and reconstruction translation/phase still need the contracts in F3 and F8.

The PCA counterexample has full-data eigenvalues 0.680066666667 and 0.666666666667. After omitting one row they are 1 and 0.340033333333. Ordered-vector cosines are both zero, yet both subspace cosines are one. This is an axis-order exchange, not disappearance of the two-dimensional variation. Near-degenerate axes require an eigenspace diagnostic and rank gating.

For column-centered blocks X and Y, specify RV = ||X'Y||_F^2 / sqrt(||X'X||_F^2 ||Y'Y||_F^2). Common covariance denominators cancel. Do not substitute column-standardized correlations without changing the declared statistic. Constant blocks make the denominator zero. The independent n=4, p=40, q=80 check gives RV=0.928931628088, permutation mean 0.937697931834 and 2.5/97.5 percentiles 0.900716863552/0.990035741167. This demonstrates the small-sample inflation the display is meant to expose. At n=4 there are only 24 possible pairings, so 999 draws repeat pairings; the range is discrete, not a precision statement about population RV.

The limitation about common superimposition is accurate. A frozen-coordinate shuffle is a conditional reference, not an independence reference that reproduces that fitting contribution. For the latter, [Klingenberg (2009), permutation section](https://onlinelibrary.wiley.com/doi/10.1111/j.1525-142X.2009.00347.x) describes refitting. The [official MorphoJ guide](https://morphometrics.uk/MorphoJ_guide/covariation/pls_intraconfig.htm) also distinguishes these cases. A descriptive conditional display is permissible under clause (c) of the shared conventions. The stronger claim that Bayesian covariance models are intractable is not demonstrated and is unnecessary to justify this design.

## Dependency and acceptance coverage

No proposed runtime feature is inherently impossible with the four packages plus custom base-R code. imager supplies morphology/labeling/contours; stats supplies prcomp, dist, hclust and cutree; standard R graphics facilities plus ggplot2 can produce the plots and PDF. Base R can implement EFA, GPA, Otsu, polygon moments, RV, CSV persistence and silhouette width. The last routine is not supplied by the stated package set: [cluster::silhouette](https://cran.r-project.org/web/packages/cluster/refman/cluster.html) is the standard additional dependency. Testthat, shinytest2, reference-generation packages and rsconnect are development/deployment dependencies, not necessarily runtime dependencies. imager is compiled and has system requirements; clean local/deployment verification remains necessary. [CRAN imager metadata](https://cran.r-project.org/web/packages/imager/index.html) confirms the stated version/date. [CRAN Momocs status](https://cran.r-project.org/web/packages/Momocs/index.html) confirms the archive date.

| Section 10 item | Objective as written? |
| --- | --- |
| Synthetic Otsu/polarity/holes/components | Partial: algorithmic objectives are clear, but fixture and tie expectations are absent. |
| Four demos | Yes for closure and warnings, but demonstrated failure; add mask identity checks. |
| Orientation | No: matching direction at ties and numerical tolerance unspecified. |
| Resampling | Partial: row count and landmark indices clear; uniqueness, endpoint and spacing metric absent. |
| Measurements/calibration | No: geometry, expected values, units and tolerances absent, with width-location ambiguity. |
| EFA | (i) numeric tolerance specified but phase missing; (ii) metric/tolerance missing; (iii) false; (iv) formula clear but tolerance unspecified. |
| GPA fixture | Partial: reference choice/settings and matching metric/tolerance missing. |
| PCA fixture | Partial: sign allowance insufficient for tied eigenspaces; numeric tolerance absent. |
| RV fixture | Partial: identical nonconstant blocks should give 1; zero denominator, centering and random fixture tolerance absent. |
| Small n, k clamp, replacement | Mostly objective; specify which panels remain available and add rank/identity cases. |
| Browser rendering and nonempty downloads | Mechanically testable but insufficient for correct data/content; no round-trip or schema checks. |
| Clean demo, all tests pass, CRAN dependency scan | Mechanically testable after repairs and environment setup. Not run on nonexistent rewrite. |

## Findings

### F1: blocking, wrong

Location: `docs/specs/2026-11-17_app_rewrite_spec.md:126-134; 317-318; 361-362`.

Claim: Automatic settings produce one closed artifact outline without warnings for every demo.

Problem: The border-majority rule reverses polarity for artifacts 1, 2 and 4. Hole filling then absorbs the light artifact into the selected dark-background component. These results touch the border and produce open contours.

Evidence: Section 6.1 calls the majority border class background. checks.log records light-border fractions 0.525119617225, 0.520742358079, 0.310846560847 and 0.539735099338. Images 1, 2 and 4 each produce two open contours. Image 3 produces one closed 725-vertex contour. All checks used the original full images, Otsu, disk opening radius 2 and 4-connectivity.

Repair suggestion: Replace or qualify border-majority polarity with a rule that handles the light image margins and top scale bar. Specify connectivity, border handling and the opening kernel. Retain manual override and reject unusable masks before analysis.

Acceptance check: Run the final specified pipeline on all four unchanged JPEGs. Assert one closed artifact contour, no border contact, no scale bar, filled enclosed ink holes, and the stated area/point-count checks. Save reference masks so closure alone cannot certify the background.

### F2: blocking, wrong

Location: `docs/specs/2026-11-17_app_rewrite_spec.md:247-253`.

Claim: An ellipse with semi-axes 3 and 1 has first-harmonic amplitudes 3 and 1 and all higher harmonics below 1e-4 under the displayed formulas.

Problem: Those are angular-parameter Fourier expectations. The specified transform parameterizes the outline by arc length, so a noncircular ellipse has higher harmonics. Correct code cannot pass this test.

Evidence: Independent R evaluation at 2000 vertices gives first-harmonic x/y amplitudes 2.62452424251 and 1.09927297350, with a third-harmonic coefficient 0.241790947734. Maximum coefficient disagreement with installed Momocs 1.4.1 is 3.47e-18 (checks.log).

Repair suggestion: Retain the arc-length equations. Replace the ellipse expectation with verified arc-length fixture values, and use a circle for the first-harmonic-only test. Keep reconstruction accuracy as a separate test.

Acceptance check: The independent implementation and Momocs agree within 1e-8 on the same polygon and phase convention; the revised ellipse and circle fixtures both pass.

### F3: major, unsupported

Location: `docs/specs/2026-11-17_app_rewrite_spec.md:224-253; 267-270`.

Claim: The specified coefficients suffice for reproducible inverse reconstruction, Momocs matching and automatic harmonic selection.

Problem: The inverse, constant terms, phase origin, evaluation positions and reconstruction error metric are absent. Landmark-centering does not force the arc-length mean to zero. Momocs starts its first segment at the last vertex, whereas the spec anchors index 1 at T. The denominator for the 99 percent power rule and the zero-higher-power case are also undefined.

Evidence: The text gives only n=1..H equations and asks for an inverse within a stated tolerance without stating one. Installed efourier.default uses Dx <- coo[,1] - coo[,1][c(nr,1:(nr-1))]. A centered triangle has nonzero Fourier DC (0.0274850704906, 0.0393039712659), recorded in checks.log.

Repair suggestion: Specify x(t)=x0+sum(a_n cos+b_n sin), the analogous y(t), the arc-length constant terms or an explicit recentering convention, t=0 and polygon closure. Reconcile the Momocs fixture phase explicitly. Define evaluation times, error units/tolerance and which original/aligned outline is overlaid. Select H against a fixed declared reference spectrum and handle zero residual power.

Acceptance check: An asymmetric polygon passes coefficient matching after the documented phase conversion and reconstructs in the declared frame. A circle, a sharp polygon, H below 16, and a zero-residual-power spectrum have deterministic results.

### F4: major, unsupported

Location: `docs/specs/2026-11-17_app_rewrite_spec.md:148-161; 363-367`.

Claim: The orientation and landmark procedure determines T, B and u for every outline.

Problem: At exactly 45 degrees both principal axes qualify if the bound is inclusive, and neither qualifies if strict. A flat or broken tip can maximize projection at multiple points or along a whole segment. The optional click override has no snapping or invalid-landmark rules. A geometric extreme is not necessarily the surviving anatomical tip.

Evidence: Step 2 says choose the axis within 45 degrees; step 4 says the outline point with the largest projection. Neither specifies ties. Step 6 permits clicking T and B but does not require distinct boundary points or define segment insertion.

Repair suggestion: Specify axis and tip tie handling, a warning/manual-resolution policy for broken or blunt ends, click projection onto segments, and rejection of coincident T/B. Define how B is inserted into the polygon. Keep the farthest-intersection rule if the authors accept its geometric meaning on concave bases.

Acceptance check: Fixtures cover exact 45 degrees, a flat tip, truncated oblique tip, near-circle, wider-than-long outline, basal notch and a midline with several crossings. Expected landmarks or required manual intervention are explicit and invariant to contour start index.

### F5: major, unsupported

Location: `docs/specs/2026-11-17_app_rewrite_spec.md:165-175; 283-287; 368`.

Claim: Two 100-point edge samples yield 200 distinct coordinates with T at 1 and B at 101.

Problem: Endpoint inclusion is unspecified. The usual inclusive sampling gives 200 rows but only 198 distinct points, duplicates B at 100/101 and T at 1/200, and creates zero-length segments in EFA. Region fractions must reverse on the second edge, but the index-to-fraction map and boundary ownership are absent.

Evidence: Section 6.3 says 100 points starting at each landmark without specifying the excluded endpoint. checks.log reproduces 200 rows and 198 distinct points for inclusive sampling. Section 7.6 defines both edges from tip 0 to base 1 although the second stored edge starts at B.

Repair suggestion: Sample each directed edge at j/100 for j=0..99, excluding its destination. Specify second-edge tip fractions as 1-j/100, the mirror relabelling, and nonoverlapping interval endpoints. Define spacing as distance along the source contour, not Euclidean chords across corners.

Acceptance check: A simple asymmetric polygon has 200 distinct vertices, T only at 1, B only at 101, all positive segment lengths, correct paired-edge fractions, and exactly one region assignment per point.

### F6: major, wrong

Location: `docs/specs/2026-11-17_app_rewrite_spec.md:204-206; 257-270; 352-353`.

Claim: Same-number leave-one-out PC cosines describe stability, and PC1 to PC3 reconstructions are available once PCA is enabled at n=3.

Problem: Near-tied eigenvectors can exchange order or rotate within an unchanged subspace. The cosine is a valid ordered-vector sensitivity measure but cannot distinguish that from loss of the represented shape variation. At n=3, leave-one-out fits have rank at most 1, and the full fit has rank at most 2. Degenerate specimens can reduce rank further.

Evidence: For rows (1.01,0), (-1.01,0), (0,1), (0,-1), omitting the first row gives PC1/PC2 cosines 0/0 while both subspace principal-angle cosines equal 1. checks.log also records ranks 2 and 1 before/after leaving one of three specimens out.

Repair suggestion: Use singular values of V_full^T V_refit to report principal angles for a declared leading eigenspace, with a rule for clusters of near-tied eigenvalues. Report individual-axis cosines only where identifiable. Freeze H and state whether GPA is refitted. Gate reconstructions and comparisons by numerical rank, and describe n-1 as an upper bound.

Acceptance check: The order-swap example reports stable subspace and ambiguous axes. Three specimens do not produce PC3 or a rank-2 leave-one-out stability claim. Identical and rank-deficient shapes yield explicit unavailable messages without arbitrary zero-eigenvalue axes.

### F7: major, unsupported

Location: `docs/specs/2026-11-17_app_rewrite_spec.md:119-121; 179-190; 347-350`.

Claim: Maximum width has a defined longitudinal position, and calibration can be applied to all images after independent downscaling.

Problem: The transverse bounding extent can use extrema at different heights, so it has no unique position along length. Rectangles also have tied maximum-width positions. A shared millimetres-per-downscaled-pixel factor is wrong when images have different resize factors. Full-resolution outline is ambiguous after mandatory downscaling.

Evidence: Width is defined as extent perpendicular to u, followed by position of maximum width. All coordinates are declared downscaled in 6.1, while 6.4 requests a full-resolution outline and a checkbox applying current calibration to all images.

Repair suggestion: Distinguish transverse extent from maximum cross-sectional width. Define cross-sections, multiple intervals at notches, tied positions, and whether base means B or minimum axial extent. Store original dimensions and resize/crop transforms. Transfer scale through original-pixel units, or prohibit batch transfer across unequal resize factors. State area units, dimensionless outputs and the meaning of full-resolution.

Acceptance check: Known rectangle, ellipse, slanted/asymmetric and notched fixtures yield declared widths/positions. Images with resize factors 1 and 0.5 give the same physical dimensions before/after crop and after batch calibration.

### F8: major, unsupported

Location: `docs/specs/2026-11-17_app_rewrite_spec.md:194-200; 210-220; 283-287; 319-336`.

Claim: The specimen CSV restores the store, and setting changes reproducibly update analysis while preserving saved specimens.

Problem: A coordinate-only CSV does not define restoration of calibration, names/groups, full contours, landmarks or processing settings. Editable names and duplicate upload basenames leave identity unresolved. The spec does not say whether changing controls mutates saved outlines, whether symmetrization occurs before or after GPA and renormalization, or which measurements remain raw. Averaging with a mirror can change centroid size.

Evidence: The store holds both 200-point and full-resolution outlines plus metadata, but its persistence format is only a CSV of outline coordinates. Results recompute when specimens or settings change. Symmetrized shapes are described only as the average of each configuration and its relabelled mirror.

Repair suggestion: Define a versioned CSV schema with stable image/specimen IDs, both outline representations and all required metadata. Separate preview processing, saved raw specimens and global analysis settings. Specify symmetrization before a fresh GPA (or another fully stated normalization sequence), measurement provenance, H/boundary invalidation and reload/merge behavior.

Acceptance check: Export, restart and reload reproduce measurements, groups, outlines and analysis. Renaming and duplicate basenames do not overwrite unrelated specimens. Toggling symmetry twice restores original results; H/region changes do not alter saved pixel coordinates or calibration.

### F9: major, unsupported

Location: `docs/specs/2026-11-17_app_rewrite_spec.md:294-311`.

Claim: The frozen-coordinate permutation range represents RV with no association, and Bayesian covariance modeling is intractable here.

Problem: The permutation is a valid descriptive reference for re-pairing the already aligned blocks, but does not preserve covariance induced by the common GPA. It is not an independence reference for anatomical parts processed through the whole pipeline. RV centering, normalization and zero-variance behavior are unspecified. No model or timing evidence establishes Bayesian intractability; the descriptive exception already permits this choice.

Evidence: The spec permutes one block after the common fit and labels the range as no association, while acknowledging common-fit covariance. Klingenberg (2009), permutation section, distinguishes frozen-coordinate permutation from refitting for a joint-fit independence reference: https://onlinelibrary.wiley.com/doi/10.1111/j.1525-142X.2009.00347.x . Shared conventions explicitly allow descriptive frequentist methods.

Repair suggestion: Give the centered-block trace formula and undefined-denominator behavior. Label the present range explicitly as conditional re-pairing of fixed Procrustes coordinates, not removal of superimposition effects; alternatively specify and justify a refitting reference. Use the descriptive-method exception, or substantiate the separate feasibility claim.

Acceptance check: RV matches the trace formula, equals 1 for a nonconstant block against itself, and is unavailable for a constant block. Fixed-seed permutations reproduce the reference. UI wording states what is randomized and does not present its range as a credible/confidence interval or as correction for GPA-induced association.

### F10: minor, unsupported

Location: `docs/specs/2026-11-17_app_rewrite_spec.md:89-110; 275-279; 357-386`.

Claim: The named runtime dependencies fully specify the functions needed for every feature.

Problem: Mean silhouette width has no supplied function in shiny, bslib, imager, ggplot2 or base R. It is implementable in base R, but unlike EFA and GPA the design does not assign or specify that implementation. Test and deployment dependencies also need separate declarations.

Evidence: Section 7.5 promises mean silhouette width. Installed package inspection identifies cluster::silhouette; the official cluster manual documents it: https://cran.r-project.org/web/packages/cluster/refman/cluster.html . Section 10 additionally requires testthat, shinytest2 and rsconnect.

Repair suggestion: Either declare cluster as a direct runtime dependency or specify a tested base-R silhouette routine, including singleton and coincident-point cases. List test/fixture/deployment packages separately. The four-package runtime design is otherwise feasible with its promised custom routines and standard R graphics/statistics facilities.

Acceptance check: The silhouette implementation matches cluster::silhouette fixtures for ordinary clusters, singleton clusters and duplicate points. A clean test environment installs all declared development dependencies, and runtime dependency scanning matches the chosen design.

### F11: major, not_demonstrated

Location: `docs/specs/2026-11-17_app_rewrite_spec.md:247-253; 355-386`.

Claim: The acceptance list provides an objective pass/fail gate.

Problem: Several entries specify topics rather than expected outputs. Orientation says within tolerance without a value; inverse reconstruction promises a stated tolerance but omits it. GPA and measurements lack error metrics/tolerances. Bimodal Otsu thresholds can tie. PCA sign-only comparison is insufficient for tied eigenvalues. Nonempty downloads do not establish valid content or round-trip restoration.

Evidence: Section 10 lists measurements on known dimensions and calibration conversion without values; 7.2(ii) says within a stated tolerance; the app test checks downloads only for non-emptiness.

Repair suggestion: Define fixture geometries, expected values, absolute/relative tolerances and tie rules. Choose the GPA reference settings and comparison metric. Require parsed CSV schemas/units and round-trip checks, valid PDF pages/content, and explicit edge-case messages. Keep numeric acceptance independent of the implementation under test.

Acceptance check: Each test entry has a deterministic assertion and specified tolerance where numerical. Incorrect units, swapped specimen labels, missing calibration, duplicate endpoints and an invalid-but-nonempty download are detected.

### F12: minor, unsupported

Location: `docs/specs/2026-11-17_app_rewrite_spec.md:54-55 (D11); bifaceAnalysis/app.R:1088-1104`.

Claim: Analysis fails with fewer than three specimens as an unconditional sample-size defect.

Problem: One specimen and k greater than n do fail. For two specimens, the reproduced PCA/clustering path succeeds with k=2, although only one PC has nonzero rank. The default k=3 explains the demonstrated two-specimen failure. A universal two-specimen runtime failure is not established by this code.

Evidence: small_n.R reproduces hclust failure at n=1 and cutree failure at n=2,k=3. At n=2,k=2 it returns two clusters and a 2-by-2 score matrix whose PC2 scores are approximately 1e-17. app.R sets k directly from input before cutree.

Repair suggestion: State the conditions precisely: missing n/k guards, default k=3 fails at n=2, and additional PCs are rank-deficient. Distinguish invalid interpretation from a runtime error.

Acceptance check: The defect description matches tests at n=1, n=2 with k=2 and k=3, and rank-deficient samples. It does not claim the complete Shiny observer was exercised by the isolated numerical test.

## Held

- Read the complete spec, all 2237 lines of app.R, all four demo JPEGs and the complete README in the requested order, after ~/.claude/CLAUDE.md. HEAD is a85cb88d0e3f966e553b23e0be5958612b7a38e7; scoped files have no working-tree diff. The project CLAUDE.md is absent.
- D1 is confirmed by installed Momocs 1.4.1: efourier.Out constructs A1..AH, B1..BH, C1..CH, D1..DH and concatenates vectors in that order. app.R:1565 groups consecutive quartets incorrectly.
- D2-D6 are supported by the cited source: no landmark-based reindexing; correlation/distance ratio instead of covariance ratio; stale curvature and same-direction blade ranges; blade clicks absent from region calculations; specimen overlays instead of thin-plate splines. D4 has a success return but no bounded recursion fallback. D5 blade clicks are used for display and completeness validation, not region construction.
- D7-D10 are supported: successful processing appends without removal; demo artifacts are light and left-pointing; contour selection uses vertex count and enclosed ink holes remain; click coordinates become CSS pixels. D12 has no anatomical tip anchoring in coo_align or normalized EFA. D14 duplicates plotting in the roughly 288-line PDF handler and promises CSV/RDS exports absent from app.R.
- CRAN confirms Momocs archival on 2026-03-18 and imager 1.0.8 publication on 2025-12-23. D13 holds for ordinary CRAN installation; archived Momocs source remains obtainable. Sources: https://cran.r-project.org/web/packages/Momocs/index.html and https://cran.r-project.org/web/packages/imager/index.html .
- The displayed EFA coefficient equations and power sum-of-squares/2 are correct. Independent coefficients agree with Momocs within 3.47e-18. Procrustes coordinates are a valid EFA input, and no first-ellipse normalization is necessary once scale, orientation, traversal and phase are consistently fixed. Raw pixel measurements can remain separate from normalized analysis.
- The farthest-midline-intersection rule resolves ordinary multiple crossings geometrically; it does not identify an anatomical base by itself. Under forced light-foreground diagnostics, all demos have one base-side crossing, a unique projection maximum and axial length/width greater than 1. None tests exact 45 degrees or wider-than-long behavior. Blunt/irregular apices and lateral concavities are visible, but breakage cannot be diagnosed confidently from these photos.
- The fixed arc-length correspondence limitation is accurate for descriptive geometric regions. It does not establish anatomical homology or prove 100 points per edge adequate. Pooling edges does not eliminate asymmetric variation caused by arbitrary photographed faces.
- RV small-sample/high-dimension inflation is reproduced: independent 4-by-40 and 4-by-80 Gaussian blocks give observed RV 0.928931628088 and permutation mean 0.937697931834, with range 0.900716863552 to 0.990035741167. RV(X,X)=1. Common-fit-induced covariation is a valid limitation, also documented in the official MorphoJ guide: https://morphometrics.uk/MorphoJ_guide/covariation/pls_intraconfig.htm .
- The descriptive permutation display is allowed by the shared research conventions. No Bayesian replacement or covariance-ratio implementation is required by this audit. Standard R statistics/graphics plus the declared packages can implement the proposed app; silhouette needs an explicit additional dependency or custom routine.

## Not examined

- The rewrite, Shiny browser behavior, shinyapps.io deployment, actual download files and rsconnect dependency scan were not run because the rewrite does not exist. shinytest2 is not installed. Existing app defects were checked statically, supplemented by isolated numerical tests, not a complete interactive run.
- Original Kuhl and Giardina (1982) and Escoufier (1973) PDFs were not obtained. Numerical formulas were checked independently and against installed Momocs 1.4.1. The original Adams paper was not fully audited; dropping the invented statistic is an accepted author decision.
- The external Momocs raw-source URL returned Cache miss; installed 1.4.1 supplied the function bodies. PMC returned reCAPTCHA for Klingenberg (2009); relevant publisher HTML and official MorphoJ documentation were available.
- The command-line script, root outputs, deployed-source identity at f8eeed2, prior audits and resolution ledgers were outside the brief reading scope. The README script path was not independently checked for existence. No existing author decision was reopened.
- No proof of anatomical correspondence, sensitivity study establishing adequacy of 100 points per edge, or Bayesian computational feasibility benchmark was performed. The segmentation run instantiated unspecified details explicitly; it is an independent prototype, not a test of future code.
- A final Python here-document intended to assert JSON finding IDs and existence of every listed evidence path was blocked: zsh could not create its temporary file (operation not permitted). It was not rerun by another mechanism. Earlier python3 -m json.tool validation succeeded. R emitted an OpenMP temporary-file warning but completed the calculations. An exploratory border-summary command had a syntax error and yielded no evidence; recorded segmentation diagnostics subsequently succeeded.

## Questions for the authors

- For broken/blunt tips and concave bases, should the app require manual landmarks, or accept purely geometric extremes with a warning?
- Should maximum-width position describe a cross-section, and which location should be reported when widths tie or a notch creates several intervals?
- Should symmetrization affect only analysis, leaving measurements and saved pixel outlines unchanged?
- Must the outline CSV restore the complete session, including full contours, calibration and settings, or only import coordinates for shape analysis?
- Should the RV range remain a conditional descriptive re-pairing reference, or model independence of parts with the superimposition step represented?

## Files written

- codex/audits/2026-11-17_spec_report.md
- codex/evidence/2026-11-17_spec/checks.R
- codex/evidence/2026-11-17_spec/checks.log
- codex/evidence/2026-11-17_spec/small_n.R
- codex/evidence/2026-11-17_spec/small_n.log
- codex/evidence/2026-11-17_spec/sources.txt
- codex/evidence/2026-11-17_spec/make_report.py
- codex/evidence/2026-11-17_spec/audit.json
- codex/evidence/2026-11-17_spec/data/segmentation1.png
- codex/evidence/2026-11-17_spec/data/segmentation2.png
- codex/evidence/2026-11-17_spec/data/segmentation3.png
- codex/evidence/2026-11-17_spec/data/segmentation4.png
