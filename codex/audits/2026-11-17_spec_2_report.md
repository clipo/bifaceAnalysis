# Collections specification audit

Date: 2026-11-17. Mode: spec. Stem: `2026-11-17_spec_2`.

Audited working-tree HEAD, from `git rev-parse HEAD`:
`5481f0b2f346604195c2d4c6ebf070c6a158eaae`. This matches the brief. The
extension's introductory reference to `8d387bf` was not used as the audited
tree. No tracked modifications appeared in `git status --short`; two
pre-existing untracked audit prompt/findings files were not read or edited.

The core partition arithmetic and fixed-data reassignment designs hold. The
specification is not ready for unconditional build acceptance: a correct
haversine implementation fails its pole-distance test, and accepted map inputs,
unavailable PCs, label normalization, and export/state behavior need explicit
rules. The findings below concern the specification, not defects in an
implementation of the extension, which does not yet exist.

## Scope and execution

Read the shared `~/.claude/CLAUDE.md`, attempted the project `CLAUDE.md`, then
read the collections specification and the specified main-specification
sections/build notes, followed by the named R files in the requested order.
The project-root `CLAUDE.md` is absent (`cat` reported no such file and
`git ls-files CLAUDE.md` returned no entry). Read the scoped R files in full,
including the surrounding app wiring and plotting/report helpers. No prior
audit ledgers were opened because the brief limits the reading scope. No
recorded author choice of a static map, descriptive analysis, or exclusion of
group-separating ordinations is reopened here.

Independent checks are in `evidence/2026-11-17_spec_2/checks.R`, with all cited
numeric results in `checks.log`. Run from `codex/`:

```
Rscript evidence/2026-11-17_spec_2/checks.R > evidence/2026-11-17_spec_2/checks.log 2>&1
```

Environment: R 4.6.2 (2026-10-29), ggplot2 4.0.3, shiny 1.14.0, testthat
3.3.2. The calculations use base/stats R and the scoped shape/EFA/plot helpers.
Shiny and testthat versions were queried, but neither app server tests nor the
project test suite were run. `maps` is not installed. No packages were
installed. An OpenMP warning about sizing a `/tmp` file appeared; the R commands
completed successfully. No command was blocked by the sandbox, and no blocked
operation was retried by another route.

## Findings

### F1. Pole-distance acceptance value contradicts the specified radius

Severity: blocking. Classification: wrong.

Location: `docs/specs/2026-11-17_collections_spec.md:353-356`, with section 5.6,
lines 192-196.

The test requires: “The equator to the pole is 10007.543 km within 0.001.”
The specified radius, 6371.0088 km, gives `pi * R / 2 = 10007.557221018` km.
The discrepancy is 0.014221018 km, over fourteen times the allowed tolerance.
10007.543 is consistent with rounding a calculation using 6371 km instead.

Repair: retain the stated radius and change the expected value to
10007.557 km, tolerance 0.001 km. State explicitly that input degrees are
converted to radians before applying the formula.

Acceptance: all four geographic fixtures pass the same function and radius.
The independent results are 111.195080234, 10007.557221018, 0, and
111.195080234 km, respectively. The meridian case also needs an explicit
tolerance rather than an exact comparison to a rounded decimal.

### F2. The map window and scale rules do not cover the accepted coordinate domain

Severity: major. Classification: unsupported.

Location: collections specification section 6, lines 223-240; accepted
latitudes at lines 91-92.

The specification accepts both poles but requires a minimum longitude span
with “the equivalent ground distance in longitude at the middle latitude.”
The local conversion is `0.5 / cos(latitude)` degrees. It is undefined at a
pole. At 89.9 degrees, a minimum 0.5-degree latitude window centered on the
position extends to 90.15 degrees and needs approximately 286.479 degrees of
longitude. The spec has no clamping, withholding, or projection rule. This can
occur with two eligible collections at the same allowed position, so the
two-collection map threshold does not prevent it. Padding order relative to
the minimum span is also unspecified.

A scale bar specified at the middle latitude is locally meaningful there.
Its placement latitude and conversion back to longitude are not specified.
For the same angular bar, ground length at 60 degrees is only 0.707107 of
that at 45 degrees. Large north-south windows therefore cannot use an
unqualified, universally accurate scale bar. `coord_quickmap` sets an aspect
ratio, not a projection that fixes this.

The warning “Outside the contiguous United States a window smaller than a
country shows no boundary” is false for a window intersecting a coastline or
international boundary. Conversely, a small US window wholly inside a county
need not show a county line. Overlap with the US box correctly permits state
and county layers in a partially overlapping window; there is no need to
discard the non-US portion or truncate the overall extent to that box.

Repair: specify padding/minimum order, supported latitude and extent limits,
and behavior when those limits are exceeded, retaining the author's static-map
choice. Define the bar's placement, distance approximation, and tie rule.
Say “may show no boundaries” for small areas. Retain full polygon geometry and
zoom with coordinate limits, rather than dropping vertices outside the window.

Acceptance: inspect extents and bar endpoint distances for two coincident
nonpolar collections, near-pole and pole positions, wide latitude ranges,
and windows crossing each side of the US box. Test actual background geometry
with `maps` installed. A synthetic full polygon with a zoomed
`coord_quickmap` window built successfully here; this does not validate the
absent databases.

### F3. Valid zero-variation and rank-one data have undefined required outputs

Severity: major. Classification: unsupported.

Location: collections specification sections 5.5, 6 and 8, lines 180-184,
218-221 and 297-298.

Section 6 requires “mean PC1 score (default), mean PC2 score” without the
rank qualification used in main 7.3. Two eligible collections can consist of
identical shapes, yielding rank zero, or vary on only one component. The
existing `shape_pca` deliberately returns zero or one score columns in those
cases. Synthetic matrices confirmed ranks 0 and 1. Collection eligibility
alone cannot make these map selections available. Section 8 also requires
the ratio of largest to smallest variance when the denominator can be zero,
including the all-zero case. Section 10 explicitly includes collections with
no within-collection variation.

Repair: state which selectors are disabled/omitted, the fallback map quantity
when the selected PC disappears, the unavailable-field CSV representation,
and the text for positive/zero and zero/zero variance ratios. Keep valid
non-PC results available. Define report section inclusion for these states.

Acceptance: exercise rank 0, 1 and 2; all-zero within variance; one zero and
one positive variance; and a selected PC disappearing after an edit. Maps,
downloads and reports must return the specified results without indexing a
missing component or displaying NaN as an interpretable ratio.

### F4. Existing group displays do not implement the new collection identity

Severity: major. Classification: wrong.

Location: collections specification lines 61-63, 269-270 and 283-285;
`bifaceAnalysis/R/store.R:88,193,221`; `bifaceAnalysis/app.R:507-510,567,586,655-658,748-752`;
`bifaceAnalysis/R/plots.R:20-24`.

The claim that group colouring “already colours the morphospace by collection”
does not hold for all accepted outline files. Step 2 saving and step 3 editing
trim group labels, but the CSV loader preserves text and `groups()` returns
that text unchanged. Thus `A` and ` A ` are different existing PCA colours
but one proposed collection. A spaces-only loaded label counts as nonempty
in existing plots but as unlabelled in the extension. The existing scores CSV
also exports the untrimmed labels. Literal `NA` is a valid text label in the
existing loader and must not silently become missing in the new reader.

The palette diagnosis is correct: nine levels produce only eight colours.
Changing the shared helper to a discrete viridis fallback can cover every
inspected group-coloured figure, including the outline grid's `fill` aesthetic,
the stack, measurements and morphospace. Merely applying the fallback to new
figures leaves the advertised existing PCA view repeating colours. Such a
shared-helper change is an exception to the promises that existing tabs are
unchanged. Per-figure level sets can also assign a different colour to the
same collection after ineligible groups are omitted.

Repair: define a canonical collection-label vector and explicitly decide
whether existing displays/score exports adopt it or are described as showing
raw groups. State the narrow exception to unchanged existing tabs if
normalization and the palette helper are changed. Cover both `colour` and
`fill`, with a shared named palette if cross-figure colour identity is intended.

Acceptance: load labels `A`, ` A `, spaces-only, empty and literal `NA`, then
edit a label and verify membership, matching, template, PCA legend and exports.
For 9 groups inspect both colour and fill scales, not just whether plots build.
Verify that all existing numerical analysis reactives remain unchanged.

### F5. Unlabelled specimens can change automatic H, not only alignment

Severity: major. Classification: wrong.

Location: collections specification section 4.1, lines 70-75;
`bifaceAnalysis/app.R:624-638`.

The required caution says an unlabelled specimen can change results “slightly
through the alignment.” The code computes power and automatic H from all
aligned specimens before computing `coefs()`. Adding an excluded specimen can
change the coefficient dimension and the PCA substantially. In the independent
synthetic example, the same helpers selected H = 7 for six elongated shapes
and H = 17 after adding one unlabelled lobed shape. PC axes also depend directly
on all rows even if alignment and H were held fixed.

Repair: retain the shared full-data analysis, but name alignment, automatic
harmonic selection and PCA as sources of dependence. Remove the bound implied
by “slightly.” Explain that exports and distances belong to the current full
analysis/settings.

Acceptance: a server test adds/removes an excluded specimen that changes
automatic H, verifies reuse of the full coefficient matrix and PCA, and checks
that the caution describes those dependencies. Do not refit only the eligible
collections.

### F6. The unbiasedness claim conflates the untruncated estimate with the displayed quantity

Severity: major. Classification: wrong.

Location: collections specification section 5.3, lines 142-151; section 5.4,
lines 171-176; section 10, lines 343-347.

The displayed formula is `sqrt(max(0, d^2 - v_g/n_g - v_h/n_h))`, followed by
“The quantity under the square root is an unbiased estimate.” Only the
untruncated quantity Q can have that property. The clipped value `max(0,Q)`
is biased, even before taking its square root. In 4000 independent normal
draws with equal population means, mean Q was 0.01076448 (MCSE 0.02620148),
whereas mean `max(0,Q)` was 0.65458671 (MCSE 0.01724860). Mean displayed
distance was 0.50103134 despite a population distance of zero.

There is also a qualification missing from the population claims. The
derivation assumes independent coefficient vectors in a fixed representation,
independence across collections, and finite second moments. Independence of
the original specimens does not by itself establish that assumption after
jointly estimated GPA rotations and data-selected H. These operations are
present in the scoped code. Fixed-X reassignment references remain valid as
descriptions, but exact repeated-sampling unbiasedness of the app's transformed
data is not established by the stated assumption alone.

Repair: name Q separately, state its idealized assumptions, and describe both
clipping and square-root transformation. Say `v_g/n_g` estimates the mean's
expected squared error, rather than equating a random estimate to that
expectation. State that dc = 0 means the estimated correction equals or exceeds
the observed squared distance; it does not establish identical population
means or rule out a real difference. Qualify the population interpretation of
variance/correction for fitted preprocessing, or demonstrate it with a separate
derivation and simulation.

Acceptance: test Q independently of clipping, verify exact hand-computed
values, and confirm that text never calls `max(0,Q)` unbiased. The Gaussian
simulation in this audit supports the idealized Q identity, not unbiasedness
of the whole shape-processing pipeline.

### F7. File and export contracts leave consequential choices unresolved

Severity: major. Classification: unsupported.

Location: collections specification section 4.2, lines 83-102; section 7,
lines 272-276; section 10, lines 364-386.

Case/space-insensitive headers allow two columns to normalize to the same
required name, for example `latitude` and ` Latitude `. Neither rejection nor
column precedence is specified. Required-column presence alone cannot decide
which coordinate to use. The spec describes preserving old locations on
failure but does not explicitly define atomic full replacement on success,
or rematching after group edits. The template implies all nonempty collection
names, whereas eligibility controls analysis; that distinction should be
explicit. The file-size cross-reference does not give a number in main 6.5;
the actual app request limit is 40 * 1024^2 bytes (`app.R:10`).

The summary CSV points to section 5.5, which does not explicitly list a
collection identifier column, exact headers, order, or how absent PCs and
unavailable mm values are represented. “Empty files with a header” is not an
implementable fixed-schema contract when these choices are unresolved. The
distance file also needs a stated missing-km representation. These choices
matter for downstream models and for the promised exact-column tests.

Repair: reject duplicate normalized required headers, specify text-preserving
name parsing and duplicate checks after trimming, state successful replacement
and live exact-name rematching, make template membership explicit, give the
byte limit, and publish the two complete CSV schemas with missing-value rules.

Acceptance: fixtures with conflicting normalized headers and literal text
labels; valid replacement omitting an old row; failed replacement; group
rename/removal and reappearance; an ineligible named collection in the template;
and exact headers/values for unavailable PCs, mixed scaling, missing locations
and header-only downloads.

### F8. Several acceptance checks can pass without correct statistics or figures

Severity: major. Classification: not_demonstrated.

Location: collections specification section 10, particularly lines 343-350,
357-360 and 374-386.

The correction test specifies 4000 draws and a three-MCSE tolerance, but not
the distribution, covariance, simulation seed, or MCSE calculation. Thus it
is not reproducible from the text. “Quantity under the square root” also
inherits F6. The leave-one-out test only requires ranges to “bracket” hand
values, which lets an excessively wide range pass. A ggplot object that
builds without warnings can have missing points, incorrect colours, wrong
extents or the wrong values. The tests do not explicitly check pair-reference
values or summary/mean-shape calculations. A counts-and-r=1 check alone does
not verify the full matrix-reference distribution. Detailed coverage is below.

Repair: fix distributions and seeds, define MCSE as sample SD divided by the
square root of the number of independent draws, give exact finite examples,
and inspect computed plot data as well as object class. Include numerical
checks for each output and branch. For example N(mu,I12), seed 118, n=5 and 8,
delta=2 and 0 gives expected raw squared distances 7.9 and 3.9 and expected Q
4 and 0. Test that an intentionally omitted correction fails.

Acceptance: exact LOO endpoints, golden enumerated references, tested summary
values and inverse reconstructions, and plot data/extent/scale assertions.
Tests should fail for deliberately wrong divisors, all-specimen rather than
within-pair reassignment, omitted correction, empty plots and repeated palettes.

### F9. The five-location justification and reference cautions need correction

Severity: minor. Classification: unsupported.

Location: collections specification section 5.6, lines 188-209; section 8,
lines 299-315.

“Three locations give 3 pairs and four give 6, which cannot show a
relationship” does not justify a statistical boundary at five. A correlation
can be computed at three or four nondegenerate locations; five does not make
the ten pair distances independent. Keep five as the chosen classroom display
threshold, but do not present it as a mathematical threshold for relationships.

The text calls the procedure a reference and withholds p-values, which is
consistent with the descriptive decision. It does not explicitly require the
new reassignment ranges to be labelled as distributions of rearranged fixed
data rather than uncertainty intervals, as main 7.6 does. Nor does it require
the essential spatial caution that many plotted pairs share the same small
set of collections. A single collection can affect many points. The map also
needs to say that a mapped mean does not display within-collection variation
or its sampling uncertainty. Rounding coordinates changes the geographic
distances, not just disclosure precision.

Repair: call five an operational threshold and add those explicit cautions,
including that zero corrected distance is not established equality (F6).

Acceptance: tab/report text explains the number of collections as well as
pairs, dependent pairs, fixed-data reference meaning, and limits of point-map
colours. No range is described as a confidence interval or used to declare
a real population or spatial difference.

## Mathematical checks that held

For fixed coefficient rows let `y_i = x_i - c` and `SST = sum |y_i|^2 > 0`.
For a simple random subset of size n from the N rows,

```
E |mean(y_subset)|^2 = (N - n) SST / [n N (N - 1)].
E SSB = sum_g (N - n_g) SST / [N (N - 1)]
      = (G - 1) SST / (N - 1).
```

SST is fixed, so `E(SSB/SST) = (G-1)/(N-1)` for unequal as well as equal
sizes, with no distributional assumption on the rows. For three collections
of five this is 0.142857142857, correctly rounded to 0.14. With 50,000
reassignments of a fixed nonnormal matrix, means were 0.142369226020 (MCSE
0.000177353200) for sizes 5,5,5 and 0.143235437440 (MCSE 0.000186478966) for
3,5,7. Exhausting all 1680 labelled assignments for sizes 3,3,3 gave exactly
0.25. Partition sums agreed to rounding error, at most 1.14e-13 in these
checks, and shares agreed with pooled multivariate linear-model residual sums
of squares. `summary(lm(...))` does not directly provide an ANOVA sums-of-squares
table; residual sums and total sums, or `anova`, give the intended oracle.

For mutually independent coefficient vectors with collection mean mu_g and
covariance Sigma_g in a fixed coordinate system,

```
E v_g = trace(Sigma_g)
E |m_g - m_h|^2 = |mu_g - mu_h|^2
                   + trace(Sigma_g)/n_g + trace(Sigma_h)/n_h.
```

Subtracting the two unbiased variance estimates divided by sample size gives
an unbiased Q. Equal covariances and normality are not required for this
identity. The expectation of the raw norm is at least the population norm by
convexity; it need not be strictly larger for every possible distribution.
Normal simulations with delta=2 yielded mean Q 3.979827006 (MCSE 0.044409805)
and raw squared distance 7.877283159 (MCSE 0.043615152), against expectations
4 and 7.9. The corresponding zero-distance figures are recorded in F6.
Omitting correction clearly fails these chosen fixtures. This does not
establish unbiasedness after jointly estimated preprocessing.

The n-1 divisor makes expected variance independent of n under the same
fixed-distribution independent-sampling assumptions. Normal simulations with
12 unit-variance variables and n=3,5,8,20 gave mean v 11.96428025,
12.01501977, 11.97607929 and 11.98792663 against expectation 12. Their MCSEs
were 0.02442061, 0.01732827, 0.01312426 and 0.00795474.

The leave-one-out range is correctly called a sensitivity range. In fact,
for n>2, the average of the n leave-one-out variances equals v, so their
minimum and maximum enclose v, apart from numerical error. For scalar rows
0,1,2,7, v=29/3 and the leave-one-out values are 31/3,13,43/3,1. The exact
range is [1,43/3]. It is not a confidence or credible interval and its width
does not supply a calibrated precision measure.

The pair reference is matched to the displayed raw distance, as required by
lines 153-156. It is a valid conditional arbitrary-partition reference.
Interpreting it as a population test would require exchangeability, stronger
than merely equal means when dispersions differ. The specification does not
require that inferential interpretation.

The haversine formula is correct on the chosen sphere with angles in radians.
That sphere approximates WGS84 positions; it is not an ellipsoidal geodesic.
The one-degree, identity and meridian-crossing values hold. F1 identifies the
single wrong numerical oracle.

Reassigning whole location rows among located collections is well defined.
There are 5! = 120 labelled assignments, including the identity once. Using
all 120 on a fixed example with shape distance equal to geographic distance
gave observed and maximum r=1, mean reference approximately 0, and quantiles
-0.511752672815 and 0.743248859451. Coincident rows are not a problem: with
two at one position and three at another, 120 assignments produce 10 distinct
coordinate arrangements, each 12 times. Keep that multiplicity or use
equivalent weights. There can be fewer distinct r values still. At least two
positions with G>=5 guarantees nonconstant geographic pair distances, but
corrected shape distances can all be zero. The latter guard is specified in
5.6 but lacks an explicit test.

The existing code can supply all required shapes, coefficients, measurements
and PCA scores without changing numerical results. `shapes()` applies the
symmetric option before `aligned()`, `power()` uses all specimens and 30
harmonics, `coefs()` uses the chosen H, and `pca()` consumes that matrix.
`groups()` does not enter these computations. Collection rows can be selected
after these reactives. `efa_vector` uses a1..aH,b1..bH,c1..cH,d1..dH;
`efa_inverse` uses the same layout and defaults x0=y0=0. Its linearity means
the reconstructed mean vector is the mean of the zero-offset reconstructions.
PCA means and measurements need to be passed to collection summaries from the
existing reactives; X, labels and locations alone are not enough to reconstruct
original measured lengths or the all-specimen PCA axes.

`region_rv` provides the required seed restoration pattern, including when
there was no `.Random.seed`. The same plotting functions can be reused for
screen and PDF. The report helper allows continuation pages, so “one page for
each figure” must be understood as one figure section with any necessary text
continuations, not a guarantee of exactly one page regardless of text length.
Section 3's six-figure count can be met by making the share output numeric;
if it is plotted separately, the inventory becomes seven. Freeze that inventory
before asserting page counts. The existing sidebar can accept two downloads;
the statement that it is “unchanged” needs that obvious explicit exception.

## Section 10 coverage assessment

| Proposed test | Assessment |
| --- | --- |
| Partition identity and pooled lm share | Correct for nonconstant X. Use an independent pooled sum-of-squares oracle, not an average of per-column R-squared values. |
| Replicated collections give 0 | Correct when the replicated row set varies. If all rows are identical, SST=0 and 5.2 requires unavailable, not zero. Specify the fixture. |
| No within variation gives 1 | Correct when means differ enough to pass the numerical-zero guard. |
| 999-draw mean within 0.02 | Target correct; fixture distribution and seed need fixing to avoid an unspecified stochastic test. No collection-effect assumption is necessary for the conditional identity. |
| Repeatability and RNG state | Useful, including absence of an initial seed. Explicitly apply to pair and geographic random references too. |
| Raw Euclidean mean distance | Correct if “row means” means means across specimen rows yielding coefficient vectors; not R `rowMeans(X)` for each specimen. Specify a numeric fixture. |
| 4000-pair correction simulation | Correct target only for Q before clipping and under the independent-vector model; distribution/seed/MCSE are missing (F6/F8). |
| Negative correction gives zero | Correct and useful. Also test a positive and exact-zero Q. |
| LOO range | Zero example correct but weak; assert exact endpoints for a specified nonconstant example. |
| Four haversine values | Three correct, pole oracle wrong (F1); supply tolerances for all rounded values. |
| Correlation of lower triangles | Correct, provided both matrices share collection order and omit the diagonal. |
| 120 assignments and identity | Correct for labelled rows. Assert each permutation or an independently enumerated distribution, not only length 120. |
| Dshape=Dgeo gives r=1 and maximum | Correct for nonconstant distances. Insufficient alone to check shuffling of both matching rows and columns or reference quantiles. |
| All positions coincide | Correct unavailable case. Distinct positions with all corrected distances zero is missing. |
| Reader acceptance/rejections | Covers all listed rejection classes in principle. Missing decisions in F7 prevent complete fixtures. Test boundary row counts and duplicates after trimming. |
| Collection/location thresholds | Correct intended outputs; also check 0/1 located and two located at a single position, all-PC-unavailable data, and eligible versus noneligible scaling. |
| Figures build for 2,3,6,9 | Smoke coverage only. Geography figure requires >=5, so function-level synthetic inputs must be distinguished from server withholding. Need data/aesthetic/extent checks (F8). |
| Report page counts | No exact fixture counts or section/text allocation specified. Existing report can add continuation pages. Check section presence/absence and text as well as counts. |
| Existing numerics unchanged after labels | Supported by current dependency graph; numerical equality is the right assertion. Labels necessarily change legends/group tables. |
| Load/replace/clear and downloads | Useful behavior tests once F7's replacement/schema decisions are fixed. Explicitly verify km for located pairs even below the five-collection correlation threshold. |
| Existing suite and manual checklist | Appropriate future requirements. Neither was run as part of this specification audit. |

Rules/statistics without explicit numerical or behavioral assertions include:
mean coefficient reconstruction and zero offsets; within variance's n-1
divisor on nonconstant input; pair-specific reassignment pools, sizes and
quantiles; each summary mean and its identity; reuse of all-specimen PCA axes,
alignment, H and symmetry; the all-eligible-specimens scale gate; numerical
zero for the partition; geographic references at the 7/8 boundary and partial
coincidence; no geographic statistic for constant shape distances; map
minimum extent, padding, polar behavior, bar length/placement and boundary
layer selection at the US box and 3-degree threshold; actual plotted values,
labels and categorical colour uniqueness; rank-dependent selectors; header-only
schemas and optional-field values; live name rematching; report cautions and
omission of the map only; and prevention of recomputation when only the mapped
quantity changes. Generic smoke tests do not verify these contracts.

## Not examined and limits

- `maps` databases, real map rendering and their actual low-resolution
  boundaries: package absent. No installation attempted. Only the installed
  `ggplot2::map_data` implementation, `CoordQuickmap` aspect method and a
  synthetic polygon zoom were checked. `map_data` delegates to `maps::map`
  and preserves groups/order, which supports the proposed architecture but
  is not a rendered test of world/state/county.
- Project-root `CLAUDE.md`: absent as described above.
- Existing project tests, Shiny browser/server execution, PDF generation,
  deployment and performance: not run. This is a pre-implementation
  specification audit; the extension is not available to exercise. Section
  11's under-two-second target remains an acceptance target, not a demonstrated
  result. Six collections imply 15 pairs and 720 exhaustive spatial
  assignments; 5040 is the maximum at seven collections.
- Prior audit reports/ledgers, image data and unlisted source files: excluded
  by the brief. The supplied main-spec build notes were read and respected.
- Flagged literature attributions: no permitted source PDFs were supplied or
  opened. The author already identifies them as unverified. This report makes
  no literature-based finding or claim.
- Sampling independence, archaeological representativeness, benefit of
  symmetrizing actual collections, and applicability of the correction to
  fitted GPA in real datasets: not established by synthetic fixed-vector
  checks. Repeated sample-image uploads exercise UI only, as the spec says.

## Questions for the authors

1. What latitude/extent domain should the static map support, and what should
   be shown outside it? Where should its scale bar be valid?
2. Should raw group labels in existing plots and score exports be canonicalized
   to collection labels, and should a collection keep one colour across views?
3. What map quantity replaces a PC that does not exist, and what exact CSV
   schemas and unavailable-value representations should downstream users rely on?
4. Does successful locations upload replace the complete table, with group
   edits rematched by current name only? How should ambiguous normalized
   headers be rejected?
5. Will the population claims be limited to independent coefficient vectors
   in a fixed representation, with fitted preprocessing described as an
   approximation?
6. Is the share a numeric output, fixing the intended six figures, and should
   report acceptance count sections while allowing text continuation pages?

## Files written

- `codex/audits/2026-11-17_spec_2_report.md`
- `codex/evidence/2026-11-17_spec_2/checks.R`
- `codex/evidence/2026-11-17_spec_2/checks.log`
