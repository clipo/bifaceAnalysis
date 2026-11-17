---
date: 2026-11-17
mode: spec
title: Rewrite of the biface morphometrics Shiny app
scope:
  - docs/specs/2026-11-17_app_rewrite_spec.md
  - bifaceAnalysis/app.R
  - bifaceAnalysis/README.md
  - bifaceAnalysis/www
note: Design spec for a rewrite; the code under audit for the defect claims is the existing bifaceAnalysis/app.R
prior:
review:
commit: a85cb88d0e3f966e553b23e0be5958612b7a38e7
model: gpt-6-astra
effort: high
codex_version: codex-cli 0.169.0
started: 2026-11-17T14:44:06-05:00
message_sha256: f9312ffd5860246cfde1cfd98309ef7fb4f78b3c1a2ea6b4fe415d5700c0d737
session_id: 01a1eb09-bcf0-7ee0-8a0f-1d5bacd77627
sandbox: workspace-write [workdir, $TMPDIR]
finished: 2026-11-17T14:54:57-05:00
exit_status: ok
changed_outside_stem:
---

Mode: spec
Stem: 2026-11-17_spec
Audited commit: a85cb88d0e3f966e553b23e0be5958612b7a38e7
Date: 2026-11-17
Repository root: the parent of your working directory (..); the scope paths in the brief are relative to it, so read src/x.md as ../src/x.md. Your own files go under the working directory.

## Project

This repository holds a Shiny app (`bifaceAnalysis/app.R`, deployed on shinyapps.io) that teaches outline-based morphometrics of bifaces and projectile points to students. The spec proposes a rewrite: it lists defects in the current app, replaces the archived Momocs package with base-R elliptical Fourier analysis, rebuilds the region analysis on tip-and-base-anchored corresponding points, and adds measurements, group labels, exports, and interpretation text. Nothing in the spec has been implemented. The repository has no project CLAUDE.md.

## Read in this order

1. `docs/specs/2026-11-17_app_rewrite_spec.md`, in full. Section 2 makes claims about the current code; sections 6 and 7 define the new methods; section 10 defines the tests.
2. `bifaceAnalysis/app.R`, in full. The functions `custom_mshape`, `find_inflection_points`, and `define_regions` (lines 19 to 352), the `crop_btn` observer (lines 697 to 817), the `run_analysis` observer (lines 1059 to 1515), and the `harmonic_plot` and `modularity_plot` outputs (lines 1531 to 1947) are the ones the spec's defect list cites.
3. `bifaceAnalysis/www/`, the four demo JPEG images the new pipeline must segment.
4. `bifaceAnalysis/README.md`, for what the current app claims to provide.

## Check

1. For each defect D1 to D14 in spec section 2, does the cited code at the cited lines do what the spec says? For D1, confirm the column order of the coefficient matrix returned by `Momocs::efourier` from the Momocs 1.4.1 source or documentation if you can obtain it, and say so under not_examined if you cannot.
2. Are the elliptical Fourier formulas in section 7.2 correct as written (coefficients, the power definition, and the claim that no first-ellipse normalization is needed once position, size, rotation, and start point are fixed by sections 6.2, 6.3, and 7.1)? Check the formulas numerically on an ellipse with semi-axes 3 and 1 if R or Python is available.
3. Section 6.2: is the orientation and landmark procedure fully determined for every input? Look for cases where steps 2 to 5 give no answer, more than one answer, or a wrong answer: a principal axis exactly 45 degrees from the declared direction, a midline that meets the outline more than once on the base side (basal notches, barbs), a tip that is broken or blunt, and a shape wider than long. Do the four demo images present any of these cases?
4. Section 6.3: does resampling each edge to 100 arc-length-spaced points give 200 distinct points with T at index 1 and B at index 101 as stated, with no duplicated or skipped landmark? Is the stated limitation of this correspondence accurate and sufficient for the region analysis in section 7.6?
5. Section 7.1 and 7.2: is it sound to compute Fourier coefficients on Procrustes-aligned coordinates, or does that choice conflict with anything else in the spec (for example the measurements in 6.4, the symmetrized option, or the reconstruction panel)?
6. Section 7.6: is the RV coefficient with a permutation reference range computed and described correctly, and are the two stated limitations (small-sample bias, covariation induced by a common superimposition) stated accurately? Is the justification for a non-Bayesian reference consistent with the conventions in `~/.claude/CLAUDE.md`?
7. Section 7.3: is the leave-one-out cosine a valid description of axis stability when PC1 and PC2 have close eigenvalues and can exchange order between refits? If not, what should replace it?
8. Section 6.1: will steps 3 to 8 produce one closed outline for each of the four demo images, given the scale bar at the top of each image and the ink numbers on the artifacts? State what you ran or why you could not run it.
9. Section 5 and section 10: can every promised behaviour be delivered with shiny, bslib, imager, and ggplot2 alone? Name any function the design needs that those packages and base R do not supply. Is each acceptance test in section 10 specific enough to pass or fail without judgement?
10. Does any part of the spec contradict another part, or leave a behaviour undefined that an implementer would have to decide (for example what happens to saved specimens when the region boundaries, the harmonic count, or the symmetrize setting change, or how calibration interacts with downscaling and cropping)?

## Note from the author

Design spec for a rewrite; the code under audit for the defect claims is the existing bifaceAnalysis/app.R

## Milestone focus

This is a specification or design written before any implementation. Audit whether it is complete enough to build from without going back to the author, whether its stated assumptions hold for the data, methods, and tools it names, whether any decision is unjustified or contradicts another, and whether anything it promises cannot be delivered with the stated tools.

## Rules

- Read every file in scope in full from the working tree. Do not rely on summaries the authors wrote about their own work when the underlying file is in scope.
- Treat every claim as unverified until you have checked it against the files or the tools.
- Classify each finding as wrong, unsupported, or plausible but not demonstrated (not_demonstrated). Give a concrete repair and an acceptance check that would show the repair holds.
- List what you checked and found correct under held. List what you did not examine, with the reason, under not_examined. List the questions you would have asked the authors under questions.
- Put scripts and logs under evidence/2026-11-17_spec/, anything derived from project data (extracts, copies, page images, large tables) under evidence/2026-11-17_spec/data/, and your prose report at audits/2026-11-17_spec_report.md, all relative to your working directory. Files under data/ are not committed, so every value your report or a finding cites must appear in the report itself or in a tracked log. List every file you wrote under evidence_files.
- Record `git rev-parse HEAD` in the report.
- If the sandbox blocks a command, record it under not_examined rather than working around it.
- Do not report style preferences. Do not pad the list. If you find nothing, say so and return an empty findings list.
- Respond only with JSON matching the provided schema.
