# Manual browser checklist

These checks need a browser and are not automated (shinytest2 cannot be
installed on the development machine; see the specification, section 10).
Run `shiny::runApp("bifaceAnalysis")` and record the date, browser, and result
of each item.

Status: NOT YET COMPLETED. During the build on 2026-11-17 the browser used for
checking could not load the app from the local server. Server tests in
`tests/testthat/test-server.R` cover the handling of crop, scale, and landmark
inputs, image identity on upload, scale clearing, and loading of outline
files, with inputs supplied by the test and not by a browser. Unchecked:
pointer behaviour, the file-upload control, updates of the text fields in the
page, rendering, and layout.

| # | Check | Expected | Result |
|---|-------|----------|--------|
| 1 | Open the app | Step 1 page appears with no error banner | |
| 2 | Click "Run the demonstration" | Step 4 opens with 4 specimens; every tab shows its figures and a "How to read this" note | |
| 3 | Step 2, sample image 1 | Orange outline follows the artifact; tip and base crosses at the two ends; the mask panel shows the artifact in white; no warning box | |
| 4 | Mode "Crop", drag a box around the artifact | Dashed yellow box remains; outline unchanged; the drag selection clears | |
| 5 | Drag a box that cuts the artifact | A warning that the outline touches the edge appears; no error | |
| 6 | "Clear crop" | Box disappears | |
| 7 | Mode "Set scale", click the two ends of 10 scale-bar divisions, enter the distance | Green line between the clicks; scale text shows pixels per mm | |
| 8 | Tick the shared-scale box, go to the next image | Scale text says the shared scale is in use | |
| 8a | Press "Clear scale for this image" there | Scale text says the scale is cleared and measurements will be in pixels; the first image keeps its scale | |
| 9 | Mode "Set tip and base", click the base end then the tip end | After the first click a yellow circle marks it and a message asks for the base; after the second the preview is turned upside down | |
| 9a | Click twice at the same place on the outline | An error notice appears and the previous tip and base remain | |
| 10 | "Use automatic tip and base" | Preview returns to the automatic orientation | |
| 11 | Change tip direction to "Right" | Tip and base markers exchange ends | |
| 12 | Untick automatic threshold and move the slider to both extremes | Outline changes or a message appears; no error | |
| 13 | Save each image | Name field updates for each image; the count increases; after the last image a notice points to step 3 | |
| 14 | Step 3: rename, label, and remove a specimen | Grid and table update; a duplicate name receives a numeric suffix | |
| 15 | Download outlines, reload the page, load the file | The same specimens and measurements return | |
| 16 | Load a CSV that is not an outline file | A message explains the problem; existing specimens remain | |
| 17 | Step 4 with 1, 2, and 3 specimens | Tabs that need more specimens show a sentence saying how many are needed | |
| 18 | Change harmonics, clusters, colouring, symmetric option, region boundaries | Figures and notes update; no error | |
| 19 | Download the PDF report and each CSV | Files open; the PDF has a cover page and one page per figure with its note | |
| 20 | Upload a PNG and a JPEG of your own, including one larger than 1200 pixels | Both load; the large one responds at normal speed | |
| 20a | Upload a file named like a sample image, after cropping that sample image | The upload shows no crop, scale, or landmarks; saving it adds a specimen and does not replace the sample's | |
| 20b | Set two region boundaries to the same value | The region tab asks for three different values; other tabs are unaffected | |
| 21 | Narrow the window to phone width | Controls remain usable; no overlapping text | |
