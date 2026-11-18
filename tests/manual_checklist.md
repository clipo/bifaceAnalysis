# Manual browser checklist

These checks need a browser and are not automated (shinytest2 cannot be
installed on the development machine; see the specification, section 10).
Run `shiny::runApp("bifaceAnalysis")`, or open a deployed copy, and record the
date, browser, and result of each item.

Status: COMPLETED on 2026-11-17 in Google Chrome on macOS, against the test
deployment https://lipolab.shinyapps.io/bifaceAnalysis-test/ (the browser used
for checking cannot load the app from the local server). The checks were made
through browser automation: pointer clicks and drags on the page, the page's
own file-upload controls, and scripted changes to form controls. The first
pass ran on commit daa5cf1 and found four defects, listed below the table.
They were repaired, the app was deployed again to the test address, and the
affected items were repeated. Results marked "after repair" were observed on
the repaired version. Limits of the run are stated in the Result column.

| # | Check | Expected | Result |
|---|-------|----------|--------|
| 1 | Open the app | Step 1 page appears with no error banner | Pass. |
| 2 | Click "Run the demonstration" | Step 4 opens with 4 specimens; every tab shows its figures and a "How to read this" note | Pass (after repair). All seven tabs showed figures and notes. |
| 3 | Step 2, sample image 1 | Orange outline follows the artifact; tip and base crosses at the two ends; the mask panel shows the artifact in white; no warning box | Pass. |
| 4 | Mode "Crop", drag a box around the artifact | Dashed yellow box remains; outline unchanged; the drag selection clears | Pass. |
| 5 | Drag a box that cuts the artifact | A warning says the artifact may be cut off and asks for the crop to be redrawn; no error | Pass (after repair). First pass: failed, see defect 2. The outline that touches the edge is not selected in this case, so the expected text was revised. |
| 6 | "Clear crop" | Box disappears | Pass. |
| 7 | Mode "Set scale", click the two ends of 10 scale-bar divisions, enter the distance | Green line between the clicks; scale text shows pixels per mm | Pass. 9.43 pixels per mm for 10 mm. |
| 8 | Tick the shared-scale box, go to the next image | Scale text says the shared scale is in use | Pass. |
| 8a | Press "Clear scale for this image" there | Scale text says the scale is cleared and measurements will be in pixels; the first image keeps its scale | Pass. |
| 9 | Mode "Set tip and base", click the base end then the tip end | After the first click a yellow circle marks it and a message asks for the base; after the second the preview is turned upside down | Pass. |
| 9a | Click twice at the same place on the outline | An error notice appears and the previous tip and base remain | Pass (after repair). First pass: not testable, see defect 1. |
| 10 | "Use automatic tip and base" | Preview returns to the automatic orientation | Pass. |
| 11 | Change tip direction to "Right" | Tip and base markers exchange ends | Pass. Changed through the page's select control by script. |
| 12 | Untick automatic threshold and move the slider to both extremes | Outline changes or a message appears; no error | Pass. 0.98 gave three warning boxes and no error. 0.02 gave a normal outline. |
| 13 | Save each image | Name field updates for each image; the count increases; after the last image a notice points to step 3 | Pass. |
| 14 | Step 3: rename, label, and remove a specimen | Grid and table update; a duplicate name receives a numeric suffix | Pass. A second "artifact1" became "artifact1 (2)". |
| 15 | Download outlines, reload the page, load the file | The same specimens and measurements return | Partly checked. The download link returned an outline file of 219,525 bytes. An outline file written by `specimens_to_csv` on the development machine was loaded with three specimens present: four specimens resulted, with the same measurements as before. The file downloaded by the browser was not itself reloaded, because the automation cannot read the browser's download folder. |
| 16 | Load a CSV that is not an outline file | A message explains the problem; existing specimens remain | Pass on two of three attempts read from the page. On the first attempt no message was found 3 seconds after the upload, and the cause was not established. Existing specimens remained in every attempt. |
| 17 | Step 4 with 1, 2, and 3 specimens | Tabs that need more specimens show a sentence saying how many are needed | Pass with 1, 2, and 3 specimens. No error output in any tab. The "How to read this" card is empty when its figure is withheld. |
| 18 | Change harmonics, clusters, colouring, symmetric option, region boundaries | Figures and notes update; no error | Pass. Harmonics set to 12 by hand, clusters set to 3, colouring by cluster, symmetric option, and boundaries were each changed. |
| 19 | Download the PDF report and each CSV | Files open; the PDF has a cover page and one page per figure with its note | Partly checked. Each download link returned a file: a PDF of 17 page objects with the `%PDF-1.4` header, and three CSV files with a header row. The files were read by script. They were not opened in a PDF viewer or spreadsheet. |
| 20 | Upload a PNG and a JPEG of your own, including one larger than 1200 pixels | Both load; the large one responds at normal speed | Pass. A 2000-pixel PNG and a JPEG loaded and were outlined within a few seconds. |
| 20a | Upload a file named like a sample image, after cropping that sample image | The upload shows no crop, scale, or landmarks; saving it adds a specimen and does not replace the sample's | Pass. The upload showed no crop box while the sample image of the same name had one, and was saved as a separate specimen. |
| 20b | Set two region boundaries to the same value | The region tab asks for three different values; other tabs are unaffected | Pass. |
| 21 | Narrow the window to phone width | Controls remain usable; no overlapping text | Pass (after repair) at a 396-pixel viewport, tested in a frame inside the page because the browser window could not be resized. The navigation collapses to a menu button, cards stack in one column, and the page has no horizontal scroll. |

## Defects found in the first pass and repaired

1. Layout shift under the pointer. Messages were shown above the image, so a
   message that appeared after the first landmark or scale click moved the
   image down by about 40 pixels before the second click. Messages are now
   shown below the image inside the same card.
2. A crop that cut through the artifact selected an inked number inside it and
   reported only that other objects were present. `segment_image` now warns
   that the artifact may be cut off when regions at least four times the size
   of the selected one touch the border in both brightness classes. A test was
   added in `tests/testthat/test-image.R`.
3. The page was fitted to the window height, which compressed the figures and
   cut off the "How to read this" text on the PCA and Regions tabs. The page
   now scrolls and figures keep their stated heights.
4. The specimen pickers in step 3 and on the Harmonics tab opened empty. They
   now select the first specimen. The Outlines note also referred to
   "coloured lines" when the lines are grey without group labels, and now
   says "thin lines".

## Not covered by this run

Browsers other than Chrome, touch input, and a real phone were not tested.
Item 15 and item 19 were checked through the download links and not through
files saved by the browser. The author's own images were not used: the uploads
were the sample images renamed and enlarged.
