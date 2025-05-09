# Biface Shape Analysis (Momocs Pipeline)
This project provides tools for morphometric analysis of biface or artifact outlines using the **Momocs** package in R.
There are two versions:
1. **Shiny App** – An interactive GUI that allows you to upload images, crop them, extract outlines, and run PCA/cluster analysis.
2. **R Script** – A command-line version that runs the full pipeline on a directory of images.

---

## 📆 Features
- Upload or use sample images
- Interactive cropping (click to crop region)
- **Adaptive thresholding with live preview**
- **Support for both light-on-dark and dark-on-light artifacts**
- **Configurable noise removal**
- **Automatic largest contour selection**
- Outline extraction and interpolation
- **Enhanced region analysis with tip, blade, shoulder, and base regions**
- **Manual key point selection for precise region identification**
- Elliptical Fourier Analysis (EFA)
- Principal Component Analysis (PCA)
- Hierarchical clustering
- **Thin-plate spline deformation visualization**
- **Integration/modularity analysis of shape regions**
- Publication-quality plots:
  - Outlines grid
  - **Stacked and aligned outlines**
  - **EFA harmonic power analysis (actual, not theoretical)**
  - PCA plots with shape morphospace
  - **Shape variation at PCA extremes**
  - **Regional variance heatmaps**
  - **Region correlation matrices**
- Downloadable outputs:
  - PCA results (`.csv`)
  - Extracted outlines (`.rds`)
  - **Comprehensive analysis report (`.pdf`)**

---

## 📁 Folder Structure
```
bifaceAnalysis/
├── app.R                 # Shiny app version
├── analysis_script.R     # Standalone version
├── www/                  # Folder for sample/demo images
│   ├── artifact1.jpg
│   ├── artifact2.jpg
│   └── ...
└── README.md
```

---

## 🚀 Running the Shiny App

### Installation
1. Clone or download this repository to your local machine
2. Install the required R packages if not already installed:
```r
install.packages(c("shiny", "Momocs", "imager", "dplyr", "ggplot2"))
```

### Launching the App
1. Open RStudio  
2. Set working directory to this folder:
```r
setwd("path/to/bifaceAnalysis")
```
3. Run the app:
```r
shiny::runApp()
```
The app will launch in your browser or in RStudio's viewer pane.

---

## 🖼️ Getting Started: Step-by-Step Guide

### Loading Images
1. **Use Demo Images**: The app automatically looks for images in the `www/` folder
2. **Upload Your Own**: Select "Upload My Own" and use the file browser to select JPEG images

### Image Processing
1. **Select an image** from the list (if multiple are available)
2. **Click two points** on the image to define the crop region (diagonal corners)
3. **Adjust settings** as needed:
   - Select "Darker than background" or "Lighter than background" based on your artifact
   - Adjust the threshold slider to get a clean binary image
   - Use the noise removal slider to eliminate small artifacts
4. **Click "Process Image"** to extract the outline

### Region Analysis Settings
1. Set **Tip direction** based on your artifact's orientation
2. Adjust the proportion sliders for:
   - **Tip proportion**: The pointed end of the projectile
   - **Blade proportion**: The main cutting edges
   - **Shoulder proportion**: Area where blade contracts to form the base (NEW)
   - **Base proportion**: The hafting element

### Running the Analysis
1. Process all images using the "Next Image" button
2. Click "Run Analysis" to perform:
   - Elliptical Fourier Analysis
   - Principal Component Analysis
   - Clustering
   - Region-based morphological analysis
3. View results in the "Results" tab
4. Download a comprehensive report using the "Download Graphics (PDF)" button

---

## 📋 Using Manual Key Points Selection

For precise region identification on projectile points:

1. Process an image to extract the outline
2. Switch to "Manual key points" mode in the sidebar
3. Click each button and then click on the corresponding point on the outline:
   - **Tip**: The pointed end of the projectile
   - **Left/Right Blade Junction**: The points where the tip transitions to blade
   - **Left/Right Shoulder Junction** (NEW): The points where the blade begins to contract
   - **Base**: The bottom/hafting end of the point
4. Run analysis to get accurate region-based calculations

### Tip for Better Results
- For the most accurate analysis, select key points in a clockwise or counter-clockwise order
- The shoulder points should be placed where there is a noticeable change in the outline curvature
- Ensure all six key points (tip, blade junctions, shoulder junctions, base) are selected before running analysis

---

## 🔍 Understanding Region Analysis

The app now divides projectile points into four distinct regions for enhanced morphological analysis:

1. **Tip Region**: The distal end, typically pointed
2. **Blade Region**: The main cutting edges 
3. **Shoulder Region** (NEW): The area where the blade contracts to form the base
4. **Base Region**: The proximal/hafting element

This four-region approach provides:
- More nuanced analysis of shape variation
- Better identification of design elements and manufacturing techniques
- Enhanced detection of functional differences between artifacts
- More precise integration/modularity analysis between regions

---

## ⚙️ Running the Standalone Script

Open `analysis_script.R` and set:
```r
threshold_value <- 0.1
force_rebuild_outlines <- TRUE  # or FALSE if using cached outlines
```
Then run:
```r
source("analysis_script.R")
```
This will:
- Extract outlines from all `.jpg` images in `/images/`
- Save outlines to `cached_outlines.rds`
- Generate plots and PCA results in the working directory

---

## 📝 Best Practices

### Image Preparation
- Use high-contrast images with clear silhouettes
- Ensure artifacts are properly oriented
- Remove background noise or other artifacts
- Consistent lighting and camera angles improve results

### Analysis Tips
- Process at least 10 specimens for meaningful PCA results
- Use harmonics value of 20 for most artifact types
- For detailed shape analysis, use the manual key points method
- When comparing different artifact types, ensure consistent orientation

### Interpreting Results
- The PCA plot shows major axes of shape variation
- Region variance bar plots indicate which areas have the most shape variation
- Correlation matrices show which regions vary together
- The modularity index indicates if regions change independently

---

## 📄 License
This project is licensed under the Creative Commons Attribution 4.0 International (CC BY 4.0).
You are free to:
Share — copy and redistribute the material in any medium or format
Adapt — remix, transform, and build upon the material for any purpose

Under the following terms:
Attribution — You must give appropriate credit, provide a link to the license, and indicate if changes were made.

---

## 🤩 Requirements
Install the following R packages if not already installed:
```r
install.packages(c("shiny", "Momocs", "imager", "tidyverse", "here", "ggplot2"))
```

---

## ✨ Credits
Developed using the `Momocs` package by Vincent Bonhomme.
Shiny UI and workflow by Carl Lipo, Binghamton University (2025).