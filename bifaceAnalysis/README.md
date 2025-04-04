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
- Elliptical Fourier Analysis (EFA)
- Principal Component Analysis (PCA)
- Hierarchical clustering
- **Thin-plate spline deformation visualization**
- **Integration/modularity analysis of shape regions**
- Publication-quality plots:
  - Outlines grid
  - **Stacked and aligned outlines**
  - **EFA harmonic power distribution**
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

### 1. Open RStudio  
### 2. Set working directory to this folder  
### 3. Run the app:
```r
shiny::runApp("bifaceAnalysis")
```

The app will launch in your browser.

---

## 🖼️ Using Sample Images

If no files are uploaded, the app will automatically use demo images found in the `www/` folder. You can also upload your own JPEG images.

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

## 🤩 Requirements

Install the following R packages if not already installed:
```r
install.packages(c("shiny", "Momocs", "imager", "tidyverse", "here", "ggplot2"))
```

---

## 📝 Notes

- All images should be **simple artifact silhouettes** on a white background.
- Thresholding assumes **dark object on light background**.
- Outlines are interpolated to 100 points by default.

---

## 📄 License


This project is licensed under the Creative Commons Attribution 4.0 International (CC BY 4.0).

You are free to:

Share — copy and redistribute the material in any medium or format

Adapt — remix, transform, and build upon the material for any purpose

Under the following terms:

Attribution — You must give appropriate credit, provide a link to the license, and indicate if changes were made.
---

## ✨ Credits

Developed using the `Momocs` package by Vincent Bonhomme.

Shiny UI and workflow by Carl Lipo, Binghamton University (2025).