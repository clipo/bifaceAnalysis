# Biface Shape Analysis (Momocs Pipeline)

This project provides tools for morphometric analysis of biface or artifact outlines using the **Momocs** package in R.

There are two versions:
1. **Shiny App** – An interactive GUI that allows you to upload images, crop them, extract outlines, and run PCA/cluster analysis.
2. **R Script** – A command-line version that runs the full pipeline on a directory of images.

---

## 📆 Features

- Upload or use sample images
- Interactive cropping (click to crop region)
- **Adaptive thresholding** with live preview
- **Support for both light-on-dark and dark-on-light artifacts**
- **Noise removal with configurable intensity**
- **Automatic selection of largest contour**
- Outline extraction and interpolation
- Elliptical Fourier Analysis (EFA)
- Principal Component Analysis (PCA)
- Clustering (hierarchical)
- Publication-quality plots:
  - PCA plots (with shape outlines or ggplot2)
  - Dendrogram
  - Outlines grid
  - **Stacked overlay of aligned outlines**
  - **Harmonic power visualization**
- Downloadable outputs:
  - PCA results (`.csv`)
  - Extracted outlines (`.rds`)
  - Combined plots (`.pdf`)

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

## 🖱️ Using the Application

1. **Select image source** - Use demo images or upload your own
2. **Configure image processing settings**:
   - Choose whether artifacts are darker or lighter than background
   - Adjust threshold value manually if needed
   - Set noise removal intensity
3. **Process images**:
   - Click two points on the image to define crop region (top-left and bottom-right)
   - Press "Process Image" button
   - Verify outline extraction in the preview
   - Use "Next Image" to proceed to the next image
4. **Run analysis**:
   - Set number of harmonics and clusters
   - Click "Run Analysis" button
   - View results in the Results tab
   - Download outputs as needed

---

## 🤩 Requirements

Install the following R packages if not already installed:
```r
install.packages(c("shiny", "Momocs", "imager", "dplyr", "ggplot2"))
```

---

## 📝 Notes

- Images can be either **dark objects on light background** or **light objects on dark background**
- The application automatically detects background/foreground contrast
- Outlines are interpolated to 100 points by default
- The results tab contains multiple visualizations including:
  - Grid of all outlines
  - Stacked and aligned outlines
  - Harmonic power distribution
  - PCA morphospace with shape deformations

---

## 📄 License

This project is licensed under the Creative Commons Attribution 4.0 International (CC BY 4.0).

You are free to:
- Share — copy and redistribute the material in any medium or format
- Adapt — remix, transform, and build upon the material for any purpose

Under the following terms:
- Attribution — You must give appropriate credit, provide a link to the license, and indicate if changes were made.

---

## ✨ Credits

Developed using the `Momocs` package by Vincent Bonhomme.
Shiny UI and workflow by Carl Lipo, Binghamton University (2025).