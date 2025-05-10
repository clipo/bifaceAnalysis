# Biface Shape Analysis (Momocs Pipeline)

This application provides a computational framework for morphometric analysis of lithic bifaces and projectile points using the **Momocs** package in R. The analytical architecture incorporates elliptical Fourier analysis, principal component analysis, and region-based morphological assessment to quantify shape variation in archaeological assemblages.

The software implementation includes:
1. **Shiny App** – An interactive graphical interface facilitating image processing, outline extraction, and multivariate analysis
2. **R Script** – A command-line implementation supporting automated batch processing

---

## 📆 Features

- Upload or utilize sample images
- Interactive cropping functionality
- **Adaptive thresholding with real-time visualization**
- **Artifact/background contrast handling (light-on-dark or dark-on-light)**
- **Configurable noise elimination parameters**
- **Contour detection with largest-area prioritization**
- Outline extraction and interpolation
- **Three distinct region detection methodologies:**
  - **Geometric detection (morphology-based)**
  - **Proportional detection (percentage-based)**
  - **Manual landmark identification**
- Elliptical Fourier Analysis (EFA)
- Principal Component Analysis (PCA)
- Hierarchical clustering
- **Thin-plate spline deformation visualization**
- **Integration/modularity analysis of morphological regions**
- Publication-quality visualizations:
  - Outline arrays
  - **Superimposed alignment visualization**
  - **EFA harmonic power spectral analysis**
  - PCA morphospace projection
  - **Extremal shape variation visualization**
  - **Regional variance heatmaps**
  - **Inter-region correlation matrices**
- Exportable analytical products:
  - PCA coordinates (`.csv`)
  - Extracted outlines (`.rds`)
  - **Comprehensive analytical report (`.pdf`)**

---

## 📁 Repository Structure
```
bifaceAnalysis/
├── app.R                 # Shiny application implementation
├── analysis_script.R     # Command-line implementation
├── www/                  # Sample image repository
│   ├── artifact1.jpg
│   ├── artifact2.jpg
│   └── ...
└── README.md
```

---

## 🚀 Application Deployment

### Installation Prerequisites
1. Clone or download this repository
2. Install required R dependencies:
```r
install.packages(c("shiny", "Momocs", "imager", "dplyr", "ggplot2"))
```

### Launching the Application
1. Open RStudio  
2. Set working directory to the repository location:
```r
setwd("path/to/bifaceAnalysis")
```
3. Execute the application:
```r
shiny::runApp()
```
The application will initialize in your browser or RStudio's integrated viewer.

---

## 🖼️ Analytical Workflow

### Image Acquisition
1. **Sample Images**: The application automatically detects images in the `www/` directory
2. **Custom Dataset**: Select "Upload My Own" to import JPEG images from local storage

### Image Processing
1. **Select an image** from the available inventory
2. **Define crop region** by selecting diagonal corners of the region of interest
3. **Configure processing parameters**:
   - Artifact contrast type (darker or lighter than background)
   - Threshold value for binary conversion
   - Noise elimination parameters
4. **Initiate processing** to extract the artifact outline

### Region Detection Methodology
The application implements three distinct approaches to morphological region identification, each with specific analytical advantages:

#### 1. Geometric Detection (Shape-based)
This method employs computational geometry algorithms to identify morphologically significant regions based on intrinsic shape characteristics. The implementation:
- Calculates curvature values along the outline
- Detects significant inflection points through adaptive thresholding
- Identifies critical morphological landmarks (tip, shoulders, base)
- Defines regions based on these natural boundaries

**Configuration Parameters:**
- **Curvature Sensitivity**: Controls threshold for identifying significant inflection points (higher values detect only more pronounced features)
- **Curvature Smoothing**: Determines the degree of noise reduction in curvature calculation (higher values produce smoother, more generalized results)

#### 2. Proportional Detection (Percentage-based)
This approach divides the outline into regions based on user-defined proportions:
- **Tip Proportion**: Percentage allocated to the distal extremity
- **Blade Proportion**: Percentage allocated to the cutting edges
- **Shoulder Proportion**: Percentage allocated to the transitional zone
- **Base Proportion**: Percentage allocated to the proximal element

The percentages are applied relative to the specified orientation.

#### 3. Manual Key Point Selection
This method provides direct user control through interactive selection of critical morphological landmarks:
- **Tip**: The distal extremity
- **Blade Junctions** (Left/Right): Transitions between tip and blade edges
- **Shoulder Junctions** (Left/Right): Locations where blade contracts to form base
- **Base**: The proximal extremity

### Region Detection Selection Guide

| Method | Optimal Applications | Advantages | Limitations |
|--------|----------------------|------------|-------------|
| **Geometric** | Typologically consistent assemblages with well-defined morphological features | • Automated identification of natural morphological transitions<br>• Adaptability to diverse artifact shapes<br>• Reproducible results | • Sensitivity to outline noise<br>• Requires parameter calibration |
| **Proportional** | Assemblages with consistent proportional relationships | • Straightforward implementation<br>• Ensures consistent region allocation<br>• Effective for standardized typologies | • Imposes arbitrary boundaries<br>• Lacks sensitivity to morphological features |
| **Manual** | Highly variable assemblages requiring expert assessment | • Maximum control over region definition<br>• Allows incorporation of expert knowledge<br>• Accommodates irregular specimens | • Time-intensive<br>• Potential subjectivity<br>• Requires domain expertise |

### Executing Analysis
1. Process all artifacts using the "Next Image" function
2. Select "Run Analysis" to perform:
   - Elliptical Fourier transformation
   - Principal Component decomposition
   - Hierarchical clustering
   - Region-specific morphological analysis
3. Examine results in the "Results" tab
4. Generate comprehensive report via "Download Graphics (PDF)"

---

## 📋 Implementing Geometric Region Detection

The geometric detection approach represents a methodological advancement in projectile point analysis, identifying regions through computational assessment of morphological characteristics rather than arbitrary proportions.

### Theoretical Framework
The methodology employs curvature analysis to detect significant inflection points along artifact outlines. These inflection points correspond to archaeologically meaningful morphological transitions between functional regions. The algorithm:

1. Calculates local curvature at each outline point using vector-based angular assessment
2. Applies adaptive thresholding to identify statistically significant inflection points
3. Classifies these points into morphological landmarks (tip, shoulders, base)
4. Defines regions based on these natural boundaries and their geometric relationships

### Implementation Guide
To utilize geometric region detection:

1. Process an image to extract the outline
2. Select "Geometric (Shape-based)" from the Region Detection Method options
3. Adjust sensitivity parameters if needed:
   - **Curvature Sensitivity**: Controls detection threshold (1.0-2.0 optimal for most projectile points)
   - **Curvature Smoothing**: Controls noise reduction (3-7 optimal for most specimens)
4. The outline display will visualize curvature values using a color gradient (blue=low, red=high)
5. Execute analysis to apply geometric region detection

### Parameter Optimization
For optimal results when analyzing diverse assemblages:
- **Well-defined points** (e.g., Clovis, Eden): Use higher sensitivity (1.8-2.0) with moderate smoothing (5-7)
- **Transitional forms** (e.g., mixed assemblages): Use moderate sensitivity (1.3-1.7) with moderate smoothing (5-7)
- **Irregular specimens**: Use lower sensitivity (1.0-1.2) with higher smoothing (7-11)

---

## 📝 Manual Key Point Selection

For precise control over region identification:

1. Process an image to extract the outline
2. Select "Manual key points" from Region Detection Method options
3. Select critical morphological landmarks in sequence:
   - **Tip**: The distal extremity
   - **Blade Junctions**: Points where the tip transitions to blade edges
   - **Shoulder Junctions**: Points where blade contracts to form base
   - **Base**: The proximal/hafting element
4. Execute analysis to apply manual region definitions

### Recommendations for Optimal Results
- Select landmarks in clockwise or counter-clockwise sequence
- Place shoulder points at positions with notable changes in outline curvature
- Ensure all six required landmarks are selected prior to analysis
- Maintain consistent landmark selection criteria across assemblages

---

## 🔍 Region Analysis Interpretation

The application divides projectile points into four morphologically and functionally distinct regions:

1. **Tip Region**: The distal penetrating element
2. **Blade Region**: The primary cutting edges 
3. **Shoulder Region**: The transitional zone where blade contracts
4. **Base Region**: The proximal hafting element

This four-region analytical framework facilitates:
- Nuanced assessment of morphological variation
- Identification of manufacturing techniques and design elements
- Detection of functional adaptations
- Precise integration/modularity analysis between regions

The analytical outputs include:
- Region-specific variance quantification
- Inter-region correlation analysis
- Modularity index calculation
- Regional contributions to principal components

---

## ⚙️ Command-Line Implementation

The standalone script provides equivalent functionality without graphical interface:

```r
# Configuration
threshold_value <- 0.1
force_rebuild_outlines <- TRUE  # Set to FALSE to use cached outlines
region_method <- "geometric"    # Options: "geometric", "proportional", "manual"

# Execution
source("analysis_script.R")
```

This implementation:
- Extracts outlines from all `.jpg` images in `/images/` directory
- Caches outlines to `cached_outlines.rds`
- Generates analytical outputs in the working directory

---

## 📝 Methodological Recommendations

### Image Preparation
- Utilize high-contrast imagery with distinct silhouettes
- Ensure consistent orientation across specimens
- Eliminate background artifacts that may affect contour detection
- Standardize lighting and camera positioning

### Analytical Considerations
- Process minimum of 10 specimens for robust PCA results
- Utilize 20 harmonics for most projectile point typologies
- Implement geometric detection for morphologically diverse assemblages
- Apply manual key point selection for irregular specimens
- Ensure consistent orientation when comparing diverse artifact types

### Results Interpretation
- PCA projections visualize primary axes of shape variation
- Region-specific variance quantifies morphological plasticity in functional elements
- Inter-region correlation matrices reveal covarying morphological elements
- Modularity indices indicate functional integration/independence between regions

---

## 📄 License
This project is licensed under the Creative Commons Attribution 4.0 International (CC BY 4.0).

---

## 🤩 Dependencies
```r
install.packages(c("shiny", "Momocs", "imager", "dplyr", "ggplot2"))
```

---

## ✨ Attribution
Developed utilizing the `Momocs` package by Vincent Bonhomme.
Shiny implementation and analytical framework by Carl Lipo, Binghamton University (2025).