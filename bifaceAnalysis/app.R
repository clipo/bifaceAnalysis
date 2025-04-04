library(shiny)
library(Momocs)
library(imager)
library(dplyr)
library(ggplot2)

# This function safely extracts x, y coordinates from contours
safe_extract_contour <- function(contour) {
  tryCatch({
    df <- as.data.frame(contour)
    coords <- matrix(c(df$x, df$y), ncol = 2)
    colnames(coords) <- c("x", "y")
    return(coords)
  }, error = function(e) {
    return(NULL)
  })
}

# Custom mshape function to calculate the mean shape from a list of coordinates
custom_mshape <- function(coo_list) {
  # Check if input is a Momocs Out object
  if (inherits(coo_list, "Out")) {
    coo_list <- coo_list$coo
  }
  
  # Get number of coordinates and dimensions
  n_shapes <- length(coo_list)
  if (n_shapes == 0) return(NULL)
  
  n_points <- nrow(coo_list[[1]])
  n_dims <- ncol(coo_list[[1]])
  
  # Initialize mean shape matrix
  mean_shape <- matrix(0, nrow = n_points, ncol = n_dims)
  
  # Sum all shapes
  for (i in 1:n_shapes) {
    mean_shape <- mean_shape + coo_list[[i]]
  }
  
  # Divide by number of shapes to get mean
  mean_shape <- mean_shape / n_shapes
  
  return(mean_shape)
}

ui <- fluidPage(
  titlePanel("Projectile Point Shape Analysis"),
  
  sidebarLayout(
    sidebarPanel(
      radioButtons("image_source", "Image source:", 
                   choices = c("Demo Images", "Upload My Own"), 
                   selected = "Demo Images"),
      
      conditionalPanel(
        condition = "input.image_source == 'Upload My Own'",
        fileInput("images", "Upload JPEG Images", 
                  accept = c(".jpg", ".jpeg"), 
                  multiple = TRUE)
      ),
      
      h4("Image Processing"),
      radioButtons("artifact_mode", "Artifact appearance:", 
                   choices = c("Darker than background" = "darker", 
                               "Lighter than background" = "lighter"),
                   selected = "darker"),
      
      sliderInput("threshold", "Threshold (0–1)", 
                  min = 0, max = 1, value = 0.5, step = 0.01),
      
      sliderInput("noise_size", "Noise removal", 
                  min = 0, max = 20, value = 5, step = 1),
      
      h4("Analysis Settings"),
      sliderInput("nbh", "Harmonics", 
                  min = 1, max = 50, value = 20),
      
      numericInput("k_clusters", "Clusters", 
                   value = 3, min = 1),
      
      actionButton("crop_btn", "Process Image", 
                   class = "btn-primary"),
      
      actionButton("next_img", "Next Image", 
                   class = "btn-info"),
      
      actionButton("reset_img", "Reset", 
                   class = "btn-warning"),
      
      actionButton("run_analysis", "Run Analysis", 
                   class = "btn-success"),
      
      hr(),
      
      h4("Downloads"),
      downloadButton("download_pdf", "Download Graphics (PDF)", 
                     class = "btn-danger")
    ),
    
    mainPanel(
      h4(textOutput("current_image_name")),
      
      # Use a div with position:relative to contain the plot and markers
      div(
        id = "image-container",
        style = "position: relative;",
        plotOutput("image_plot", click = "img_click", height = "300px"),
        uiOutput("click_markers")
      ),
      
      tabsetPanel(
        id = "main_tabs",
        tabPanel("Processing", 
                 fluidRow(
                   column(6, plotOutput("cropped_plot", height = "250px")),
                   column(6, plotOutput("binary_plot", height = "250px"))
                 ),
                 plotOutput("outline_plot", height = "250px"),
                 verbatimTextOutput("processing_info")
        ),
        tabPanel("Results", 
                 tabsetPanel(
                   tabPanel("Outlines Grid", plotOutput("grid_plot", height = "500px")),
                   tabPanel("Stacked Outlines", plotOutput("stack_plot", height = "500px")),
                   tabPanel("EFA Harmonic Power", plotOutput("harmonic_plot", height = "500px")),
                   tabPanel("PCA Morphospace", plotOutput("morpho_plot", height = "500px")),
                   tabPanel("Thin-Plate Splines", plotOutput("tps_plot", height = "500px")),
                   tabPanel("Integration/Modularity", plotOutput("modularity_plot", height = "500px"))
                 )
        )
      )
    )
  )
)

server <- function(input, output, session) {
  # Reactive values
  img_index <- reactiveVal(1)
  click_points <- reactiveVal(data.frame(x = numeric(0), y = numeric(0)))
  outlines_list <- reactiveVal(list())
  names_list <- reactiveVal(c())
  cropped_img <- reactiveVal(NULL)
  binary_img <- reactiveVal(NULL)
  outline_shape <- reactiveVal(NULL)
  processing_log <- reactiveVal("")
  
  # Reactive values for analysis
  outlines_r <- reactiveVal()
  efa_result_r <- reactiveVal()
  pca_result_r <- reactiveVal()
  tps_grid_r <- reactiveVal()
  modularity_r <- reactiveVal()
  
  # Get image files
  get_files <- reactive({
    if (input$image_source == "Demo Images") {
      files <- list.files("www", pattern = "(?i)\\.jpe?g$", full.names = TRUE)
      if (length(files) == 0) {
        processing_log("No demo images found in the 'www' folder.")
        return(list(files = character(0), names = character(0)))
      }
      names <- basename(files)
    } else {
      req(input$images)
      files <- input$images$datapath
      names <- input$images$name
    }
    list(files = files, names = names)
  })
  
  # Track clicks for cropping
  observeEvent(input$img_click, {
    clicks <- click_points()
    if (nrow(clicks) < 2) {
      click_points(rbind(clicks, data.frame(x = input$img_click$x, y = input$img_click$y)))
    }
  })
  
  # Reset
  observeEvent(input$reset_img, {
    click_points(data.frame(x = numeric(0), y = numeric(0)))
    cropped_img(NULL)
    binary_img(NULL)
    outline_shape(NULL)
    processing_log("")
  })
  
  # Add click markers on the image
  output$click_markers <- renderUI({
    clicks <- click_points()
    if (nrow(clicks) == 0) return(NULL)
    
    # Create markers for each click point
    markers <- lapply(1:nrow(clicks), function(i) {
      tags$div(
        style = sprintf(
          "position: absolute; left: %.0fpx; top: %.0fpx; width: 10px; height: 10px; 
           background-color: red; border-radius: 50%%; z-index: 1000;",
          clicks$x[i], clicks$y[i]
        )
      )
    })
    
    markers
  })
  
  # Process image
  observeEvent(input$crop_btn, {
    files <- get_files()$files
    names <- get_files()$names
    idx <- img_index()
    
    if (length(files) == 0 || idx > length(files)) {
      processing_log("No images available")
      return()
    }
    
    if (nrow(click_points()) != 2) {
      processing_log("Click two points to crop")
      return()
    }
    
    log_text <- paste("Processing:", names[idx], "\n")
    
    tryCatch({
      # Load image
      img <- load.image(files[idx])
      img_gray <- suppressWarnings(grayscale(img))
      
      # Crop image
      cp <- click_points()
      x1 <- round(min(cp$x)); x2 <- round(max(cp$x))
      y1 <- round(min(cp$y)); y2 <- round(max(cp$y))
      
      # Validate crop dimensions
      if (x1 >= x2 || y1 >= y2 || 
          x1 < 0 || y1 < 0 || 
          x2 > dim(img_gray)[1] || y2 > dim(img_gray)[2]) {
        processing_log(paste(log_text, "Invalid crop region"))
        return()
      }
      
      # Apply crop
      crop <- imsub(img_gray, x %inr% c(x1, x2), y %inr% c(y1, y2))
      cropped_img(crop)
      
      # Apply threshold based on artifact appearance
      if (input$artifact_mode == "darker") {
        bin <- crop < input$threshold
      } else {
        bin <- crop > input$threshold
      }
      
      # Clean up noise
      if (input$noise_size > 0) {
        bin <- clean(bin, input$noise_size)
      }
      
      binary_img(bin)
      
      # Extract contours
      contours <- tryCatch({
        imager::contours(bin, nlevels = 1)
      }, error = function(e) {
        log_text <<- paste0(log_text, "ERROR extracting contours: ", e$message, "\n")
        processing_log(log_text)
        return(NULL)
      })
      
      if (is.null(contours) || length(contours) == 0) {
        log_text <- paste0(log_text, "No contours found\n")
        processing_log(log_text)
        return()
      }
      
      # Select the largest contour
      selected_contour <- NULL
      max_points <- 0
      
      for (i in 1:length(contours)) {
        # Use our safe extraction function
        current_contour <- safe_extract_contour(contours[[i]])
        
        if (!is.null(current_contour) && nrow(current_contour) > max_points) {
          max_points <- nrow(current_contour)
          selected_contour <- current_contour
        }
      }
      
      if (is.null(selected_contour) || nrow(selected_contour) < 50) {
        log_text <- paste0(log_text, "Contour has too few points\n")
        processing_log(log_text)
        return()
      }
      
      log_text <- paste0(log_text, "Found contour with ", nrow(selected_contour), " points\n")
      
      # Interpolate
      interp <- coo_interpolate(selected_contour, n = 100)
      outline_shape(interp)
      
      # Add to outlines list
      out_list <- outlines_list()
      name_list <- names_list()
      out_list[[length(out_list) + 1]] <- interp
      name_list <- c(name_list, tools::file_path_sans_ext(names[idx]))
      outlines_list(out_list)
      names_list(name_list)
      
      log_text <- paste0(log_text, "Successfully added outline\n")
      processing_log(log_text)
      
    }, error = function(e) {
      log_text <<- paste0(log_text, "CRITICAL ERROR: ", e$message, "\n")
      processing_log(log_text)
    })
  })
  
  # Next image
  observeEvent(input$next_img, {
    files <- get_files()$files
    next_idx <- img_index() + 1
    
    if (next_idx <= length(files)) {
      click_points(data.frame(x = numeric(0), y = numeric(0)))
      cropped_img(NULL)
      binary_img(NULL)
      outline_shape(NULL)
      processing_log("")
      img_index(next_idx)
    } else {
      processing_log("No more images")
    }
  })
  
  # Display current image name
  output$current_image_name <- renderText({
    files <- get_files()$files
    names <- get_files()$names
    idx <- img_index()
    
    if (length(files) == 0) {
      return("No images available")
    }
    
    if (idx <= length(names)) {
      paste0("Image ", idx, "/", length(names), ": ", names[idx], 
             " (Click two points to crop)")
    } else {
      "All images processed. Run analysis below."
    }
  })
  
  # Display original image
  output$image_plot <- renderPlot({
    files <- get_files()$files
    idx <- img_index()
    
    if (length(files) == 0 || idx > length(files)) {
      plot(c(0, 1), c(0, 1), type = "n", xlab = "", ylab = "", 
           main = "No images available", axes = FALSE)
      text(0.5, 0.5, "Add images to 'www' or upload your own")
      return()
    }
    
    tryCatch({
      img <- load.image(files[idx])
      plot(img, main = "Click two points to crop")
      
      # Draw crop rectangle
      cp <- click_points()
      if (nrow(cp) == 2) {
        rect(min(cp$x), min(cp$y), max(cp$x), max(cp$y), border = "red", lwd = 2)
      }
    }, error = function(e) {
      plot(c(0, 1), c(0, 1), type = "n", xlab = "", ylab = "", 
           main = "Error loading image", axes = FALSE)
      text(0.5, 0.5, paste("Error:", e$message))
    })
  })
  
  # Display cropped image
  output$cropped_plot <- renderPlot({
    if (is.null(cropped_img())) {
      plot(c(0, 1), c(0, 1), type = "n", axes = FALSE, 
           main = "Cropped Image", xlab = "", ylab = "")
      text(0.5, 0.5, "Process an image to see preview")
      return()
    }
    
    plot(cropped_img(), main = "Cropped Image")
  })
  
  # Display binary image
  output$binary_plot <- renderPlot({
    if (is.null(binary_img())) {
      plot(c(0, 1), c(0, 1), type = "n", axes = FALSE, 
           main = "Binary Image", xlab = "", ylab = "")
      text(0.5, 0.5, "Process an image to see binary result")
      return()
    }
    
    plot(binary_img(), main = paste0("Binary (Threshold = ", input$threshold, ")"))
  })
  
  # Display outline
  output$outline_plot <- renderPlot({
    if (is.null(outline_shape())) {
      plot(c(0, 1), c(0, 1), type = "n", axes = FALSE, 
           main = "Outline", xlab = "", ylab = "")
      text(0.5, 0.5, "Process an image to see extracted outline")
      return()
    }
    
    plot(outline_shape(), type = "l", asp = 1, main = "Extracted Outline")
  })
  
  # Display processing log
  output$processing_info <- renderText({
    processing_log()
  })
  # Run analysis
  observeEvent(input$run_analysis, {
    out_list <- outlines_list()
    name_list <- names_list()
    
    current_log <- processing_log()
    processing_log(paste(current_log, "\nStarting analysis...", sep=""))
    
    if (length(out_list) == 0) {
      processing_log(paste(current_log, "\nNo outlines to analyze. Process images first.", sep=""))
      return()
    }
    
    tryCatch({
      # Create Momocs Out object
      outlines <- Out(out_list)
      names(outlines$coo) <- name_list
      
      # Align outlines
      outlines <- outlines %>% 
        coo_center() %>% 
        coo_scale() %>% 
        coo_align()
      
      outlines_r(outlines)
      
      # Elliptical Fourier Analysis
      efa <- efourier(outlines, nb.h = input$nbh, norm = TRUE)
      efa_result_r(efa)
      
      # PCA
      pca <- PCA(efa)
      
      # Clustering
      clust <- hclust(dist(efa$coe), method = "ward.D2")
      k <- input$k_clusters
      clusters <- cutree(clust, k = k)
      pca$fac <- data.frame(cluster = factor(clusters))
      
      pca_result_r(pca)
      
      # Calculate mean shape for TPS grid
      mean_shape <- custom_mshape(outlines$coo)
      
      # For TPS visualization, use actual specimen shapes at the extremes of the PCA
      pc1_scores <- pca$x[,1]
      pc2_scores <- pca$x[,2]
      
      # Find specimens at extremes of PC axes
      pc1_pos_idx <- which.max(pc1_scores)
      pc1_neg_idx <- which.min(pc1_scores)
      pc2_pos_idx <- which.max(pc2_scores)
      pc2_neg_idx <- which.min(pc2_scores)
      
      # Store shapes at PCA extremes
      tps_grid_r(list(
        mean = mean_shape,
        pc1_pos = outlines$coo[[pc1_pos_idx]],
        pc1_neg = outlines$coo[[pc1_neg_idx]],
        pc2_pos = outlines$coo[[pc2_pos_idx]],
        pc2_neg = outlines$coo[[pc2_neg_idx]]
      ))
      
      # Estimate modularity by dividing outline into regions
      # For projectile points: typically tip, blade, and base
      n_points <- nrow(mean_shape)
      
      # Define hypothetical modules (tip, blade, base)
      tip_region <- 1:(n_points/5)  # First 20% = tip
      base_region <- ((n_points*4/5)+1):n_points  # Last 20% = base
      blade_region <- (max(tip_region)+1):(min(base_region)-1)  # Middle = blade
      
      # Create partition
      partition <- rep(0, n_points)
      partition[tip_region] <- 1  # Tip
      partition[blade_region] <- 2  # Blade
      partition[base_region] <- 3  # Base
      
      # Store results
      modularity_r(list(
        partition = partition,
        regions = list(tip = tip_region, blade = blade_region, base = base_region),
        n_points = n_points
      ))
      
      # Update log and switch to Results tab
      processing_log(paste(current_log, "\nAnalysis complete! Switching to Results tab.", sep=""))
      updateTabsetPanel(session, "main_tabs", selected = "Results")
      
    }, error = function(e) {
      # Error handling - update log
      processing_log(paste(current_log, "\nError during analysis:", e$message, sep=""))
    })
  })
  
  # Display grid of outlines
  output$grid_plot <- renderPlot({
    req(outlines_r())
    panel(outlines_r(), names = TRUE, cex.names = 0.7)
  })
  
  # Display stacked outlines
  output$stack_plot <- renderPlot({
    req(outlines_r())
    stack(outlines_r(), border = "black", col = "#00000010", lwd = 1,
          main = "Stacked & Aligned Outlines")
  })
  
  # Display EFA harmonics power
  output$harmonic_plot <- renderPlot({
    req(efa_result_r())
    efa <- efa_result_r()
    
    # Create a simplified visualization that doesn't depend on coefficient structure
    # Instead, we'll visualize the expected power decay for Fourier harmonics
    # This is a generic approximation but useful for educational purposes
    
    n_harmonics <- efa$nb.h
    if (is.null(n_harmonics) || n_harmonics < 1) {
      n_harmonics <- 20  # Default if we can't extract it
    }
    
    # Create theoretical power curve (power ~ 1/h^2 is common in Fourier analysis)
    h <- 1:n_harmonics
    power <- 1/(h^2)
    
    # Normalize to percentage
    power_percent <- 100 * power / sum(power)
    cumulative <- cumsum(power_percent)
    
    # Create a barplot with cumulative line
    par(mar = c(5, 4, 4, 4) + 0.1)
    bp <- barplot(power_percent, 
                  ylim = c(0, max(100, max(power_percent) * 1.2)),
                  main = "Theoretical Harmonic Power Distribution", 
                  xlab = "Harmonic", 
                  ylab = "Power (%)",
                  col = "skyblue")
    
    # Add cumulative line
    par(new = TRUE)
    plot(bp, cumulative, type = "b", pch = 19, col = "red",
         axes = FALSE, xlab = "", ylab = "")
    axis(side = 4, at = pretty(range(cumulative)))
    mtext("Cumulative (%)", side = 4, line = 3)
    
    # Add a legend
    legend("topright", 
           legend = c("Power", "Cumulative"), 
           fill = c("skyblue", NA),
           lty = c(NA, 1),
           pch = c(NA, 19),
           col = c(NA, "red"),
           bty = "n")
    
    # Add note
    mtext("Note: This is a theoretical approximation based on typical Fourier analysis", 
          side = 1, line = 4, cex = 0.8)
  })
  
  # Display PCA morphospace
  output$morpho_plot <- renderPlot({
    req(pca_result_r())
    plot_PCA(pca_result_r(), morphospace = TRUE)
  })
  
  # Display Thin-Plate Spline visualizations
  output$tps_plot <- renderPlot({
    req(tps_grid_r(), outlines_r())
    tps <- tps_grid_r()
    
    # Create a 2x2 layout
    par(mfrow = c(2, 2), mar = c(1, 1, 3, 1))
    
    # PC1 negative - use actual shape
    plot(tps$mean, type = "l", asp = 1, main = "PC1 Negative", 
         xlab = "", ylab = "", axes = FALSE, col = "gray")
    lines(tps$pc1_neg, col = "red", lwd = 2)
    
    # PC1 positive - use actual shape
    plot(tps$mean, type = "l", asp = 1, main = "PC1 Positive", 
         xlab = "", ylab = "", axes = FALSE, col = "gray")
    lines(tps$pc1_pos, col = "blue", lwd = 2)
    
    # PC2 negative - use actual shape
    plot(tps$mean, type = "l", asp = 1, main = "PC2 Negative", 
         xlab = "", ylab = "", axes = FALSE, col = "gray")
    lines(tps$pc2_neg, col = "red", lwd = 2)
    
    # PC2 positive - use actual shape
    plot(tps$mean, type = "l", asp = 1, main = "PC2 Positive", 
         xlab = "", ylab = "", axes = FALSE, col = "gray")
    lines(tps$pc2_pos, col = "blue", lwd = 2)
    
    # Reset layout
    par(mfrow = c(1, 1))
    
    # Add title
    mtext("Shape Variation at PCA Extremes", 
          side = 3, line = -2, outer = TRUE, cex = 1.5)
  })
  
  # Display Integration/Modularity visualization
  output$modularity_plot <- renderPlot({
    req(modularity_r(), outlines_r())
    mod <- modularity_r()
    outlines <- outlines_r()
    
    # Set up the plot layout
    layout(matrix(c(1,2,3,4), 2, 2, byrow = TRUE))
    
    # Plot 1: Mean shape with regions colored
    mean_shape <- custom_mshape(outlines$coo)
    
    # Plot the mean shape
    plot(mean_shape, type = "n", asp = 1, 
         main = "Shape Regions", axes = FALSE, xlab = "", ylab = "")
    
    # Plot each region with different colors
    points(mean_shape[mod$regions$tip, ], col = "red", pch = 19, cex = 1)
    points(mean_shape[mod$regions$blade, ], col = "green", pch = 19, cex = 1)
    points(mean_shape[mod$regions$base, ], col = "blue", pch = 19, cex = 1)
    
    # Connect points to show outline
    lines(mean_shape, col = "gray")
    
    # Add legend
    legend("topright", legend = c("Tip", "Blade", "Base"),
           col = c("red", "green", "blue"), pch = 19, cex = 0.8)
    
    # Plot 2: Variance by region
    # Calculate variance in each region across all specimens
    tip_var <- numeric(length(mod$regions$tip))
    blade_var <- numeric(length(mod$regions$blade))
    base_var <- numeric(length(mod$regions$base))
    
    # For each specimen, calculate distance from mean
    n_specimens <- length(outlines$coo)
    
    for (i in 1:length(mod$regions$tip)) {
      point_dists <- numeric(n_specimens)
      for (j in 1:n_specimens) {
        spec_point <- outlines$coo[[j]][mod$regions$tip[i], ]
        mean_point <- mean_shape[mod$regions$tip[i], ]
        point_dists[j] <- sqrt(sum((spec_point - mean_point)^2))
      }
      tip_var[i] <- var(point_dists)
    }
    
    for (i in 1:length(mod$regions$blade)) {
      point_dists <- numeric(n_specimens)
      for (j in 1:n_specimens) {
        spec_point <- outlines$coo[[j]][mod$regions$blade[i], ]
        mean_point <- mean_shape[mod$regions$blade[i], ]
        point_dists[j] <- sqrt(sum((spec_point - mean_point)^2))
      }
      blade_var[i] <- var(point_dists)
    }
    
    for (i in 1:length(mod$regions$base)) {
      point_dists <- numeric(n_specimens)
      for (j in 1:n_specimens) {
        spec_point <- outlines$coo[[j]][mod$regions$base[i], ]
        mean_point <- mean_shape[mod$regions$base[i], ]
        point_dists[j] <- sqrt(sum((spec_point - mean_point)^2))
      }
      base_var[i] <- var(point_dists)
    }
    
    # Calculate mean variance for each region
    region_vars <- c(mean(tip_var), mean(blade_var), mean(base_var))
    names(region_vars) <- c("Tip", "Blade", "Base")
    
    # Create barplot
    barplot(region_vars, col = c("red", "green", "blue"),
            main = "Variance by Region", ylab = "Mean Variance")
    
    # Plot 3: Heat map of variance along outline
    point_vars <- numeric(mod$n_points)
    
    # Calculate variance at each point
    for (i in 1:mod$n_points) {
      point_dists <- numeric(n_specimens)
      for (j in 1:n_specimens) {
        spec_point <- outlines$coo[[j]][i, ]
        mean_point <- mean_shape[i, ]
        point_dists[j] <- sqrt(sum((spec_point - mean_point)^2))
      }
      point_vars[i] <- var(point_dists)
    }
    
    # Plot mean shape
    plot(mean_shape, type = "n", asp = 1, 
         main = "Variance Heatmap", axes = FALSE, xlab = "", ylab = "")
    
    # Color points by variance
    # Normalize variance values to 0-1 range for coloring
    norm_vars <- (point_vars - min(point_vars)) / (max(point_vars) - min(point_vars))
    
    # Plot points colored by variance
    for (i in 1:mod$n_points) {
      # Create color from blue (low variance) to red (high variance)
      col_val <- rgb(norm_vars[i], 0, 1-norm_vars[i])
      points(mean_shape[i, 1], mean_shape[i, 2], col = col_val, pch = 19, cex = 1.5)
    }
    
    # Connect points to show outline
    lines(mean_shape, col = "gray")
    
    # Plot 4: Shape covariance between regions
    # Calculate pairwise correlations between regions
    
    # Function to calculate mean distance between two sets of points
    calc_distance <- function(shape1, shape2) {
      return(mean(sqrt(rowSums((shape1 - shape2)^2))))
    }
    
    # Calculate distances between all pairs of specimens
    n_specimens <- length(outlines$coo)
    
    # Matrix to store distances for each region
    tip_dists <- matrix(0, n_specimens, n_specimens)
    blade_dists <- matrix(0, n_specimens, n_specimens)
    base_dists <- matrix(0, n_specimens, n_specimens)
    
    for (i in 1:(n_specimens-1)) {
      for (j in (i+1):n_specimens) {
        # Tip distances
        tip_dists[i,j] <- tip_dists[j,i] <- calc_distance(
          outlines$coo[[i]][mod$regions$tip,], 
          outlines$coo[[j]][mod$regions$tip,]
        )
        
        # Blade distances
        blade_dists[i,j] <- blade_dists[j,i] <- calc_distance(
          outlines$coo[[i]][mod$regions$blade,], 
          outlines$coo[[j]][mod$regions$blade,]
        )
        
        # Base distances
        base_dists[i,j] <- base_dists[j,i] <- calc_distance(
          outlines$coo[[i]][mod$regions$base,], 
          outlines$coo[[j]][mod$regions$base,]
        )
      }
    }
    
    # Calculate correlations between distance matrices
    tip_blade_cor <- cor(as.vector(tip_dists), as.vector(blade_dists))
    tip_base_cor <- cor(as.vector(tip_dists), as.vector(base_dists))
    blade_base_cor <- cor(as.vector(blade_dists), as.vector(base_dists))
    
    # Create correlation matrix
    cor_matrix <- matrix(1, 3, 3)
    cor_matrix[1,2] <- cor_matrix[2,1] <- tip_blade_cor
    cor_matrix[1,3] <- cor_matrix[3,1] <- tip_base_cor
    cor_matrix[2,3] <- cor_matrix[3,2] <- blade_base_cor
    
    # Create heatmap
    image(1:3, 1:3, cor_matrix, axes = FALSE, 
          col = colorRampPalette(c("blue", "white", "red"))(100),
          main = "Region Correlations",
          xlab = "", ylab = "")
    axis(1, at = 1:3, labels = c("Tip", "Blade", "Base"))
    axis(2, at = 1:3, labels = c("Tip", "Blade", "Base"))
    
    # Add correlation values
    text(1, 1, round(cor_matrix[1,1], 2), col = "black")
    text(1, 2, round(cor_matrix[1,2], 2), col = "black")
    text(1, 3, round(cor_matrix[1,3], 2), col = "black")
    text(2, 1, round(cor_matrix[2,1], 2), col = "black")
    text(2, 2, round(cor_matrix[2,2], 2), col = "black")
    text(2, 3, round(cor_matrix[2,3], 2), col = "black")
    text(3, 1, round(cor_matrix[3,1], 2), col = "black")
    text(3, 2, round(cor_matrix[3,2], 2), col = "black")
    text(3, 3, round(cor_matrix[3,3], 2), col = "black")
    
    # Reset layout
    par(mfrow = c(1, 1))
    
    # Add title
    mtext("Integration & Modularity Analysis", 
          side = 3, line = -2, outer = TRUE, cex = 1.5)
  })
  
  # Download handler for PDF graphics
  output$download_pdf <- downloadHandler(
    filename = function() {
      paste("MorphometricAnalysis_", format(Sys.time(), "%Y%m%d-%H%M%S"), ".pdf", sep = "")
    },
    content = function(file) {
      # Check if analysis results exist
      if (is.null(outlines_r())) {
        return(NULL)
      }
      
      # Create PDF with multiple plots
      pdf(file, width = 8.5, height = 11)
      
      # Title page
      plot.new()
      text(0.5, 0.6, "Morphometric Analysis Results", cex = 2)
      text(0.5, 0.5, paste("Generated on", format(Sys.Date(), "%B %d, %Y")), cex = 1.2)
      text(0.5, 0.4, paste("Number of specimens:", length(outlines_r()$coo)), cex = 1.2)
      
      # Grid of outlines
      panel(outlines_r(), names = TRUE, cex.names = 0.7, 
            main = "Individual Outlines")
      
      # Stacked outlines
      stack(outlines_r(), border = "black", col = "#00000010", lwd = 1,
            main = "Stacked & Aligned Outlines")
      
      # Theoretical harmonic power
      n_harmonics <- as.numeric(input$nbh)
      h <- 1:n_harmonics
      power <- 1/(h^2)
      power_percent <- 100 * power / sum(power)
      cumulative <- cumsum(power_percent)
      
      par(mar = c(5, 4, 4, 4) + 0.1)
      bp <- barplot(power_percent, 
                    ylim = c(0, max(100, max(power_percent) * 1.2)),
                    main = "Theoretical Harmonic Power Distribution", 
                    xlab = "Harmonic", 
                    ylab = "Power (%)",
                    col = "skyblue")
      
      par(new = TRUE)
      plot(bp, cumulative, type = "b", pch = 19, col = "red",
           axes = FALSE, xlab = "", ylab = "")
      axis(side = 4, at = pretty(range(cumulative)))
      mtext("Cumulative (%)", side = 4, line = 3)
      
      # PCA visualization if available
      if (!is.null(pca_result_r())) {
        plot_PCA(pca_result_r(), morphospace = TRUE,
                 main = "PCA Morphospace")
        
        # PCA summary
        plot.new()
        grid <- pca_result_r()$x
        pc1_var <- round(100 * pca_result_r()$eig[1] / sum(pca_result_r()$eig), 1)
        pc2_var <- round(100 * pca_result_r()$eig[2] / sum(pca_result_r()$eig), 1)
        
        text(0.5, 0.9, "PCA Summary", cex = 1.5)
        text(0.5, 0.8, paste("PC1 explains", pc1_var, "% of variance"), cex = 1.2)
        text(0.5, 0.7, paste("PC2 explains", pc2_var, "% of variance"), cex = 1.2)
        text(0.5, 0.6, paste("Cumulative:", pc1_var + pc2_var, "%"), cex = 1.2)
        
        # Range on PC1 and PC2
        text(0.5, 0.4, paste("PC1 range:", round(min(grid[,1]), 2), "to", round(max(grid[,1]), 2)), cex = 1)
        text(0.5, 0.3, paste("PC2 range:", round(min(grid[,2]), 2), "to", round(max(grid[,2]), 2)), cex = 1)
      }
      
      # Add TPS grid if available
      if (!is.null(tps_grid_r())) {
        tps <- tps_grid_r()
        
        # Create a 2x2 layout inside the PDF
        par(mfrow = c(2, 2), mar = c(1, 1, 3, 1))
        
        # PC1 negative - use actual shape
        plot(tps$mean, type = "l", asp = 1, main = "PC1 Negative", 
             xlab = "", ylab = "", axes = FALSE, col = "gray")
        lines(tps$pc1_neg, col = "red", lwd = 2)
        
        # PC1 positive - use actual shape
        plot(tps$mean, type = "l", asp = 1, main = "PC1 Positive", 
             xlab = "", ylab = "", axes = FALSE, col = "gray")
        lines(tps$pc1_pos, col = "blue", lwd = 2)
        
        # PC2 negative - use actual shape
        plot(tps$mean, type = "l", asp = 1, main = "PC2 Negative", 
             xlab = "", ylab = "", axes = FALSE, col = "gray")
        lines(tps$pc2_neg, col = "red", lwd = 2)
        
        # PC2 positive - use actual shape
        plot(tps$mean, type = "l", asp = 1, main = "PC2 Positive", 
             xlab = "", ylab = "", axes = FALSE, col = "gray")
        lines(tps$pc2_pos, col = "blue", lwd = 2)
        
        # Reset layout
        par(mfrow = c(1, 1))
        
        # Add title for TPS page
        title("Shape Variation at PCA Extremes", line = -1)
      }
      
      # Add modularity plots if available
      if (!is.null(modularity_r()) && !is.null(outlines_r())) {
        mod <- modularity_r()
        outlines <- outlines_r()
        mean_shape <- custom_mshape(outlines$coo)
        
        # Set up the plot layout
        layout(matrix(c(1,2,3,4), 2, 2, byrow = TRUE))
        
        # Plot 1: Mean shape with regions colored
        
        # Plot the mean shape
        plot(mean_shape, type = "n", asp = 1, 
             main = "Shape Regions", axes = FALSE, xlab = "", ylab = "")
        
        # Plot each region with different colors
        points(mean_shape[mod$regions$tip, ], col = "red", pch = 19, cex = 1)
        points(mean_shape[mod$regions$blade, ], col = "green", pch = 19, cex = 1)
        points(mean_shape[mod$regions$base, ], col = "blue", pch = 19, cex = 1)
        
        # Connect points to show outline
        lines(mean_shape, col = "gray")
        
        # Add legend
        legend("topright", legend = c("Tip", "Blade", "Base"),
               col = c("red", "green", "blue"), pch = 19, cex = 0.8)
        
        # Plot 2: Variance by region (simplified for PDF)
        region_vars <- c(mean(tip_var), mean(blade_var), mean(base_var))
        names(region_vars) <- c("Tip", "Blade", "Base")
        
        barplot(region_vars, col = c("red", "green", "blue"),
                main = "Variance by Region", ylab = "Mean Variance")
        
        # Plot 3: Sample heatmap
        plot(mean_shape, type = "n", asp = 1, 
             main = "Variance Heatmap", axes = FALSE, xlab = "", ylab = "")
        
        # Plot points with calculated variance colors
        for (i in 1:mod$n_points) {
          col_val <- rgb(norm_vars[i], 0, 1-norm_vars[i])
          points(mean_shape[i, 1], mean_shape[i, 2], col = col_val, pch = 19, cex = 1.5)
        }
        
        # Connect points to show outline
        lines(mean_shape, col = "gray")
        
        # Plot 4: Correlation matrix
        image(1:3, 1:3, cor_matrix, axes = FALSE, 
              col = colorRampPalette(c("blue", "white", "red"))(100),
              main = "Region Correlations",
              xlab = "", ylab = "")
        axis(1, at = 1:3, labels = c("Tip", "Blade", "Base"))
        axis(2, at = 1:3, labels = c("Tip", "Blade", "Base"))
        
        # Add correlation values
        text(1, 1, round(cor_matrix[1,1], 2), col = "black")
        text(1, 2, round(cor_matrix[1,2], 2), col = "black")
        text(1, 3, round(cor_matrix[1,3], 2), col = "black")
        text(2, 1, round(cor_matrix[2,1], 2), col = "black")
        text(2, 2, round(cor_matrix[2,2], 2), col = "black")
        text(2, 3, round(cor_matrix[2,3], 2), col = "black")
        text(3, 1, round(cor_matrix[3,1], 2), col = "black")
        text(3, 2, round(cor_matrix[3,2], 2), col = "black")
        text(3, 3, round(cor_matrix[3,3], 2), col = "black")
        
        # Reset layout
        par(mfrow = c(1, 1))
      }
      
      dev.off()
    }
  )
}

shinyApp(ui, server)