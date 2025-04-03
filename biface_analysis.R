# Load libraries
library(Momocs)
library(magick)
library(imager)
library(tidyverse)
library(here)
library(ggplot2)

#--------------------------------------------------
# 1. Configurations
#--------------------------------------------------
force_rebuild_outlines <- FALSE     # Set TRUE to force regeneration
image_folder <- here("images")
outlines_cache <- here("cached_outlines.rds")

#--------------------------------------------------
# 2. Enhanced Image Processor (Adaptive Thresholding)
#--------------------------------------------------
process_image <- function(image_path, num_points = 100) {
  cat("🔍 Processing:", image_path, "\n")
  img_orig <- load.image(image_path)
  img_gray <- suppressWarnings(grayscale(img_orig))
  
  # Cropping
  plot(img_gray, main = "Click TOP-LEFT and BOTTOM-RIGHT to crop")
  coords <- locator(2)
  if (length(coords$x) < 2) {
    warning("❌ Two clicks not detected. Skipping:", image_path)
    return(NULL)
  }
  x1 <- round(min(coords$x)); x2 <- round(max(coords$x))
  y1 <- round(min(coords$y)); y2 <- round(max(coords$y))
  img_crop <- imsub(img_gray, x %inr% c(x1, x2), y %inr% c(y1, y2))
  
  # Background Analysis - Determine if artifact is darker or lighter than background
  # Sample border pixels to estimate background color
  width <- dim(img_crop)[1]
  height <- dim(img_crop)[2]
  
  # Get border pixels (10% inset from edge to avoid potential artifacts at the very edge)
  border_inset <- 0.1
  top <- img_crop[, round(height * border_inset):round(height * (1-border_inset)), 1, drop=FALSE]
  bottom <- img_crop[, round(height * border_inset):round(height * (1-border_inset)), height, drop=FALSE]
  left <- img_crop[1, round(height * border_inset):round(height * (1-border_inset)), , drop=FALSE]
  right <- img_crop[width, round(height * border_inset):round(height * (1-border_inset)), , drop=FALSE]
  
  # Convert to vectors and combine
  border_values <- c(as.vector(top), as.vector(bottom), as.vector(left), as.vector(right))
  background_mean <- mean(border_values, na.rm = TRUE)
  
  # Get central area (potential artifact region)
  center_width <- width * 0.5
  center_height <- height * 0.5
  center_x_range <- round((width - center_width)/2):round((width + center_width)/2)
  center_y_range <- round((height - center_height)/2):round((height + center_height)/2)
  center_region <- img_crop[center_x_range, center_y_range, ,drop=FALSE]
  center_mean <- mean(as.vector(center_region), na.rm = TRUE)
  
  # Determine if artifact is darker or lighter than background
  artifact_is_darker <- center_mean < background_mean
  
  # Try otsu thresholding first
  # Calculate histogram
  hist_data <- hist(as.vector(img_crop), plot = FALSE, breaks = seq(0, 1, by = 0.01))
  
  # Otsu's method to find optimal threshold
  otsu_threshold <- function(hist_obj) {
    total <- sum(hist_obj$counts)
    sum_total <- sum(hist_obj$mids * hist_obj$counts)
    
    weight_background <- cumsum(hist_obj$counts)
    weight_foreground <- total - weight_background
    
    sum_background <- cumsum(hist_obj$mids * hist_obj$counts)
    mean_background <- sum_background / weight_background
    mean_foreground <- (sum_total - sum_background) / weight_foreground
    
    # Calculate between-class variance
    between_var <- weight_background * weight_foreground * (mean_background - mean_foreground)^2
    
    # Find the threshold that maximizes between-class variance
    max_idx <- which.max(between_var)
    return(hist_obj$mids[max_idx])
  }
  
  # Get Otsu threshold
  threshold <- otsu_threshold(hist_data)
  
  # Apply thresholding based on whether artifact is darker or lighter
  if (artifact_is_darker) {
    img_bin <- img_crop < threshold
    cat("📊 Artifact appears darker than background. Threshold:", threshold, "\n")
  } else {
    img_bin <- img_crop > threshold
    cat("📊 Artifact appears lighter than background. Threshold:", threshold, "\n")
  }
  
  # Show user the binary image
  par(mfrow=c(1,2))
  plot(img_crop, main = "Cropped Grayscale")
  plot(img_bin, main = paste("Binary (Threshold =", round(threshold, 3), ")"))
  par(mfrow=c(1,1))
  
  # Ask user if the threshold is acceptable
  cat("Is this thresholding acceptable? (y/n): ")
  user_response <- tolower(readline())
  
  # If not acceptable, allow manual threshold adjustment
  if (user_response != "y") {
    cat("Enter new threshold value (0-1): ")
    new_threshold <- as.numeric(readline())
    
    if (!is.na(new_threshold) && new_threshold >= 0 && new_threshold <= 1) {
      threshold <- new_threshold
      
      if (artifact_is_darker) {
        img_bin <- img_crop < threshold
      } else {
        img_bin <- img_crop > threshold
      }
      
      # Show updated binary image
      plot(img_bin, main = paste("Updated Binary (Threshold =", threshold, ")"))
    } else {
      cat("❌ Invalid threshold value. Using original threshold.\n")
    }
  }
  
  # Noise cleanup - remove small connected components
  img_bin <- clean(img_bin, 5)
  
  # Outline extraction
  contours <- imager::contours(img_bin, nlevels = 1)
  if (length(contours) == 0) {
    warning("❌ No contour found:", image_path)
    return(NULL)
  }
  
  # Check for multiple contours and select the largest one
  selected_contour <- NULL
  max_points <- 0
  
  for (i in 1:length(contours)) {
    current_contour <- as.data.frame(contours[[i]]) %>% select(x, y) %>% as.matrix()
    if (nrow(current_contour) > max_points) {
      max_points <- nrow(current_contour)
      selected_contour <- current_contour
    }
  }
  
  if (is.null(selected_contour) || nrow(selected_contour) < num_points) {
    warning("❌ Contour has too few points:", image_path)
    return(NULL)
  }
  
  # Interpolate to get uniform number of points
  coords_interp <- coo_interpolate(selected_contour, n = num_points)
  
  # Final preview
  plot(coords_interp, type = "l", asp = 1, main = paste("✅ Outline:", basename(image_path)))
  
  return(coords_interp)
}

#--------------------------------------------------
# 3. Batch Processing
#--------------------------------------------------
read_images_as_out <- function(folder_path, num_points = 100) {
  files <- list.files(folder_path, pattern = "(?i)\\.jpe?g$", full.names = TRUE)
  if (length(files) == 0) stop("❌ No JPEG images found.")
  
  outlines_list <- list()
  names_list <- c()
  
  for (file in files) {
    coords <- process_image(file, num_points)
    if (!is.null(coords)) {
      outlines_list[[length(outlines_list) + 1]] <- coords
      names_list <- c(names_list, basename(file))
    } else {
      warning("⚠️ Skipping:", file)
    }
  }
  
  if (length(outlines_list) == 0) stop("❌ No valid outlines extracted.")
  
  outlines <- Out(outlines_list)
  names(outlines$coo) <- names_list
  return(outlines)
}

#--------------------------------------------------
# 4. Load or Generate Outlines
#--------------------------------------------------
if (file.exists(outlines_cache) && !force_rebuild_outlines) {
  cat("💾 Using cached outlines from", outlines_cache, "\n")
  outlines <- readRDS(outlines_cache)
  outlines <- outlines %>%
    coo_center() %>%
    coo_scale() %>%
    coo_align()
  stack(outlines, border = "black", col = "#00000010", lwd = 1,
        main = "Overlay of Cached Outlines")
} else {
  cat("🔄 Generating outlines...\n")
  outlines <- read_images_as_out(image_folder, num_points = 100)
  saveRDS(outlines, outlines_cache)
  cat("✅ Cached outlines saved to", outlines_cache, "\n")
}

#--------------------------------------------------
# 5. Align & Shape Analysis
#--------------------------------------------------
outlines <- outlines %>%
  coo_center() %>%
  coo_scale() %>%
  coo_align()

efa_result <- efourier(outlines, nb.h = 20, norm = TRUE)
pca_result <- PCA(efa_result)

#--------------------------------------------------
# 6. Clustering
#--------------------------------------------------
dist_matrix <- dist(efa_result$coe)
clust_result <- hclust(dist_matrix, method = "ward.D2")
k <- 3
cluster_assignments <- cutree(clust_result, k = k)

# Attach cluster info to PCA object
pca_result$fac <- data.frame(cluster = factor(cluster_assignments))

# Save scores with clusters
pca_df <- as.data.frame(pca_result$x)
pca_df$file <- names(outlines$coo)
pca_df$cluster <- pca_result$fac$cluster
write.csv(pca_df, "PCA_results.csv", row.names = FALSE)

#--------------------------------------------------
# 7. Standard Visualizations (Console + PDF)
#--------------------------------------------------

# 7.1 Grid of all outlines
panel(outlines, names = TRUE, cex.names = 0.5)
pdf("all_outlines_grid.pdf", width = 10, height = 10)
panel(outlines, names = TRUE, cex.names = 0.5)
dev.off()

# 7.2 Stacked overlay
stack(outlines, border = "black", col = "#00000010", lwd = 1)
pdf("stacked_outlines_overlay.pdf", width = 8, height = 8)
stack(outlines, border = "black", col = "#00000010", lwd = 1)
dev.off()

# 7.3 PCA deformation grid
PCcontrib(pca_result, nax = 1:2, sd.r = 3, amp.shp = 1, amp = 1)
pdf("pca_deformation_grids.pdf", width = 10, height = 5)
PCcontrib(pca_result, nax = 1:2, sd.r = 3, amp.shp = 1, amp = 1)
dev.off()

# 7.4 Dendrogram
plot(clust_result, labels = names(outlines$coo), main = "Hierarchical Clustering")
pdf("cluster_dendrogram.pdf", width = 8, height = 6)
plot(clust_result, labels = names(outlines$coo), main = "Hierarchical Clustering")
dev.off()

#--------------------------------------------------
# 8. ggplot2 PCA Plot Colored by Cluster
#--------------------------------------------------

pca_plot <- ggplot(pca_df, aes(x = PC1, y = PC2, color = cluster, label = file)) +
  geom_point(size = 3) +
  geom_text(size = 2, vjust = -0.5) +
  theme_minimal(base_size = 14) +
  scale_color_brewer(palette = "Set2") +
  labs(
    title = "PCA of Artifact Shapes (Momocs)",
    subtitle = paste("k =", k, "clusters"),
    x = "Principal Component 1",
    y = "Principal Component 2",
    color = "Cluster"
  )

print(pca_plot)
ggsave("ggplot2_pca_clusters.pdf", pca_plot, width = 10, height = 8)

#--------------------------------------------------
# 9. Momocs PCA Morphospace with Shape Outlines
#--------------------------------------------------

# Console preview
plot_PCA(pca_result, morphospace = TRUE)

# Save as PDF
pdf("pca_morphospace_shapes.pdf", width = 10, height = 8)
plot_PCA(pca_result, morphospace = TRUE)
dev.off()

#--------------------------------------------------
# 10. Summary to Console
#--------------------------------------------------
summary(pca_result)