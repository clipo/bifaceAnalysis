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
threshold_value <- 0.1              # Set global threshold
force_rebuild_outlines <- FALSE     # Set TRUE to force regeneration
image_folder <- here("images")
outlines_cache <- here("cached_outlines.rds")

#--------------------------------------------------
# 2. Image Processor (Interactive Crop, No Threshold Prompt)
#--------------------------------------------------
process_image <- function(image_path, num_points = 100, threshold_value = 0.1) {
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
  
  # Threshold + Binary
  img_bin <- img_crop < threshold_value
  plot(img_bin, main = paste("Threshold =", threshold_value))
  
  # Outline extraction
  contours <- imager::contours(img_bin, nlevels = 1)
  if (length(contours) == 0) {
    warning("❌ No contour found:", image_path)
    return(NULL)
  }
  coords <- as.data.frame(contours[[1]]) %>% select(x, y) %>% as.matrix()
  if (nrow(coords) < num_points) {
    warning("❌ Too few points:", image_path)
    return(NULL)
  }
  coords_interp <- coo_interpolate(coords, n = num_points)
  
  # Final preview
  plot(coords_interp, type = "l", asp = 1, main = paste("✅ Outline:", basename(image_path)))
  
  return(coords_interp)
}

#--------------------------------------------------
# 3. Batch Processing
#--------------------------------------------------
read_images_as_out <- function(folder_path, num_points = 100, threshold_value = 0.1) {
  files <- list.files(folder_path, pattern = "(?i)\\.jpe?g$", full.names = TRUE)
  if (length(files) == 0) stop("❌ No JPEG images found.")
  
  outlines_list <- list()
  names_list <- c()
  
  for (file in files) {
    coords <- process_image(file, num_points, threshold_value)
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
  outlines <- read_images_as_out(image_folder, num_points = 100, threshold_value = threshold_value)
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