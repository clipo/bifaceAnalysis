library(shiny)
library(Momocs)
library(imager)
library(tidyverse)
library(ggplot2)

ui <- fluidPage(
  titlePanel("Momocs Shape Analysis with Interactive Cropping"),
  
  sidebarLayout(
    sidebarPanel(
      radioButtons("image_source", "Select image source:", choices = c("Demo Images", "Upload My Own"), selected = "Demo Images"),
      conditionalPanel(
        condition = "input.image_source == 'Upload My Own'",
        fileInput("images", "Upload JPEG Images", accept = c(".jpg", ".jpeg"), multiple = TRUE)
      ),
      sliderInput("threshold", "Threshold (0–1)", min = 0, max = 1, value = 0.1, step = 0.01),
      sliderInput("nbh", "Number of Harmonics", min = 1, max = 50, value = 20),
      numericInput("k_clusters", "Number of Clusters", value = 3, min = 1),
      checkboxInput("show_threshold", "Show thresholded image preview", value = TRUE),
      actionButton("crop_btn", "Crop This Image"),
      actionButton("next_img", "Next Image"),
      actionButton("run_analysis", "Run EFA + PCA"),
      br(),
      downloadButton("download_pca", "Download PCA Results"),
      downloadButton("download_outlines", "Download Outlines RDS"),
      downloadButton("download_plots", "Download All Plots (PDF)")
    ),
    
    mainPanel(
      h4(textOutput("current_image_name")),
      plotOutput("image_plot", click = "img_click"),
      conditionalPanel("input.show_threshold == true", plotOutput("cropped_plot")),
      plotOutput("outline_plot"),
      tabsetPanel(
        tabPanel("Outlines Grid", plotOutput("grid_plot")),
        tabPanel("Stacked Outlines", plotOutput("stack_plot")),
        tabPanel("PCA (ggplot2)", plotOutput("gg_pca_plot")),
        tabPanel("PCA (Outlines)", plotOutput("morpho_plot")),
        tabPanel("Dendrogram", plotOutput("dendro_plot"))
      )
    )
  )
)

server <- function(input, output, session) {
  img_index <- reactiveVal(1)
  click_points <- reactiveVal(data.frame(x = numeric(0), y = numeric(0)))
  outlines_list <- reactiveVal(list())
  names_list <- reactiveVal(c())
  cropped_img <- reactiveVal(NULL)
  outline_shape <- reactiveVal(NULL)
  
  get_files <- reactive({
    if (input$image_source == "Demo Images") {
      files <- list.files("www", pattern = "\\.jpe?g$", full.names = TRUE)
      names <- basename(files)
    } else {
      req(input$images)
      files <- input$images$datapath
      names <- input$images$name
    }
    list(files = files, names = names)
  })
  
  observeEvent(input$img_click, {
    clicks <- click_points()
    if (nrow(clicks) < 2) {
      click_points(rbind(clicks, data.frame(x = input$img_click$x, y = input$img_click$y)))
    }
  })
  
  observeEvent(input$crop_btn, {
    files <- get_files()$files
    names <- get_files()$names
    idx <- img_index()
    req(idx <= length(files))
    req(nrow(click_points()) == 2)
    
    img <- load.image(files[idx]) %>% grayscale()
    cp <- click_points()
    x1 <- round(min(cp$x)); x2 <- round(max(cp$x))
    y1 <- round(min(cp$y)); y2 <- round(max(cp$y))
    crop <- imsub(img, x %inr% c(x1, x2), y %inr% c(y1, y2))
    cropped_img(crop)
    
    bin <- crop < input$threshold
    contours <- contours(bin, nlevels = 1)
    
    if (length(contours) > 0) {
      pts <- as.data.frame(contours[[1]]) %>% select(x, y)
      if (nrow(pts) >= 100) {
        interp <- coo_interpolate(as.matrix(pts), n = 100)
        outline_shape(interp)
        
        out_list <- outlines_list()
        name_list <- names_list()
        out_list[[length(out_list) + 1]] <- interp
        name_list <- c(name_list, tools::file_path_sans_ext(names[idx]))
        outlines_list(out_list)
        names_list(name_list)
      }
    }
  })
  
  observeEvent(input$next_img, {
    click_points(data.frame(x = numeric(0), y = numeric(0)))
    img_index(img_index() + 1)
  })
  
  output$current_image_name <- renderText({
    names <- get_files()$names
    idx <- img_index()
    if (idx <= length(names)) {
      paste("Click TWO points to crop:", names[idx])
    } else {
      "✅ All images processed. Run analysis below."
    }
  })
  
  output$image_plot <- renderPlot({
    files <- get_files()$files
    idx <- img_index()
    req(idx <= length(files))
    img <- load.image(files[idx]) %>% grayscale()
    plot(img, main = "Click TWO points to crop")
  })
  
  output$cropped_plot <- renderPlot({
    req(cropped_img())
    plot(cropped_img(), main = "Cropped & Thresholded")
  })
  
  output$outline_plot <- renderPlot({
    req(outline_shape())
    plot(outline_shape(), type = "l", asp = 1, main = "Extracted Outline")
  })
  
  outlines_r <- reactiveVal()
  pca_result_r <- reactiveVal()
  pca_df_r <- reactiveVal()
  clust_result_r <- reactiveVal()
  
  observeEvent(input$run_analysis, {
    req(outlines_list())
    outlines <- Out(outlines_list())
    names(outlines$coo) <- names_list()
    outlines <- outlines %>% coo_center() %>% coo_scale() %>% coo_align()
    outlines_r(outlines)
    
    efa <- efourier(outlines, nb.h = input$nbh, norm = TRUE)
    pca <- PCA(efa)
    clust <- hclust(dist(efa$coe), method = "ward.D2")
    k <- input$k_clusters
    clusters <- cutree(clust, k = k)
    pca$fac <- data.frame(cluster = factor(clusters))
    
    df <- as.data.frame(pca$x)
    df$file <- names(outlines$coo)
    df$cluster <- pca$fac$cluster
    pca_result_r(pca)
    pca_df_r(df)
    clust_result_r(clust)
  })
  
  output$grid_plot <- renderPlot({
    req(outlines_r())
    panel(outlines_r(), names = TRUE, cex.names = 0.5)
  })
  
  output$stack_plot <- renderPlot({
    req(outlines_r())
    stack(outlines_r(), border = "black", col = "#00000010", lwd = 1)
  })
  
  output$gg_pca_plot <- renderPlot({
    req(pca_df_r())
    ggplot(pca_df_r(), aes(x = PC1, y = PC2, color = cluster, label = file)) +
      geom_point(size = 3) +
      geom_text(size = 2, vjust = -0.5) +
      theme_minimal() +
      scale_color_brewer(palette = "Set2") +
      labs(title = "PCA (ggplot2)", x = "PC1", y = "PC2")
  })
  
  output$morpho_plot <- renderPlot({
    req(pca_result_r())
    plot_PCA(pca_result_r(), morphospace = TRUE)
  })
  
  output$dendro_plot <- renderPlot({
    req(clust_result_r())
    plot(clust_result_r(), main = "Hierarchical Clustering")
  })
  
  output$download_pca <- downloadHandler(
    filename = function() "PCA_results.csv",
    content = function(file) {
      write.csv(pca_df_r(), file, row.names = FALSE)
    }
  )
  
  output$download_outlines <- downloadHandler(
    filename = function() "Outlines.rds",
    content = function(file) {
      saveRDS(Out(outlines_list()), file)
    }
  )
  
  output$download_plots <- downloadHandler(
    filename = function() "Momocs_Analysis_Plots.pdf",
    content = function(file) {
      pdf(file, width = 10, height = 8)
      if (!is.null(outlines_r())) {
        panel(outlines_r(), names = TRUE, cex.names = 0.5)
        stack(outlines_r(), border = "black", col = "#00000010", lwd = 1)
      }
      if (!is.null(pca_result_r())) {
        plot_PCA(pca_result_r(), morphospace = TRUE)
      }
      if (!is.null(clust_result_r())) {
        plot(clust_result_r(), main = "Hierarchical Clustering")
      }
      if (!is.null(pca_df_r())) {
        print(
          ggplot(pca_df_r(), aes(x = PC1, y = PC2, color = cluster, label = file)) +
            geom_point(size = 3) +
            geom_text(size = 2, vjust = -0.5) +
            theme_minimal() +
            scale_color_brewer(palette = "Set2") +
            labs(title = "PCA (ggplot2)", x = "PC1", y = "PC2")
        )
      }
      dev.off()
    }
  )
}

shinyApp(ui, server)
