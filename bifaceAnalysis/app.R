# Biface outline analysis: a Shiny app for teaching outline-based morphometrics.
# Analysis code is in R/ (sourced automatically by Shiny). This file holds the
# interface and the reactive wiring only.
# Specification: docs/specs/2026-11-17_app_rewrite_spec.md.

library(shiny)
library(bslib)
library(ggplot2)

options(shiny.maxRequestSize = 40 * 1024^2)

MAX_HARMONICS <- 30L
DEMO_DIR <- "www"
DEMO_DIRECTION <- "left"
IMAGE_PATTERN <- "(?i)\\.(jpe?g|png)$"

how_to_read <- function(id) {
  card(card_header("How to read this"), card_body(uiOutput(id), class = "small"))
}

paragraphs <- function(x) tagList(lapply(x[nzchar(x)], tags$p))

result_row <- function(plot_id, text_id, height = "440px") {
  layout_columns(col_widths = c(8, 4),
                 card(card_body(plotOutput(plot_id, height = height))),
                 how_to_read(text_id))
}

ui <- page_navbar(
  title = "Biface outline analysis",
  id = "nav",
  theme = bs_theme(version = 5, primary = "#0072B2"),
  header = tags$head(tags$style(HTML(
    ".small p { margin-bottom: 0.6rem; }
     .warn-box { border-left: 4px solid #D55E00; padding: 0.4rem 0.8rem; background: #FFF4EC; margin-bottom: 0.5rem; }"))),

  nav_panel(
    "1. Images", value = "images",
    layout_columns(
      col_widths = c(6, 6),
      card(
        card_header("Start here"),
        p("This app measures and compares the outline shapes of bifaces and projectile points from photographs. Work through the numbered steps at the top of the page."),
        tags$ol(
          tags$li("Choose images on this page."),
          tags$li("Check the outline the app finds for each image and save it."),
          tags$li("Review the saved specimens and label them."),
          tags$li("Read the results. Each result has a note on how to read it.")
        ),
        p("To see what the results look like before using your own images, run the demonstration."),
        actionButton("run_demo", "Run the demonstration with four sample points", class = "btn-primary"),
        hr(),
        radioButtons("image_source", "Images to work on",
                     choices = c("Sample images" = "demo", "My own images" = "upload"), selected = "demo"),
        conditionalPanel(
          "input.image_source == 'upload'",
          fileInput("images", "Upload JPEG or PNG images", accept = c(".jpg", ".jpeg", ".png"), multiple = TRUE)
        ),
        textOutput("image_count"),
        actionButton("go_outline", "Continue to step 2", class = "btn-outline-primary")
      ),
      card(
        card_header("Photographs that work well"),
        tags$ul(
          tags$li("One artifact per image, lying flat, photographed from directly above."),
          tags$li("A plain background that contrasts with the artifact: dark cloth for light stone, white paper for dark stone."),
          tags$li("Background visible on every side of the artifact. Nothing else should touch it."),
          tags$li("Even lighting without hard shadows. A shadow that touches the artifact becomes part of its outline."),
          tags$li("A scale bar in the frame, not touching the artifact, if you want measurements in millimetres."),
          tags$li("All artifacts pointing the same way saves time, because you state the tip direction once.")
        ),
        p(class = "text-muted", "Uploaded images are used only during your session and are not stored.")
      )
    )
  ),

  nav_panel(
    "2. Outline", value = "outline",
    layout_sidebar(
      sidebar = sidebar(
        width = 340,
        strong(textOutput("image_title")),
        div(actionButton("prev_img", "Previous"), actionButton("next_img", "Next")),
        radioButtons("mode", "Clicking on the image will",
                     choices = c("Crop: drag a box around the artifact" = "crop",
                                 "Set scale: click the two ends of a known distance" = "scale",
                                 "Set tip and base: click the tip, then the base" = "landmarks"),
                     selected = "crop"),
        selectInput("direction", "The tip points",
                    choices = c("Left" = "left", "Right" = "right", "Up" = "up", "Down" = "down"),
                    selected = DEMO_DIRECTION),
        accordion(
          open = FALSE,
          accordion_panel(
            "Detection settings",
            checkboxInput("auto_threshold", "Choose the threshold automatically", TRUE),
            conditionalPanel("!input.auto_threshold",
                             sliderInput("threshold", "Threshold", min = 0.02, max = 0.98, value = 0.5, step = 0.01)),
            selectInput("polarity", "The artifact is",
                        choices = c("Detected automatically" = "auto", "Lighter than the background" = "light",
                                    "Darker than the background" = "dark")),
            sliderInput("speck", "Remove specks narrower than (pixels)", min = 0, max = 10, value = 3, step = 1),
            actionButton("clear_crop", "Clear crop", class = "btn-sm"),
            actionButton("clear_landmarks", "Use automatic tip and base", class = "btn-sm")
          ),
          accordion_panel(
            "Scale",
            numericInput("scale_mm", "Distance between the two clicked points (mm)", value = 10, min = 0.1),
            checkboxInput("scale_all", "Use this scale for images that have none of their own", FALSE),
            actionButton("clear_scale", "Clear scale", class = "btn-sm"),
            textOutput("scale_text")
          )
        ),
        textInput("spec_name", "Specimen name"),
        textInput("spec_group", "Group label (optional, for example a type or site)"),
        actionButton("save_spec", "Save outline and go to next image", class = "btn-primary"),
        textOutput("saved_text")
      ),
      uiOutput("outline_warnings"),
      layout_columns(
        col_widths = c(8, 4),
        card(card_header("Image with detected outline"),
             plotOutput("image_plot", height = "480px", click = "img_click",
                        brush = brushOpts("img_brush", resetOnNew = TRUE))),
        card(card_header("Standardized outline"),
             plotOutput("std_plot", height = "380px"),
             p(class = "small text-muted", "The outline turned tip up. Crosses mark the tip and base. Colours show the four regions used in step 4."))
      )
    )
  ),

  nav_panel(
    "3. Specimens", value = "specimens",
    layout_columns(
      col_widths = c(8, 4),
      card(card_header("Saved specimens"),
           plotOutput("spec_grid", height = "360px"),
           tableOutput("spec_table")),
      card(
        card_header("Edit"),
        selectInput("edit_which", "Specimen", choices = character(0)),
        textInput("edit_name", "Name"),
        textInput("edit_group", "Group label"),
        div(actionButton("edit_apply", "Apply changes", class = "btn-primary btn-sm"),
            actionButton("edit_remove", "Remove this specimen", class = "btn-outline-danger btn-sm")),
        hr(),
        p("Sessions on the server end after a period without activity. Download your outlines to keep them, and load the file later to continue."),
        downloadButton("dl_outlines", "Download outlines (CSV)", class = "btn-sm"),
        fileInput("load_outlines", "Load a saved outline file", accept = ".csv"),
        actionButton("go_results", "Continue to step 4", class = "btn-outline-primary")
      )
    )
  ),

  nav_panel(
    "4. Results", value = "results",
    layout_sidebar(
      sidebar = sidebar(
        width = 300,
        textOutput("n_text"),
        selectInput("colour_by", "Colour specimens by",
                    choices = c("Group label" = "group", "Cluster" = "cluster", "Nothing" = "none")),
        checkboxInput("symmetrize", "Analyse symmetric shapes (average each outline with its mirror image)", FALSE),
        checkboxInput("auto_harmonics", "Choose the number of harmonics automatically", TRUE),
        conditionalPanel("!input.auto_harmonics",
                         sliderInput("n_harmonics", "Harmonics", min = 2, max = MAX_HARMONICS, value = 12, step = 1)),
        numericInput("k", "Number of clusters", value = 2, min = 2, step = 1),
        accordion(open = FALSE, accordion_panel(
          "Region boundaries",
          p(class = "small", "Position along each edge, as a percentage of the distance from tip to base."),
          sliderInput("b_tip", "Tip ends at", min = 5, max = 95, value = 15, step = 5, post = "%"),
          sliderInput("b_blade", "Blade ends at", min = 5, max = 95, value = 60, step = 5, post = "%"),
          sliderInput("b_shoulder", "Shoulder ends at", min = 5, max = 95, value = 80, step = 5, post = "%")
        )),
        hr(),
        downloadButton("dl_report", "Report (PDF)", class = "btn-primary btn-sm"),
        downloadButton("dl_meas", "Measurements (CSV)", class = "btn-sm"),
        downloadButton("dl_scores", "PC scores (CSV)", class = "btn-sm"),
        downloadButton("dl_coef", "Fourier coefficients (CSV)", class = "btn-sm")
      ),
      navset_card_tab(
        id = "results_tabs",
        nav_panel("Outlines", result_row("stack_plot", "outlines_text")),
        nav_panel("Measurements", result_row("meas_plot", "meas_text", "360px"), card(tableOutput("meas_table"))),
        nav_panel("Harmonics", result_row("power_plot", "harm_text", "340px"),
                  card(selectInput("recon_which", "Specimen to rebuild", choices = character(0)),
                       plotOutput("recon_plot", height = "240px"))),
        nav_panel("PCA", result_row("morpho_plot", "pca_text"),
                  layout_columns(col_widths = c(4, 8), card(plotOutput("scree_plot", height = "280px")),
                                 result_row("stab_plot", "stab_text", "280px"))),
        nav_panel("Shape change", result_row("pcs_plot", "pcs_text")),
        nav_panel("Clusters", result_row("dend_plot", "clust_text"),
                  layout_columns(col_widths = c(6, 6), card(plotOutput("sil_plot", height = "260px")),
                                 card(card_header("Clusters against group labels"), tableOutput("cross_table")))),
        nav_panel("Regions",
                  layout_columns(col_widths = c(4, 4, 4),
                                 card(plotOutput("regions_plot", height = "320px")),
                                 card(plotOutput("pvar_plot", height = "320px")),
                                 card(plotOutput("rvar_plot", height = "320px"))),
                  result_row("rv_plot", "regions_text", "340px"))
      )
    )
  ),

  nav_panel(
    "Help", value = "help",
    card(
      card_header("What the app does"),
      tags$dl(
        tags$dt("Finding the outline"),
        tags$dd("The image is converted to grey levels and divided into light and dark at a threshold chosen by Otsu's method, which picks the level that best separates the two. The artifact is taken to be the largest light or dark region that does not touch the edge of the image or crop. Holes inside it, such as inked catalogue numbers, are filled. Images larger than 1200 pixels on their longer side are reduced to that size."),
        tags$dt("Orientation, tip, and base"),
        tags$dd("You state which way the tip points. The app computes the outline's axis of elongation nearest that direction and turns the outline so the axis is vertical with the tip up. The tip is the outline point farthest along the axis. The base is the point where the midline (the line through the centre of area, parallel to the axis) meets the far end of the outline. You can set both by clicking if the automatic choice is wrong."),
        tags$dt("Corresponding points"),
        tags$dd("The tip and base divide the outline into two edges. Each edge is represented by 100 points spaced evenly along its length. Point 30 on an edge is 30 percent of the way from tip to base on every specimen. This allows specimens to be compared point by point, but an equal fraction of edge length does not always fall on the same feature of different artifacts."),
        tags$dt("Alignment"),
        tags$dd("Generalized Procrustes alignment centres every outline, scales it to the same size, and rotates it to best fit the mean of all outlines. After this step outlines differ only in shape."),
        tags$dt("Elliptical Fourier analysis and principal components"),
        tags$dd("Each aligned outline is described by Fourier coefficients. Principal component analysis of those coefficients finds the main ways in which the outlines differ. Because orientation, size, and the starting point are already fixed by the steps above, the coefficients are used as computed, without further normalization."),
        tags$dt("Clusters"),
        tags$dd("Ward's hierarchical clustering groups specimens by the distances between their coefficient sets."),
        tags$dt("Regions"),
        tags$dd("Each edge is divided into tip, blade, shoulder, and base by position along its length. The app reports how variable each region is and how strongly pairs of regions vary together.")
      ),
      h6("Limits to keep in mind"),
      tags$ul(
        tags$li("The outline records the plan shape only. Thickness, cross-section, flaking, and raw material are not measured."),
        tags$li("Breaks and resharpening change outlines. The analysis cannot tell these apart from differences in the original form."),
        tags$li("Which face was photographed determines left and right. Use the symmetric option if this matters for your question."),
        tags$li("Results from a few specimens can change substantially when one specimen is added or removed. The PCA tab shows how much."),
        tags$li("A camera stores an orientation flag that some programs apply and others ignore, so a photograph can appear turned here. Set the tip direction to match what you see on the page.")
      )
    )
  ),
  nav_spacer(),
  nav_item(tags$a("Source code", href = "https://github.com/clipo/bifaceAnalysis", target = "_blank"))
)

server <- function(input, output, session) {
  rv <- reactiveValues(
    idx = 1L,
    state = list(),          # per image: crop, scale (2 x 2), scale_mm, tip, base
    specimens = list(),      # saved specimens, named by id
    shared_ppm = NULL,       # shared scale in ORIGINAL pixels per mm
    last_saved = ""
  )

  notify <- function(msg, type = "message") showNotification(msg, type = type, duration = 6)
  plural <- function(n) if (n == 1) "" else "s"

  # ---- images ---------------------------------------------------------------
  image_files <- reactive({
    if (identical(input$image_source, "upload")) {
      if (is.null(input$images)) return(data.frame(name = character(0), path = character(0)))
      data.frame(name = input$images$name, path = input$images$datapath, stringsAsFactors = FALSE)
    } else {
      f <- list.files(DEMO_DIR, pattern = IMAGE_PATTERN, full.names = TRUE)
      data.frame(name = basename(f), path = f, stringsAsFactors = FALSE)
    }
  })
  observeEvent(image_files(), { rv$idx <- 1L })
  output$image_count <- renderText({
    n <- nrow(image_files())
    if (n == 0) "No images yet. Upload one or more images." else sprintf("%d image%s ready.", n, plural(n))
  })
  observeEvent(input$go_outline, nav_select("nav", "outline"))
  observeEvent(input$go_results, nav_select("nav", "results"))

  current_id <- reactive({
    f <- image_files()
    if (nrow(f) == 0) return(NULL)
    f$name[min(rv$idx, nrow(f))]
  })
  current_image <- reactive({
    f <- image_files()
    req(nrow(f) > 0)
    path <- f$path[min(rv$idx, nrow(f))]
    tryCatch(load_working_image(path), error = function(e) {
      notify(paste("This image could not be read:", conditionMessage(e)), "error")
      NULL
    })
  })
  state <- reactive({
    id <- current_id()
    if (is.null(id)) return(list())
    s <- rv$state[[id]]
    if (is.null(s)) list() else s
  })
  set_state <- function(field, value) {
    id <- isolate(current_id())
    if (is.null(id)) return()
    s <- isolate(rv$state[[id]])
    if (is.null(s)) s <- list()
    s[field] <- list(value)
    rv$state[[id]] <- s
  }

  step_image <- function(by) {
    n <- nrow(image_files())
    if (n > 0) rv$idx <- min(max(rv$idx + by, 1L), n)
  }
  observeEvent(input$prev_img, step_image(-1L))
  observeEvent(input$next_img, step_image(1L))
  output$image_title <- renderText({
    f <- image_files()
    if (nrow(f) == 0) return("No images. Go to step 1.")
    sprintf("Image %d of %d: %s", min(rv$idx, nrow(f)), nrow(f), current_id())
  })
  observeEvent(current_id(), {
    id <- current_id()
    saved <- rv$specimens[[id]]
    updateTextInput(session, "spec_name", value = if (is.null(saved)) tools::file_path_sans_ext(id) else saved$name)
    if (!is.null(saved)) updateTextInput(session, "spec_group", value = saved$group)
  })

  # ---- interaction with the image ------------------------------------------
  inside <- function(p, img) p[1] >= 0 && p[1] <= img$width && p[2] >= 0 && p[2] <= img$height

  observeEvent(input$img_brush, {
    b <- input$img_brush
    img <- current_image()
    if (!is.null(b) && !is.null(img) && identical(input$mode, "crop")) {
      crop <- c(max(0, b$xmin), min(img$width, b$xmax), max(0, b$ymin), min(img$height, b$ymax))
      if (crop[2] - crop[1] >= 20 && crop[4] - crop[3] >= 20) {
        set_state("crop", crop)
        set_state("tip", NULL); set_state("base", NULL)
      }
    }
    session$resetBrush("img_brush")
  })
  observeEvent(input$img_click, {
    img <- current_image()
    p <- c(input$img_click$x, input$img_click$y)
    if (is.null(img) || length(p) != 2 || !inside(p, img)) return()
    if (identical(input$mode, "scale")) {
      sc <- state()$scale
      set_state("scale", if (is.null(sc) || nrow(sc) >= 2) matrix(p, 1) else rbind(sc, p))
      set_state("scale_mm", input$scale_mm)
    } else if (identical(input$mode, "landmarks")) {
      s <- state()
      if (is.null(s$tip) || !is.null(s$base)) { set_state("tip", p); set_state("base", NULL) }
      else set_state("base", p)
    }
  })
  observeEvent(input$scale_mm, {
    if (!is.null(isolate(state()$scale))) set_state("scale_mm", input$scale_mm)
  }, ignoreInit = TRUE)
  observeEvent(input$clear_crop, { set_state("crop", NULL); set_state("tip", NULL); set_state("base", NULL) })
  observeEvent(input$clear_landmarks, { set_state("tip", NULL); set_state("base", NULL) })
  observeEvent(input$clear_scale, set_state("scale", NULL))

  # Scale of the current image in working pixels per mm, from its own clicks.
  own_ppm <- reactive({
    s <- state()
    mm <- s$scale_mm
    if (is.null(s$scale) || nrow(s$scale) < 2 || is.null(mm) || !is.finite(mm) || mm <= 0) return(NA_real_)
    d <- sqrt(sum((s$scale[1, ] - s$scale[2, ])^2))
    if (d < 2) NA_real_ else d / mm
  })
  # The shared scale is kept in original pixels so it transfers between images
  # that were resized by different factors.
  observeEvent(own_ppm(), {
    img <- current_image()
    if (is.finite(own_ppm()) && !is.null(img)) rv$shared_ppm <- ppm_to_original(own_ppm(), img$factor)
  })
  current_ppm <- reactive({
    if (is.finite(own_ppm())) return(own_ppm())
    img <- current_image()
    if (isTRUE(input$scale_all) && !is.null(rv$shared_ppm) && !is.null(img)) return(ppm_to_working(rv$shared_ppm, img$factor))
    NA_real_
  })
  output$scale_text <- renderText({
    s <- state()
    if (is.finite(own_ppm())) sprintf("Scale for this image: %.2f pixels per mm.", own_ppm())
    else if (!is.null(s$scale) && nrow(s$scale) == 1) "One point set. Click the other end of the known distance."
    else if (is.finite(current_ppm())) sprintf("Using the shared scale: %.2f pixels per mm.", current_ppm())
    else "No scale set. Measurements will be in pixels."
  })

  # ---- processing (preview only; nothing here changes a saved specimen) -----
  settings <- debounce(reactive({
    list(auto = isTRUE(input$auto_threshold), threshold = input$threshold, polarity = input$polarity,
         speck = input$speck, direction = input$direction)
  }), 250)

  processed <- reactive({
    img <- current_image()
    req(img)
    s <- state(); st <- settings()
    both <- !is.null(s$tip) && !is.null(s$base)
    process_image(img, crop = s$crop, threshold = if (st$auto) NULL else st$threshold,
                  polarity = st$polarity, speck_side = st$speck, direction = st$direction,
                  tip = if (both) s$tip else NULL, base = if (both) s$base else NULL)
  })
  observeEvent(processed(), {
    p <- processed()
    if (isTRUE(input$auto_threshold) && !is.null(p$threshold))
      updateSliderInput(session, "threshold", value = round(p$threshold, 2))
  })

  output$outline_warnings <- renderUI({
    if (nrow(image_files()) == 0)
      return(div(class = "warn-box", "No images are loaded. Go to step 1 to choose images."))
    p <- processed()
    s <- state()
    msgs <- p$warnings
    if (!is.null(s$tip) && is.null(s$base)) msgs <- c(msgs, "Tip set. Now click the base to finish placing the landmarks.")
    if (length(msgs) == 0) return(NULL)
    tagList(lapply(msgs, function(m) div(class = "warn-box", m)))
  })

  output$image_plot <- renderPlot({
    img <- current_image()
    req(img)
    p <- processed(); s <- state()
    op <- par(mar = c(0, 0, 0, 0)); on.exit(par(op))
    plot(NA, xlim = c(0, img$width), ylim = c(img$height, 0), asp = 1, xaxs = "i", yaxs = "i",
         axes = FALSE, xlab = "", ylab = "")
    rasterImage(img$raster, 0, img$height, img$width, 0, interpolate = FALSE)
    if (!is.null(s$crop)) rect(s$crop[1], s$crop[4], s$crop[2], s$crop[3], border = "#F0E442", lty = 2, lwd = 2)
    if (!is.null(p$outline_plot)) polygon(p$outline_plot, border = "#D55E00", lwd = 2.5)
    if (isTRUE(p$ok)) {
      pts <- rbind(p$tip_plot, p$base_plot)
      points(pts, pch = 4, cex = 2.2, lwd = 3, col = "#56B4E9")
      text(pts, labels = c("tip", "base"), pos = 3, col = "#56B4E9", font = 2, cex = 1.1)
    } else if (!is.null(s$tip)) points(s$tip[1], s$tip[2], pch = 4, cex = 2.2, lwd = 3, col = "#56B4E9")
    if (!is.null(s$scale)) {
      points(s$scale, pch = 3, cex = 2, lwd = 3, col = "#009E73")
      if (nrow(s$scale) == 2) lines(s$scale, col = "#009E73", lwd = 3)
    }
  })

  output$std_plot <- renderPlot({
    p <- processed()
    validate(need(isTRUE(p$ok), "No outline to show yet."))
    reg <- assign_regions(point_positions(), boundaries())
    op <- par(mar = c(0, 0, 0, 0)); on.exit(par(op))
    plot(p$std$outline, type = "n", asp = 1, axes = FALSE, xlab = "", ylab = "")
    polygon(p$std$outline, col = "grey93", border = "grey50")
    points(p$std$points, pch = 19, cex = 0.6, col = REGION_COLOURS[as.character(reg)])
    points(rbind(p$std$tip, p$std$base), pch = 4, cex = 2, lwd = 3)
  })

  # ---- specimens ------------------------------------------------------------
  add_specimen <- function(id, name, group, p, ppm) {
    others <- rv$specimens[setdiff(names(rv$specimens), id)]
    rv$specimens[[id]] <- new_specimen(id, unique_specimen_name(name, others), trimws(group),
                                       p$std$outline, p$std$base_row, ppm, p$settings)
  }
  observeEvent(input$save_spec, {
    id <- current_id()
    p <- if (is.null(id)) NULL else processed()
    if (is.null(id) || !isTRUE(p$ok)) {
      notify("There is no usable outline to save for this image. Read the message above the image.", "error")
      return()
    }
    replaced <- !is.null(rv$specimens[[id]])
    add_specimen(id, input$spec_name, input$spec_group, p, current_ppm())
    n <- length(rv$specimens)
    rv$last_saved <- sprintf("%s %s. %d specimen%s saved in total.", if (replaced) "Replaced" else "Saved",
                             rv$specimens[[id]]$name, n, plural(n))
    if (rv$idx < nrow(image_files())) rv$idx <- rv$idx + 1L
    else notify("That was the last image. Continue to step 3 to review the specimens.")
  })
  output$saved_text <- renderText(rv$last_saved)

  observeEvent(input$run_demo, {
    f <- list.files(DEMO_DIR, pattern = IMAGE_PATTERN, full.names = TRUE)
    if (length(f) == 0) { notify("The sample images are missing.", "error"); return() }
    withProgress(message = "Processing sample images", value = 0, {
      for (path in f) {
        p <- process_image(load_working_image(path), direction = DEMO_DIRECTION)
        id <- basename(path)
        if (isTRUE(p$ok)) add_specimen(id, tools::file_path_sans_ext(id), "", p, NA_real_)
        incProgress(1 / length(f))
      }
    })
    updateRadioButtons(session, "image_source", selected = "demo")
    nav_select("nav", "results")
  })

  observe({
    ids <- names(rv$specimens)
    nm <- vapply(unname(rv$specimens), function(s) s$name, character(1))
    ch <- if (length(ids)) stats::setNames(ids, nm) else character(0)
    updateSelectInput(session, "edit_which", choices = ch, selected = isolate(input$edit_which))
    updateSelectInput(session, "recon_which", choices = unname(nm), selected = isolate(input$recon_which))
  })
  observeEvent(input$edit_which, {
    s <- rv$specimens[[input$edit_which]]
    req(s)
    updateTextInput(session, "edit_name", value = s$name)
    updateTextInput(session, "edit_group", value = s$group)
  })
  observeEvent(input$edit_apply, {
    id <- input$edit_which
    req(id, rv$specimens[[id]])
    rv$specimens[[id]]$name <- unique_specimen_name(input$edit_name, rv$specimens[setdiff(names(rv$specimens), id)])
    rv$specimens[[id]]$group <- trimws(input$edit_group)
  })
  observeEvent(input$edit_remove, {
    id <- input$edit_which
    req(id)
    rv$specimens[[id]] <- NULL
  })
  observeEvent(input$load_outlines, {
    loaded <- tryCatch(specimens_from_csv(input$load_outlines$datapath), error = function(e) e)
    if (inherits(loaded, "error")) { notify(conditionMessage(loaded), "error"); return() }
    rv$specimens <- merge_specimens(rv$specimens, loaded)
    notify(sprintf("Loaded %d specimen%s.", length(loaded), plural(length(loaded))))
  })

  specimens <- reactive({
    s <- rv$specimens
    stats::setNames(s, vapply(unname(s), function(z) z$name, character(1)))
  })
  n_spec <- reactive(length(specimens()))
  groups <- reactive(vapply(specimens(), function(s) s$group, character(1)))
  measurements <- reactive(specimen_measurements(specimens()))

  measurement_table <- function(m) {
    any_mm <- any(is.finite(m$length_mm)); all_mm <- all(is.finite(m$length_mm))
    out <- data.frame(Specimen = m$specimen, Group = m$group, check.names = FALSE, stringsAsFactors = FALSE)
    if (any_mm) {
      out[["Length (mm)"]] <- m$length_mm; out[["Width (mm)"]] <- m$width_mm; out[["Area (mm2)"]] <- m$area_mm2
    }
    if (!all_mm) { out[["Length (px)"]] <- m$length_px; out[["Width (px)"]] <- m$width_px }
    out[["Length / width"]] <- m$length_width_ratio
    out[["Widest at (fraction of length above base)"]] <- m$widest_position
    out[["Asymmetry"]] <- m$asymmetry
    out
  }
  output$spec_grid <- renderPlot({
    validate(need(n_spec() > 0, "No specimens saved yet. Save outlines in step 2, or load an outline file."))
    plot_outline_grid(lapply(specimens(), function(s) s$outline), groups())
  })
  output$spec_table <- renderTable({
    req(n_spec() > 0)
    measurement_table(measurements())
  }, digits = 3, na = "")

  # ---- analysis (settings here never alter a saved specimen) ----------------
  need_n <- function(n, what) {
    validate(need(n_spec() >= n, sprintf("%s needs at least %d specimens. You have %d.", what, n, n_spec())))
  }
  need_pca <- function(min_rank = 1L) {
    need_n(3, "Principal component analysis")
    validate(need(pca()$rank >= 1L, "The saved outlines are identical in shape, so there is no shape variation to analyse."),
             need(pca()$rank >= min_rank, sprintf("This result needs at least %d components of shape variation. These specimens have %d.", min_rank, pca()$rank)))
  }
  boundaries <- reactive(sort(c(input$b_tip, input$b_blade, input$b_shoulder)) / 100)
  shapes <- reactive({
    pts <- lapply(specimens(), function(s) s$points)
    if (isTRUE(input$symmetrize)) lapply(pts, symmetrize_shape, pair = mirror_index()) else pts
  })
  aligned <- reactive({ req(n_spec() >= 2); gpa(shapes()) })
  power <- reactive(harmonic_power(efa_matrix(aligned()$aligned, MAX_HARMONICS)))
  n_harm <- reactive({
    if (isTRUE(input$auto_harmonics)) return(choose_harmonics(power()))
    h <- suppressWarnings(as.integer(input$n_harmonics))
    if (is.na(h)) h <- 12L
    min(max(h, 2L), MAX_HARMONICS)
  })
  coefs <- reactive(efa_matrix(aligned()$aligned, n_harm()))
  pca <- reactive({ req(n_spec() >= 3); shape_pca(coefs()) })
  tree <- reactive({ req(n_spec() >= 4); shape_clusters(coefs()) })
  k_used <- reactive({
    k <- suppressWarnings(as.integer(input$k))
    if (length(k) == 0 || is.na(k)) k <- 2L
    min(max(k, 2L), n_spec() - 1L)
  })
  clusters <- reactive(stats::cutree(tree(), k_used()))
  sil <- reactive(silhouette_by_k(coefs(), tree()))
  regions <- reactive(assign_regions(point_positions(), boundaries()))
  pvar <- reactive(point_variance(aligned()$aligned))
  region_var <- reactive(tapply(pvar(), regions(), mean))
  rv_table <- reactive({
    req(n_spec() >= 4, nlevels(droplevels(regions())) >= 2)
    region_rv(aligned()$aligned, regions())
  })
  stability <- reactive({ req(n_spec() >= 4); pca_stability(coefs()) })
  colouring <- reactive({
    if (identical(input$colour_by, "group") && any(nzchar(groups()))) list(values = groups(), name = "Group")
    else if (identical(input$colour_by, "cluster") && n_spec() >= 4) list(values = as.character(clusters()), name = "Cluster")
    else list(values = NULL, name = NULL)
  })

  output$n_text <- renderText(sprintf("%d specimen%s in the analysis.", n_spec(), plural(n_spec())))

  figures <- list(
    stack = reactive({ need_n(2, "The overlay"); plot_outline_stack(aligned()$aligned, aligned()$mean, groups()) }),
    meas = reactive({ need_n(1, "The measurement plot"); plot_measurements(measurements(), groups()) }),
    power = reactive({ need_n(2, "Harmonic analysis"); plot_harmonic_power(power(), n_harm()) }),
    recon = reactive({
      need_n(2, "Harmonic analysis")
      al <- aligned()$aligned
      which <- if (!is.null(input$recon_which) && input$recon_which %in% names(al)) input$recon_which else names(al)[1]
      plot_reconstruction(al[[which]], c(1, 2, 4, 8, n_harm()))
    }),
    morpho = reactive({ need_pca(2L); plot_morphospace(pca(), colouring()$values, colouring()$name) }),
    scree = reactive({ need_pca(1L); plot_scree(pca()) }),
    stab = reactive({
      need_n(4, "The leave-one-out check"); need_pca(2L)
      plot_stability(stability())
    }),
    pcs = reactive({ need_pca(1L); plot_pc_shapes(pca()) }),
    dend = reactive({ need_n(4, "Clustering"); plot_dendrogram(tree(), k_used()) }),
    sil = reactive({ need_n(4, "Clustering"); plot_silhouette(sil(), k_used()) }),
    regions = reactive({ need_n(2, "Region analysis"); plot_regions(aligned()$mean, regions()) }),
    pvar = reactive({ need_n(2, "Region analysis"); plot_point_variance(aligned()$mean, pvar()) }),
    rvar = reactive({ need_n(2, "Region analysis"); plot_region_variance(region_var()) }),
    rv = reactive({
      need_n(4, "Covariation between regions")
      validate(need(any(is.finite(rv_table()$rv)), "Covariation between regions is not available because at least one region does not vary between these specimens."))
      plot_region_rv(rv_table())
    })
  )
  texts <- list(
    outlines = reactive({ req(n_spec() >= 2); text_outlines(n_spec(), isTRUE(input$symmetrize)) }),
    meas = reactive({ req(n_spec() >= 1); text_measurements(measurements()) }),
    harm = reactive({ req(n_spec() >= 2); text_harmonics(power(), n_harm(), isTRUE(input$auto_harmonics)) }),
    pca = reactive({ req(n_spec() >= 3, pca()$rank >= 2L); text_pca(pca(), n_spec()) }),
    stab = reactive({ req(n_spec() >= 4, !is.null(stability())); text_stability(stability()) }),
    pcs = reactive({ req(n_spec() >= 3, pca()$rank >= 1L); text_pc_shapes(pca()) }),
    clust = reactive({ req(n_spec() >= 4); text_clusters(tree(), k_used(), sil(), groups()) }),
    regions = reactive({
      req(n_spec() >= 2)
      text_regions(region_var(), if (n_spec() >= 4) rv_table() else NULL, n_spec(), boundaries())
    })
  )

  output$stack_plot <- renderPlot(figures$stack())
  output$meas_plot <- renderPlot(figures$meas())
  output$power_plot <- renderPlot(figures$power())
  output$recon_plot <- renderPlot(figures$recon())
  output$morpho_plot <- renderPlot(figures$morpho())
  output$scree_plot <- renderPlot(figures$scree())
  output$stab_plot <- renderPlot(figures$stab())
  output$pcs_plot <- renderPlot(figures$pcs())
  output$dend_plot <- renderPlot(figures$dend())
  output$sil_plot <- renderPlot(figures$sil())
  output$regions_plot <- renderPlot(figures$regions())
  output$pvar_plot <- renderPlot(figures$pvar())
  output$rvar_plot <- renderPlot(figures$rvar())
  output$rv_plot <- renderPlot(figures$rv())

  output$outlines_text <- renderUI(paragraphs(texts$outlines()))
  output$meas_text <- renderUI(paragraphs(texts$meas()))
  output$harm_text <- renderUI(paragraphs(texts$harm()))
  output$pca_text <- renderUI(paragraphs(texts$pca()))
  output$stab_text <- renderUI(paragraphs(texts$stab()))
  output$pcs_text <- renderUI(paragraphs(texts$pcs()))
  output$clust_text <- renderUI(paragraphs(texts$clust()))
  output$regions_text <- renderUI(paragraphs(texts$regions()))

  output$meas_table <- renderTable({ req(n_spec() > 0); measurement_table(measurements()) }, digits = 3, na = "")
  output$cross_table <- renderTable({
    need_n(4, "Clustering")
    validate(need(any(nzchar(groups())), "Enter group labels in step 3 to compare them with the clusters."))
    g <- ifelse(nzchar(groups()), groups(), "(no label)")
    tab <- as.data.frame.matrix(table(Cluster = clusters(), Group = g))
    cbind(Cluster = rownames(tab), tab)
  })

  # ---- downloads ------------------------------------------------------------
  stamp <- function(prefix, ext) function() sprintf("%s_%s.%s", prefix, format(Sys.time(), "%Y%m%d-%H%M%S"), ext)
  empty_csv <- function(file) utils::write.csv(data.frame(message = "No results are available yet."), file, row.names = FALSE)

  output$dl_outlines <- downloadHandler(stamp("biface_outlines", "csv"), function(file) {
    if (n_spec() == 0) empty_csv(file) else utils::write.csv(specimens_to_csv(rv$specimens), file, row.names = FALSE)
  })
  output$dl_meas <- downloadHandler(stamp("biface_measurements", "csv"), function(file) {
    if (n_spec() == 0) empty_csv(file) else utils::write.csv(measurements(), file, row.names = FALSE)
  })
  output$dl_scores <- downloadHandler(stamp("biface_pc_scores", "csv"), function(file) {
    if (n_spec() < 3 || pca()$rank < 1L) return(empty_csv(file))
    out <- data.frame(specimen = rownames(pca()$scores), group = groups(), stringsAsFactors = FALSE)
    if (n_spec() >= 4) out$cluster <- clusters()
    utils::write.csv(cbind(out, pca()$scores), file, row.names = FALSE)
  })
  output$dl_coef <- downloadHandler(stamp("biface_fourier_coefficients", "csv"), function(file) {
    if (n_spec() < 2) return(empty_csv(file))
    utils::write.csv(cbind(data.frame(specimen = rownames(coefs()), stringsAsFactors = FALSE), coefs()), file, row.names = FALSE)
  })

  # A figure or text that is unavailable (too few specimens) is left out.
  available <- function(expr) tryCatch(expr, error = function(e) NULL)
  specimen_line <- function(s) {
    st <- s$settings
    scale <- if (is.finite(s$px_per_mm)) sprintf("%.2f px per mm", s$px_per_mm) else "no scale"
    if (is.null(st$threshold)) return(sprintf("%s: processing settings not recorded, %s.", s$name, scale))
    sprintf("%s: threshold %.2f, artifact %s than background, speck removal %d px, tip %s, %s, %s, %s.",
            s$name, st$threshold, if (identical(st$polarity, "light")) "lighter" else "darker",
            as.integer(st$speck), st$direction,
            if (isTRUE(st$manual_landmarks)) "tip and base set by hand" else "automatic tip and base",
            if (isTRUE(st$cropped)) "cropped" else "whole image", scale)
  }
  output$dl_report <- downloadHandler(stamp("biface_report", "pdf"), function(file) {
    n <- n_spec()
    cover <- sprintf("Generated %s. %d specimen%s.", format(Sys.time(), "%Y-%m-%d %H:%M"), n, plural(n))
    sections <- list()
    if (n == 0) {
      cover <- c(cover, "No specimens have been saved, so there are no results to report.")
    } else {
      cover <- c(cover,
                 sprintf("Settings: %s; %s; region boundaries at %s of edge length from the tip.",
                         if (n >= 2) sprintf("%d harmonics (%s)", n_harm(), if (isTRUE(input$auto_harmonics)) "automatic" else "set by hand") else "harmonics not applicable",
                         if (isTRUE(input$symmetrize)) "symmetric shapes" else "shapes as photographed",
                         paste(fmt_pct(boundaries()), collapse = ", ")),
                 "Processing of each specimen:",
                 vapply(specimens(), specimen_line, character(1)),
                 sprintf("Software: R %s, shiny %s, imager %s, ggplot2 %s.", getRversion(), utils::packageVersion("shiny"),
                         utils::packageVersion("imager"), utils::packageVersion("ggplot2")),
                 "Methods are described on the Help page of the app.")
      add <- function(heading, fig, txt = NULL) {
        p <- available(figures[[fig]]())
        if (is.null(p)) return()
        tx <- if (is.null(txt)) NULL else available(texts[[txt]]())
        sections[[length(sections) + 1L]] <<- list(heading = heading, plot = p, text = if (is.null(tx)) character(0) else tx)
      }
      sections[[1]] <- list(heading = "Specimens", text = character(0),
                            plot = plot_outline_grid(lapply(specimens(), function(s) s$outline), groups()))
      add("Aligned outlines", "stack", "outlines")
      add("Measurements", "meas", "meas")
      add("Harmonic power", "power", "harm")
      add("Reconstruction from harmonics", "recon")
      add("Principal components", "morpho", "pca")
      add("Share of variance by component", "scree")
      add("Dependence of the components on single specimens", "stab", "stab")
      add("Shape change along components", "pcs", "pcs")
      add("Clusters", "dend", "clust")
      add("Silhouette width by number of clusters", "sil")
      add("Regions", "regions", "regions")
      add("Variance at each point", "pvar")
      add("Variance by region", "rvar")
      add("Covariation between regions", "rv")
    }
    write_report(file, cover, sections)
  })
}

shinyApp(ui, server)
