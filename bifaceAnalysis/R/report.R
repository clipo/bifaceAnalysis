# PDF report: the figures shown on screen with their interpretation text.
# Text is wrapped to the measured page width and continued on further pages
# when it does not fit, so no line is dropped.

REPORT <- list(page_w = 8.5, page_h = 11, body_w = 7.3, body_h = 10, heading_h = 0.5,
               line_h = 0.165, pad = 0.2, fontsize = 9.5, min_plot_h = 4.5)

report_gpar <- function() grid::gpar(fontsize = REPORT$fontsize, fontfamily = "sans", lineheight = 1.1)

# Width in inches of each string in the report's text font. Needs an open device.
text_width <- function(x) {
  if (length(x) == 0L) return(numeric(0))
  grid::pushViewport(grid::viewport(gp = report_gpar()))
  on.exit(grid::popViewport())
  vapply(x, function(s) grid::convertWidth(grid::stringWidth(s), "inches", valueOnly = TRUE), numeric(1), USE.NAMES = FALSE)
}

# Wrap paragraphs into lines no wider than `max_width` inches, with an empty
# line after each paragraph. A word wider than the line is split.
wrap_report_text <- function(paragraphs, max_width = REPORT$body_w - 0.1) {
  out <- character(0)
  space <- text_width(" ")
  for (p in paragraphs) {
    words <- strsplit(trimws(p), "[[:space:]]+")[[1]]
    words <- words[nzchar(words)]
    pieces <- character(0)
    for (w in words) {                       # split words that cannot fit on a line
      while (text_width(w) > max_width) {
        k <- nchar(w) - 1L
        while (k > 1L && text_width(substr(w, 1L, k)) > max_width) k <- k - 1L
        pieces <- c(pieces, substr(w, 1L, k)); w <- substr(w, k + 1L, nchar(w))
      }
      pieces <- c(pieces, w)
    }
    widths <- text_width(pieces)
    line <- character(0); used <- 0
    for (i in seq_along(pieces)) {
      add <- widths[i] + if (length(line)) space else 0
      if (length(line) && used + add > max_width) {
        out <- c(out, paste(line, collapse = " ")); line <- character(0); used <- 0; add <- widths[i]
      }
      line <- c(line, pieces[i]); used <- used + add
    }
    out <- c(out, if (length(line)) paste(line, collapse = " "), "")
  }
  out
}

# Number of text lines that fit on a page with and without a figure.
report_capacity <- function(with_plot) {
  space <- REPORT$body_h - REPORT$heading_h - REPORT$pad - if (with_plot) REPORT$min_plot_h else 0
  as.integer(floor(space / REPORT$line_h))
}

# Draw one page: heading, optional figure, and lines of text that fit.
report_page <- function(heading, plot = NULL, lines = character(0)) {
  grid::grid.newpage()
  text_height <- grid::unit(if (length(lines)) length(lines) * REPORT$line_h + REPORT$pad else 0, "inches")
  lay <- grid::grid.layout(3, 1, heights = grid::unit.c(grid::unit(REPORT$heading_h, "inches"), grid::unit(1, "null"), text_height))
  grid::pushViewport(grid::viewport(width = grid::unit(REPORT$body_w, "inches"), height = grid::unit(REPORT$body_h, "inches"), layout = lay))
  grid::pushViewport(grid::viewport(layout.pos.row = 1))
  grid::grid.text(heading, x = 0, just = "left", gp = grid::gpar(fontsize = 15, fontface = "bold", fontfamily = "sans"))
  grid::popViewport()
  if (!is.null(plot)) {
    grid::pushViewport(grid::viewport(layout.pos.row = 2))
    print(plot, newpage = FALSE)
    grid::popViewport()
  }
  if (length(lines)) {
    grid::pushViewport(grid::viewport(layout.pos.row = 3))
    grid::grid.text(paste(lines, collapse = "\n"), x = 0, y = 1, just = c("left", "top"), gp = report_gpar())
    grid::popViewport()
  }
  grid::popViewport()
}

# Draw a section on as many pages as its text needs. Returns the page count.
report_section <- function(heading, plot = NULL, paragraphs = character(0)) {
  lines <- wrap_report_text(paragraphs)
  if (length(lines) && !nzchar(lines[length(lines)])) lines <- lines[-length(lines)]
  first <- report_capacity(!is.null(plot))
  report_page(heading, plot, utils::head(lines, first))
  rest <- lines[-seq_len(min(first, length(lines)))]
  pages <- 1L
  per_page <- report_capacity(FALSE)
  while (length(rest)) {
    if (!nzchar(rest[1])) { rest <- rest[-1]; next }
    report_page(paste(heading, "(continued)"), NULL, utils::head(rest, per_page))
    rest <- rest[-seq_len(min(per_page, length(rest)))]
    pages <- pages + 1L
  }
  pages
}

# cover: character vector of paragraphs. sections: list of
# list(heading, plot, text). Writes the cover and one or more pages per section.
write_report <- function(file, cover, sections) {
  grDevices::pdf(file, width = REPORT$page_w, height = REPORT$page_h, family = "Helvetica")
  on.exit(grDevices::dev.off())
  report_section("Biface outline analysis", NULL, cover)
  for (s in sections) report_section(s$heading, s$plot, s$text)
  invisible(file)
}
