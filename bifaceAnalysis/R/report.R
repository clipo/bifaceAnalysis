# PDF report: the figures shown on screen with their interpretation text.

report_page <- function(heading, plot = NULL, paragraphs = character(0)) {
  grid::grid.newpage()
  lines <- unlist(lapply(paragraphs, function(p) c(strwrap(p, width = 112), "")))
  text_height <- grid::unit(if (length(lines)) length(lines) * 0.165 + 0.2 else 0, "inches")
  lay <- grid::grid.layout(3, 1, heights = grid::unit.c(grid::unit(0.5, "inches"), grid::unit(1, "null"), text_height))
  grid::pushViewport(grid::viewport(width = grid::unit(7.3, "inches"), height = grid::unit(10, "inches"), layout = lay))
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
    grid::grid.text(paste(lines, collapse = "\n"), x = 0, y = 1, just = c("left", "top"),
                    gp = grid::gpar(fontsize = 9.5, fontfamily = "sans", lineheight = 1.1))
    grid::popViewport()
  }
  grid::popViewport()
}

# cover: character vector of paragraphs. sections: list of
# list(heading, plot, text). Writes one cover page and one page per section.
write_report <- function(file, cover, sections) {
  grDevices::pdf(file, width = 8.5, height = 11, family = "Helvetica")
  on.exit(grDevices::dev.off())
  report_page("Biface outline analysis", NULL, cover)
  for (s in sections) report_page(s$heading, s$plot, s$text)
  invisible(file)
}
