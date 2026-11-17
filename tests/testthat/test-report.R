test_that("report text is wrapped to the page width and never dropped", {
  f <- tempfile(fileext = ".pdf"); grDevices::pdf(f, 8.5, 11, family = "Helvetica"); on.exit({ grDevices::dev.off(); unlink(f) })
  para <- c(paste(rep("WWWW MMMM wide words", 40), collapse = " "), paste(rep("x", 400), collapse = ""), "short")
  lines <- wrap_report_text(para)
  expect_lte(max(text_width(lines)), REPORT$body_w)
  expect_identical(gsub(" ", "", paste(lines, collapse = "")), gsub(" ", "", paste(para, collapse = "")))   # nothing lost
  expect_identical(sum(!nzchar(lines)), 3L)
  expect_lte(report_capacity(FALSE) * REPORT$line_h + REPORT$pad + REPORT$heading_h, REPORT$body_h)
  expect_lte(report_capacity(TRUE) * REPORT$line_h + REPORT$pad + REPORT$heading_h + REPORT$min_plot_h, REPORT$body_h)
})

test_that("a long report keeps every record on continuation pages", {
  n <- 500L
  cover <- c("FIRSTMARK", sprintf("specimen%03d: threshold 0.52, artifact lighter than background, speck removal 3 px, tip left, automatic tip and base, whole image, no scale.", seq_len(n)), "LASTMARK")
  text <- sprintf("paragraph%02d %s", 1:30, paste(rep("interpretation text that fills the line", 6), collapse = " "))
  f <- tempfile(fileext = ".pdf")
  write_report(f, cover, list(list(heading = "Figure page", plot = ggplot2::ggplot(), text = text)))
  pages <- length(grepRaw("/Type /Page /Parent", readBin(f, "raw", file.size(f)), fixed = TRUE, all = TRUE))
  expect_gte(pages, ceiling(2 * n / report_capacity(FALSE)) + 2L)
  skip_if(Sys.which("pdftotext") == "", "pdftotext is not installed")
  txt <- paste(system2("pdftotext", c("-layout", f, "-"), stdout = TRUE), collapse = "\n")
  expect_true(grepl("FIRSTMARK", txt)); expect_true(grepl("LASTMARK", txt))
  found <- unique(unlist(regmatches(txt, gregexpr("specimen[0-9]{3}", txt))))
  expect_identical(sort(found), sprintf("specimen%03d", seq_len(n)))
  expect_identical(sort(unique(unlist(regmatches(txt, gregexpr("paragraph[0-9]{2}", txt))))), sprintf("paragraph%02d", 1:30))
  expect_true(grepl("Figure page (continued)", txt, fixed = TRUE))
})
