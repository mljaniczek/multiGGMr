test_that("package sources are ASCII-only", {
  pkg <- test_path("..", "..")
  files <- c(list.files(file.path(pkg, "src"), pattern = "\\.(cpp|h)$", full.names = TRUE),
             list.files(file.path(pkg, "R"), pattern = "\\.R$", full.names = TRUE))
  skip_if(length(files) == 0, "package sources not available (installed tests)")
  bad <- Filter(function(f) any(grepl("[^\x01-\x7F]", readLines(f, warn = FALSE), useBytes = TRUE)), files)
  expect_identical(basename(bad), character(0))
})
