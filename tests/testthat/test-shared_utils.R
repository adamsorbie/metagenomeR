# Tests for shared_utils.R functions

test_that("t_df transposes data frame correctly", {
  df <- data.frame(
    a = 1:3,
    b = 4:6,
    row.names = c("row1", "row2", "row3")
  )

  result <- metagenomeR:::t_df(df)

  expect_s3_class(result, "data.frame")
  expect_equal(nrow(result), 2)  # original columns become rows

expect_equal(ncol(result), 3)  # original rows become columns
})

test_that("read_tab_delim_metag reads tab-delimited files", {
  # Create a temporary file
  tmp_file <- tempfile(fileext = ".tsv")
  test_data <- data.frame(
    sample1 = c(10, 20, 30),
    sample2 = c(15, 25, 35),
    sample3 = c(12, 22, 32),
    row.names = c("taxon_a", "taxon_b", "taxon_c")
  )

  # Write with row names as first column (tab-delimited format)
  write.table(
    test_data,
    file = tmp_file,
    sep = "\t",
    quote = FALSE,
    row.names = TRUE,
    col.names = NA
  )

  result <- read_tab_delim_metag(tmp_file)

  expect_s3_class(result, "data.frame")
  expect_equal(nrow(result), 3)
  expect_equal(ncol(result), 3)
  expect_equal(rownames(result), c("taxon_a", "taxon_b", "taxon_c"))

  # Clean up
  unlink(tmp_file)
})

test_that("filter_rownames filters data frame by row names", {
  df <- data.frame(
    value = c(1, 2, 3, 4, 5),
    row.names = c("a", "b", "c", "d", "e")
  )

  result <- metagenomeR:::filter_rownames(df, c("a", "c", "e"))

  expect_s3_class(result, "data.frame")
  expect_equal(nrow(result), 3)
  expect_equal(rownames(result), c("a", "c", "e"))
  expect_equal(result$value, c(1, 3, 5))
})

test_that("filter_rownames handles empty filter vector", {
  df <- data.frame(
    value = c(1, 2, 3),
    row.names = c("a", "b", "c")
  )

  result <- metagenomeR:::filter_rownames(df, character(0))

  expect_s3_class(result, "data.frame")
  expect_equal(nrow(result), 0)
})

test_that("filter_rownames handles non-matching names", {
  df <- data.frame(
    value = c(1, 2, 3),
    row.names = c("a", "b", "c")
  )

  result <- metagenomeR:::filter_rownames(df, c("x", "y", "z"))

  expect_s3_class(result, "data.frame")
  expect_equal(nrow(result), 0)
})
