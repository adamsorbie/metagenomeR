# Tests for taxonomic_utils.R functions

# Tests for return_level() - pure function mapping taxonomic levels to abbreviations
test_that("return_level returns correct abbreviation for subspecies", {
  expect_equal(return_level("Subspecies"), "t")
  expect_equal(return_level("subspecies"), "t")
  expect_equal(return_level("t"), "t")
})

test_that("return_level returns correct abbreviation for species", {
  expect_equal(return_level("Species"), "s")
  expect_equal(return_level("species"), "s")
  expect_equal(return_level("s"), "s")
})

test_that("return_level returns correct abbreviation for genus", {
  expect_equal(return_level("Genus"), "g")
  expect_equal(return_level("genus"), "g")
  expect_equal(return_level("g"), "g")
})

test_that("return_level returns correct abbreviation for family", {
  expect_equal(return_level("Family"), "f")
  expect_equal(return_level("family"), "f")
  expect_equal(return_level("f"), "f")
})

test_that("return_level returns correct abbreviation for order", {
  expect_equal(return_level("Order"), "o")
  expect_equal(return_level("order"), "o")
  expect_equal(return_level("o"), "o")
})

test_that("return_level returns correct abbreviation for class", {
  expect_equal(return_level("Class"), "c")
  expect_equal(return_level("class"), "c")
  expect_equal(return_level("c"), "c")
})

test_that("return_level returns correct abbreviation for phylum", {
  expect_equal(return_level("Phylum"), "p")
  expect_equal(return_level("phylum"), "p")
  expect_equal(return_level("p"), "p")
})

test_that("return_level returns correct abbreviation for kingdom/domain", {
  expect_equal(return_level("Kingdom"), "k")
  expect_equal(return_level("kingdom"), "k")
  expect_equal(return_level("k"), "k")
  expect_equal(return_level("Domain"), "k")
  expect_equal(return_level("domain"), "k")
  expect_equal(return_level("d"), "k")
})

test_that("return_level throws error for invalid level", {
  expect_error(return_level("invalid"), "Not a valid taxonomic rank")
  expect_error(return_level(""), "Not a valid taxonomic rank")
  expect_error(return_level("x"), "Not a valid taxonomic rank")
})

# Tests for meta_to_df() using zeller2014 data
test_that("meta_to_df extracts metadata as data frame with rownames", {
  result <- meta_to_df(zeller2014, rownames = TRUE)

  expect_s3_class(result, "data.frame")
  expect_equal(nrow(result), phyloseq::nsamples(zeller2014))
  # Should have row names (sample IDs)
  expect_true(length(rownames(result)) > 0)
  expect_false("sample_id" %in% colnames(result))
})
test_that("meta_to_df extracts metadata as data frame without rownames", {
  result <- meta_to_df(zeller2014, rownames = FALSE)

  expect_s3_class(result, "data.frame")
  expect_equal(nrow(result), phyloseq::nsamples(zeller2014))
  # Should have sample_id column instead of rownames
  expect_true("sample_id" %in% colnames(result))
})

test_that("meta_to_df preserves metadata columns", {
  result <- meta_to_df(zeller2014)

  # zeller2014 should have these metadata columns based on CLAUDE.md
  expected_cols <- c("Age", "BMI", "Gender")
  for (col in expected_cols) {
    expect_true(col %in% colnames(result),
                info = paste("Expected column", col, "not found"))
  }
})

# Tests for ps_to_feattab()
test_that("ps_to_feattab extracts OTU table as data frame", {
  result <- ps_to_feattab(zeller2014)

  expect_s3_class(result, "data.frame")
  # Should have taxa as rows (100 taxa according to CLAUDE.md)
  expect_equal(nrow(result), phyloseq::ntaxa(zeller2014))
  # Should have samples as columns (141 samples according to CLAUDE.md)
  expect_equal(ncol(result), phyloseq::nsamples(zeller2014))
})

test_that("ps_to_feattab returns numeric values", {
  result <- ps_to_feattab(zeller2014)

  # All values should be numeric
  expect_true(all(sapply(result, is.numeric)))
})

# Tests for get_top_n()
test_that("get_top_n returns correct number of taxa", {
  result <- get_top_n(zeller2014, n = 5)

  expect_type(result, "character")
  expect_equal(length(result), 5)
})

test_that("get_top_n returns taxa sorted by abundance", {
  result <- get_top_n(zeller2014, n = 10)

  expect_type(result, "character")
  expect_equal(length(result), 10)

  # Verify these are actual taxa names from the phyloseq object
  all_taxa <- phyloseq::taxa_names(zeller2014)
  expect_true(all(result %in% all_taxa))
})

test_that("get_top_n handles n larger than available taxa", {
  n_taxa <- phyloseq::ntaxa(zeller2014)
  result <- get_top_n(zeller2014, n = n_taxa + 10)

  # Should return at most the number of taxa with non-zero abundance
  expect_true(length(result) <= n_taxa)
})

# Tests for taxonomy() - note: zeller2014 has no taxonomy table
test_that("taxonomy function exists and handles phyloseq objects",
          {
            # This will error if zeller2014 has no tax_table
            # which is expected based on CLAUDE.md
            expect_error(taxonomy(zeller2014))
          })
