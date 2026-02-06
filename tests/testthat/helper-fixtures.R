# Test fixtures for metagenomeR tests
# This file is automatically loaded before tests run

# Load required packages for testing
# These are needed because the package code uses unqualified function calls
suppressPackageStartupMessages({
  library(phyloseq)
  library(dplyr)
  library(vegan)
  library(ggplot2)
  library(data.table)
  library(tibble)
  if (requireNamespace("ggpubr", quietly = TRUE)) {
    library(ggpubr)
  }
})

# Load the zeller_2014 sample data
data("zeller_2014", package = "metagenomeR")

# Create a simple mock functional profile for testing
create_mock_func_profile <- function() {
  # Create sample data
  n_samples <- 10
  n_features <- 50
  sample_names <- paste0("Sample", 1:n_samples)

  # Create UNIREF matrix (samples as rows, features as columns)
  uniref_mat <- matrix(
    runif(n_samples * n_features, 0, 100),
    nrow = n_samples,
    ncol = n_features,
    dimnames = list(
      sample_names,
      paste0("UniRef90_", 1:n_features)
    )
  )
  uniref_df <- as.data.frame(uniref_mat)

  # Create EC matrix
  n_ec <- 20
  ec_mat <- matrix(
    runif(n_samples * n_ec, 0, 50),
    nrow = n_samples,
    ncol = n_ec,
    dimnames = list(
      sample_names,
      paste0("EC", 1:n_ec)
    )
  )
  ec_df <- as.data.frame(ec_mat)

  # Create KO matrix
  n_ko <- 30
  ko_mat <- matrix(
    runif(n_samples * n_ko, 0, 75),
    nrow = n_samples,
    ncol = n_ko,
    dimnames = list(
      sample_names,
      paste0("K", sprintf("%05d", 1:n_ko))
    )
  )
  ko_df <- as.data.frame(ko_mat)

  # Create metadata
  metadata <- data.frame(
    row.names = sample_names,
    Condition = rep(c("Control", "Treatment"), each = n_samples / 2),
    Age = sample(25:65, n_samples, replace = TRUE),
    BMI = runif(n_samples, 18.5, 35)
  )

  list(
    UNIREF = uniref_df,
    EC = ec_df,
    KO = ko_df,
    Metadata = metadata
  )
}

# Create a simple test data frame for basic tests
create_test_df <- function() {
  data.frame(
    group = rep(c("A", "B"), each = 5),
    value = c(1:5, 6:10),
    row.names = paste0("sample", 1:10)
  )
}

# Create unequal groups data frame for testing
create_unequal_groups_df <- function() {
  data.frame(
    group = c(rep("A", 3), rep("B", 5)),
    value = 1:8,
    row.names = paste0("sample", 1:8)
  )
}
