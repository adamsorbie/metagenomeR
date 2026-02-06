#' Perform PERMANOVA (Adonis) Analysis on a Phyloseq Object
#'
#' This function performs a PERMANOVA (Adonis) test on a phyloseq object using a distance matrix and a grouping variable.
#' It returns the R-squared and p-values for the relationship between the distance matrix and the grouping variable.
#'
#' @param ps A phyloseq object containing both the feature table and metadata.
#' @param dist_matrix A distance matrix calculated using a method like Bray-Curtis or Jaccard.
#' @param group_variable A character string specifying the metadata column to be used as the grouping variable.
#' @param ... Additional arguments to be passed to the `adonis2` function.
#'
#' @return A result object from the `adonis2` function containing the R-squared and p-values.
#'
#' @note If the metadata contains missing values (NAs) in the `group_variable`, the function will print an error message and return `NULL`.
#'
#' @examples
#' # Example usage:
#' dist_matrix <- phyloseq::distance(ps, method = "bray")
#' adonis_result <- phyloseq_adonis(ps, dist_matrix, group_variable = "Treatment")
#'
#' @export
phyloseq_adonis <- function(ps, dist_matrix, group_variable, ...) {
  meta_df <- meta_to_df(ps)
  if (sum(is.na(meta_df[[group_variable]])) > 0) {
    print("metadata contains NAs, remove these samples with subset_samples
          before continuing")
    return (NULL)
  } else {
    # convert distance matix object to string
    dist_str <- deparse(substitute(dist_matrix))
    # define formula
    form <- as.formula(paste(dist_str, group_variable, sep="~"))
    ps_ad <- adonis2(form, data = meta_df, ...)
    return(ps_ad)
  }
}

#' Perform Beta Dispersion Analysis on a Phyloseq Object
#'
#' This function calculates the homogeneity of group dispersions (beta diversity) using a distance matrix and a grouping variable.
#' It returns the ANOVA results of the dispersion analysis.
#'
#' @param ps A phyloseq object containing both the feature table and metadata.
#' @param dist_matrix A distance matrix calculated using methods like Bray-Curtis or Jaccard.
#' @param group_variable A character string specifying the metadata column to be used as the grouping variable.
#' @param ... Additional arguments to be passed to the `betadisper` function.
#'
#' @return An ANOVA result object from the `betadisper` function showing the dispersion among groups.
#'
#' @note If the metadata contains missing values (NAs) in the `group_variable`, the function will print an error message and return `NULL`.
#'
#' @examples
#' # Example usage:
#' dist_matrix <- phyloseq::distance(ps, method = "bray")
#' betadisper_result <- phyloseq_betadisper(ps, dist_matrix, group_variable = "Treatment")
#'
#' @export
phyloseq_betadisper <- function(ps, dist_matrix, group_variable, ...) {
  meta_df <- meta_to_df(ps)
  if (sum(is.na(meta_df[[group_variable]])) > 0) {
    print("metadata contains NAs, remove these samples with subset_samples
          before continuing")
    return (NULL)
  } else {

    bd <- betadisper(dist_matrix, meta_df[[group_variable]])
    anova_res <- anova(bd)
    return(anova_res)
  }
}

#' Calculate Beta Diversity and Ordination on a Phyloseq Object
#'
#' This function calculates a beta diversity distance matrix and ordination using methods such as NMDS, MDS, or PCoA.
#' It supports distance metrics such as Bray-Curtis and Jaccard.
#'
#' @param ps A phyloseq object containing both the feature table and metadata.
#' @param dist A character string specifying the distance metric to use. Supported metrics are `"bray"` and `"jaccard"`.
#' @param ord_method A character string specifying the ordination method to use. Supported methods are `"NMDS"`, `"MDS"`, and `"PCoA"`. Default is `"NMDS"`.
#'
#' @return A list containing the distance matrix and ordination object.
#'
#' @note If an unsupported ordination method or distance metric is provided, the function will print an error message and stop execution.
#'
#' @examples
#' # Example usage:
#' beta_div_result <- calc_betadiv(ps, dist = "bray", ord_method = "NMDS")
#'
#' @export
calc_betadiv <- function(ps, dist, ord_method = "NMDS") {
  if (ord_method %in% c("NMDS", "MDS", "PCoA")) {

    if (dist %in% c("bray", "jaccard")) {
      dist_mat <- distance(ps, dist)
      if (ord_method == "NMDS"){
        ord <- ordinate(ps, ord_method, dist_mat, trace=FALSE)
      } else {
        ord <- ordinate(ps, ord_method, dist_mat)
      }

      return_list <- list("Distance_Matrix" = dist_mat,
                          "Ordination" = ord)
      return(return_list)
    }
    else {
      print(
        "Distance metric not supported, supported metrics are; bray, jaccard"
      )
    }
  }
  else {
    print("Ordination method not supported, supported methods are: NMDS, MDS, PCoA")
    stop()
  }

}
