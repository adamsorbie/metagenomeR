#' Perform MaAsLin2 Analysis for Taxonomic Profile
#'
#' This function performs a MaAsLin2 analysis on a taxonomic profile dataset, allowing for the inclusion
#' of fixed effects and the specification of abundance and prevalence thresholds.
#' It supports different taxonomic classifications such as Metaphlan or GTDB.
#'
#' @param tax_profile A phyloseq object or similar containing the taxonomic profiles and metadata.
#' @param out A character string specifying the output directory for MaAsLin2 results.
#' @param fixed A character vector of fixed effects (i.e., covariates) to include in the model.
#' @param abun_thresh A numeric value specifying the minimum abundance threshold for filtering taxa. Default is `0`.
#' @param prev_thresh A numeric value specifying the minimum prevalence threshold for filtering taxa. Default is `0`.
#' @param taxonomy A character string specifying the taxonomy database used (e.g., `"metaphlan"` or `"gtdb"`). Default is `"metaphlan"`.
#' @param ... Additional arguments passed to the `Maaslin2` function.
#'
#' @return The output is saved to the directory specified by the `out` parameter, containing MaAsLin2 results.
#'
#' @examples
#' # Example usage:
#' maaslin2_tax(tax_profile, out = "maaslin2_output", fixed = c("age", "sex"))
#'
#' @export
maaslin2_tax <- function(tax_profile, out,
                         fixed, abun_thresh=0, prev_thresh=0, taxonomy="metaphlan",...){

  mat_in <- ps_to_feattab(tax_profile) %>%
    t_df()
  # this only works for species level
  if (taxonomy == "metaphlan"){
    colnames(mat_in) <- sapply(strsplit(colnames(mat_in), split= "\\|", fixed = TRUE), tail, 1L)
  } else if (taxonomy == "gtdb"){
    colnames(mat_in) <- sapply(strsplit(colnames(mat_in), split= ";", fixed = TRUE), tail, 1L)
  } else {
    stop("unsupported taxonomy")
  }


  metadata_in <- meta_to_df(tax_profile)

  Maaslin2(input_data = mat_in,
           input_metadata = metadata_in,
           output = out,
           normalization = "NONE",
           transform = "NONE",
           min_abundance = abun_thresh,
           min_prevalence = prev_thresh,
           fixed_effects = fixed,
           cores=6,
           ...
  )
}

#' Perform MaAsLin2 Analysis for Functional Profiles
#'
#' This function performs a MaAsLin2 analysis on functional profiles such as UniRef, EC, KO, or pathways.
#' The function supports multiple feature types and includes filtering options for abundance and prevalence thresholds.
#'
#' @param func_profiles A list containing functional profiles and corresponding metadata.
#' @param feattype A character string specifying the type of feature (e.g., `"UNIREF"`, `"EC"`, `"KO"`, `"pathways"`).
#' @param out A character string specifying the output directory for MaAsLin2 results.
#' @param fixed A character vector of fixed effects (i.e., covariates) to include in the model.
#' @param abun_thresh A numeric value specifying the minimum abundance threshold for filtering features. Default is `0`.
#' @param prev_thresh A numeric value specifying the minimum prevalence threshold for filtering features. Default is `0`.
#' @param cores An integer specifying the number of cores to use for parallel computation. Default is `6`.
#' @param ... Additional arguments passed to the `Maaslin2` function.
#'
#' @return The output is saved to the directory specified by the `out` parameter, containing MaAsLin2 results.
#'
#' @examples
#' # Example usage:
#' maaslin2_func(func_profiles, feattype = "UNIREF", out = "maaslin2_output", fixed = c("age", "sex"))
#'
#' @export
maaslin2_func <- function(func_profiles, feattype, out,
                          fixed, abun_thresh=0, prev_thresh=0, cores=6, ...){
  # write check for filtering and AST transform
  if (feattype %in% c("UNIREF", "EC", "KO", "pathways")){
    mat_in <- func_profiles[[feattype]]
    metadata_in <- func_profiles$Metadata
  } else {
    stop("Feature type not supported, choose from: 'UNIREF', 'EC', 'KO', 'pathways'")
  }

  Maaslin2(input_data = mat_in,
           input_metadata = metadata_in,
           output = out,
           normalization = "NONE",
           transform = "NONE",
           min_abundance = abun_thresh,
           min_prevalence = prev_thresh,
           fixed_effects = fixed,
           cores=cores,
           ...
  )
}


