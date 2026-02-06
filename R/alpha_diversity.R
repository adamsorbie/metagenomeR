# Calculate shannon effective
#' @keywords internal
shannon_e <- function(x) {
  summed <- sum(x)
  sample_shannon <-
    round(exp(-sum(x[x > 0] / summed * log(x[x > 0] / summed))), digits = 2)
  return(sample_shannon)
}

# Calculate richness
#' @keywords internal
richness <- function(x, detection = 1e-5) {
  sample_richness <- sum(x > detection)
  return(sample_richness)
}



#' Calculate alpha-diversity indices from a phyloseq object
#'
#' @param ps phyloseq object
#' @param indices Alpha-diversity indices to calculate. Supported metrics are: \cr
#' `"Richness"`,`"Shannon.Effective"`, `"Shannon"` and `"Inverse.Simpson"`
#' @return Returns specified alpha-diversity metrics as a dataframe
#'
#' @examples
#' data(zeller_2014)
#' alpha_div <- calc_alpha(zeller_2014)
#' # Only richness
#' alpha_div <- calc_alpha(zeller_2014, indices = c("Richness"))
#'
#' @export
calc_alpha <- function(ps,
                       indices = c("Richness", "Shannon.Effective", "Shannon", "Inverse.Simpson")) {
  mat_in <- ps_to_feattab(ps) %>%
    t()
  if (length(indices < 4)) {
    if (indices %in% c("Richness",
                       "Shannon.Effective",
                       "Shannon",
                       "Inverse.Simpson")) {

    }
  } else {

  }
  diversity <-
    setNames(data.frame(matrix(
      ncol = length(indices), nrow = nsamples(ps)
    )), indices)
  rownames(diversity) <- rownames(meta_to_df(ps))

  diversity$Richness <- apply(mat_in, 1, Richness)
  diversity$Shannon.Effective <- apply(mat_in, 1, Shannon.E)
  diversity$Shannon <- apply(mat_in, 1, function(x)
    diversity(x, index = "shannon"))
  diversity$InvSimpson <- apply(mat_in, 1, function(x)
    diversity(x, index = "invsimpson"))

  return(diversity)
}
