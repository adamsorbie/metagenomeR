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
#' @return Returns specified alpha-diversity metrics as a dataframe
#'
#' @examples
#' data(zeller_2014)
#' alpha_div <- calc_alpha(zeller_2014)
#' # Only richness
#' alpha_div <- calc_alpha(zeller_2014, indices = c("Richness"))
#'
#' @export
calc_alpha <- function(ps, ...) {
  mat_in <- ps_to_feattab(ps) %>%
    t()

  diversity <-
    setNames(data.frame(matrix(ncol = 4, nrow = nsamples(ps))),
             c("Richness", "Shannon.Effective", "Shannon", "Inverse.Simpson"))
  rownames(diversity) <- rownames(meta_to_df(ps))

  diversity$Richness <- apply(mat_in, 1, richness, ...)
  diversity$Shannon.Effective <- apply(mat_in, 1, shannon_e)
  diversity$Shannon <- apply(mat_in, 1, function(x) diversity(x))
  diversity$Inverse.Simpson <- apply(mat_in, 1, function(x) diversity(x, index = "invsimpson"))

  diversity_out <- merge(diversity, meta_to_df(ps), by=0)
  return(diversity_out)
}
