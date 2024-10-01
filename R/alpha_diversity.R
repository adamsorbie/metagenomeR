# calculate shannon effective
#' @keywords internal
Shannon.E <- function(x) {
  summed <- sum(x)
  shannon.e <-
    round(exp(-sum(x[x > 0] / summed * log(x[x > 0] / summed))), digits = 2)
  return(shannon.e)
}

# calculate richness
#' @keywords internal
Richness <- function(x, detection = 1e-5) {
  observed <- sum(x > detection)
  return(observed)
}



#' Calculate alpha-diversity indices from a phyloseq object
#'
#' @param ps phyloseq object
#' @param indices Alpha-diversity indices to calculate. Supported metrics are `"Richness"`,`"Shannon.Effective"`, `"Shannon"` and `"Inverse.Simpson"`
#' @return Returns the specified alpha diversity metrics for each sample as a dataframe
#'
#' @examples
#'
#' data(zeller_2014)
#' alpha_div <- alpha_div(zeller_2014)
#' Only richness
#' alpha_div <- alpha_div(zeller_2014, indices=c("Richness"))
#' @export
calc_alpha <- function(ps, indices=c("Richness",
                                     "Shannon.Effective",
                                     "Shannon",
                                     "Inverse.Simpson")) {
  mat_in <- ps_to_feattab(ps) %>%
    t()
  if (length(indices < 4)){
    if (indices %in% c("Richness", "Shannon.Effective", "Shannon", "Inverse.Simpson")){

    }
  } else {

  }
  diversity <-
    setNames(data.frame(matrix(ncol = length(indices), nrow = nsamples(ps))),
             indices)
  rownames(diversity) <- rownames(meta_to_df(ps))

  diversity$Richness <- apply(mat_in, 1, Richness)
  diversity$Shannon.Effective <- apply(mat_in, 1, Shannon.E)
  diversity$Shannon <- apply(mat_in, 1, function(x) diversity(x, index="shannon"))
  diversity$InvSimpson <- apply(mat_in, 1, function(x) diversity(x, index = "invsimpson"))

  return(diversity)
}


