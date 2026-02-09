#' Transform taxonomic abundance data
#'
#' Apply transformations to taxonomic abundance data in a phyloseq object.
#'
#' @param ps A phyloseq object containing taxonomic abundance data.
#' @param transform Character string specifying the transformation to apply.
#'   One of "relative", "arcsin", or "log10". Default is "relative".
#' @param offset Numeric offset to add before log10 transformation to avoid
#'   log of zero. Default is 1.
#'
#' @return A phyloseq object with transformed abundance data.
#'
#' @export
transform_tax <- function(ps,
                          transform = "relative",
                          offset = 1) {
  if (length(is(ps)) == 1 && class(ps) == "phyloseq") {
    x <- ps_to_feattab(ps)
    mean_sum <- mean(colSums(x))
  }
  else {
    print("not a phyloseq object, exiting")
    stop()
  }

  if (transform %in% c("relative", "arcsin", "log10")) {
    if (transform == "relative") {
      ps_t <- t(t(x) / colSums(x))

    } else if (transform == "arcsin") {
      if (mean_sum == 1) {
        ps_t <- asin(sqrt(x))
      } else if (mean_sum == 100) {
        x <- x / 100
        ps_t <- asin(sqrt(x))
      } else {
        stop("not proportion data, exiting")
      }
    } else if (transform == "log10") {
      if (max_x == 100 | max_x == 1) {
        warning("data are relative abundances", call. = F)
        ps_t <- log10(x + offset)
      }   else {
        ps_t <- log10(x + offset)
      }
    }
    otu_table(ps)@.Data <- as.matrix(ps_t)

    return(ps)

  } else {
    print("Not a valid transform, exiting")
    stop()
  }

}

filter_ps <- function(ps, abun, prev = NULL) {
  x <- ps_to_feattab(ps) %>%
    as.matrix()


  if (is.null(prev)) {
    ps_filt <- x[rowSums(x)  >= abun, ]
  } else if (!is.null(prev)) {
    ps_filt <- x[rowSums(x >= abun) >= ncol(x) * prev, ]
  }

  otu_table(ps)@.Data <- ps_filt

  return(ps)
}

subset_samples_func <- function(func_profile, filter_statement) {
  keep <- func_profile$Metadata %>%
    dplyr::filter(eval(rlang::parse_expr(filter_statement))) %>%
    rownames()

  func_profile_filt <- map(func_profile, function(x)
    filter_rownames(x, keep))

  return(func_profile_filt)

}


prune_samples_func  <- function(func_profile, keep) {
  func_profile_filt <- map(func_profile, function(x)
    filter_rownames(x, keep))
  return(func_profile_filt)

}


asin_sqrt <- function(x) {
  mean_sum <- mean(colSums(x))
  if (mean_sum == 1) {
    ps_t <- asin(sqrt(x))
  } else if (mean_sum == 100) {
    x <- x / 100
    ps_t <- asin(sqrt(x))
  } else {
    stop("not proportion data, exiting")
  }
}

transform_func <- function(func_profile,
                           transform,
                           features_are_rows = TRUE,
                           offset = 1e-6) {
  if ("Metadata" %in% names(func_profile)) {
    feat_tables <- get_feat_tables(func_profile)
  } else {
    feat_tables <- func_profile
  }

  if (features_are_rows == FALSE) {
    feat_tables <- map(feat_tables, function(x)
      t_df(x))
  }
  if (transform %in% c("relative", "arcsin", "log10")) {
    if (transform == "relative") {
      feat_tables <- map(feat_tables, function(x)
        t_df(100 * t_df(x) / colSums(x)))
      feat_tables$Metadata <- func_profile$Metadata
    }
    else if (transform == "arcsin") {
      feat_tables <- map(feat_tables, function(x)
        asin_sqrt(x))
      feat_tables$Metadata <- func_profile$Metadata
    } else if (transform == "log10") {
      feat_tables <- map(feat_tables, function(x)
        log10(x + offset))
      feat_tables$Metadata <- func_profile$Metadata
    }
  } else {
    stop("Not a valid transform")
  }



  if (features_are_rows == TRUE) {
    return(feat_tables)
  } else {
    feat_tables <- map(feat_tables, function(x)
      t_df(x))
    feat_tables$Metadata <- func_profile$Metadata
    return(feat_tables)
  }
}

filt_prev_abund  <- function(x, abund, prev, orientation = "cols") {
  if (orientation %in% c("cols", "columns", "c", "r", "rows")) {
    if (orientation %in% c("cols", "columns", "c")) {
      filt <- x[which(rowSums(x > abund)
                      >= prev * ncol(x)), ]
      return(filt)
    } else {
      filt <- x[, which(colSums(x > abund)
                        >= prev * nrow(x))]
      return(filt)
    }

  }

  return(filt)
}

filter_func <- function(func_profile, abund, prev, renorm = T) {
  feat_tables <- get_feat_tables(func_profile)
  func_profile_filt <- map(feat_tables, function(x)
    filt_prev_abund(x, abund, prev, orientation = "rows"))

  func_profile_filt$Metadata <- func_profile$Metadata

  if (renorm == T) {
    func_profile_filt <- transform_func(func_profile_filt, "relative", features_are_rows = F)
  }
  return(func_profile_filt)
}
