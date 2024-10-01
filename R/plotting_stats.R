# XY formula
#' @param ps phyloseq object
#' @return Returns the taxonomy table slot as a dataframe
#' @export
#'
#' @examples
#' data(dietswap)
#' taxa <- taxonomy(dietswap)
#' @keywords internal
xyform <- function (y_var, x_vars) {
  # y_var: a length-one character vector
  # x_vars: a character vector of object names
  as.formula(sprintf("%s ~ %s", y_var, paste(x_vars, collapse = " + ")))
}

# check sample sizes for paired test
#' @param ps phyloseq object
#' @return Returns the taxonomy table slot as a dataframe
#' @export
#'
#' @examples
#' data(dietswap)
#' taxa <- taxonomy(dietswap)
#' @keywords internal
check_equal_sample_sizes <- function(df, variable_col) {
  # Count the number of occurrences for each group
  sample_sizes <- df %>%
    group_by_at(variable_col) %>%
    summarise(count = n()) %>%
    pull(count)

  # Check if all counts are equal
  if (length(unique(sample_sizes)) == 1) {
    return(TRUE)  # All groups have equal sample sizes
  } else {
    return(FALSE) # Sample sizes are not equal across groups
  }
}

#
stat_plot <- function(df, formula, variable_col, assume_normality=F, multiple_groups=F, paired=F) {

  # Check paired
  if (paired == TRUE){
    if (check_equal_sample_sizes(df, variable_col) == F){
      stop("Sample sizes are not equal for paired statistical test")
    }
  }

  if (multiple_groups == TRUE) {
    if (paired == TRUE) {

      stat_variance <- df %>%
        friedman_test(formula)
      stat_test <- df %>%
        pairwise_wilcox_test(
          formula,
          comparisons = comparisons_list,
          p.adjust.method = "BH",
          paired = TRUE
        ) %>%
        add_significance() %>%
        add_xy_position(x = variable_col) %>%
        filter(p.adj < 0.05)
    }
    else {
      stat_variance <- df %>%
        kruskal_test(formula)
      stat_test <- df %>%
        pairwise_wilcox_test(formula,
                             comparisons = comparisons_list,
                             p.adjust.method = "BH") %>%
        add_significance() %>%
        add_xy_position(x = variable_col) %>%
        filter(p.adj < 0.05)
    }
  }
  else if (multiple_groups == FALSE) {
    if (paired == TRUE) {
      stat_test <- df %>%
        wilcox_test(formula, paired = TRUE) %>%
        add_significance() %>%
        add_xy_position(x = variable_col) %>%
        filter(p < 0.05)
    }
    else {
      stat_test <- df %>%
        wilcox_test(formula) %>%
        add_significance() %>%
        add_xy_position(x = variable_col) %>%
        filter(p < 0.05)
    }

  }
}

# plot boxplot with stats
plot_boxplot <- function(df,
                         variable_col,
                         value_col,
                         comparisons_list,
                         fill_var = variable_col,
                         xlab = variable_col,
                         ylab = value_col,
                         p_title = NULL,
                         multiple_groups = FALSE,
                         cols = NULL,
                         group.order = NULL,
                         paired = FALSE,
                         normal=FALSE,
                         ...) {
  # extend color palette with transparent value - required due to way we are
  # layering plot
  if (is.null(cols)) {
    cols <- pal_npg()(length(unique(df[, variable_col])))
  }
  cols <- c(cols, "transparent")

  if (!is.null(group.order)) {
    df[, variable_col] <-
      factor(df[, variable_col], levels = group.order)
  }

  formula <- xyform(value_col, variable_col)

  if (multiple_groups == T) {
    stat_variance <- stat_plot(
      df,
      formula,
      variable_col,
      assume_normality = normal,
      multiple_groups = multiple_groups,
      paired = paired
    )
  }
  stat_test <- stat_plot(
    df,
    formula,
    variable_col,
    assume_normality = normal,
    multiple_groups = multiple_groups,
    paired = paired
  )

  # aes string accepts strings as column names, this code plots boxplot and adds error bars
  plot <- ggplot(
    df,
    aes_string(
      x = variable_col,
      y = value_col,
      fill = variable_col,
      color = variable_col
    )
  ) +
    geom_boxplot(
      color = "black",
      alpha = 0.8,
      outlier.shape = 5,
      outlier.size = 1
    ) +
    geom_point(size = 1.5, position = position_jitterdodge()) +
    labs(x = xlab, y = ylab) +
    stat_boxplot(color = "black",
                 geom = "errorbar",
                 width = 0.2)
  # creates new 'finalised plot' and adds statistical significance, labels and adjusts theme and title
  final_plot <- plot +
    theme_classic2() +
    ggtitle(p_title) +
    theme(
      axis.text.x = element_text(size = 14),
      axis.text.y = element_text(size = 14),
      axis.title.x = element_text(size = 18),
      axis.title.y = element_text(size = 18),
      axis.line = element_line(colour = "black"),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_blank(),
      panel.border = element_blank(),
      plot.title = element_text(hjust = 0.5),
      legend.position = "None"
    ) +
    scale_fill_manual(values = cols) +
    scale_color_manual(values = cols) +
    rotate_x_text(angle = 45) +
    lemon::coord_capped_cart()

  if (dim(stat_test)[1] == 0) {
    plot_out <- final_plot
  }
  else {
    if (multiple_groups == T) {
      label <- "p.adj.signif"
    } else {
      label <- "p.signif"
    }
    plot_out <- final_plot +
      stat_pvalue_manual(
        stat_test,
        label = label,
        hide.ns = T,
        inherit.aes = FALSE,
        ...
      )
  }

  return(plot_out)
}



# plot scatter plot with correlation if desired
plot_scatter <- function(df,
                         x,
                         y,
                         point_color,
                         line_color,
                         fill_color,
                         xlab,
                         ylab,
                         corr.method = NULL,
                         ...) {
  p <-
    ggplot(data = df, mapping = aes(x = .data[[x]], y = .data[[y]])) +
    geom_point(aes(color = point_color), size = 2.5) +
    geom_smooth(method = "lm",
                color = line_color,
                fill = fill_color) +
    theme_bw() +
    theme(
      legend.position = "None",
      axis.title.x = element_text(size = 14),
      axis.text.x = element_text(size = 12),
      axis.title.y = element_text(size = 14),
      axis.text.y = element_text(size = 12)
    ) +
    xlab(xlab) +
    ylab(ylab)

  if (!is.null(corr.method)) {
    p <- p + stat_cor(method = corr.method, ...)
    return(p)
  }
  else {
    return(p)
  }
}



plot_pca <- function(func_profile,what, var, colours, add_ellipse=F) {

  # calc PCA
  metadata_filt <- func_profile[["Metadata"]] %>% select(all_of(var))
  dat <- merge(func_profile[[what]], metadata_filt, by=0) %>%
    column_to_rownames("Row.names")
  # this will select numerical metadata
  pca <- prcomp(dat %>% select(where(is.numeric)), scale. = T)


  # Plot
  scree_plot <- fviz_eig(pca)
  pca_plot <- autoplot(pca, data = dat,
                       colour = var) +
    geom_point(aes(color= .data[[ var ]] ), size=3, alpha=0.75) +
    scale_color_manual(values = colours) +
    scale_fill_manual(values = colours) +
    theme_cowplot()
  pca_plot$layers[[1]] <- NULL

  if (add_ellipse == TRUE){
    pca_plot <- pca_plot +
      geom_polygon(stat = "ellipse", aes(fill = .data [[ var ]] ), alpha = 0.3)
  }

  return(pca_plot)
}

plot_pcoa <- function(func_profile, what, var, colours, add_ellipse=F, size=3, plot="MDS"){


  if (plot == "MDS"){
    dist <- calc_dist(func_profile, feat_type = what)
    mds <- cmdscale(dist, eig = T)
    x_label <- "PCo1"
    y_label <- "PCo2"
    colnames(mds$points) <- c(x_label, y_label)


  } else if (plot == "NMDS"){
    # NMDS
    mds <- metaMDS(comm = func_profile[[what]], distance = "bray", trace = FALSE, autotransform = FALSE)
    x_label <- "NMDS1"
    y_label <- "NMDS2"
    colnames(mds$points) <- c("PCo1", "PCo2")
  }


  mds_dat <- mds$points %>%
    merge(func_profile[["Metadata"]], by=0)
  expvar <- (eigenvals(mds)/sum(eigenvals(mds)))[1:2]

  mds_plot <- ggplot(mds_dat, aes(PCo1, PCo2, color = .data[[ var ]])) +
    geom_point(size=size, alpha=0.75) +
    scale_color_manual(values = colours) +
    scale_fill_manual(values = colours) +
    theme_cowplot() +
    xlab(paste(x_label, round(expvar[1]*100, digits = 1), "%")) +
    ylab(paste(y_label, round(expvar[2]*100, digits = 1), "%" ))


  if (add_ellipse == TRUE){
    mds_plot <- mds_plot +
      geom_polygon(stat = "ellipse", aes(fill = .data [[ var ]] ), alpha = 0.3)
  }

  return(mds_plot)

}


plot_volcano <- function(results_df, pthresh=0.05, FCthresh=0.58, ...) {
  res <- results_df$res
  # replace rownames with highest classified level
  rownames(res) <- as.character(lapply(strsplit(as.character(rownames(res)), split="\\|"),
                                       tail, n=1))
  groups <- results_df$groups
  EnhancedVolcano(res, x="log2FC", y="p.adj",
                  lab=rownames(res),
                  labSize = 4,
                  pCutoff = pthresh,
                  FCcutoff = FCthresh,
                  title = paste(groups[1],
                                "versus",
                                groups[2],
                                sep=" "),
                  subtitle = NULL,
                  max.overlaps = 30,
                  ...)
}
