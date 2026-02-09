#' Create a Formula for Linear Models
#'
#' This function generates a formula to be used in linear models or statistical tests. It constructs the formula in the format `y ~ x1 + x2 + ...`.
#'
#' @param y_var A character string representing the dependent variable.
#' @param x_vars A character vector of independent variables.
#'
#' @return A formula object of the form `y_var ~ x_vars`.
#'
#' @keywords internal
xyform <- function (y_var, x_vars) {
  # y_var: a length-one character vector
  # x_vars: a character vector of object names
  as.formula(sprintf("%s ~ %s", y_var, paste(x_vars, collapse = " + ")))
}

#' Check if Sample Sizes are Equal Across Groups
#'
#' This function checks whether the sample sizes for each group in a specified column of a data frame are equal.
#'
#' @param df A data frame containing the data.
#' @param variable_col A character string representing the column with groupings to check for equal sample sizes.
#'
#' @return A logical value: `TRUE` if all groups have equal sample sizes, and `FALSE` otherwise.
#'
#' @examples
#' # Example usage:
#' equal_sizes <- check_equal_sample_sizes(df, "group_column")
#'
#' @export
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

#' Perform Statistical Tests and create dataframe to add to plot results
#'
#' This function performs statistical tests (such as Wilcoxon, Kruskal-Wallis, or Friedman) on a data frame based on the provided formula. It also checks for paired tests and multiple groups.
#'
#' @param df A data frame containing the data.
#' @param formula A formula specifying the relationship between dependent and independent variables.
#' @param variable_col A character string representing the column with groupings.
#' @param assume_normality A logical value indicating whether to assume normality for parametric tests. Default is `FALSE`.
#' @param multiple_groups A logical value indicating whether multiple groups are present. Default is `FALSE`.
#' @param paired A logical value indicating whether the test is paired. Default is `FALSE`.
#'
#' @return A data frame with the statistical test results, including p-values and significance.
#'
#' @note The function will stop if the sample sizes for paired tests are not equal.
#'
#' @examples
#' # Example usage:
#' stat_results <- stat_plot(df, formula = y ~ x, variable_col = "group")
#'
#' @importFrom rstatix wilcox_test kruskal_test friedman_test pairwise_wilcox_test add_significance add_xy_position
#' @export
stat_plot <- function(df,
                      formula,
                      variable_col,
                      assume_normality = F,
                      multiple_groups = F,
                      paired = F) {
  # Check paired
  if (paired == TRUE) {
    if (check_equal_sample_sizes(df, variable_col) == F) {
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
        dplyr::filter(p.adj < 0.05)
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
        dplyr::filter(p.adj < 0.05)
    }
  }
  else if (multiple_groups == FALSE) {
    if (paired == TRUE) {
      stat_test <- df %>%
        wilcox_test(formula, paired = TRUE) %>%
        add_significance() %>%
        add_xy_position(x = variable_col) %>%
        dplyr::filter(p < 0.05)
    }
    else {
      stat_test <- df %>%
        wilcox_test(formula) %>%
        add_significance() %>%
        add_xy_position(x = variable_col) %>%
        dplyr::filter(p < 0.05)
    }

  }
}

#' Plot Boxplot with Statistical Significance
#'
#' This function generates a boxplot for a given data frame and performs statistical tests to visualize significance. It supports plotting for paired or unpaired tests with options for multiple groups.
#'
#' @param df A data frame containing the data to be plotted.
#' @param variable_col A character string representing the independent variable column in the data frame.
#' @param value_col A character string representing the dependent variable column in the data frame.
#' @param comparisons_list A list of comparisons to be made for significance testing.
#' @param fill_var A character string specifying the variable to use for filling the boxplot. Default is the same as `variable_col`.
#' @param xlab A character string for the x-axis label. Default is the same as `variable_col`.
#' @param ylab A character string for the y-axis label. Default is the same as `value_col`.
#' @param p_title A character string for the plot title. Default is `NULL`.
#' @param multiple_groups A logical value indicating if the plot contains multiple groups. Default is `FALSE`.
#' @param cols A vector of colors for the boxplot. Default is `NULL`.
#' @param group.order A character vector specifying the order of the groups. Default is `NULL`.
#' @param paired A logical value indicating if the test is paired. Default is `FALSE`.
#' @param normal A logical value indicating if normality should be assumed. Default is `FALSE`.
#' @param ... Additional arguments passed to the plotting functions.
#'
#' @return A ggplot object representing the boxplot with optional statistical significance annotations.
#'
#' @examples
#' # Example usage:
#' plot <- plot_boxplot(df, "group", "value", comparisons_list = list(c("A", "B")))
#'
#' @export
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
                         normal = FALSE,
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
    # could maybe just add if clause here and convert this function into general purpose plotting
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


#' Plot Scatter Plot with Optional Correlation
#'
#' This function generates a scatter plot with an optional linear regression line and calculates correlation statistics if specified.
#'
#' @param df A data frame containing the data to be plotted.
#' @param x A character string representing the column to be used for the x-axis.
#' @param y A character string representing the column to be used for the y-axis.
#' @param point_color A character string specifying the color of the points.
#' @param line_color A character string specifying the color of the regression line.
#' @param fill_color A character string specifying the fill color of the regression line.
#' @param xlab A character string for the x-axis label.
#' @param ylab A character string for the y-axis label.
#' @param corr.method A character string specifying the correlation method (e.g., `"pearson"`, `"spearman"`, `"kendall"`). Default is `NULL`.
#' @param ... Additional arguments passed to the `stat_cor` function for correlation statistics.
#'
#' @return A ggplot object representing the scatter plot, with optional correlation results.
#'
#' @examples
#' # Example usage:
#' plot <- plot_scatter(df, x = "height", y = "weight", point_color = "blue", line_color = "red", fill_color = "lightblue", xlab = "Height", ylab = "Weight", corr.method = "pearson")
#'
#' @export
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


#' Plot PCA (Principal Component Analysis) for Functional Profile
#'
#' This function calculates and plots a PCA from a functional profile with metadata, optionally adding ellipses for grouping.
#'
#' @param func_profile A list containing functional profiles and corresponding metadata.
#' @param what A character string specifying which part of the functional profile to use for PCA.
#' @param var A character string specifying the grouping variable in the metadata.
#' @param colours A vector of colors to use for the groups.
#' @param add_ellipse A logical value indicating whether to add ellipses around the groups. Default is `FALSE`.
#'
#' @return A ggplot object representing the PCA plot.
#'
#' @examples
#' # Example usage:
#' pca_plot <- plot_pca(func_profile, what = "UNIREF", var = "Condition", colours = c("red", "blue"), add_ellipse = TRUE)
#'
#' @export
plot_pca <- function(func_profile,
                     what,
                     var,
                     colours,
                     add_ellipse = F) {
  # calc PCA
  metadata_filt <- func_profile[["Metadata"]] %>% dplyr::select(all_of(var))
  dat <- merge(func_profile[[what]], metadata_filt, by = 0) %>%
    column_to_rownames("Row.names")
  # this will select numerical metadata
  pca <- prcomp(dat %>% dplyr::select(where(is.numeric)), scale. = T)


  # Plot
  scree_plot <- fviz_eig(pca)
  pca_plot <- autoplot(pca, data = dat, colour = var) +
    geom_point(aes(color = .data[[var]]), size = 3, alpha = 0.75) +
    scale_color_manual(values = colours) +
    scale_fill_manual(values = colours) +
    theme_cowplot()
  pca_plot$layers[[1]] <- NULL

  if (add_ellipse == TRUE) {
    pca_plot <- pca_plot +
      geom_polygon(stat = "ellipse", aes(fill = .data [[var]]), alpha = 0.3)
  }

  return(pca_plot)
}

#' Plot PCoA (Principal Coordinates Analysis) or NMDS (Non-metric Multidimensional Scaling)
#'
#' This function performs and plots PCoA or NMDS from a functional profile with metadata. It supports adding ellipses and specifying plot dimensions.
#'
#' @param func_profile A list containing functional profiles and corresponding metadata.
#' @param what A character string specifying which part of the functional profile to use for PCoA or NMDS.
#' @param var A character string specifying the grouping variable in the metadata.
#' @param colours A vector of colors to use for the groups.
#' @param add_ellipse A logical value indicating whether to add ellipses around the groups. Default is `FALSE`.
#' @param size A numeric value specifying the size of the points. Default is `3`.
#' @param plot A character string specifying the type of plot, either `"MDS"` (default) or `"NMDS"`.
#'
#' @return A ggplot object representing the PCoA or NMDS plot.
#'
#' @examples
#' # Example usage:
#' pcoa_plot <- plot_pcoa(func_profile, what = "KO", var = "Treatment", colours = c("red", "green"), add_ellipse = TRUE)
#'
plot_pcoa <- function(func_profile,
                      what,
                      var,
                      colours,
                      add_ellipse = F,
                      size = 3,
                      plot = "MDS") {
  if (plot == "MDS") {
    dist <- calc_dist(func_profile, feat_type = what)
    mds <- cmdscale(dist, eig = T)
    x_label <- "PCo1"
    y_label <- "PCo2"
    colnames(mds$points) <- c(x_label, y_label)


  } else if (plot == "NMDS") {
    # NMDS
    mds <- metaMDS(
      comm = func_profile[[what]],
      distance = "bray",
      trace = FALSE,
      autotransform = FALSE
    )
    x_label <- "NMDS1"
    y_label <- "NMDS2"
    colnames(mds$points) <- c("PCo1", "PCo2")
  }


  mds_dat <- mds$points %>%
    merge(func_profile[["Metadata"]], by = 0)
  expvar <- (eigenvals(mds) / sum(eigenvals(mds)))[1:2]

  mds_plot <- ggplot(mds_dat, aes(PCo1, PCo2, color = .data[[var]])) +
    geom_point(size = size, alpha = 0.75) +
    scale_color_manual(values = colours) +
    scale_fill_manual(values = colours) +
    theme_cowplot() +
    xlab(paste(x_label, round(expvar[1] * 100, digits = 1), "%")) +
    ylab(paste(y_label, round(expvar[2] * 100, digits = 1), "%"))


  if (add_ellipse == TRUE) {
    mds_plot <- mds_plot +
      geom_polygon(stat = "ellipse", aes(fill = .data [[var]]), alpha = 0.3)
  }

  return(mds_plot)

}

#' Plot Beta Diversity Ordination with Statistical Annotations
#'
#' This function creates a beta diversity ordination plot (PCoA/NMDS) from a phyloseq
#' object with PERMANOVA (adonis) results displayed as a caption. It also checks for
#' homogeneity of group dispersions and warns if significant.
#'
#' @param ps A phyloseq object containing the microbiome data and sample metadata.
#' @param dist_matrix A distance matrix (class "dist") calculated from the phyloseq object.
#' @param ordination An ordination object (e.g., from `ordinate()` or `calc_betadiv()`).
#' @param group_variable A character string specifying the metadata variable to use for
#'   grouping samples.
#' @param add_ellipse A logical value indicating whether to add confidence ellipses around
#'   groups. Default is `FALSE`.
#' @param cols A vector of colors to use for the groups. Default is `NULL` (auto-generated).
#' @param shape A character string specifying an optional metadata variable to use for
#'   point shapes. Default is `NULL`.
#'
#' @return A ggplot object representing the ordination plot with PERMANOVA statistics
#'   in the caption.
#'
#' @examples
#' # Example usage:
#' # beta <- calc_betadiv(ps, dist = "bray", ord_method = "PCoA")
#' # plot_beta_div(ps, beta$Distance_Matrix, beta$Ordination, "Treatment")
#'
#' @export
plot_beta_div <- function(ps,
                          dist_matrix,
                          ordination,
                          group_variable,
                          add_ellipse = FALSE,
                          cols = NULL,
                          shape = NULL) {


  if (!inherits(ps, "phyloseq")) {
    stop("Input 'ps' must be a phyloseq object")
  }

  if (is.null(cols)) {
    cols <- calc_pal(ps, group_variable)
  }
  # significance
  ad <- phyloseq_adonis(ps, dist_matrix, group_variable)
  betadisp <- phyloseq_betadisper(ps, dist_matrix, group_variable)

  # check sample data dimensions >1 otherwise phyloseq
  # fails to colour samples due to bug:https://github.com/joey711/phyloseq/issues/541
  if (dim(sample_data(ps))[2] < 2) {
    # add repeat of first column as dummy
    sample_data(ps)[, 2] <- sample_data(ps)[, 1]
  }
  plot <- plot_ordination(ps, ordination , color = group_variable,  shape = shape)
  plot$layers[[1]] <- NULL

  plot_out <- plot + geom_point(size = 3, alpha = 0.75) +
    theme_bw() +
    scale_fill_manual(values = cols) +
    scale_color_manual(values = cols) +
    labs(caption = bquote(Adonis ~ R ^ 2 ~ .(round(ad$R2[1], 2)) ~
                            ~ p - value ~ .(ad$`Pr(>F)`[1])))
  if (add_ellipse == TRUE){
    plot_out <- plot_out +
      geom_polygon(stat = "ellipse", aes(fill = .data [[ group_variable ]] ), alpha = 0.3)
  }

  # ideally this needs to be added to the main plot eventually underneath adonis
  if (betadisp$`Pr(>F)`[[1]] < 0.05){
    warning("Group dispersion is not homogenous, interpret results carefully",
            call. = F)
  }

  return(plot_out)
}

#' Plot Volcano Plot for Differential Analysis
#'
#' This function generates a volcano plot from a results data frame, showing fold changes and adjusted p-values. It highlights significant features based on specified thresholds.
#'
#' @param results_df A list containing the results of a differential analysis, including log2 fold changes and adjusted p-values.
#' @param pthresh A numeric value specifying the p-value cutoff for significance. Default is `0.05`.
#' @param FCthresh A numeric value specifying the log2 fold change cutoff for significance. Default is `0.58`.
#' @param ... Additional arguments passed to the `EnhancedVolcano` function.
#'
#' @return A volcano plot highlighting significant features based on fold change and p-value thresholds.
#'
#' @examples
#' # Example usage:
#' volcano_plot <- plot_volcano(results_df, pthresh = 0.05, FCthresh = 0.58)
#'
#' @export

plot_volcano <- function(results_df,
                         pthresh = 0.05,
                         FCthresh = 0.58,
                         ...) {
  # this should be modified to plot results from maaslin
  res <- results_df$res
  # replace rownames with highest classified level
  rownames(res) <- as.character(lapply(strsplit(as.character(rownames(
    res
  )), split = "\\|"), tail, n = 1))
  groups <- results_df$groups
  EnhancedVolcano(
    res,
    x = "log2FC",
    y = "p.adj",
    lab = rownames(res),
    labSize = 4,
    pCutoff = pthresh,
    FCcutoff = FCthresh,
    title = paste(groups[1], "versus", groups[2], sep = " "),
    subtitle = NULL,
    max.overlaps = 30,
    ...
  )
}

#' Plot Taxonomic Composition as Stacked Bar Chart
#'
#' This function creates a stacked bar chart showing the relative abundance of taxa
#' at a specified taxonomic level, faceted by a grouping variable.
#'
#' @param ps A phyloseq object containing the taxonomic data.
#' @param tax_level A character string specifying the taxonomic rank to plot.
#'   Must be one of: "Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species".
#' @param var A character string or symbol specifying the metadata variable to facet by.
#' @param ord A character vector specifying the order of factor levels for the grouping
#'   variable. Default is `NULL` (uses default ordering).
#' @param n_taxa An integer specifying the number of top taxa to display. Default is `10`.
#' @param per_group A logical value indicating whether to select top taxa per group
#'   rather than overall. Default is `FALSE`.
#' @param groups A character vector of group names to use when `per_group = TRUE`.
#'   Default is `NULL`.
#' @param agg A logical value indicating whether to aggregate taxa at the specified
#'   taxonomic level before plotting. Default is `FALSE`.
#'
#' @return A ggplot object representing the stacked bar chart of taxonomic composition.
#'
#' @examples
#' # Example usage:
#' # comp_plot <- plot_taxonomic_comp(ps, "Phylum", "Treatment", n_taxa = 8)
#' # comp_plot <- plot_taxonomic_comp(ps, "Genus", "Group", ord = c("Control", "Treatment"))
#'
#' @export
plot_taxonomic_comp  <- function(ps, tax_level, var, ord=NULL, n_taxa=10,
                                 per_group=F, groups=NULL, agg=F) {
  ranks <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
  if (!tax_level %in% ranks){
    stop("Provide a valid taxonomic rank")
  }
  var <- enquo(var)

  if (per_group == T){
    top_tax <- map(groups, function(x) ps %>%
                     get_top_n_group(n = n_taxa, level = tax_level,var= !!var,
                                     group =x, agg=agg))
    top_tax <- c(unique(unlist(top_tax)), "other")
    # redefine n_taxa as in this case it reflects n across groups
    n_taxa <- length(top_tax)
  } else {
    top_tax <- c(get_top_n(ps, n=n_taxa, level = tax_level), "other")
  }

  if (!is.null(ord)) {
    ps <- ps %>%
      ps_mutate(plot_var = factor(.data[[var]], levels = ord))
  }

  if (agg == T){
    melt_ps <- ps %>%
      aggregate_taxa(level = tax_level) %>%
      psmelt()
  } else {
    melt_ps <- ps %>%
      psmelt()
  }

  taxa_plot <- melt_ps %>% dplyr::filter(OTU %in% top_tax)

  # group rest into other
  repl_cols <- c("OTU", ranks)

  other_plot <- melt_ps %>% dplyr::filter(!OTU %in% top_tax) %>%
    group_by(Sample) %>%
    mutate(Abundance = sum(Abundance)) %>%
    ungroup %>%
    distinct(Sample, .keep_all = T) %>%
    mutate_at(vars(repl_cols), ~ paste("other"))

  # bind dfs
  plot_df <- bind_rows(taxa_plot, other_plot)

  taxa_pal <- c(stacked_bar.palette[1:length(top_tax)-1], "#DCDCDC", "white")
  names(taxa_pal) <- plot_df %>% dplyr::filter(OTU %in% top_tax) %>%
    pull(.data[[ tax_level ]]) %>%
    unique()

  comp_fig <- plot_df %>%
    ggplot(aes(fill=.data[[ tax_level ]], y=Abundance, x=Sample)) +
    geom_bar(position="fill", stat="identity") +
    scale_fill_manual(values = taxa_pal) +
    facet_grid(. ~ factor(plot_var),
               scales = "free", space = "free") +
    theme_cowplot() +
    scale_y_continuous(labels = scales::percent, expand = c(0, 0)) +
    theme(
      panel.background = ggplot2::element_blank(),
      panel.grid = ggplot2::element_blank(),
      axis.text.x = element_blank(),
      axis.title.x = element_blank(),
      axis.title.y = element_text(size = 22),
      axis.text.y = element_text(size = 18),
      axis.ticks = element_blank()
    )
  #return_list <- list("Figure" = comp_fig, "Other" = other_plot)

  return(comp_fig)
}
