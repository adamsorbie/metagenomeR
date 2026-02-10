#' Transpose Data Frame
#'
#' This function transposes a data frame, converting rows to columns and vice versa.
#'
#' @param x A data frame to be transposed.
#'
#' @return A transposed version of the input data frame.
#'
#' @examples
#' # Example usage:
#' transposed_df <- t_df(original_df)
#'
#' @keywords internal
t_df <- function(x) {
  return(as.data.frame(t(x)))
}

#' Read a Tab-Delimited Metagenomics File
#'
#' This function reads a tab-delimited metagenomics file into a data frame using specific parameters such as row names and custom comment characters.
#'
#' @param df A character string representing the file path of the tab-delimited file.
#' @param comment_char A character string specifying the comment character to use. Default is `"#"`.
#'
#' @return A data frame containing the data from the tab-delimited file.
#'
#' @examples
#' # Example usage:
#' metag_data <- read_tab_delim_metag("data/tab_delimited_file.txt")
#'
#' @importFrom data.table fread
#' @export
read_tab_delim_metag <- function(df, comment_char) {
  # read all tab delimited files using these params
  df_out <-
    data.frame(
      fread(
        df,
        header = T,
        check.names = F,
        quote = ""
      ),
      row.names = 1,
      check.names = F
    )

  return(df_out)
}

#' Filter Data Frame by Row Names
#'
#' This function filters the rows of a data frame based on a specified vector of row names.
#'
#' @param df A data frame to be filtered.
#' @param filt_vector A vector of row names to keep in the data frame.
#'
#' @return A data frame filtered to include only the specified row names.
#'
#' @examples
#' # Example usage:
#' filtered_df <- filter_rownames(df, c("row1", "row2"))
#'
#' @keywords internal
filter_rownames <- function(df, filt_vector) {
  # wrapper script around filter to handle rownames
  df_filt <- df %>%
    rownames_to_column(var = "id") %>%
    dplyr::filter(id %in% filt_vector) %>%
    column_to_rownames(var = "id")
  return(df_filt)
}

#' Generate Color Palette from Phyloseq Object
#'
#' This internal function generates a color palette based on the unique levels
#' of a grouping variable in a phyloseq object's sample data.
#'
#' @param ps A phyloseq object containing sample metadata.
#' @param group_variable A character string specifying the metadata variable
#'   to use for determining the number of colors.
#'
#' @return A vector of colors from the NPG palette, or NULL with a message
#'   if the number of groups exceeds the palette limit.
#'
#' @importFrom ggsci pal_npg
#' @keywords internal
calc_pal <- function(ps, group_variable) {
  # update function to accept non-ps objects and add support for continuous palettes
  meta <- meta_to_df(ps)
  groups <- unique(meta[, group_variable])

  if (length(groups) < 10) {
    pal <- pal_npg()(length(groups))
  } else {
    message("Exceeded colour palette limit, use custom palette")
    pal <- NULL
  }
  pal
}
#' Custom colour palette for stacked barplots
#'
#' A character vector of 30 hex colour codes designed for visualizing
#' taxonomic composition in stacked barplots. The first 10 colours are
#' from the NPG (Nature Publishing Group) palette, with additional
#' distinguishable colours for datasets with many taxa.
#'
#' @format A character vector of length 30 containing hex colour codes.
#'
#' @examples
#' # Use with ggplot2 scale_fill_manual
#' # scale_fill_manual(values = stacked_bar_palette)
#'
#' # Preview the palette
#' scales::show_col(stacked_bar_palette)
#'
#' @export
stacked_bar_palette <- c(
  "#E64B35FF",
  "#4DBBD5FF",
  "#00A087FF",
  "#3C5488FF",
  "#F39B7FFF",
  "#8491B4FF",
  "#91D1C2FF",
  "#DC0000FF",
  "#7E6148FF",
  "#B09C85FF",
  "#E4E9B2",
  "#F9A620",
  "#054A29",
  "#52414C",
  "#D81E5B",
  "#331832",
  "#27474E",
  "#573D1C",
  "#404E4D",
  "#DAD4EF",
  "#E86A92",
  "#044389",
  "#6C4B5E",
  "#4E6E58",
  "#826AED",
  "#FF0054",
  "#9E0059",
  "#387D7A",
  "#395E66",
  "#1BE7FF"
)
