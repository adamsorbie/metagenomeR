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
