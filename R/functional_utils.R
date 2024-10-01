#' Read Tab-Delimited Metagenomics File
#'
#' This function reads a tab-delimited file into a data frame using specific parameters, such as no quotes and not checking column names.
#'
#' @param df A character string representing the file path of the tab-delimited file.
#'


#' Import Functional Profile
#'
#' This function imports functional profiles (such as UniRef, EC, KO, and pathways), processes them, and applies optional filtering and renormalization.
#'
#' @param meta A data frame or file path representing the metadata.
#' @param un A data frame or file path representing the UniRef data.
#' @param ec A data frame or file path representing the EC (Enzyme Commission) data.
#' @param ko A data frame or file path representing the KO (KEGG Orthology) data.
#' @param py A data frame or file path representing the pathway data.
#' @param filter_ungrouped A logical value indicating whether to filter out ungrouped, unmapped, and unintegrated data. Default is `TRUE`.
#' @param renorm A logical value indicating whether to renormalize the data after filtering. Default is `TRUE`.
#' @param read A logical value indicating whether to read the input data from files. Default is `FALSE`.
#' @param fix_names A logical value indicating whether to clean the sample names by removing suffixes. Default is `TRUE`.
#'
#' @return A list containing the processed UniRef, EC, KO, and pathway functional profiles.
#'
#' @examples
#' # Example usage:
#' func_profile <- import_func_profile(meta, un, ec, ko, py, read = TRUE)
#'
#' @export
import_func_profile <- function(meta, un, ec, ko, py, filter_ungrouped=T, renorm=T,
                                read=F, fix_names=T){

  if (read == TRUE){

    un <- read_tab_delim_metag(un)
    ec <- read_tab_delim_metag(ec)
    ko <- read_tab_delim_metag(ko)
    py <- read_tab_delim_metag(py)
    meta <- read_tab_delim_metag(meta)
  }

  feat <- list(un, ec, ko, py)
  names(feat) <- c("UNIREF", "EC", "KO", "pathways")

  if (fix_names == T){
    colnames(feat$UNIREF) <- gsub("_merged_clean_Abundance-RPKs", "", colnames(feat$UNIREF))
    colnames(feat$EC) <- gsub("_merged_clean_Abundance-RPKs", "", colnames(feat$EC))
    colnames(feat$KO) <- gsub("_merged_clean_Abundance-RPKs", "", colnames(feat$KO))
    colnames(feat$pathways) <- gsub("_merged_clean_Abundance", "", colnames(feat$pathways))

  }

  sample_names <- meta %>% rownames()

  feat <- map(feat, function(x) x %>%
                select(all_of(sample_names)) )

  if (filter_ungrouped == TRUE){
    feat <- map(feat, function(x) x %>%
                  rownames_to_column("tmp") %>%
                  filter(! tmp %in% c("UNGROUPED", "UNMAPPED", "UNINTEGRATED")) %>%
                  column_to_rownames("tmp"))
    if (renorm == T){
      feat <- map(feat, function(x) t_df(100 * t_df(x) / colSums(x)))
    }
  }
  # annotate
  feat <- map(c("UNIREF", "EC","KO", "pathways"), function(x) annotate_func(func_profile = feat,
                                                                            feat_type = x))

  names(feat) <- c("UNIREF", "EC", "KO", "pathways")

  # add EC prefix
  colnames(feat[["EC"]]) <- paste0("EC:", colnames(feat[["EC"]]))

  return_list <- list("UNIREF" = feat[["UNIREF"]],
                      "EC" = feat[["EC"]],
                      "KO" = feat[["KO"]],
                      "pathways" = feat[["pathways"]],
                      "Metadata" = meta)

  return(return_list)

}

#' Extract Feature Tables from Functional Profile
#'
#' This function extracts and returns the feature tables (excluding metadata) from a functional profile object.
#'
#' @param func_profile A list containing functional profiles and metadata.
#'
#' @return A list containing only the feature tables (e.g., UniRef, EC, KO, pathways) without metadata.
#'
#' @examples
#' # Example usage:
#' feat_tables <- get_feat_tables(func_profile)
#'
#' @keywords internal
get_feat_tables <- function(func_profile) {
  feat_tables <- within(func_profile, rm("Metadata"))
  return(feat_tables)
}



#' Get Database Path for Functional Annotations
#'
#' This function returns the path to a specific functional annotation database file depending on the operating system and the type of database requested.
#'
#' @param db A character string specifying the type of functional database. Supported options are `"KO"`, `"EC"`, and `"pathways"`.
#'
#' @return A character string representing the path to the functional database file.
#'
#' @examples
#' # Example usage:
#' db_path <- get_db_paths("KO")
#'
#' @keywords internal
get_db_paths <- function(db) {
  if (db == "KO") {
    # needs to be stored online
    return("D:/Users/adam-/Data/Metagenomics/data/ko_info.tsv")
  } else if (db == "EC") {
    return("D:/Users/adam-/Data/Metagenomics/data/ec_level4_info.tsv")
  } else if (db == "pathways") {
    return("D:/Users/adam-/Data/Metagenomics/data/metacyc_pathways_info.txt")
  }

}

#' Annotate Functional Profiles with Descriptive Labels
#'
#' This function annotates functional profiles (e.g., UniRef, EC, KO, pathways) with descriptive names using external mapping files.
#'
#' @param func_profile A list containing functional profiles.
#' @param feat_type A character string specifying the feature type to annotate. Supported types are `"UNIREF"`, `"EC"`, `"KO"`, and `"pathways"`.
#'
#' @return The functional profile annotated with descriptive labels for each feature.
#'
#' @examples
#' # Example usage:
#' annotated_profile <- annotate_func(func_profile, feat_type = "KO")
#'
#' @export
annotate_func <- function(func_profile, feat_type=NULL) {



  if (feat_type %in% c("UNIREF", "EC", "KO", "pathways")){

    if (feat_type == "KO"){
      mapfile <- read_tsv(get_db_paths("KO"),
                          col_names = c("feat", "name"))
    } else if (feat_type == "EC"){
      mapfile <- read_tsv(get_db_paths("EC"),
                          col_names = c("feat", "name"))
    } else {
      mapfile <- read_tsv(get_db_paths("pathways"),
                          col_names = c("feat", "name"))
    }

  }

  if (feat_type == "UNIREF"){
    return(func_profile[["UNIREF"]] %>% t_df())
  } else {
    lookup <- mapfile$feat
    names(lookup) <- mapfile$name

    feat <- func_profile[[feat_type]]

    feat$name <- lapply(rownames(feat), function(x) names(lookup)[match(x, lookup)])

    feat <- feat %>% rownames_to_column("id") %>%
      unite(feat_type, all_of(c("id", "name")), sep = "_", remove = T) %>%
      column_to_rownames(var="feat_type") %>%
      t_df()

    func_profile[feat_type] <- feat
  }

}


