#' Extract Taxonomy Table from a Phyloseq Object
#'
#' This function extracts the taxonomy table from a phyloseq object and converts it into a data frame.
#'
#' @param ps A phyloseq object containing the taxonomy data.
#'
#' @return A data frame representing the taxonomy table.
#'
#' @examples
#' # Example usage:
#' data(zeller_2014)
#' taxa <- taxonomy(zeller_2014)
#'
#' @export
taxonomy <- function (ps) {
  return(as.data.frame(tax_table(ps)))
}

#' Convert Sample Metadata to Data Frame
#'
#' This function converts the sample metadata in a phyloseq object into a data frame.
#' Optionally, it includes row names as a column.
#'
#' @param ps A phyloseq object containing the sample metadata.
#' @param rownames A logical value indicating whether to include row names as a column. Default is `TRUE`.
#'
#' @return A data frame containing the sample metadata.
#'
#' @examples
#' # Example usage:
#' metadata <- meta_to_df(ps, rownames = FALSE)
#'
#' @export
meta_to_df <- function(ps, rownames=T) {
  if (rownames == T){
    return(as(sample_data(ps), "data.frame"))
  } else {
    return(as(sample_data(ps), "data.frame") %>%
             rownames_to_column("sample_id"))
  }

}

#' Convert OTU Table to Data Frame
#'
#' This function converts the OTU (operational taxonomic unit) table of a phyloseq object into a data frame.
#'
#' @param ps A phyloseq object containing the OTU table.
#'
#' @return A data frame representing the OTU table.
#'
#' @examples
#' # Example usage:
#' otu_table <- ps_to_feattab(ps)
#'
#' @export
ps_to_feattab <- function(ps) {
  return(as.data.frame(ps@otu_table))
}

#' Get Top N Taxa by Mean Abundance
#'
#' This function extracts the top N taxa from a phyloseq object based on mean relative abundance.
#' Optionally, taxa can be aggregated to a specified taxonomic level.
#'
#' @param ps A phyloseq object containing the taxa data.
#' @param n An integer specifying the number of top taxa to return.
#' @param level A character string specifying the taxonomic level to aggregate taxa. Default is `"species"`.
#' @param agg A logical value indicating whether to aggregate taxa at the specified level. Default is `FALSE`.
#'
#' @return A character vector of the top N taxa by mean abundance.
#'
#' @examples
#' # Example usage:
#' top_taxa <- get_top_n(ps, n = 10, level = "species", agg = TRUE)
#'
#' @export
get_top_n <- function(ps, n, level = "species", agg=F) {
  if (level != "species" & agg == T) {
    ps <- ps %>% tax_fix %>%
      tax_agg(rank = level)
    ps <- ps_get(ps)
  }

  topn <- ps %>%
    transform(transform = "relative") %>%
    psmelt() %>%
    group_by(OTU) %>%
    summarise(Mean_abund = mean(Abundance)) %>%
    # remove those with zero abundance
    filter(Mean_abund > 0) %>%
    slice_max(Mean_abund, n = n) %>%
    pull(OTU)
  return(topn)
}

#' Get Top N Taxa by Mean Abundance per Group
#'
#' This function extracts the top N taxa from a phyloseq object based on mean relative abundance for a specific group.
#' Optionally, taxa can be aggregated to a specified taxonomic level.
#'
#' @param ps A phyloseq object containing the taxa data.
#' @param n An integer specifying the number of top taxa to return.
#' @param level A character string specifying the taxonomic level to aggregate taxa. Default is `"species"`.
#' @param var A character string specifying the grouping variable in the phyloseq object.
#' @param group A character string specifying the group to filter on.
#' @param agg A logical value indicating whether to aggregate taxa at the specified level. Default is `FALSE`.
#'
#' @return A character vector of the top N taxa by mean abundance for the specified group.
#'
#' @examples
#' # Example usage:
#' top_taxa_group <- get_top_n_group(ps, n = 10, level = "species", var = "Treatment", group = "Control")
#'
#' @export
get_top_n_group <- function(ps, n, level = "species",var,
                            group=NULL, agg=F) {
  if (level != "species" & agg == T) {
    ps <- ps %>% tax_fix %>%
      tax_agg(rank = level)
    ps <- ps_get(ps)
  }
  topn <- ps %>%
    transform(transform = "relative") %>%
    psmelt() %>%
    filter({{ var }} == {{ group }}) %>%
    group_by(OTU) %>%
    summarise(Mean_abund = mean(Abundance)) %>%
    # remove those with zero abundance
    filter(Mean_abund > 0) %>%
    slice_max(Mean_abund, n = n) %>%
    pull(OTU)

  return(topn)
}



#' Load Phyloseq Object from ASV Table, Taxonomy, and Mapping
#'
#' This function converts ASV (amplicon sequence variant) table, taxonomy, and sample mapping data into a phyloseq object. It also supports adding a phylogenetic tree.
#'
#' @param asvtab A data frame representing the ASV table.
#' @param taxa A data frame representing the taxonomy data.
#' @param mapping A data frame representing the sample metadata.
#' @param tree A character string specifying the file path of the phylogenetic tree (optional).
#'
#' @return A phyloseq object containing the ASV table, taxonomy, metadata, and optionally the phylogenetic tree.
#'
#' @examples
#' # Example usage:
#' phyloseq_obj <- load_phylo(asvtab, taxa, mapping, tree = "tree.nwk")
#'
#' @export
load_phylo <- function(asvtab, taxa, mapping, tree = NULL) {
  # convert to phyloseq and return list
  phylo_asv <- otu_table(asvtab, taxa_are_rows = T)

  phylo_tax <- tax_table(as.matrix(taxa))

  phylo_map <- sample_data(mapping)

  if (exists("tree")) {
    phylo_tree <- read_tree(tree)
    return(merge_phyloseq(phylo_asv, phylo_tax, phylo_tree, phylo_map))
  }
  else {
    return(merge_phyloseq(phylo_asv, phylo_tax, phylo_map))
  }
}

#' Map Taxonomic Level to Abbreviation
#'
#' This function maps a specified taxonomic level to its corresponding abbreviation.
#'
#' @param level A character string representing the taxonomic level (e.g., "Species", "Genus", "Family").
#'
#' @return A character string representing the abbreviation of the specified taxonomic level.
#'
#' @examples
#' # Example usage:
#' level_abbreviation <- return_level("Genus")
#'
#' @export
return_level <- function(level) {
  if (level %in% c("Subspecies", "subspecies", "t")) {
    return("t")
  } else if (level %in% c("Species", "species", "s")) {
    return("s")
  } else if (level %in% c("Genus", "genus", "g")) {
    return("g")
  } else if (level %in% c("Family", "family", "f")) {
    return("f")
  } else if (level %in% c("Order", "order", "o")) {
    return("o")
  } else if (level %in% c("Class", "class", "c")) {
    return("c")
  } else if (level %in% c("Phylum", "phylum", "p")) {
    return("p")
  } else if (level %in% c("Kingdom", "kingdom", "k",
                          "Domain", "domain", "d")) {
    return("k")
  } else {
    stop("Not a valid taxonomic rank")
  }
}

#' Select Taxonomic Rank from Merged Table
#'
#' This function extracts the specified taxonomic rank from a merged table that contains taxonomic classifications. It supports different taxonomic formats such as Metaphlan and GTDB.
#'
#' @param merged_table A data frame containing the merged count and taxonomic information.
#' @param level A character string specifying the taxonomic level to extract (e.g., "Genus", "Species").
#' @param tax_format A character string specifying the format of the taxonomic classification. Options are `"metaphlan"` or `"gtdb"`. Default is `"metaphlan"`.
#'
#' @return A list containing two data frames:
#' - `counts`: The count data filtered to include only the specified taxonomic level.
#' - `tax_table`: The taxonomy table filtered to include only the specified taxonomic level.
#'
#' @examples
#' # Example usage:
#' result <- select_rank(merged_table, level = "Genus", tax_format = "metaphlan")
#' counts <- result$counts
#' tax_table <- result$tax_table
#'
select_rank <- function(merged_table, level, tax_format = "metaphlan") {
  level <- return_level(level)

  tax_table <-  merged_table %>%
    rownames_to_column("taxa") %>%
    dplyr::select(taxa)

  if (tax_format == "metaphlan") {
    tax_table_sep <- tax_table %>%
      rowwise() %>%
      separate(
        taxa,
        into = c(
          "Kingdom",
          "Phylum",
          "Class",
          "Order",
          "Family",
          "Genus",
          "Species",
          "Subspecies"
        ),
        remove = F,
        sep = "\\|"
      )
  } else if (tax_format == "gtdb") {
    tax_table_sep <- tax_table %>%
      rowwise() %>%
      separate(
        taxa,
        into = c(
          "Kingdom",
          "Phylum",
          "Class",
          "Order",
          "Family",
          "Genus",
          "Species",
          "Subspecies"
        ),
        remove = F,
        sep = ";"
      )
  }

  tax_table_fill <- t(zoo::na.locf(t(tax_table_sep))) %>%
    as.data.frame()

  tax_table_out <- tax_table_fill %>%
    rowwise() %>%
    separate(Subspecies, into = c("Level", "tmp"), sep = "_") %>%
    select(-tmp) %>%
    filter(Level == level) %>%
    column_to_rownames("taxa") %>%
    select(-Level)

  merged_table_out <- merged_table %>%
    rownames_to_column("taxa") %>%
    dplyr::filter(taxa %in% rownames(tax_table_out)) %>%
    column_to_rownames("taxa")

  return_list <-
    list(counts = merged_table_out, tax_table = tax_table_out)
  return(return_list)
}

#' Import Phyloseq Object from Merged Table and Metadata
#'
#' This function reads in a merged taxonomic count table and metadata, selects the specified taxonomic level, and converts the data into a phyloseq object. It supports different taxonomic formats such as Metaphlan and GTDB.
#'
#' @param merged_table_path A character string representing the file path of the merged count table.
#' @param metapath A character string representing the file path of the metadata.
#' @param level A character string specifying the taxonomic level to extract (e.g., "Genus", "Species").
#' @param tax_format A character string specifying the format of the taxonomic classification. Options are `"metaphlan"` or `"gtdb"`. Default is `"metaphlan"`.
#' @param table_comment_char A character string specifying the comment character in the merged table file. Default is `"#"`.
#' @param meta_comment_char A character string specifying the comment character in the metadata file. Default is `"#"`.
#'
#' @return A phyloseq object containing the filtered count data, taxonomic table, and metadata.
#'
#' @examples
#' # Example usage:
#' phyloseq_obj <- import_pseq_metag("path/to/merged_table.txt", "path/to/metadata.txt", level = "Genus")
#'
#' @export
import_pseq_metag <- function(merged_table_path, metapath, level,
                              tax_format="metaphlan") {
  # read files

  merged_table <- read_tab_delim_metag(merged_table_path)
  metadata <- read_tab_delim_metag(metapath)

  merged_table_rank <- select_rank(merged_table, level=level, tax_format=tax_format)

  out <- load_phylo(merged_table_rank$counts, merged_table_rank$tax_table, metadata)

  return(out)
}
