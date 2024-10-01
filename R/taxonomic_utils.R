#' @param ps phyloseq object
#' @return Returns the taxonomy table slot as a dataframe
#' @export
#'
#' @examples
#' data(zeller_2014)
#' taxa <- taxonomy(zeller_2014)
#'
taxonomy <- function (ps) {
  return(as.data.frame(tax_table(ps)))
}


meta_to_df <- function(ps, rownames=T) {
  if (rownames == T){
    return(as(sample_data(ps), "data.frame"))
  } else {
    return(as(sample_data(ps), "data.frame") %>%
             rownames_to_column("sample_id"))
  }

}

ps_to_feattab <- function(ps) {
  return(as.data.frame(ps@otu_table))
}

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

filter_rownames <- function(df, filt_vector) {
  # wrapper script around filter to handle rownames
  df_filt <- df %>%
    rownames_to_column(var="id") %>%
    filter(id %in% filt_vector) %>%
    column_to_rownames(var="id")
  return(df_filt)
}

common_cols <- function(list_df, n=length(list_df)) {
  col_list <- map(list_df,names)

  if (n < length(list_df)){
    n_cols <- table(unlist(col_list))
    return(names(n_cols[n_cols >= n]))
  }

  return(Reduce(intersect,col_list))
}


combine_counts <- function(counts_list, taxa_are_cols=T, prev_thresh=0.3) {
  if (taxa_are_cols == F){
    counts_list <- map(counts_list, function(x) t(x) %>% as.data.frame)
  }
  # purely presence absence
  shared_taxa <- common_cols(counts_list, n=round(length(counts_list) * prev_thresh))

  combined_counts <- bind_rows(counts_list) %>%
    select(all_of(shared_taxa))
  return(combined_counts)
}

combine_meta <- function(meta_list, col_list) {

  combined_meta <- bind_rows(meta_list) %>%
    select(all_of(col_list))
}

t_df <- function(x) {
  return(as.data.frame(t(x)))
}

read_tab_delim_metag <- function(df, comment_char="#") {
  # read all tab delimited files using these params
  df_out <-
    read.table(
      df,
      row.names = 1,
      header = 1,
      sep = "\t",
      check.names = F,
      quote = "",
      comment.char = comment_char
    )
  return(df_out)
}

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

import_pseq_metag <- function(merged_table_path, metapath, level,
                              tax_format="metaphlan", table_comment_char="#",
                              meta_comment_char="#") {
  # read files

  merged_table <- read_tab_delim_metag(merged_table_path, comment_char = table_comment_char)
  metadata <- read_tab_delim_metag(metapath, comment_char = meta_comment_char)

  merged_table_rank <- select_rank(merged_table, level=level, tax_format=tax_format)

  out <- load_phylo(merged_table_rank$counts, merged_table_rank$tax_table, metadata)

  return(out)
}
