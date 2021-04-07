
## DEDICATED UTILITY FUNCTIONS OPERATING ON MAIN DATASET ##
library(dplyr)
library(tidyr)

#' Returns a vector of all Conditions present in the dataset
#' 
#' Given a vector of celltypes, will return all uniquely present condition values, 
#' as character vector and sorted alphabetically.
#' 
#' @param celltypes a character . Defaults to all celltypes.
#' @return vector with conditions
#' @examples 
#' get_conditions_vector("dicty")
#' get_conditions_vector()
get_conditions_vector <- function(celltypes = get_celltypes()) {
    tmp <- all_data %>%
        filter(source %in% celltypes) %>%
        select(condition) %>% 
        unique() %>% 
        pull(condition) #vectorized instead of tibble
#    names(tmp) <- tmp
    sort(tmp)
}

#' Returns a vector of all Baits present in the dataset
#' 
#' Given a vector of celltypes, will return all uniquely present bait values, 
#' as character vector and sorted alphabetically.
#' 
#' @param celltypes a character. Defaults to all celltypes.
#' @return vector with celltypes
#' @examples 
#' get_baits_vector("dicty")
#' get_baits_vector()
get_baits_vector <- function(celltypes = get_celltypes()) {
    tmp <- all_data %>%
        filter(source %in% celltypes) %>%
        select(bait) %>% 
        unique() %>% 
        pull(bait) #vectorized instead of tibble
#    names(tmp) <- tmp
    sort(tmp)
    
}

#' Returns a vector of all genes present in the dataset
#' 
#' Given a celltype, will return all uniquely present genes, 
#' as character vector, with uniprot IDs as value and long ID as names.
#' 
#' @param celltypes the celltypes. Defaults to all available types
#' @return named vector with genes 
#' @examples 
#' get_genes_vector()
#' get_genes_vector("dicty")
get_genes_vector <- function(celltypes = get_celltypes()) {
    tmp <- all_data %>%
        filter(source %in% celltypes) %>%
        mutate("id_plus_description" = paste(uniprot, long_id, sep = "; "))
    tmp_v <- pull(tmp, uniprot)
    names(tmp_v) <- tmp$id_plus_description
    c("Search" = "", tmp_v)
}


#' Returns a vector of all source organisms found for captured proteins.
#' 
#' Given (a) celltype(s), will return all uniquely present source
#' organisms, as character vector.
#' 
#' @param celltypes the celltypes. Defaults to all available types
#' @return named vector with source organisms 
#' @examples 
#' get_source_organism_vector()
#' get_source_organism_vector("dicty")
get_source_organism_vector <- function(celltypes = get_celltypes()) {
    tmp <- all_data %>%
        filter(source %in% celltypes) %>%
        group_by(organism) %>%
        summarise(count = n(), .groups = "drop") %>%
        mutate(organism_count = paste0(organism, " [", count, "]")) %>%
        pull(organism_count)
    # tmp_v <- pull(tmp, organism_count)
    # names(tmp_v) <- tmp$organism
    tmp
}

#get_source_organism_vector()

#' Fetch all measurements for a gene
#' 
#' Will return a tibble with measurements for this gene
#' 
#' @param gene_id the uniprot identifier
#' @return tibble with all data for the gene
#' @examples 
#' get_gene_conditions_table("sp|P0A6M8|EFG_ECOLI")
#' get_gene_conditions_table("sp|P02887|DIS1B_DICDI")
get_gene_conditions_table <- function(gene_ids) {
    if (! (class(gene_ids) == "character")) {
        warning("class of gene ID is not character. Selection returns emb=pty")
        return(NULL)
    }
    #message(paste0("[get_gene_conditions_table]: selected uniprot: __", gene_ids, "__"))
    all_data %>% 
        filter(uniprot %in% gene_ids) %>%
        select(uniprot, name, condition, bait, spectral_count, z_score_spectral) %>%
        arrange(uniprot, name, condition, bait)
}


#' Creates an HTML element displaying gene info
#' 
#' Creates an html info element for a gene, to be inserted within a shiny app.
#' 
#' @param gene_id the ID of the gene to get info on.
#' @return html div fragment with several h4 header elements.
#' @examples 
#' get_gene_info_header("sp|P0A6M8|EFG_ECOLI")
#' get_gene_info_header("sp|P02887|DIS1B_DICDI")
get_gene_info_header <- function(gene_id) {
    # print(paste0("B: selected uniprot: ", gene_id))

    tmp <- all_data %>%
        filter(uniprot == gene_id) %>%
        head(1)
    used_tag <- tags$h4
    #list that results from tags$div() needs to be converted to character for
    #the downstream cat() function (within renderText() function)
    as.character(
        tags$div(used_tag(paste0("Description: ", tmp$long_id[1])), 
        used_tag(paste0("MW (kDa): ", tmp$mw[1])),
        used_tag(paste0("Source: ", tmp$source[1]))))
}



get_matching_genes_for_long_id_or_uniprot <- function(pattern) {
    if (is.null(pattern) && nchar(pattern < 3)) {
        return(tibble(long_id = character(), uniprot = character()))
    } else {
        return(all_data %>%
            filter(str_detect(long_id, regex(pattern, ignore_case = TRUE)) |
                str_detect(uniprot, regex(pattern, ignore_case = TRUE))) %>%
            select(long_id, uniprot) %>%
            distinct())
    }
}

get_matching_genes_for_uniprot <- function(uniprot_query) {
    if (is.null(uniprot_query) && length(uniprot_query == 0)) {
        return(tibble(long_id = character(), uniprot = character()))
    } else {
        return(all_data %>%
            filter(uniprot %in% uniprot_query) %>%
            group_by(uniprot) %>%
            summarize(long_id = long_id, 
                      uniprot = uniprot, 
                      conditions = n(),
                      .groups = "drop") %>%
            distinct())
    }
}



