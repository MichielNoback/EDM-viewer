
## DEDICATED UTILITY FUNCTIONS OPERATING ON MAIN DATASET ##
suppressMessages(library(dplyr))
suppressMessages(library(tidyr))


#' Returns a tibble of genes for the given search string
#' 
#' Will look in both long-id and uniprot columns.
#' Only starts looking when pattern is at least 3 characters long.
#' @param pattern the pattern to look for.
#' @return tibble with rows matching the ID pattern 
#' @examples 
#' get_matching_genes_for_long_id_or_uniprot("chap")
get_matching_genes_for_long_id_or_uniprot <- function(pattern, only_dictyostelium = TRUE) {
    if (is.null(pattern) || nchar(pattern) < 3) {
        return(tibble(long_id = character(), uniprot = character()))
    } else {
        tmp <- all_data %>%
            filter(str_detect(long_id, regex(pattern, ignore_case = TRUE)) |
                str_detect(uniprot, regex(pattern, ignore_case = TRUE))) 

        if (only_dictyostelium) {
            tmp <- tmp %>% filter(str_detect(organism, "Dictyostelium"))
        }
        tmp <- tmp %>% 
            select(long_id, uniprot) %>%
            distinct()
        return(tmp)
    }
}

#' Returns a summary table for the given input, collapsed on conditions
#' 
#' @param genes_tibble the input of selected genes.
#' @return tibble with summary of condition counts 
#' @examples 
#' get_matching_genes_for_long_id_or_uniprot("chap")
get_genes_summary_table <- function(genes_tibble) {
    return(genes_tibble %>%
        group_by(uniprot) %>%
        summarize(long_id = long_id, 
                  uniprot = uniprot, 
                  conditions = n(),
                  .groups = "drop") %>%
        distinct())
}


#' Returns a vector of experiments represented in the given genes
#' 
#' An experiment is defined as a condition + bait combination
#' 
#' @param genes_tibble the genes to select from
#' @return vector with represented experiments
get_represented_experiments <- function(genes_tibble) {
    genes_tibble %>%
        unite(col = "experiment", condition, bait, sep = " / ") %>%
        pull(experiment) %>% 
        unique() %>% 
        sort()
}

#' returns a vector of experiments that should be preselected in selection checkbox group
get_preselected_experiments <- function() {
    c("background / Roco4", "control / RapA", "normal / Galpha2", "starv / Galpha2")
}

#' Fetch all data for a selection of genes
#' 
#' Will return a tibble all data for the requested genes
#' 
#' @param uniprot_ids the uniprot identifiers
#' @return tibble with all data for the gene
#' @examples 
#' get_gene_data_for_uniprots("sp|P0A6M8|EFG_ECOLI")
#' get_gene_data_for_uniprots(c("sp|P0A6Y8|DNAK_ECOLI", "sp|P02887|DIS1B_DICDI"))
get_gene_data_for_uniprots <- function(uniprot_ids) {
    all_data %>% 
    filter(uniprot %in% uniprot_ids) 
} 

#' Creates a summary tibble 
#' 
#' Creates a summary tibble given a larger tibble of gene data
#' 
#' @param gene_tibble the rows of gene data ro be summarized
#' @return tibble with summary data for the genes
get_gene_summary_table <- function(gene_tibble) {
        return(gene_tibble %>%
            group_by(uniprot) %>%
            summarize(long_id = long_id,
                      uniprot = uniprot,
                      conditions = n(),
                      .groups = "drop") %>%
            distinct())
}




#' [TO BE DEPRECATED] Fetch all measurements for a gene
#' 
#' Will return a tibble with measurements for this gene
#' 
#' @param gene_id the uniprot identifier
#' @return tibble with all data for the gene
#' @examples 
#' get_gene_conditions_table("sp|P0A6M8|EFG_ECOLI")
#' get_gene_conditions_table("sp|P02887|DIS1B_DICDI")
# get_gene_conditions_table <- function(gene_ids) {
#     if (! (class(gene_ids) == "character")) {
#         warning("class of gene ID is not character. Selection returns empty")
#         return(NULL)
#     }
#     #message(paste0("[get_gene_conditions_table]: selected uniprot: __", gene_ids, "__"))
#     all_data %>% 
#         filter(uniprot %in% gene_ids) %>%
#         select(uniprot, name, condition, bait, spectral_count, z_score_spectral) %>%
#         arrange(uniprot, name, condition, bait)
# }

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
# get_conditions_vector <- function(celltypes = get_celltypes()) {
#     tmp <- all_data %>%
#         filter(source %in% celltypes) %>%
#         select(condition) %>% 
#         unique() %>% 
#         pull(condition) #vectorized instead of tibble
# #    names(tmp) <- tmp
#     sort(tmp)
# }

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
# get_baits_vector <- function(celltypes = get_celltypes()) {
#     tmp <- all_data %>%
#         filter(source %in% celltypes) %>%
#         select(bait) %>% 
#         unique() %>% 
#         pull(bait) #vectorized instead of tibble
# #    names(tmp) <- tmp
#     sort(tmp)
# }

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
# get_source_organism_vector <- function(celltypes = get_celltypes()) {
#     tmp <- all_data %>%
#         filter(source %in% celltypes) %>%
#         group_by(organism) %>%
#         summarise(count = n(), .groups = "drop") %>%
#         mutate(organism_count = paste0(organism, " [", count, "]")) %>%
#         pull(organism_count)
#     # tmp_v <- pull(tmp, organism_count)
#     # names(tmp_v) <- tmp$organism
#     tmp
# }


#'
#'Function in development and not used right now
# filter_genes <- function(celltypes, 
#              conditions = NULL, 
#              baits = NULL, 
#              mw = c(-Inf, Inf),
#              spectral_count_range = c(-Inf, Inf),
#              z_score_range = c(-Inf, Inf),
#              organisms = NULL) {
#     if (is.null(celltypes) || length(celltypes) == 0) {
#         return(all_data[0, ]) #returns empty row
#     } 
#     all_data %>%
#         filter(if (!is.null(baits)) bait %in% baits else TRUE) %>%
#         filter(if (!is.null(organisms)) organism %in% organisms else TRUE )
#         #summarize(count = n())
#         #head()
#     
# }



#### DEPRECATED FUNCTIONS ####


#' Returns a tibble of genes for the given uniprot IDs
#' 
#' @param uniprot_query the IDs of the genes to get info on.
#' @return tibble with info on the requested IDs
#' @examples 
#' get_matching_genes_for_uniprot("sp|P0A6Y8|DNAK_ECOLI")
#' get_matching_genes_for_uniprot(c("sp|P0A6Y8|DNAK_ECOLI", "sp|P02887|DIS1B_DICDI")
# get_matching_genes_for_uniprot <- function(uniprot_query) {
#     if (is.null(uniprot_query) || length(uniprot_query) == 0) {
#         return(tibble(long_id = character(), uniprot = character()))
#     } else {
#         return(all_data %>%
#             filter(uniprot %in% uniprot_query) %>%
#             group_by(uniprot) %>%
#             summarize(long_id = long_id,
#                       uniprot = uniprot,
#                       conditions = n(),
#                       .groups = "drop") %>%
#             distinct())
#     }
# }

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
# get_genes_vector <- function(celltypes = get_celltypes()) {
#     tmp <- all_data %>%
#         filter(source %in% celltypes) %>%
#         mutate("id_plus_description" = paste(uniprot, long_id, sep = "; "))
#     tmp_v <- pull(tmp, uniprot)
#     names(tmp_v) <- tmp$id_plus_description
#     c("Search" = "", tmp_v)
# }

#' [Deprecated] Creates an HTML element displaying gene info
#' 
#' Creates an html info element for a gene, to be inserted within a shiny app.
#' 
#' @param gene_id the ID of the gene to get info on.
#' @return html div fragment with several h4 header elements.
#' @examples 
#' get_gene_info_header("sp|P0A6M8|EFG_ECOLI")
#' get_gene_info_header("sp|P02887|DIS1B_DICDI")
# get_gene_info_header <- function(gene_id) {
#     # print(paste0("B: selected uniprot: ", gene_id))
# 
#     tmp <- all_data %>%
#         filter(uniprot == gene_id) %>%
#         head(1)
#     used_tag <- tags$h4
#     #list that results from tags$div() needs to be converted to character for
#     #the downstream cat() function (within renderText() function)
#     as.character(
#         tags$div(used_tag(paste0("Description: ", tmp$long_id[1])),
#         used_tag(paste0("MW (kDa): ", tmp$mw[1])),
#         used_tag(paste0("Source: ", tmp$source[1]))))
# }
