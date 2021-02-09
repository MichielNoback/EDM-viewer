library(shiny)
library(dplyr)
library(tidyr)
library(ggplot2)
library(readr)
library(stringr)
library(here)

## DATA LOADING FUNCTIONS##

#' Returns the celltypes currently contained with the app
#'
#' @return A character array of celltypes
#' @examples
#' get_celltypes()
get_celltypes <- function() {
    #c("dicty", "neutro")
    c("dicty")
}


#' Reads the primary data, given a celltype. 
#' 
#' Serves as helper function to 
#' \code{prepare_data()}
#' @param celltype The base name of the files to be loaded. File names are #' assumed to be named as follows: \code{all_data_<celltype>.txt}
#' @return A tibble holding the data for the given celltype.
#' @examples
#' read_primary("dicty)
read_primary <- function(celltype) {
    read_delim(file = paste0(here(), "/data/all_data_", celltype, ".txt"),
               delim = "\t+",
               col_types = c("ccdliffd"))
}

#' Prepares all data for the app.
#' 
#' Will fetch the available cell types from \code{get_celltypes()} and 
#' return all data merged into a single tibble.
#' @return All data
#' @examples
#' prepare_data()
prepare_data <- function() {
    all_data <- tibble(
        long_id = character(),
        uniprot = character(),
        mw = double(),
        is_grouping = logical(),
        spectral_count = integer(),
        bait = factor(),
        condition = factor(),
        source = factor(),
        z_score_spectral = double()
    )
    for (celltype in get_celltypes()) {
        tmp <- read_primary(celltype)
        tmp$source <- celltype
        all_data <- bind_rows(all_data, tmp)
    }
    
    ## create readable legend labels
    names <- str_match(all_data$uniprot, "[^\\|]+\\|([^\\|]+)$")[, 2]
    names <- ifelse(test = is.na(names), 
                     yes = all_data$uniprot,
                     no = names)
    all_data <- all_data %>% mutate(name = names)
    all_data
}

## DATA LOADING ##

#' All data is loaded automatically 
all_data <- prepare_data()


## UTILITY FUNCTIONS ##


## DEDICATED UTILITY FUNCTIONS OPERATING ON MAIN DATASET ##

#' Returns a vector of all Conditions present in the dataset
#' 
#' Given a celltype, will return all uniquely present condition values, 
#' as character vector and sorted alphabetically.
#' 
#' @param celltypes a character 
#' @return named vector with celltypes
#' @examples 
#' get_conditions_vector("dicty")
get_conditions_vector <- function(celltypes) {
    tmp <- all_data %>%
        filter(source %in% celltypes) %>%
        select(condition) %>% 
        unique() %>% 
        pull(condition) #vectorized instead of tibble
#    names(tmp) <- tmp
    sort(tmp)
}

#' Returns a vector of all genes present in the dataset
#' 
#' Given a celltype, will return all uniquely present genes, 
#' as character vector, with uniprot IDs as value and long ID as names.
#' 
#' @param celltypes the celltypes
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
            filter(str_detect(long_id, pattern) |
                           str_detect(uniprot, pattern)) %>%
            select(long_id, uniprot) %>%
            distinct())
    }
}
