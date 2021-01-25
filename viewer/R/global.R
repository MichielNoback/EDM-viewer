library(shiny)
library(dplyr)
library(ggplot2)
library(readr)
library(here)

## DATA LOADING ##

get_celltypes <- function() {
    #c("dicty", "neutro")
    c("dicty")
}

## Reads the primary data, given a celltype
##long_id	uniprot	mw	is_grouping	spectral_count	bait	condition	z_score_spectral
read_primary <- function(celltype) {
    read_delim(file = paste0(here(), "/data/all_data_", celltype, ".txt"),
               delim = "\t+",
               col_types = c("ccdliffd"))
}

## Prepares all data for the app
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
    all_data
}

## LOAD THE DATA ##
all_data <- prepare_data()


## UTILITY FUNCTIONS ##


## DEDICATED UTILITY FUNCTIONS OPERATING ON MAIN DATASET ##

## Returns a vector of all Conditions present in the dataset
get_conditions_vector <- function(celltypes) {
    tmp <- all_data %>%
        filter(source %in% celltypes) %>%
        select(condition) %>% 
        unique() %>% 
        pull(condition) #vectorized instead of tibble
#    names(tmp) <- tmp
    sort(tmp)
}

## Returns a vector of all genes present in the dataset
get_genes_vector <- function(celltypes = get_celltypes()) {
    tmp <- all_data %>%
        filter(source %in% celltypes) %>%
        select(long_id, uniprot)
    tmp_v <- pull(tmp, long_id)
    names(tmp_v) <- tmp$uniprot
    tmp_v
}

