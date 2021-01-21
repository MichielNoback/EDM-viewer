library(shiny)
library(dplyr)
library(ggplot2)
library(readr)
library(here)

## UTILITY FUNCTIONS ##

get_celltypes <- function() {
    c("dicty", "neutro")
}

## Reads the primary data, given a celltype
read_primary <- function(celltype) {
    read_delim(file = paste0(here(), "/data/all_data_", celltype, ".txt"),
               delim = "\t+",
               col_types = c("ccdliff"))
}
# read_primary <- function(celltype) {
#     read_delim(file = paste0(here(), "/data/all_data_", celltype, ".txt"),
#                delim = "\t+",
#                col_names = c("long_ID", "short_ID", "weight", "flag", "signal", "Condition"))
# }


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
        source = factor()
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


## DEDICATED UTILITY FUNCTIONS OPERATING ON MAIN DATASET ##

## Returns a vector of all Conditions present in the dataset
get_conditions_vector <- function(celltypes) {
    tmp <- all_data %>%
        filter(source %in% celltypes) %>%
        select(Condition) %>% 
        unique() %>% 
        pull(Condition) #vectorized instead of tibble
    names(tmp) <- tmp
    sort(tmp)
}

## Returns a vector of all genes present in the dataset
get_genes_vector <- function(celltypes) {
    tmp_tbl <- all_data %>%
        filter(source %in% celltypes) %>%
        select(long_ID, short_ID)
    tmp <- pull(tmp_tbl, long_ID)
    names(tmp) <- pull(tmp_tbl, short_ID)
    tmp
}

