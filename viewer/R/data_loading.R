library(readr)
library(here)
suppressMessages(library(dplyr))
suppressMessages(library(stringr))

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
    tmp <- readr::read_delim(file = paste0(here(), "/viewer/data/all_data_", celltype, ".txt"),
               delim = "\t",
               col_types = c("ccclffidcfccccc"))
    tmp %>% 
        mutate(mw = stringr::str_replace(mw, "(\\d) kDa", "\\1")) %>%
        #suppress "NAs introduced by coercion" warning
        mutate(mw = suppressWarnings(as.integer(mw)))
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
        bait = factor(),
        condition = factor(),
        spectral_count = integer(),
        z_score_spectral = double(),
        uniprot_short = character(),
        organism = factor(),
        protein_names = character(),
        gene_names = character(),
        go_biological_process = character(),
        go_molecular_function = character(),
        go_cellular_component = character(),
        source = factor()
    )
    for (celltype in get_celltypes()) {
        tmp <- read_primary(celltype)
        names(tmp) <- names(all_data[1:15])
        tmp$source <- celltype
        all_data <- bind_rows(all_data, tmp)
    }
    all_data
}

#' Creates a column that holds a short name for the protein
#' 
#' Adds a new column to the globalEnv dataset all_data with the primary accession number, if extractable.
#' @examples
#' create_name_column()
create_name_column <- function(all_data) {
    names <- str_match(all_data$uniprot, "[^\\|]+\\|([^\\|]+)$")[, 2]
    names <- ifelse(test = is.na(names), 
                     yes = all_data$uniprot,
                     no = names)
    all_data %>% mutate(name = names)
}

# create_source_organism_column <- function() {
#     organisms <- str_match(all_data$long_id, "OS=(.+) (GN)|(PE)=")[, 2]
#     organisms <- ifelse(test = is.na(organisms),
#                         yes = "NA",
#                         no = organisms)
#     all_data <<- all_data %>% mutate(organism = organisms)
# }

