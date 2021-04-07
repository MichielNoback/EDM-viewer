library(readr)
library(here)
library(stringr)

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
    all_data
}

create_name_column <- function() {
    names <- str_match(all_data$uniprot, "[^\\|]+\\|([^\\|]+)$")[, 2]
    names <- ifelse(test = is.na(names), 
                     yes = all_data$uniprot,
                     no = names)
    all_data <<- all_data %>% mutate(name = names)
    all_data
}

create_source_organism_column <- function() {
    organisms <- str_match(all_data$long_id, "OS=(.+) (GN)|(PE)=")[, 2]
    organisms <- ifelse(test = is.na(organisms),
                        yes = "NA",
                        no = organisms)
    all_data <<- all_data %>% mutate(organism = organisms)
}

## DATA LOADING ##

#' All data is loaded automatically 
all_data <- prepare_data()

#' create readable legend labels
create_name_column()

#' create source organism column
create_source_organism_column()

#report
message(paste0("[global.R]: loaded ", nrow(all_data), " data rows"))

