library(dplyr)
library(readr)
library(tidyr)
library(stringr) 

data_dir <- '../data/'

my_files_neutro <- c(
  'Neutrophils_GST_GppNHp.txt',
  'Neutrophils_Galphai_AluF.txt',
  'Neutrophils_Galphai_AluFl_GST.txt',
  'Neutrophils_Galphai_GDP.txt',
  'Neutrophils_Galphai_GDP_GST.txt',
  'Neutrophils_Galphai2_GDP.txt',
  'Neutrophils_Galphai_GppNHp.txt',
  'Neutrophils_Rac_GppNHp.txt'
)

my_files_dicty <- c(
  'Chemotaxis_Dicty_RasC.txt',
  'Chemotaxis_Dicty_Ric8.txt',
  'Chemotaxis_Dicty_Roco4_BACKGROUND.txt',
  'Chemotaxis_Dicty_Galpha8_GppNHp.txt',
  'Chemotaxis_Dicty_Gbeta1.txt',
  'Chemotaxis_Dicty_Galpha8_GDP.txt',
  'Chemotaxis_Dicty_RapA.txt',
  'Chemotaxis_Dicty_Galpha2_starv.txt',
  'Chemotaxis_Dicty_Gbeta2.txt',
  'Chemotaxis_Dicty_Galpha4.txt',
  'Chemotaxis_Dicty_Galpha8.txt',
  'Chemotaxis_Dicty_RasB.txt',
  'Chemotaxis_Dicty_Rac1.txt',
  'Chemotaxis_Dicty_RasG1.txt',
  'Chemotaxis_Dicty_Galpha2_veg.txt'
)

# function to read exported unique protein count file from scaffold
# exported with the exported all option, otherwise the same file will contain 
# many meta information lines
get_data <- function(){
  
  spdat <- read_delim(paste0(data_dir, "overview_unique_peptide_counts_MISSING_GALPHA_GDP_and_GTP.csv"), col_names=T, delim="\t", skip=2)
  
  spdat <- spdat %>% select(!c(1,2,3,6,9)) %>%
    rename('long_id' = `Identified Proteins (2302)`, 
           'uniprot' = `Accession Number`,
           'mw' = `Molecular Weight`,
           'is_grouping' = `Protein Grouping Ambiguity`) %>%
    pivot_longer(!c(long_id, uniprot, mw, is_grouping), 
                 names_to = c("bait", "condition"), 
                 names_sep = "_", 
                 values_to = "spectral_count") %>%
    mutate(across(c(bait, condition), as.factor)) %>%
    mutate(z_score_spectral = scale(spectral_count))
}



# obsolete funtion!, used for the separate txt files with counts for each experiment
# read raw data, pass: 'dicty' or 'neutro' as argument
get_data_old <- function(what){
  data_list <- list()
  
  for (file_name in get(paste0("my_files_", what))){
    spdat <- read_delim(paste0(data_dir, file_name), col_names=F, delim="\t", skip=1)
    names(spdat) <- c('long_id','uniprot','mw','is_grouping','spectral_count')
    
    if (what == 'dicty'){
      exp_name <- gsub(file_name, pattern="Chemotaxis_Dicty_", replace="") %>% 
        gsub(pattern="\\.txt", replace="")
    }
    else if (what == 'neutro'){
      exp_name <- gsub(file_name, pattern="Neutrophils_", replace="") %>% 
        gsub(pattern="\\.txt", replace="")
    }
    
    split_expname <- strsplit(exp_name, split = "_")
    
    # will introduce na's as there are ? for kDa
    spdat <- spdat %>% 
      mutate(mw = as.numeric(gsub("kDa", "", mw)))
    
    # split bait and condition, i.e. Galpha2_veg
    spdat <- spdat %>% 
      mutate(bait = split_expname[[1]][1])
    
    # when no condition could be split, the condition was normal
    spdat <- spdat %>% 
      mutate(condition = split_expname[[1]][2]) %>% 
      replace_na(list(condition = "Normal"))
    
    data_list[[file_name]]<- spdat %>% 
      select(long_id, uniprot, mw, is_grouping, spectral_count, bait, condition)
  }
  
  all_data <- bind_rows(data_list)
  all_data <- all_data %>% 
    mutate(across(c(bait, condition), as.factor))
    
    # Z-score normalisation, by group and total
    all_data <- all_data %>% 
      group_by(bait, condition) %>% 
      mutate(grouped_z_score_spectral = scale(spectral_count))
    
    all_data <- all_data %>% 
      mutate(z_score_spectral = scale(spectral_count))
}


# get the data
all_data_dicty <- get_data()
# all_data_neutro <- get_data_old('neutro')


filter_data <- function(celltype, experiment_bool=FALSE, bait_bool=FALSE, condition_bool=FALSE, condition_in=NULL, bait_in=NULL, uniprot_in=NULL, normalized=FALSE){
  
  my_data <- NULL
  my_title <- NULL
  my_x <- NULL
  my_y <- if_else(normalized, "z_score_spectral","spectral_count")

  # get summed count for a given uniprot, group per bait_condition (experiment)
  if (!is.null(uniprot_in) & is.null(condition_in) & experiment_bool){
    
    my_x <- 'bait_condition'
    my_title <- uniprot_in
    
    my_data <- filter(get(paste0("all_data_", celltype)), grepl(uniprot_in, uniprot)) %>%
      group_by(bait, condition) %>%
      dplyr::summarise(.groups = 'drop', across(my_y, sum)) %>%
      tidyr::unite(bait_condition, c(bait, condition), sep = "_")
  }
  
  # get count for a given condition
  else if (!is.null(condition_in) & bait_bool){
    
    my_x <- 'bait'
    my_title <- condition_in
    
    my_data <- filter(get(paste0("all_data_", celltype)), grepl(condition_in, condition))
  }
  
  # get count for a given bait
  else if (!is.null(bait_in) & !experiment_bool){

    my_x <- 'condition'
    my_title <- bait_in
    
    my_data <- filter(get(paste0("all_data_", celltype)), grepl(bait_in, bait))
    
    # get count for a given bait group by bait_condition (experiment)
  }else if (experiment_bool){
    
    my_title <- bait_in
    my_x <- 'bait_condition'
    
    my_data <- filter(get(paste0("all_data_", celltype)), grepl(bait_in, bait)) %>%
      group_by(bait, condition) %>% 
      tidyr::unite(bait_condition, c(bait, condition), sep = "_")
  }

  # store labels that can be used to plot
  my_data_fig_labels <- list()
  my_data_fig_labels$title <- my_title
  my_data_fig_labels$x_label <- my_x
  my_data_fig_labels$y_label <- my_y
  
  return(list(count_data = my_data, fig_labels = my_data_fig_labels))
  
}


#write.table(all_data_dicty , file = "all_data_dicty.txt", row.names = FALSE, quote = FALSE, sep = '\t')
#write.table(all_data_neutro , file = "all_data_neutro.txt", row.names = FALSE, quote = FALSE, sep = '\t')
