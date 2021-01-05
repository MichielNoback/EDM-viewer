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


# read raw data, pass: 'dicty' or 'neutro' as argument
get_data <- function(what){
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
  
  if (what == 'dicty'){

    all_data <- do_dicty_background_correction(all_data)
    
    # Z-score normalisation, by group and total, background should be excluded?
    all_data <- all_data %>% 
      group_by(bait, condition) %>% 
      mutate(grouped_z_score_spectral = scale(corr_spectral_count))
    
    all_data <- all_data %>% 
      mutate(z_score_spectral = scale(corr_spectral_count))
    
  } else if (what == 'neutro'){
    all_data <- all_data %>% 
      group_by(bait, condition) %>% 
      mutate(grouped_z_score_spectral = scale(spectral_count))
    
    all_data <- all_data %>% 
      mutate(z_score_spectral = scale(spectral_count))
  }
}

do_dicty_background_correction <- function(my_dicty_data){
  
  # can we subtract raw counts or do we need to subtract the normalized score?
  
  all_data_dicty_corrected <- my_dicty_data %>% 
    inner_join(my_dicty_data, by = 'uniprot')  %>%
    mutate(corr_spectral_count = if_else(condition.y != 'BACKGROUND', spectral_count.x, spectral_count.x - spectral_count.y))  %>% 
    filter(condition.x != 'BACKGROUND' & condition.y == 'BACKGROUND') %>% 
    select(-ends_with(".y")) %>% 
    rename_with(~ gsub(".x", "", .x))
  
  all_data_dicty_non_corrected <- my_dicty_data %>% 
    inner_join(my_dicty_data, by = 'uniprot')  %>%
    mutate(corr_spectral_count = if_else(condition.y != 'BACKGROUND', spectral_count.x, spectral_count.x - spectral_count.y))  %>% 
    filter(condition.x == condition.y | condition.x != 'BACKGROUND' & condition.y == 'BACKGROUND') %>%
    filter(!uniprot %in% unique(my_dicty_data$uniprot[my_dicty_data$condition == 'BACKGROUND'])) %>%
    select(-ends_with(".y")) %>% 
    rename_with(~ gsub(".x", "", .x))
  
  # background needs work to match columns
  all_data_dicty_background <- my_dicty_data %>%
    filter(condition == 'BACKGROUND') %>% 
    mutate(corr_spectral_count = spectral_count)
  
  all_data_dicty_corrected_combined <-rbind(all_data_dicty_corrected, all_data_dicty_non_corrected, all_data_dicty_background)
  
  # TODO, negative counts?
  #all_data_dicty$spectral_count[all_data_dicty$uniprot == "sp|Q556G3|GSTA2_DICDI" & all_data_dicty$condition == 'BACKGROUND']
  #[1] 397
  
  #all_data_dicty_corrected_combined$spectral_count[all_data_dicty_corrected_combined$uniprot == 'sp|Q556G3|GSTA2_DICDI']
  #[1] 217  74  47 423 112  24 456 319  47 194 922  63 743 397
  
  #all_data_dicty_corrected_combined$corr_spectral_count[all_data_dicty_corrected_combined$uniprot == 'sp|Q556G3|GSTA2_DICDI']
  #[1] -180 -323 -350   26 -285 -373   59  -78 -350 -203  525 -334  346  397
  
  #length(which(all_data_dicty_corrected_combined$corr_spectral_count <0))
  #[1] 6240
  
  
  # TODO, different size in -> out; -1760 in original
  # all_data_dicty = 17609
  # all_data_dicty_corrected_combined = 19369
  
  #sum(duplicated(select(all_data_dicty_corrected_combined, -(corr_spectral_count))))
  #[1] 1760
  
  # get unique lines
  all_data_dicty_corrected_combined <- all_data_dicty_corrected_combined %>% 
    distinct()
}


# get the data
all_data_dicty <- get_data('dicty')
all_data_neutro <- get_data('neutro')



#write.table(all_data_dicty , file = "all_data_dicty.txt", row.names = FALSE, quote = FALSE, sep = '\t')
#write.table(all_data_neutro , file = "all_data_neutro.txt", row.names = FALSE, quote = FALSE, sep = '\t')
