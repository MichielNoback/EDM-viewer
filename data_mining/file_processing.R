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
    
    # Z-score normalisation, by group and total, background should be excluded?
    all_data <- all_data %>% 
      group_by(bait, condition) %>% 
      mutate(grouped_z_score_spectral = scale(corr_spectral_count))
    
    all_data <- all_data %>% 
      mutate(z_score_spectral = scale(corr_spectral_count))
}


# get the data
all_data_dicty <- get_data('dicty')
all_data_neutro <- get_data('neutro')


#write.table(all_data_dicty , file = "all_data_dicty.txt", row.names = FALSE, quote = FALSE, sep = '\t')
#write.table(all_data_neutro , file = "all_data_neutro.txt", row.names = FALSE, quote = FALSE, sep = '\t')
