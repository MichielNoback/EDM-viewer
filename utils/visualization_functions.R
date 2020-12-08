source("../data_mining/global.R")
library(ggplot2) 

# Given a uniprot id and 'dicty' or 'neutro', plot the spectral count vs the bait + condition
plot_spectral_bait_condition_using_uniprot<- function(uniprot_id, what){
  
  # filter data on substring using grepl, leaves the original uniprot annotation intact 
  my_data <- filter(get(paste0("all_data_", what)), grepl(uniprot_id, uniprot)) %>% 
    group_by(bait, condition) %>% 
    # .groups = drop is the default, will get warnings if you not specify explicitly 
    summarise(.groups = 'drop', across(spectral_count, sum)) %>% 
    unite(bait_condition, c(bait, condition), sep = "_")
  
  p <- ggplot(my_data, aes(x = bait_condition, y = spectral_count))
  p <- p + geom_col() + labs(title = uniprot_id)
  p + theme(axis.text.x = element_text(angle = 45, hjust = 1))
}


plot_spectral_using_condition<- function(condition_in, what){
  
  my_data <- filter(get(paste0("all_data_", what)), grepl(condition_in, condition))
    
  p <- ggplot(my_data, aes(x = bait, y = spectral_count))
  p <- p + geom_col() + labs(title = condition_in)
  p + theme(axis.text.x = element_text(angle = 45, hjust = 1))
}


plot_spectral_using_bait<- function(bait_in, what){
  
  my_data <- filter(get(paste0("all_data_", what)), grepl(bait_in, bait))
  
  p <- ggplot(my_data, aes(x = condition, y = spectral_count))
  p <- p + geom_col() + labs(title = bait_in)
  p + theme(axis.text.x = element_text(angle = 45, hjust = 1))
}


boxplot_spectral_using_bait<- function(bait_in, what){
  
  my_data <- filter(get(paste0("all_data_", what)), grepl(bait_in, bait))
  
  p <- ggplot(my_data, aes(x = condition, y = log10(spectral_count)))
  p <- p + geom_boxplot() + labs(title = bait_in) 

  p + theme(axis.text.x = element_text(angle = 45, hjust = 1))
}



