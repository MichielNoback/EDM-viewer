source("../data_mining/file_processing.R")
library(ggplot2) 

filter_data_and_plot <- function(celltype, experiment=FALSE, condition_in=NULL, bait_in=NULL, uniprot_id=NULL, normalized=FALSE, viz_type='bar', log_scale=NULL){
  
  my_data <- NULL
  my_title <- NULL
  my_x <- NULL
  
  my_y <- if_else(normalized, "z_score_spectral","spectral_count")
  my_y_label <- if_else(!is.null(log_scale), paste0("log",log_scale, "(", my_y, ")"), my_y)
  
  viz_type <- if_else(viz_type=="bar", "geom_col", "geom_boxplot")
  viz_function <- get(viz_type, envir=environment(ggplot))
  
  # get count for a given uniprot
  if (!is.null(uniprot_id) & is.null(condition_in)){
    my_x <- 'bait_condition'
    my_title <- uniprot_id
    
    my_data <- filter(get(paste0("all_data_", celltype)), grepl(uniprot_id, uniprot)) %>%
      group_by(bait, condition) %>%
      dplyr::summarise(.groups = 'drop', across(my_y, sum)) %>%
      unite(bait_condition, c(bait, condition), sep = "_")
  }
  # get count for a given condition
  else if (!is.null(condition_in)){
    my_x <- 'bait'
    my_title <- condition_in
    
    my_data <- filter(get(paste0("all_data_", celltype)), grepl(condition_in, condition))
  }
  # get count for a given bait
  else if (!is.null(bait_in) & !experiment){
    my_x <- 'condition'
    my_title <- bait_in
    
    my_data <- filter(get(paste0("all_data_", celltype)), grepl(bait_in, bait))
  # get count for a given bait per experiment (bait_condition)
  }else if (experiment){
    my_title <- bait_in
    my_x <- 'bait_condition'
    
    my_data <- filter(get(paste0("all_data_", celltype)), grepl(bait_in, bait)) %>%
      group_by(bait, condition) %>% 
      tidyr::unite(bait_condition, c(bait, condition), sep = "_")
  }
  
  # should the data be log transformed
  my_data[my_y] <- ifelse(!is.null(log_scale), log(my_data[my_y], base = log_scale), my_data[my_y])
  
  p <- ggplot(my_data, aes_string(x = my_x, y =my_y)) + 
    viz_function() + labs(title = my_title) + 
    ylab(my_y_label) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  print(p)
  
  my_data
}
