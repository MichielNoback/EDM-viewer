source("../data_mining/data_processing.R")

library(ggplot2) 

# function to plot (bar or box), for the filtered data
EDM_plot <- function(my_data, viz_type='bar', log_scale=NULL){
  
  counts <- my_data$count_data
  
  title <- my_data$fig_labels$title
  x_label <- my_data$fig_labels$x_label
  y_label <- my_data$fig_labels$y_label
  my_y_label <- if_else(!is.null(log_scale), paste0("log",log_scale, "(", y_label, ")"), y_label)
  
  viz_type <- if_else(viz_type=="bar", "geom_col", "geom_boxplot")
  viz_function <- get(viz_type, envir=environment(ggplot))
  
  counts[y_label] <- ifelse(!is.null(log_scale), log(counts[y_label], base = log_scale), counts[y_label])
  
  p <- ggplot(counts, aes_string(x = x_label, y = y_label)) + 
    viz_function() + labs(title = title) + 
    ylab(my_y_label) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  print(p)
}
