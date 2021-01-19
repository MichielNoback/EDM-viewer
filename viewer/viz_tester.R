source('../utils/visualization_functions.R')


# barplot over experiments based on summed counts for a given uniprot
vis_testing <- filter_data(celltype = 'dicty', uniprot_id = 'GST', experiment_bool = TRUE)
EDM_plot(my_data = vis_testing)

# barplot over baits based on summed counts for a given condition
vis_testing <- filter_data(celltype = 'dicty', condition_in = 'Normal', bait_bool = TRUE)
EDM_plot(my_data = vis_testing)

# barplot over condition based on summed counts for a given bait
vis_testing <- filter_data(celltype = 'dicty', bait_in = 'Galpha8', condition_bool = TRUE)
EDM_plot(my_data = vis_testing)

# boxplot over condition based on summed counts for a given bait
EDM_plot(my_data = vis_testing, viz_type = 'box')

# boxplot over condition based on summed counts for a given bait + log scaling
EDM_plot(my_data = vis_testing, viz_type = 'box', log_scale = 10)

# barplot over condition based on summed counted for a given bait
vis_testing <- filter_data(celltype = 'dicty', bait_in = 'Galpha8', normalized = TRUE, condition_bool = TRUE)
EDM_plot(my_data = vis_testing)

# boxplot counts over experiments for a given bait + log scaling + normalized
vis_testing <- filter_data(celltype = 'dicty', bait_in = 'Galpha8', normalized = FALSE, experiment = TRUE)
EDM_plot(my_data = vis_testing, viz_type = 'box', log_scale = 10)

