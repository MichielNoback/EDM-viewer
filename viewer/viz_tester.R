source('../utils/visualization_functions.R')


# test uniprot + non normal
vis_testing <- filter_data_and_plot(celltype = 'dicty', uniprot_id = 'GST', normalized = FALSE)
# test uniprot + normal
vis_testing <- filter_data_and_plot(celltype = 'dicty', uniprot_id = 'GST', normalized = TRUE)

# test condition + normal
vis_testing <- filter_data_and_plot(celltype = 'dicty', condition = 'Normal', normalized = TRUE)

# test condition + non normal
vis_testing <- filter_data_and_plot(celltype = 'dicty', condition = 'Normal', normalized = FALSE)

# test bait + normal
vis_testing <- filter_data_and_plot(celltype = 'dicty', bait_in = 'Galpha8', normalized = TRUE)

# test bait + non normal
vis_testing <- filter_data_and_plot(celltype = 'dicty', bait_in = 'Galpha8', normalized = FALSE)

# test condition + boxplot  
vis_testing <- filter_data_and_plot(celltype = 'dicty', bait_in = 'Galpha8', normalized = FALSE, viz_type = 'box')

# test condition + boxplot + log scaling
vis_testing <- filter_data_and_plot(celltype = 'dicty', bait = 'Galpha8', normalized = FALSE, viz_type = 'box', log_scale = 10)

# test condition + boxplot + log scaling + normalized
vis_testing <- filter_data_and_plot(celltype = 'dicty', bait = 'Galpha8', normalized = TRUE, viz_type = 'box', log_scale = 10)

# test condition + boxplot + log scaling + normalized
vis_testing <- filter_data_and_plot(celltype = 'dicty', experiment = TRUE, bait = 'Galpha8', normalized = FALSE, viz_type = 'box', log_scale = 10)

