source('../utils/visualization_functions.R')

plot_spectral_bait_condition_using_uniprot('P36415', 'dicty')
plot_spectral_bait_condition_using_uniprot('GST', 'dicty')
plot_spectral_bait_condition_using_uniprot('GST', 'neutro')
plot_spectral_bait_condition_using_uniprot('Q54YT4', 'dicty')
plot_spectral_using_condition("Normal", 'dicty')
plot_spectral_using_bait("Galpha8", 'dicty')

condition_boxplot_spectral_using_bait("Galpha8", 'dicty')
experiment_boxplot_spectral_using_bait("Galpha8", 'dicty')


plot_normalized_spectral_bait_condition_using_uniprot('P36415', 'dicty')
condition_boxplot_normalized_spectral_using_bait("Galpha8", 'dicty')
experiment_boxplot_normalized_spectral_using_bait("Galpha8", 'dicty')
