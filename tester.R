source("./global.R")

conditions <- get_conditions_vector("dicty")


all_data %>% 
    filter(Condition %in% conditions[sample.int(length(conditions), 3)]) %>%
    mutate(signal_log2 = log2(signal)) %>%
    ggplot(mapping = aes(x = Condition, y = signal_log2, color = Condition)) +
        geom_boxplot() +
        ylab("Log2 of the signal") +
        theme(legend.position = "none")



