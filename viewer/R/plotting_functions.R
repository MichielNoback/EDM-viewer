library(ggplot2)
library(gridExtra)
library(dplyr)
library(tidyr)

#' Create bar graphs for spectral count and its z-score
#' 
#' Creates a bar graph of both the spectral_count column and the z_score_spectral
#' column for a single gene over all conditions, show side by side. The z-score
#' bar graph also displays a "95% confidence interval" (z-score = 2).
#' @param gene_data a tibble containing all measurements for a single gene
#' @return a ggplot2 object
#' @examples 
#' print(plot_bargraph_for_gene(get_gene_conditions_table("sp|P0A6M8|EFG_ECOLI")))
#' print(plot_bargraph_for_gene(get_gene_conditions_table("sp|P02887|DIS1B_DICDI")))
plot_bargraph_for_gene <- function(gene_data) {
    if(is.null(gene_data) || nrow(gene_data) == 0) {
        message("returning NULL")
        return(NULL)
    }
    #print(gene_data)
    ##define local vars
    plot_margins <- c(0.5, 0.5, 0.5, 1)
    caption <- paste(strwrap("Shaded area indicates z-score=2 boundary of 95%", 40),
                     collapse = "\n")
    
    ##process input data
    tmp <- gene_data %>% 
        #arrange(condition, bait) %>%
        unite(cond_bait, condition, bait, sep = " / ")
    
    ##create barplot of count data
    count_plt <- ggplot(data = tmp,
                        mapping = aes(x = cond_bait, 
                                      y = spectral_count, 
                                      fill = name)) +
        geom_bar(stat = "identity", 
                 position = position_dodge2(width = 0.9, preserve = "single")) +
        ylab("spectral count") +
        theme_minimal() +
        theme(axis.text.x = element_text(angle = 60, hjust = 1),
              axis.title.x = element_blank(),
              legend.position = "none",
              plot.margin = unit(plot_margins, "cm"))

    ##create barplot of z-score data
    zscore_plt <- ggplot(data = tmp,
                        mapping = aes(x = cond_bait, 
                                      y = z_score_spectral, 
                                      fill = name)) +
        geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = -2.0, ymax = 2.0),
            alpha = 0.01,
            fill = "red") +
        geom_bar(stat = "identity", 
                 position = position_dodge2(width = 0.9, preserve = "single")) +
        ylab("z-score") +
        theme_minimal() +
        theme(axis.text.x = element_text(angle = 60, hjust = 1),
              axis.title.x = element_blank(),
              plot.margin = unit(plot_margins, "cm")) +
        annotate(geom = "text", 
                 x = -Inf, 
                 y = range(tmp$z_score_spectral)[2], 
                 label = caption, 
                 hjust = 0, 
                 vjust = 1, 
                 size = 3) +
        labs(fill = "short ID")
    
    ## arrange and return
    plot_layout <- rbind(c(1, 1, 1, 2, 2, 2, 2))
    
    return(grid.arrange(count_plt, 
                       zscore_plt, 
                       layout_matrix = plot_layout))
}

#test <- get_gene_conditions_table(c("sp|P15064|RASG_DICDI", "sp|P32253|RASC_DICDI"))
#plot_bargraph_for_gene(test)


