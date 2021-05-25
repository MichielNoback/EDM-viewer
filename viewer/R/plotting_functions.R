suppressMessages(library(ggplot2))
suppressMessages(library(gridExtra))
suppressMessages(library(dplyr))
suppressMessages(library(tidyr))

#' Create bar graphs for spectral count and its z-score
#' 
#' Creates a bar graph of both the spectral_count column and the z_score_spectral
#' column for a single gene over all conditions, show side by side. The z-score
#' bar graph also displays a "95% confidence interval" (z-score = 2).
#' @param gene_data a tibble containing all measurements for a single gene
#' @return a ggplot2 object
#' @examples 
#' print(plot_bargraph_for_gene(get_gene_data_for_uniprots("sp|P0A6M8|EFG_ECOLI")))
#' print(plot_bargraph_for_gene(get_gene_data_for_uniprots("sp|P02887|DIS1B_DICDI")))
plot_bargraph_for_genes <- function(gene_data, 
                                   experiments = "all",
                                   normalized_input_values = TRUE,
                                   faceted = FALSE) {
    if(is.null(gene_data) || nrow(gene_data) == 0) {
        message("returning NULL")
        return(NULL)
    }

    ##process input data
    tmp <- gene_data %>% 
        unite(cond_bait, condition, bait, sep = " / ")
    
    if(! experiments[1] == "all") {
        #print("filtering")
        tmp <- tmp %>% filter(cond_bait %in% experiments)
    }

    if (faceted) {
        plot_faceted_bargraph_for_genes(tmp, normalized_input_values)
    } else {
        plot_side_by_side_barplot_for_genes(tmp, normalized_input_values)
    }
    

    ## arrange and return
    # plot_layout <- rbind(c(1, 1, 1, 2, 2, 2, 2))
    # 
    # return(grid.arrange(count_plt, 
    #                    zscore_plt, 
    #                    layout_matrix = plot_layout))
}

plot_side_by_side_barplot_for_genes <- function(gene_data, normalized_input_values) {
    if (normalized_input_values) {
        ##create barplot of z-score data
        caption <- paste(strwrap("Shaded area indicates z-score=2 boundary of 95%", 40),
                     collapse = "\n")
        zscore_plt <- ggplot(data = gene_data,
                            mapping = aes(x = cond_bait, 
                                      y = z_score_spectral, 
                                      fill = name)) +
        annotate("rect", xmin = -Inf, xmax = Inf, ymin = -2.0, ymax = 2.0, 
                 alpha=0.2, fill="red") +
        
        geom_bar(stat = "identity", 
                 position = position_dodge2(width = 0.9, preserve = "single")) +
        ylab("z-score") +
        theme_minimal() +
        theme(axis.text.x = element_text(angle = 60, hjust = 1),
              axis.title.x = element_blank()) +
        annotate(geom = "text", 
                 x = -Inf, 
                 y = range(gene_data$z_score_spectral)[2], 
                 label = caption, 
                 hjust = 0, 
                 vjust = 1, 
                 size = 3) +
        labs(fill = "short ID")
        return(zscore_plt)
    } else {
        ##create barplot of count data
        count_plt <- ggplot(data = gene_data,
                            mapping = aes(x = cond_bait, 
                                          y = spectral_count, 
                                          fill = name)) +
            geom_bar(stat = "identity", 
                     position = position_dodge2(width = 0.9, preserve = "single")) +
            ylab("spectral count") +
            theme_minimal() +
            theme(axis.text.x = element_text(angle = 60, hjust = 1),
                  axis.title.x = element_blank())   
        return(count_plt)
    }
}

plot_faceted_bargraph_for_genes <- function(gene_data, normalized_input_values) {
    if (normalized_input_values) {
        ##create barplot of z-score data
        caption <- paste(strwrap("Shaded area indicates z-score=2 boundary of 95%", 40),
                     collapse = "\n")
        zscore_plt <- ggplot(data = gene_data,
                            mapping = aes(x = cond_bait, 
                                      y = z_score_spectral, 
                                      fill = name)) +
        annotate("rect", xmin = -Inf, xmax = Inf, ymin = -2.0, ymax = 2.0, 
                 alpha=0.2, fill="red") +
        
        geom_bar(stat = "identity") +
        annotate(geom = "text", 
                 x = -Inf, 
                 y = range(gene_data$z_score_spectral)[2], 
                 label = caption, 
                 hjust = 0, 
                 vjust = 1, 
                 size = 3) +
        facet_wrap(. ~ name) +
        ylab("z-score") +
        theme_minimal() +
        theme(axis.text.x = element_text(angle = 60, hjust = 1),
              axis.title.x = element_blank(),
              legend.position = "none") + 
        labs(fill = "short ID")
        return(zscore_plt)
    } else {
        ##create barplot of count data
        # plot_margins <- c(0.5, 0.5, 0.5, 1)
        count_plt <- ggplot(data = gene_data,
                            mapping = aes(x = cond_bait, 
                                          y = spectral_count, 
                                          fill = name)) +
            geom_bar(stat = "identity") +
            facet_wrap(. ~ name) +
            ylab("spectral count") +
            theme_minimal() +
            theme(axis.text.x = element_text(angle = 60, hjust = 1),
                  axis.title.x = element_blank(),
                  legend.position = "none")
        return(count_plt)
    }
}


#q <- c("sp|P0A6F5|CH60_ECOLI","sp|P08622|DNAJ_ECOLI","sp|Q869S8|PSMG2_DICDI", "sp|P0A6M8|EFG_ECOLI")
#test <- get_gene_conditions_table(q)
#print(plot_bargraph_for_gene(get_gene_conditions_table(q[1:3])))


