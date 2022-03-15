suppressMessages(library(ggplot2))
suppressMessages(library(gridExtra))
suppressMessages(library(dplyr))
suppressMessages(library(tidyr))
library(igraph)
library(visNetwork)

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
        geom_hline(yintercept = 2, linetype = "dashed", color = "red") +
        geom_hline(yintercept = -2, linetype = "dashed", color = "red") +
        # annotate("rect", xmin = -Inf, xmax = Inf, ymin = -2.0, ymax = 2.0, 
        #          alpha=0.2, fill="red") +
        
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
        geom_hline(yintercept = 2, linetype = "dashed", color = "red") +
        geom_hline(yintercept = -2, linetype = "dashed", color = "red") +

        # annotate("rect", xmin = -Inf, xmax = Inf, ymin = -2.0, ymax = 2.0, 
        #          alpha=0.2, fill="red") +
        
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

#' Create an interactive igraph network plot. 
#' Interacting proteins are represented as circles and baits are depicted as triangles
#' Edges between nodes are scaled by spectral count (z-score)
#' Packages required:
#' library(igraph)
#' library(visNetwork)
#' @param my_data all_dicty or subset (i.e. background corrected version) of this dataframe
#' @param filter_cutof used cutoff value for filtering on minimal spectral count, default 30
#' @example show_igraph_network_plot(all_data_dicty, filter_cutof = 20)
#' show_igraph_network_plot(my_data = all_data_dicty)
show_igraph_network_plot <- function(
    my_data = all_data_dicty, 
    filter_cutof = 30){
  my_data <- my_data %>%
    group_by(bait, condition) %>% 
    tidyr::unite(bait_condition, c(bait, condition), sep = "_")
  
  # split the uniprot column in the indivdual parts
  my_data_splitted_uniprot <- my_data %>% 
    separate(uniprot, c("a", "b", "c", "protid", "genid", "source"), extra = "merge", fill = "left", sep = '[_|]')
  
  # remove "garbage" columns and keep the protid, geneid and source columns
  my_data_splitted_uniprot <- my_data_splitted_uniprot %>% 
    dplyr::select(-one_of(c("a", "b", "c")))
  
  # resulting in 17k+ interactions, hairball!, need to filter
  my_data_splitted_uniprot_threshold_filtered <- my_data_splitted_uniprot %>% 
    filter(spectral_count > 30)
  
  nodes <- as.data.frame(unique(my_data_splitted_uniprot_threshold_filtered$bait_condition))
  nodes$shape = 'triangle'
  colnames(nodes) <- c('id', 'shape')
  
  # create the nodes
  nodes <- nodes %>% 
    add_row(id = unique(my_data_splitted_uniprot_threshold_filtered$genid), shape = 'circle')
  nodes$label <- nodes$id
  
  # define edges
  my_edges <- data.frame(my_data_splitted_uniprot_threshold_filtered$bait_condition, 
                         my_data_splitted_uniprot_threshold_filtered$genid, 
                         my_data_splitted_uniprot_threshold_filtered$z_score_spectral)
  colnames(my_edges) <- c('from', 'to', 'value')
  
  # create the graph
  graph <- graph_from_data_frame(my_edges, directed = FALSE)
  
  # show the network
  visNetwork(nodes = nodes, edges = my_edges) %>%
    visOptions(highlightNearest = list(enabled = TRUE, 
                                       algorithm = "hierarchical", 
                                       degree = list(from = 1, to = 1)), 
                                       nodesIdSelection = TRUE) %>%
                                          visPhysics(stabilization = FALSE) %>%
                                          visLayout(improvedLayout = TRUE)
}

#q <- c("sp|P0A6F5|CH60_ECOLI","sp|P08622|DNAJ_ECOLI","sp|Q869S8|PSMG2_DICDI", "sp|P0A6M8|EFG_ECOLI")
#test <- get_gene_conditions_table(q)
#print(plot_bargraph_for_gene(get_gene_conditions_table(q[1:3])))


