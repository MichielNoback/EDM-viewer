library(ggplot2)
library(UniprotR)
library(ggpubr)
library(gridExtra)
library(topGO)


#' Plot (bar or box) using filtered data
#' 
#' @param my_data filtered data from the filter_data() function
#' @param viz_type one of bar (default) or box
#' @param log_scale log scale to be used
#' @examples 
#' viz_testing <- data_processing::filter_data(celltype = 'dicty', bait_in = 'Galpha8', experiment_bool = TRUE)
#' EDM_plot(my_data = viz_testing)
#' EDM_plot(EDM_plot(my_data = vis_testing, viz_type = 'box')
#' EDM_plot(my_data = vis_testing, viz_type = 'box', log_scale = 2)
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


#' Show overlapping proteins in a Venn diagram
#'
#' @param z_treshold z_value cutoff to be used
#' @param bait1, bait2 and bait3 specifies which baits should be used to calculate overlap
#' @example
#' show_overlapping_peptides(z_threshold = 1.5, bait1 = 'Galpha2', bait2 = 'Galpha4', bait3 = 'Galpha8')
show_overlapping_peptides <- function(my_data = all_data_dicty_roco4_corrected, z_threshold = 1.5, bait1 = 'Galpha2', bait2 = 'Galpha4', bait3 = 'Galpha8'){
  
  bait <- my_data %>% 
    ungroup() %>%
    dplyr::select(uniprot, bait, z_score_spectral) %>% 
    filter(z_score_spectral > z_threshold | z_score_spectral < -z_threshold)  %>%
    dplyr::select(!z_score_spectral)
  
  list_of_uniprot_per_bait <- split(bait$uniprot, bait$bait)
  
  toy_bait = Venn(list_of_uniprot_per_bait)
  ggvenn(toy_bait, slice = c(bait1, bait2, bait3))
}


#' Show overlapping proteins in a Venn diagram, correct spectral counts using ROCO4 as background
#'
#' @param my_dicty_data all_data_dicty_roco4_corrected dataframe from the data_processing::get_data() function
#' @param bait1, bait2 and bait3 specifies which baits should be used to calculate overlap
#' @param scaled boolean (default TRUE), use z_score cutoff to select proteins for correction. When set to FALSE all proteins that have an interaction with ROCO4 will be used.
#' @param z_treshold z_value cutoff to be used
#' @return venn_overlapping_proteins list of lists with protein identifiers for each part of the Venn diagram
#' @example
#' all_data_dict <- data_processing::get_data()
#' show_overlapping_peptides_ROCO4_Background_corrected(my_dicty_data = all_data_dicty_roco4_corrected, bait1 = 'Galpha2', bait2 = 'Galpha4', bait3 = 'Galpha8', z_threshold = 1, scaled = TRUE)
show_overlapping_peptides_ROCO4_Background_corrected <- function(my_dicty_data, bait1='Galpha2', bait2='Galpha4', bait3='Galpha8', scaled=TRUE, z_threshold=1, show_venn=FALSE){
  # TODO: also include condition?, in over-representation and GSEA 2 conditions cannot be separeted
  
  my_dicty_data$bait <- as.character(my_dicty_data$bait)
  
  # if you don't scale, all uniprot id's will be selected, resulting in no significant id's left over
  roco4_uniprot <- my_dicty_data %>% 
    filter(bait=='Roco4') %>% 
    dplyr::select(uniprot, z_score_spectral)
  
  if (scaled) {
    roco4_uniprot <- roco4_uniprot %>% 
      filter(z_score_spectral > z_threshold | z_score_spectral < -z_threshold)
  }
  
  bait <- my_dicty_data %>% 
    ungroup() %>%
    dplyr::select(uniprot, bait, z_score_spectral, spectral_count) %>% 
    filter(z_score_spectral > z_threshold | z_score_spectral < -z_threshold)  %>%
    #select(!z_score_spectral) %>%
    filter(!(uniprot %in% roco4_uniprot$uniprot))
  
  list_of_uniprot_per_bait <- split(bait$uniprot, bait$bait)
  list_of_uniprot_per_bait <- list_of_uniprot_per_bait[c(bait1, bait2, bait3)]
  
  toy_bait = Venn(list_of_uniprot_per_bait)
  if (show_venn){
    print(ggvenn(toy_bait, slice = c(bait1, bait2, bait3)))
  }
  
  e <- enrichment_test(toy_bait, bait1, bait2)
  print("Significance of overlap")
  print(e$Significance)
  
  venn_overlapping_proteins <- list()
  venn_overlapping_proteins$shared_among_all <- overlap(toy_bait, slice = c(bait1, bait2, bait3))
  
  # counts uniq for the first baits vs the other
  # i.e. set1 = a,b,c and set2 = b,c,d,e -> set1 .. set2 = a
  venn_overlapping_proteins$uniq_pairs <- discern_pairs(toy_bait)
  
  #counts shared among baits
  venn_overlapping_proteins$overlap_pairs <- overlap_pairs(toy_bait)
  
  venn_overlapping_proteins$bait1_uniq <- discern(toy_bait, slice1 = bait1, slice2 = c(bait2, bait3))
  venn_overlapping_proteins$bait2_uniq <- discern(toy_bait, slice1 = bait2, slice2 = c(bait1, bait3))
  venn_overlapping_proteins$bait3_uniq <- discern(toy_bait, slice1 = bait3, slice2 = c(bait1, bait2))
  names(venn_overlapping_proteins)[4] <- bait1   
  names(venn_overlapping_proteins)[5] <- bait2
  names(venn_overlapping_proteins)[6] <- bait3
  
  venn_overlapping_proteins$bait1_2_uniq <- discern(toy_bait, slice1 = c(bait1, bait2), slice2 = bait3)
  venn_overlapping_proteins$bait2_3_uniq <- discern(toy_bait, slice1 = c(bait2, bait3), slice2 = bait1)
  venn_overlapping_proteins$bait1_3_uniq <- discern(toy_bait, slice1 = c(bait1, bait3), slice2 = bait2)
  names(venn_overlapping_proteins)[7] <- paste(bait1, bait2, sep = "_")
  names(venn_overlapping_proteins)[8] <- paste(bait2, bait3, sep = "_")
  names(venn_overlapping_proteins)[9] <- paste(bait3, bait1, sep = "_")
  
  #return (list_of_uniprot_per_bait) 
  return (venn_overlapping_proteins) 
  #return (bait)
  #return (e)
}


GO_bar <- function(DF, my_fill, my_xlab){
  DF <- DF[1:10, ]
  DF <- na.omit(DF)
  
  GO_Plot <- ggplot(data = DF, 
                           aes(
                             x = reorder(Goterm, Count), 
                             y = Count)) + 
                          geom_bar(
                            stat = "identity", 
                            fill = my_fill) + 
                          xlab(my_xlab) + 
                          ylab("Protein count") + 
                          theme_bw() + 
                          theme(
                            text = element_text(
                              size = 12, 
                              face = "bold", 
                              colour = "black"), 
                            axis.text.x = element_text(vjust = 2)) +
                            coord_flip()
  GO_Plot
}

combine_GO <- function(accessions) 
{
  #get them from uniprot site, slow!
  #GeneOntologyObj <- GetProteinGOInfo(accessions) 
  
  #get GO annotation from saved image, allready downloaded from uniprot
  load("dicty_annotated.Rdata")
  
  #filter for the accession we have in the overlap
  GeneOntologyObj <- dicty_annotated %>% 
    filter(row.names(dicty_annotated) %in% accessions)
  print("Done filtering")
  
  BiologicalDF <- Goparse(GeneOntologyObj, 4)
  MolecularDF <- Goparse(GeneOntologyObj, 5)
  CellularDF <- Goparse(GeneOntologyObj, 6)
  print("parsed GO data")
  
  BiologicalPlot <- GO_bar(DF = BiologicalDF, my_fill = "darkred", my_xlab = "Biological Process")
  MolecularPlot <- GO_bar(DF = MolecularDF, my_fill = "darkgreen", my_xlab = "Molecular Function")
  CellularPlot <- GO_bar(DF = CellularDF, my_fill = "darkblue", my_xlab = "Cellular component")
  print("Barplots created")

  Goplots <- ggarrange(BiologicalPlot, MolecularPlot, CellularPlot, nrow = 3, ncol = 1, align = "hv")
  plot(Goplots)
}



#' Do GO term analysis, plots counts/frequency of observed GO terms on the 3 levels:
#' Biological process, Molecular function, Cellular component
#' 
#' @param accessions character vector with uniprot accessions numbers
#' @return TaxaObj annotation list containing:
#'  [1] "Entry.name"                   "Gene.names"                   "Gene.names...primary.."      
#' [4] "Gene.names...synonym.."       "Gene.names...ordered.locus.." "Gene.names...ORF.."          
#' [7] "Organism"                     "Organism.ID"                  "Protein.names"               
#' [10] "Proteomes"                    "Taxonomic.lineage..ALL."      "Virus.hosts"  
#' @example
#' taxobj <- go_term_analysis(accessions = overlap_accessions)
go_term_analysis <- function(accessions){

  #print(overlapping_accessions)
  # convertID, see https://www.uniprot.org/help/api_idmapping
  #ConvertID(Accessions, ID_from = "ACC+ID", ID_to = "GENENAME")
  
  # taxa info
  #TaxaObj <- GetNamesTaxa(accessions) 
  #PlotChromosomeInfo(TaxaObj)
  #PlotGenesNetwork(TaxaObj)
  
  
  # Get GO terms usign Uniprot accessions
  GeneOntologyObj <- GetProteinGOInfo(accessions) 
  #PlotGOBiological(GeneOntologyObj, Top = 10) 
  #Plot.GOSubCellular(GeneOntologyObj) 
  
  # Use this one
  PlotGoInfo(GeneOntologyObj)
  
  # enrichment analysis
  #Enrichment.MF(overlapping_accessions, significant = FALSE)
  
  #GetproteinNetwork(Accessions,  "/Users/ronald/Desktop/") 
  #return (TaxaObj)
  #return (GeneOntologyObj)
}


#' Do GO term analysis using Fisher exact testing, show ordered table with enriched terms, plot GO Ontology tree
#' @param accession character vector with short uniprot ids
#' @param ontology one of: "MF", "BP" or "CC", which is short for: Molecular function, Biological process and Cellular component
#' @param tapleTopNodes int number of GO terms to show in reported (ordered) table with Fisher p-values
#' @param figNodeSize int cutoff value for each node, Nodes with number of found significant GO terms above threshold are shown
#' @param figSigNodes int the number of top scoreing GO terms which should be shown in the plot, shown as rectangles
#' 
#' @example 
#' get significant overlapping peptides
#' overlapping_ids <- show_overlapping_peptides_ROCO4_Background_corrected(my_dicty_data = all_data_dicty, bait1 = 'Galpha2', bait2 = 'Galpha4', bait3 = 'Galpha8', z_threshold = 1, scaled = TRUE)
#' convert long uniprot id to short id, selecting part of the overlap
#' accessions <- get_uniprot(overlapping_ids$Galpha2_Galpha4)
#' do the go enrichment analysis
#' go_enrichment_analysis(accessions = accessions, ontology = "MF", tableTopNodes = 10, figNodeSize = 5, figSigNodes = 5)
go_enrichment_analysis <- function(accessions = accessions, ontology = "MF", pvalCutOff = NULL, tableTopNodes = 10, figNodeSize = 5, figSigNodes = 5, save_file = FALSE, bait_condition = "sampleFile", show_fig = TRUE){
  # ontology: MF, BP or CC
  load(file = "UniverseGo.RData")
  
  #create annotation object uniprot -> GO
  geneID2GO <- as.list(universeGO$Gene.ontology.IDs)
  geneID2GO <- str_split(geneID2GO, pattern = "; ")
  names(geneID2GO) <- rownames(universeGO)
  
  #get universe
  # TODO: Do we need to limit this to the uniprot id's in the Venn diagram, now taking id's on organism level
  geneNames <- unique(all_data_dicty_roco4_corrected$uniprot_short)
  #annotate which genes we are interested in
  geneList <- factor(as.integer(geneNames %in% accessions))
  names(geneList) <- geneNames

  #create topGO object
  GOdata <- new("topGOdata", 
                ontology = ontology,
                allGenes = geneList,
                annot = annFUN.gene2GO, 
                nodeSize = figNodeSize,
                gene2GO = geneID2GO)
  
  # Run fisher test
  resultFisher <- runTest(GOdata, algorithm = "classic", statistic = "fisher")
  
  # Build top scoring Go term table
  allRes <- GenTable(GOdata, classicFisher = resultFisher,
                       orderBy = "classicFisher", ranksOf = "classicFisher", topNodes = tableTopNodes)
  
  if (save_file){
    print("Saving file")
    printGraph(GOdata, resultFisher, firstSigNodes = figSigNodes, useInfo = 'all', fn.prefix = paste(bait_condition, ontology, sep = "_"))
  }
  # plot the GO term graph
  if (show_fig){
    showSigOfNodes(GOdata, score(resultFisher), firstSigNodes = figSigNodes, useInfo = 'all')
  }
  return(allRes)
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
show_igraph_network_plot <- function(my_data = all_data_dicty, filter_cutof = 30){
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
