# build top-level tab for genes query
geneSearchTab <- sidebarLayout(
    sidebarPanel(
        textInput(inputId = "gene_search",
                  label = "Find gene(s)"),
    
        dataTableOutput(outputId = "suggested_genes")
    ),

    mainPanel(
        dataTableOutput(outputId = "selected_genes_for_barplot"),
        h2("Spectral count (left) and Z-score normalized"),
        plotOutput("gene_conditions_bar"),
    )
)


# build top-level tab for gene filtering

spectral_count_range <- c(0, 7500)
z_score_limit <- 50
mw_range <- c(1, 750)

geneFilterTab <- sidebarLayout(
    sidebarPanel(
        checkboxGroupInput(inputId = "celltypes_for_genes",
                        label = "Cell type",
                        choices = get_celltypes()),
        
        checkboxGroupInput(inputId = "source_organism_for_genes",
                           label = "Source organism for captured proteins",
                           inline = TRUE,
                           choices = NULL),
        
        checkboxGroupInput(inputId = "condition_filter_for_genes",
                           label = "Conditions",
                           inline = TRUE,
                           choices = NULL),

        checkboxGroupInput(inputId = "bait_filter_for_genes",
                           label = "Baits",
                           inline = TRUE,
                           choices = NULL),        

        ## TODO max should be set dynamically
        sliderInput(inputId = "spectral_count_filter", 
                    label = "Spectral count",
                    value = spectral_count_range,
                    min = spectral_count_range[1], max = spectral_count_range[2]),
        
        ## TODO min and max should be set dynamically
        sliderInput(inputId = "spectral_z_score_filter", 
                    label = "Z-score spectral count",
                    value = c(-z_score_limit, z_score_limit),
                    min = -z_score_limit, max = z_score_limit), 
        
        ## TODO min and max should be set dynamically
        sliderInput(inputId = "mw_filter", 
                    label = "Molecular Weight (kDa)",
                    value = mw_range,
                    min = mw_range[1], max = mw_range[2])
    ),

    mainPanel(
        h2("Matching genes"),
        dataTableOutput(outputId = "matching_genes")
    )
)
# ui <- function(request) {
#     #ui code
# }

# main UI
ui <- fluidPage(
    titlePanel("EDM data browser"),
    tabsetPanel(id = "main-tabs", 
                type = "tabs",
                tabPanel("Select Genes", geneSearchTab),
                tabPanel("Filter Genes", geneFilterTab),
                tabPanel("Investigate Conditions", "Conditions in interactive Venn diagram, heatmap and PCA plot, Ontology"))
)