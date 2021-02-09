#build top-level tab for genes query
genesTab <- sidebarLayout(
    sidebarPanel(
        # checkboxGroupInput(inputId = "celltypes_for_genes",
        #                 label = "Cell type",
        #                 choices = get_celltypes()),
        textInput(inputId = "gene_search",
                  label = "Find gene(s)"),
    
        dataTableOutput(outputId = "suggested_genes")
    ),

    mainPanel(
        # htmlOutput(outputId = "gene_name", container = tags$h2),
        # htmlOutput(outputId = "gene_info"),
        h2("Absolute count (left) and Z-score"),
        plotOutput("gene_conditions_bar"),
        # dataTableOutput("gene_conditions_table")
    )
)

ui <- function(request) {
    #ui code
}

# main UI
ui <- fluidPage(
    titlePanel("EDM data browser"),
    tabsetPanel(id = "main-tabs", 
                type = "tabs",
                tabPanel("Single Genes", genesTab),
                tabPanel("Filter Genes", "Filter genes on count, z-score, bait, condition"),
                tabPanel("Conditions", "Conditions in heatmap and PCA plot"))
)