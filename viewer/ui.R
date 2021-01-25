#build main tab for genes query
genesTab <- sidebarLayout(
    sidebarPanel(
        checkboxGroupInput(inputId = "celltypes_for_genes",
                        label = "Cell type",
                        choices = get_celltypes()),      
        selectizeInput("gene_select",
                       label = "Select a gene by name or identifier",
                       choices = NULL),
    ),

    mainPanel(
        textOutput(outputId = "gene_name"),
        textOutput(outputId = "gene_info"),
        plotOutput("gene_conditions_bar")
    )
)



ui <- fluidPage(
    titlePanel("EDM data browser"),
    tabsetPanel(id = "main-tabs", 
                type = "tabs",
                tabPanel("Genes", genesTab),
                tabPanel("Conditions", "Conditions"))
)