#build main tab for genes query
genesTab <- sidebarLayout(
    sidebarPanel(
        selectizeInput(inputId = "gene_id", 
                    label = "Select a gene by identifier", 
                    choices = c("", "tr|B4E2X5|B4E2X5_HUMAN", "tr|B1AH78|B1AH78_HUMAN")),
    ),

    # Show a plot of the generated distribution
    mainPanel(
        textOutput(outputId = "gene_name"),
        textOutput(outputId = "gene_info"),
        plotOutput("gene_conditions_bar")
    )
)



ui <- fluidPage(
    # Application title
    titlePanel("EDM data browser"),
    tabsetPanel(id = "main-tabs", 
                type = "tabs",
                tabPanel("Genes", genesTab),
                tabPanel("Conditions", "Conditions"))
)