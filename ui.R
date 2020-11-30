
ui <- fluidPage(

    # Application title
    titlePanel("EDM data browser"),

    # Sidebar with a slider input for number of bins 
    sidebarLayout(
        sidebarPanel(
            checkboxGroupInput(inputId = "cell_types",
                        label = "Cell type",
                        choices = get_celltypes()),
            uiOutput("condition_select"),
            uiOutput("gene_select"),
        ),

        # Show a plot of the generated distribution
        mainPanel(
            fluidRow(column(10, verbatimTextOutput("selected_condition"))),
            fluidRow(column(10, verbatimTextOutput("selected_gene"))),
            plotOutput("conditions_box")
        )
    )
)
