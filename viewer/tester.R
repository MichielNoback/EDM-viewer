
library(DT)
source("R/plotting_functions.R")


ui <- fluidPage(
    titlePanel("Hello Shiny!"),
    sidebarLayout(
        sidebarPanel(
            textInput(inputId = "gene_search",
                      label = "Find gene(s)",
                      ),
            dataTableOutput(outputId = "suggested_genes")

        ),
        mainPanel(
            dataTableOutput(outputId = "selected_genes"),
            plotOutput(outputId = "gene_conditions_bar")
        )
    )
)

server <- function(input, output) {
    #the gene being searched, as it is typed
    gene_search <- reactive({
        req(input$gene_search)
        input$gene_search
    })
    
    # the tibble with matching entries from the all_data dataset
    matching_genes <- reactive({
        get_matching_genes_for_long_id_or_uniprot(gene_search())
    })
    
    # observes selected entries in the matching genes
    observe({
        req(input$suggested_genes_rows_selected)
        selection <- matching_genes()[input$suggested_genes_rows_selected[1], "uniprot"]
        selection <- pull(selection, uniprot)
        print(paste0("[server] selection: ", selection))
        selected_table <- get_gene_conditions_table(selection)
        #print("selected table:")
        #print(selected_table)

        output$gene_conditions_bar <- renderPlot({
            plot_bargraph_for_gene(selected_table)
        })
        
        output$selected_genes <- DT::renderDataTable({
            DT::datatable(matching_genes()[input$suggested_genes_rows_selected, ])
        })
    })
    
    # shows the matching_genes tibble to select from
    output$suggested_genes <- DT::renderDataTable({
        DT::datatable(matching_genes(),
                      options = list(dom = 't'))
    })

}

shinyApp(ui, server)
