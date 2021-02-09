# op basis van celltype -- condities laden
# dan datatable laten zien
# hier selectie in maken
# in tabs verschillende visualisaties
# gebruiker genen laten zoeken die voldoen aan criteria; bv z-score filtering of spectral-count en bait

# studenten: fair data viewer voor data collectie

library(DT)
source("R/global.R")
source("R/plotting_functions.R")


server <- function(input, output, session) { 
    #the gene being searched, as it is typed
    gene_search <- reactive({
        req(input$gene_search)
        message(paste0("[gene_search]: ", input$gene_search))
        input$gene_search
    })
    
    # the tibble with matching entries from the all_data dataset
    matching_genes <- reactive({
        message(paste0("looking for genes matching ", gene_search()))
        get_matching_genes_for_long_id_or_uniprot(gene_search())
    })
    
    # observes selected entries in the matching genes
    observe({
        req(input$suggested_genes_rows_selected)
        selection <- matching_genes()[input$suggested_genes_rows_selected, "uniprot"] %>% 
            pull(uniprot)
        
        message(paste0("[server] selection: ", paste(selection, collapse = ";")))
        selected_table <- get_gene_conditions_table(selection)

        output$gene_conditions_bar <- renderPlot({
            plot_bargraph_for_gene(selected_table)
        })

        # ## Fills the "gene_name" text box
        # output$gene_name <- renderText(gene_search())
        # output$gene_info <- renderText(get_gene_info_header(gene_search()))

        ## Create a datatable with values for the selected gene
        # output$gene_conditions_table <- renderDataTable( 
        #     ##? DT::renderDatatable en DT::datatable gebruiken
        #     datatable(get_gene_conditions_table(gene_search())) %>% 
        #         formatRound(columns = "z_score_spectral", digits = 2))

        # output$selected_genes <- DT::renderDataTable({
        #     DT::datatable(matching_genes()[input$suggested_genes_rows_selected, ])
        # })
    })
    
    ## shows the matching_genes tibble to select from
    output$suggested_genes <- DT::renderDataTable({
        DT::datatable(matching_genes(),
                      options = list(dom = 't'))
    })

    
    ## Populates the gene_select selectize; dependent on "celltypes_for_genes"
    # observe({
    #     req(input$celltypes_for_genes)
    #         updateSelectizeInput(session = session, 
    #                      inputId = "gene_select", 
    #                      choices = get_genes_vector(input$celltypes_for_genes), 
    #                      server = TRUE)
    # })

    ## Reactive depending on the "gene_select" input
    # gene_id <- reactive({
    #     req(input$gene_select)
    #     input$gene_select
    # })

    ##
    # output$gene_conditions_bar <- renderPlot(
    #     plot_bargraph_for_gene(get_gene_conditions_table(gene_id()))
    # )

}