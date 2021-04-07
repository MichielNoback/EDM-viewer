library(DT)

source("R/data_functions.R")
source("R/plotting_functions.R")


server <- function(input, output, session) {
        

    ###### Gene filter (gf_xxx) related -- TODO: modularize ######
    # gf_cell_type <- reactive({
    #     req(input$celltypes_for_genes)
    #     message(paste0("[celltypes_for_genes]", input$celltypes_for_genes))
    #     input$celltypes_for_genes
    # })
    # 
    # observe({
    #     message(paste0("[updating condition checkboxes]", gf_cell_type()))
    #     
    #     updateCheckboxGroupInput(session = session,
    #                              inputId = "source_organism_for_genes",
    #                              choices = get_source_organism_vector(gf_cell_type()),
    #                              inline = FALSE)
    #     
    #     updateCheckboxGroupInput(session = session,
    #                              inputId = "condition_filter_for_genes",
    #                              choices = get_conditions_vector(gf_cell_type()),
    #                              inline = TRUE)
    #     
    #             
    #     updateCheckboxGroupInput(session = session,
    #                              inputId = "bait_filter_for_genes",
    #                              choices = get_baits_vector(gf_cell_type()),
    #                              inline = TRUE)
    # })
    
    
    ###### Gene search related -- TODO: modularize ######
    
    # the gene being searched, as it is typed
    # returns the tibble with matching entries from the all_data dataset
    suggested_genes <- reactive({
        req(input$gene_search)
        message(paste0("[suggested_genes]: ", input$gene_search))
        message(paste0("[suggested_genes]: ", input$gene_search), 
                "; filtering ", nrow(all_data), " data rows")
        get_matching_genes_for_long_id_or_uniprot(input$gene_search)
    })
    
    ## shows the suggested_genes tibble to select from
    output$suggested_genes <- DT::renderDataTable({
        DT::datatable(suggested_genes(),
                      options = list(dom = 'tp'))
    })
    
    ## a reactive that holds the selected rows from suggested_genes
    selected_rows_from_suggested <- reactive({
        req(input$suggested_genes_rows_selected)
        suggested_genes()[input$suggested_genes_rows_selected, ]
    })
    
    ## a reactive that holds the selected uniprot_ids from the suggested_genes
    ## simply pulls the uniprot column
    selected_uniprot_ids_from_suggested <- reactive({
        message(paste0("[selected_uniprot_ids_from_suggested]: ",
                       paste(selected_rows_from_suggested()$uniprot, 
                             collapse = ";")))
        selected_rows_from_suggested()$uniprot
    })
    
    #STORE SELECTED GENES AS USER DATA
    user_data <- reactiveValues(uniprot_ids = character(0))
    
    #observes changes in selected genes from suggested genes
    observeEvent(selected_uniprot_ids_from_suggested(), {
        user_data$uniprot_ids <- 
            unique(c(user_data$uniprot_ids, 
                     selected_uniprot_ids_from_suggested()))
        message(paste0("[observeEvent:selected_uniprot_...] new user selection: \n\t",
                        paste(user_data$uniprot_ids,
                              collapse = ";")))
    })
    
    ## creates a datatable view of the user-selected genes
    output$selected_genes_for_barplot <- DT::renderDataTable({
        selected_genes <- get_matching_genes_for_uniprot(user_data$uniprot_ids)
            DT::datatable(selected_genes,
                      options = list(dom = 't'))
    })

    ## creates the barplot of selected genes
    output$gene_conditions_bar <- renderPlot({
        req(user_data$uniprot_ids)
        message(paste0("[output$gene_conditions_bar..] plotting from: \n\t",
                        paste(user_data$uniprot_id,
                              collapse = ";")))
        selected_genes_full_table <- 
            get_gene_conditions_table(user_data$uniprot_ids)
        plot_bargraph_for_gene(selected_genes_full_table)
    })
    
}