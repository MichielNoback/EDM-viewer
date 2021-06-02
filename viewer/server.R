suppressMessages(library(DT))

## Tab gene selection: ID search en conditions
## Tab voor visualisatie met subtabs voor bar/heat etx

source("R/data_functions.R")
source("R/plotting_functions.R")
source("R/data_loading.R")

## DATA LOADING ##
all_data <<- prepare_data()
# create readable legend labels
all_data <<- create_name_column(all_data)
# create source organism column
#create_source_organism_column()

message(paste0("[global.R]: loaded ", nrow(all_data), " data rows"))



server <- function(input, output, session) {

    #STORE SELECTED GENES AS USER DATA (both uniprot IDs and tibble)
    user_data <- reactiveValues(
        ## used for gene selection
        uniprot_ids = character(0),
        ## used for gene filtering; [0, ] will select 
        ## an empty tibble with matching column names
        selected_genes = all_data[0,] 
    )
    
    ###### Gene search related -- TODO: modularize ######
    
    # the gene being searched, as it is typed
    # returns the tibble with matching entries from the all_data dataset
    suggested_genes <- reactive({
        #req(input$gene_search)
        message(paste0("[suggested_genes]: ", input$gene_search), 
                "; filtering ", nrow(all_data), " data rows")
        ##resets the row selection when a new search is being performed
        selected_rows_indices(character(0))
        get_matching_genes_for_long_id_or_uniprot(
            input$gene_search, 
            input$only_dicdi_genes)
    })

    ## holds the selected row indices of the suggestion list
    selected_rows_indices <- reactiveVal(
        character(0)
    )
    
    ## observes the event that the the row selection changes
    observeEvent(input$suggested_genes_rows_selected, {
        message(paste0("[observe suggested_genes selection]: selected=",
                       input$suggested_genes_rows_selected), 
                " \n\tcurrent selected_rows_indices()=", 
                selected_rows_indices())
        selected_rows_indices(input$suggested_genes_rows_selected)
    })
    
    ## shows the suggested_genes tibble to select from
    output$suggested_genes <- DT::renderDataTable({
            DT::datatable(suggested_genes(),
                          options = list(dom = 'tp'))
    })

    ## a reactive that holds the selected rows from suggested_genes
    selected_rows_from_suggested <- reactive({
        req(selected_rows_indices())#input$suggested_genes_rows_selected)
        message(paste0("[selected_rows_from_suggested] new rows selected: ", 
                       selected_rows_indices()))
        suggested_genes()[selected_rows_indices(), ]#input$suggested_genes_rows_selected, ]
    })
    
    ## a reactive that holds the selected uniprot_ids from the suggested_genes
    ## simply pulls the uniprot column
    selected_uniprot_ids_from_suggested <- reactive({
        message(paste0("[selected_uniprot_ids_from_suggested]: ",
                       paste(selected_rows_from_suggested()$uniprot, 
                             collapse = ";")))
        selected_rows_from_suggested()$uniprot
    })
    
    #observes changes in selected genes from suggested genes
    observeEvent(selected_uniprot_ids_from_suggested(), {
        new_selection <- unique(c(user_data$uniprot_ids, 
                     selected_uniprot_ids_from_suggested()))
        user_data$uniprot_ids <- new_selection
            
        message(paste0(
            "[observeEvent:selected_uniprot_...] new user selection: \n\t",
            paste(user_data$uniprot_ids, collapse = ";")))
        #store all gene data as a subset for this user
        user_data$selected_genes <- get_gene_data_for_uniprots(new_selection)
    })
    
    ## creates a DataTable view of the user-selected genes
    output$selected_genes_for_barplot <- DT::renderDataTable({
        selected_genes_summary <- get_gene_summary_table(user_data$selected_genes)
            DT::datatable(selected_genes_summary,
                      options = list(dom = 't'))
    })
    
    ## observes selected genes to update conditions checkbox group
    observeEvent(user_data$selected_genes, {
        represented_experiments <- get_represented_experiments(user_data$selected_genes)
        message("Updating conditions")

        updateCheckboxGroupButtons(session = session,
                                 inputId = "experiment_checkboxes",
                                 choices = represented_experiments,
                                 selected = get_preselected_experiments(),
                                 status = "primary",
                                checkIcon = list(
                                      yes = icon("ok", 
                                            lib = "glyphicon"),
                                      no = icon("remove",
                                            lib = "glyphicon")))
    })
    
    #observes both input$experiment_checkboxes and barplot_layout_switch
    observe({
        message("experiment checkboxes updated")

        req(user_data$uniprot_ids)
        message(paste0("[output$gene_conditions_bar..] plotting from: \n\t",
                        paste(user_data$uniprot_ids, collapse = ";")))
        output$gene_conditions_bar <- renderPlot({
            plot_bargraph_for_genes(user_data$selected_genes, 
                                   experiments = input$experiment_checkboxes,
                                   normalized_input_values = input$barplot_input_switch,
                                   faceted = input$barplot_layout_switch)
        })
        
    })
    
    output$download_selected_genes <- downloadHandler(
        filename = "selected_genes.csv",
        content = function(file) {
            write.csv(user_data$selected_genes, file, row.names = FALSE)
        }
    )
    
    ## creates the barplot of selected genes
    # output$gene_conditions_bar <- renderPlot({
    #     req(user_data$uniprot_ids)
    #     message(paste0("[output$gene_conditions_bar..] plotting from: \n\t",
    #                     paste(user_data$uniprot_ids, collapse = ";")))
    #     plot_bargraph_for_gene(user_data$selected_genes)
    # })
    
    
    ###### Gene filter (gf_xxx) related -- TODO: modularize ######
    
    # gf_cell_type <- reactive({
    #     req(input$celltypes_for_genes)
    #     message(paste0("[celltypes_for_genes]", input$celltypes_for_genes))
    #     input$celltypes_for_genes
    # })
    # 
    # observeEvent(gf_cell_type(), {
    #     message(paste0("[update gene filter checkboxes] for ", gf_cell_type()))
    # 
    #     user_data$filtered_genes <- filter_genes(celltypes = gf_cell_type())
    #     
    #     updateCheckboxGroupInput(session = session,
    #                              inputId = "source_organism_for_genes",
    #                              choices = get_source_organism_vector(gf_cell_type()),
    #                              inline = FALSE)
    # 
    # 
    #     updateCheckboxGroupInput(session = session,
    #                              inputId = "bait_filter_for_genes",
    #                              choices = get_baits_vector(gf_cell_type()),
    #                              inline = TRUE)
    # })
    # 

}