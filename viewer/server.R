# op basis van celltype -- laadt condities
# dan datatable laten zien
# hier selectie in maken
# in tabs verschillende visualisaties

# studenten: fair data viewer voor data collectie
#

server <- function(input, output) {
    output$condition_select <- renderUI({
        cell_types <- input$cell_types
        if (is.null(cell_types)) {
            conditions <- ""
            names(conditions) <- "choose celltype(s)"
        } else {
            conditions <- get_conditions_vector(cell_types)
        }
        selectInput(inputId = "condition_select", 
                    multiple = TRUE,
                    label = "Condition",
                    
                    choices = conditions)
    })
    
    output$selected_condition <- renderPrint({
        input$condition
    })
    
#use package dt for selections
    
    output$gene_select <- renderUI({
        cell_types <- input$cell_types
        if (is.null(cell_types)) {
            genes <- ""
            names(genes) <- "select celltype"
        } else {
            genes <- get_genes_vector(cell_types)
        }
        selectInput(inputId = "gene", 
                    multiple = FALSE,
                    label = "Gene or identifier",
                    choices = genes)
    })
    
    output$selected_gene <- renderPrint({
        input$gene
    })
    
    output$conditions_box <- renderPlot({
        
        # plt <- all_data %>% 
        #     filter(Condition %in% conditions[sample.int(length(conditions), 3)]) %>%
        #     mutate(signal_log2 = log2(signal)) %>%
        #     ggplot(mapping = aes(x = Condition, y = signal_log2, color = Condition)) +
        #         geom_boxplot() +
        #         ylab("Log2 of the signal") +
        #         theme(legend.position = "none")
        # return(plt)
                
        print(paste0("[1]", input$cell_types, "__|", input$condition_select, "|__"))
        #beter!
        #req(input$cell_types)
        if (is.null(input$cell_types) 
            || is.null(input$condition_select) 
            || input$condition_select == "") {
            print("NOT OK")
        } else {
            print("OK")
            plt <- all_data %>% 
                filter(Condition %in% input$condition_select) %>%
                mutate(signal_log2 = log2(signal)) %>%
                ggplot(mapping = aes(x = Condition, y = signal_log2, color = Condition)) +
                    geom_boxplot() +
                    ylab("Log2 of the signal") +
                    theme(legend.position = "none")
            return(plt)
        }
    })
}

