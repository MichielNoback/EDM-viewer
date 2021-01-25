# op basis van celltype -- condities laden
# dan datatable laten zien
# hier selectie in maken
# in tabs verschillende visualisaties

# studenten: fair data viewer voor data collectie
#

server <- function(input, output, session) { 
    gene_name <- reactive({
        req(input$gene_select)
        input$gene_select
    })
    
    gene_celltype <- reactive({
        req(input$celltypes_for_genes)
        input$celltypes_for_genes
    })
    
    
    observe({
        req(input$celltypes_for_genes)
            updateSelectizeInput(session = session, 
                         inputId = "gene_select", 
                         choices = get_genes_vector(input$celltypes_for_genes), 
                         server = TRUE)
    })
    ##genes reactive; depends on choice of cell type
    # genes <- reactive({
    #     #cell_types <- input$cell_types
    #     cell_types <- c("dicty")
    #     if (is.null(input$cell_types)) {
    #         genes_vector <- ""
    #         names(genes_vector) <- "select celltype"
    #     } else {
    #         genes_vector <- get_genes_vector(cell_types)
    #     }
    #     return(genes_vector)
    # })
    
    
    output$gene_name <- renderText(gene_name())
}