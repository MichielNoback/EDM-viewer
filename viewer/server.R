# op basis van celltype -- condities laden
# dan datatable laten zien
# hier selectie in maken
# in tabs verschillende visualisaties

# studenten: fair data viewer voor data collectie
#

server <- function(input, output) { 
    gene_name <- reactive({
        req(input$gene_id)
        substr(input$gene_id, 1, 9)
        #"tr|B4E2X5|B4E2X5_HUMAN"
    })
    output$gene_name <- renderText(gene_name())
}