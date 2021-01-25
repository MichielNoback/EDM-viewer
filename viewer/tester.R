histogramUI <- function(id) {
  tagList(
    selectInput(NS(id, "var"), "Variable", names(mtcars)),
    numericInput(NS(id, "bins"), "bins", 10, min = 1),
    plotOutput(NS(id, "hist"))
  )
}

histogramServer <- function(id) {
  moduleServer(id, function(input, output, session) {
    data <- reactive(mtcars[[input$var]])
    output$hist <- renderPlot({
      hist(data(), breaks = input$bins, main = input$var)
    }, res = 96)
  })
}

histogramApp <- function() {
  ui <- fluidPage(
    histogramUI("hist1")
  )
  server <- function(input, output, session) {
    histogramServer("hist1")
  }
  shinyApp(ui, server)  
}

#ff checke

#all_data %>% ggplot(mapping = aes(x = z_score_spectral)) + geom_histogram(aes(y =..density..), binwidth = 2)

# cutoff <- 2
# all_data %>% 
#     mutate(z_label = cut(z_score_spectral, 
#                          breaks = c(-Inf, -cutoff, cutoff, Inf), 
#                          labels = c("low_z", "unchanged", "high_z"))) %>%
#     select(z_label) %>%
#     table()



shinyApp(ui, server)

