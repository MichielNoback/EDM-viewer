ui <- fluidPage(
  checkboxInput("error", "error?"),
  textOutput("result")
)

server <- function(input, output, session) {
  a <- reactiveVal({
      tryCatch({
              req(!input$error) 
              return(1)
          }, 
          error = function(cond) {
              message(cond)
              return(NA)
           })
  })
  b <- reactive(a() + 1)
  c <- reactive(b() + 1)
  output$result <- renderText(c())
}

shinyApp(ui, server)