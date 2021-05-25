library(shinyWidgets)

# build top-level tab for genes query
geneSearchTab <- sidebarLayout(
    sidebarPanel(
        textInput(inputId = "gene_search",
                  label = "Find gene(s)"),
        dataTableOutput(outputId = "suggested_genes"),
        width = 5
    ),

    mainPanel(
        h4("Selected genes"),
        dataTableOutput(outputId = "selected_genes_for_barplot"),
        br(""),
        downloadButton("download_selected_genes", "Download selection"),
        width = 7
    )
)

# build top-level tab for genes visualization
geneVisualizationTab <- sidebarLayout(
    sidebarPanel(
        h5("Plot layout"),
        switchInput(
            inputId = "barplot_layout_switch",
            label = "facet plot",
            labelWidth = "80px",
            value = FALSE),
        h5("Input values"),
        switchInput(
            inputId = "barplot_input_switch",
            label = "normalized",
            labelWidth = "80px",
            value = TRUE),        
        checkboxGroupButtons(
            inputId = "experiment_checkboxes",
            label = "Represented experiments", 
            choices = c("Not Available"),
            status = "primary",
            checkIcon = list(
                  yes = icon("ok", 
                        lib = "glyphicon"),
                  no = icon("remove",
                        lib = "glyphicon")))
    ),

    mainPanel(
        plotOutput("gene_conditions_bar"),
    )
)


# main UI
ui <- fluidPage(
    titlePanel("EDM data browser"),
    tabsetPanel(id = "main-tabs", 
                type = "tabs",
                tabPanel("Select Genes", geneSearchTab),
                tabPanel("Visualize", geneVisualizationTab)))
#                tabPanel("Filter Genes", geneFilterTab),
#)


# build top-level tab for gene filtering

# spectral_count_range <- c(0, 7500)
# z_score_limit <- 50
# mw_range <- c(1, 750)
# 
# geneFilterTab <- sidebarLayout(
#     sidebarPanel(
#         checkboxGroupInput(inputId = "celltypes_for_genes",
#                         label = "Select one or more Cell types",
#                         choices = get_celltypes()),
#         
#         checkboxGroupInput(inputId = "source_organism_for_genes",
#                            label = "Source organism for captured proteins",
#                            inline = TRUE,
#                            choices = NULL),
#         
# 
#         checkboxGroupInput(inputId = "bait_filter_for_genes",
#                            label = "Baits",
#                            inline = TRUE,
#                            choices = NULL),        
#         ## TODO min and max should be set dynamically
#         sliderInput(inputId = "mw_filter", 
#                     label = "Molecular Weight (kDa)",
#                     value = mw_range,
#                     min = mw_range[1], max = mw_range[2]),
# 
#         ## TODO max should be set dynamically
#         sliderInput(inputId = "spectral_count_filter", 
#                     label = "Spectral count",
#                     value = spectral_count_range,
#                     min = spectral_count_range[1], max = spectral_count_range[2]),
#         
#         ## TODO min and max should be set dynamically
#         sliderInput(inputId = "spectral_z_score_filter", 
#                     label = "Z-score spectral count",
#                     value = c(-z_score_limit, z_score_limit),
#                     min = -z_score_limit, max = z_score_limit), 
#         
#     ),
#     mainPanel(
#         h2("Genes matching the filters"),
#         dataTableOutput(outputId = "matching_genes")
#     )
# )




