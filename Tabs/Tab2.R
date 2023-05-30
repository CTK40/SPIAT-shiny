library(shinyFiles)
Tab2 <- tabPanel(
    "Plot cell categories",
    fluidPage(
        # Show part of the table
        # Select columns for genes/markers                
        fluidRow(
            # This box is for selecting columns
            column(3, uiOutput('categories')),
            column(3, uiOutput('feature_colname')),
            column(4, uiOutput('colour'))
        ),
        plotOutput("plot_cells", width = "600px", height = "600px")
    )
)
