library(shinyFiles)
Tab2 <- tabPanel(
    "Second tab",
    fluidPage(
        sidebarLayout(
            sidebarPanel(
                # load object
                # fileInput("spe_object", "Choose spe object file",
                #           accept = c(".Rda", ".rda", ".RData", ".RDS"))
                shinyFilesButton("Btn_GetFile", "Choose a file" ,
                                 title = "Please select a file:", multiple = FALSE,
                                 buttonType = "default", class = NULL)
            ),
            mainPanel(
                # Show part of the table
                # Select columns for genes/markers
                plotOutput("plot_cells", width = "600px", height = "600px"),
                fluidRow(
                    # This box is for selecting columns
                    column(3, uiOutput('feature_colname')),
                    column(3, uiOutput('categories_of_interest')),
                    column(4, uiOutput('colour_vector')),
                )
            )
        )
    )
)
