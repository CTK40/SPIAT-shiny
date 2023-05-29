library(shinyFiles)
Tab3 <- tabPanel(
    "Plot marker intensities",
    fluidPage(
        sidebarLayout(
            sidebarPanel(
                # load object
                # fileInput("spe_object", "Choose spe object file",
                #           accept = c(".Rda", ".rda", ".RData", ".RDS"))

            ),
            mainPanel(
                # Show part of the table
                # Select columns for genes/markers    
                uiOutput('marker_level'),
                plotOutput("plot_marker", width = "600px", height = "600px")
            )
        )
    )
)
