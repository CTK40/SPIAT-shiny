Tab3 <- tabPanel(
    "Plot marker intensities",
    fluidPage(
        mainPanel(
            # Select columns for genes/markers    
            uiOutput('marker_level'),
            plotOutput("plot_marker", width = "600px", height = "600px")
        )
    )
)
