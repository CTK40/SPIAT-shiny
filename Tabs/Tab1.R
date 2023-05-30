Tab1 <- tabPanel(
    "Read, format and save",
    fluidPage(
        sidebarLayout(
            sidebarPanel(
                # Select file format
                selectInput("format", "Select file format",
                            choices = c("general", "inForm", "HALO", "Xenium", "Visium", 
                                        "MERSCOPE", "CosMX", "cellprofiler", "CODEX"),
                            selected = "inForm"),
                uiOutput("format"),
                # checkboxInput("header", "Header", TRUE)
            ),
            mainPanel(
                # Show part of the table
                # Select columns for genes/markers
                fluidRow(
                    # This box is for selecting columns
                    conditionalPanel(
                        condition = "input.format != 'Visium'",
                        column(3, uiOutput('cellID_gene')),
                        column(3, uiOutput('fov_gene')),
                        column(4, uiOutput('var_gene_select')),
                        column(4,  uiOutput('var_gene_ignore'))
                    ),
                    # Partial example
                    conditionalPanel(
                        condition = "input.format == 'inForm'",
                        numericInput("n_markers", "The number of markers:", value = 1),
                        uiOutput("marker")
                    )),
                    
                # Select the range to show rows in the data (fast)
                fluidRow(
                    conditionalPanel(
                        condition = "input.format != 'Visium' && input.format != 'Xenium'",
                        column(3, numericInput("row1", "Show rows from", value = 1, min = 1)),
                        column(3,  numericInput("row2", "Show rows to", value = 10, min = 1))
                    ),
                    # Add a button to save the df object for marker intensity/gene expression
                    tags$head(
                        tags$style(HTML('#do{background-color:orange}'))
                    ),
                    column(3, actionButton("do", "Save the SpatialExperiment object"))),
                tableOutput("markerOrGene"),
                fluidRow(
                    column(3, uiOutput('cellID_metadata')),
                    column(3, uiOutput("sample")),
                    column(3, uiOutput("phenotype")),
                    column(3, uiOutput("coord_x")),
                    column(3, uiOutput("coord_y"))
                ),
                tableOutput("metadata")
            )
        )
    )
)