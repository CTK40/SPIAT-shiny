Tab1 <- tabPanel(
    "Read, format and save",
    fluidPage(
        sidebarLayout(
            sidebarPanel(
                # Select file format
                selectInput("format", "Select file format",
                            choices = c("general", "inForm", "HALO", "Xenium", "Visium", 
                                        "MERSCOPE", "CosMX"),
                            selected = "inForm"),
                uiOutput("format"),
            ),
            mainPanel(
                # Show part of the table
                # Select columns for genes/markers
                conditionalPanel(
                    condition = "input.format !== 'Visium' && input.format !== 'Xenium' 
                                 && input.format !== 'inForm' && input.format !== 'HALO'",
                    fluidRow(
                    # This box is for selecting columns
                        column(2, uiOutput('cellID_gene')),
                        column(2, uiOutput('fov_gene')),
                        column(2, uiOutput('var_gene_select')),
                        column(2, uiOutput('var_gene_ignore'))
                )),
                    # Partial example
                conditionalPanel(
                    condition = "input.format == 'inForm' || input.format == 'HALO'",
                    fluidRow(
                        column(3, numericInput("n_markers", "The number of markers:", value = 1)),
                        column(3, uiOutput("marker")),
                        column(3, uiOutput('intensity_columns')),
                        conditionalPanel(
                            condition = "input.format == 'HALO'",
                            uiOutput('dye_columns')
                        )
                    )
                ),
                    
                # Select the range to show rows in the data (fast)
                fluidRow(
                    # Add a button to save the df object for marker intensity/gene expression
                    tags$head(
                        tags$style(HTML('#do{background-color:orange}'))
                    ),
                    column(3, actionButton("do", "Save the SpatialExperiment object"))),
                tableOutput("markerOrGene"),
                fluidRow(
                    conditionalPanel(
                        condition = "input.format !== 'Visium' && input.format !== 'Xenium'",
                        column(2, numericInput("row1", "Show rows from", value = 1, min = 1)),
                        column(2,  numericInput("row2", "Show rows to", value = 10, min = 1))
                    )
                ),
                conditionalPanel(
                    condition = "input.format == 'MERSCOPE' || input.format == 'CosMX' || input.format == 'general'",
                    fluidRow(
                        column(2, uiOutput('cellID_metadata')),
                        column(2, uiOutput("sample")),
                        column(2, uiOutput("phenotype")),
                        column(2, uiOutput("coord_x")),
                        column(2, uiOutput("coord_y"))
                    ),
                    tableOutput("metadata")
                )
            )
        )
    )
)
