library(shinyFiles)
load_spe_object <- tabPanel(
    "Load object",
    fluidPage(
        # fileInput("spe_object", "Choose spe object file",
        #           accept = c(".Rda", ".rda", ".RData", ".RDS"))
        shinyFilesButton("Btn_GetFile", "Choose a file" ,
                         title = "Please select a file:", multiple = FALSE,
                         buttonType = "default", class = NULL)
    )
)
