library(shinyFiles)
load_spe_object <- tabPanel(
    "Load object",
    fluidPage(
        shinyFilesButton("Btn_GetFile", "Choose a file" ,
                         title = "Please select a file:", multiple = FALSE)
    )
)
