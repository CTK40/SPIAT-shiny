# Made it read csv file and enable column selection.
library(shiny)
source("functions.R")
source("ui.R")
source("server.R")

options(shiny.maxRequestSize=10000*1024^2)
shinyApp(ui, server)
