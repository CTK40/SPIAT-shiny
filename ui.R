source("Tabs/Tab1.R")
source("Tabs/load_object.R")
source("Tabs/Tab2.R")
source("Tabs/Tab3.R")
ui <- fluidPage(
    tabsetPanel(
        Tab1, # read file
        load_spe_object,
        Tab2, # plot cell categories
        Tab3  # plot marker intensities
    )
)
