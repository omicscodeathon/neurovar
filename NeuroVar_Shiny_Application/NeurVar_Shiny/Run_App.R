#  NeuroVar
## Run the following  code in order
source("Dependencies_NeuroVar.R")
source("UI.R",local = TRUE)

source("Server.R",local = TRUE)

# Run the application
shinyApp(ui = UI, server = Server)
