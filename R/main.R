# You can learn more about package authoring with RStudio at:
#
#   http://r-pkgs.had.co.nz/
#
# Some useful keyboard shortcuts for package authoring:
#
#   Build and Reload Package:  'Ctrl + Shift + B'
#   Check Package:             'Ctrl + Shift + E'
#   Test Package:              'Ctrl + Shift + T'

#' @export
AnalysisTool.run <- function()
{
    library(flowCore)
    library(microbenchmark)
    library(ncdfFlow)
    library(shiny)
    library(shinydashboard)
    library(shinyjs)
    library(DT)
    library(RColorBrewer)
    library(ggplot2)


    appDir <- system.file("shinyApp", "app", package = "AnalysisTool")
    if (appDir == "")
    {
        stop("Could not find app directory. Try re-installing `AnalysisTool`.", call. = FALSE)
    }
    
    if(!dir.exists(paste0(appDir,"/Projects/")))
    {
        dir.create(paste0(appDir,"/Projects/"))
    }
    shiny::runApp(appDir, display.mode = "normal", launch.browser = T)
}
