#' Launch Shiny App for ProtConservR
#'
#' A function that launches the Shiny app for ProtConservR.
#' This app allows for the upload of a protein family multifasta file of choice,
#' perform alignment of choice (ClustalW or MUSCLE), compute conservation scores
#' and plot them as either a line plot or bar plot, and map these sequence based
#' conservation scores onto protein structure (with a specified style). Note: for
#' the 3D protein structure visualization feature, selection cartoon results in
#' the best visualization, as the other options include other colours for the
#' atoms-- it can get confusing. The code has been placed in \code{./inst/shiny-scripts}.
#'
#' @return No return value but open up a Shiny page.
#'
#' @examples
#' \dontrun{
#'
#' ProtConservR::runProtConservR()
#' }
#'
#' @references
#' Grolemund, G. (2015). Learn Shiny - Video Tutorials. \href{https://shiny.rstudio.com/tutorial/}{Link}
#'
#' @export
#' @importFrom shiny runApp

runProtConservR <- function() {
  appDir <- system.file("shiny-scripts",
                        package = "ProtConservR")
  actionShiny <- shiny::runApp(appDir, display.mode = "normal")
  return(actionShiny)
}
# [END]
