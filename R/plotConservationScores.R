#' Plot Conservation Scores
#'
#' This function creates a plot of conservation scores for each position in a
#' multiple sequence alignment.
#'
#' @param conservation_data List containing:
#'   \itemize{
#'     \item scores: Numeric vector of conservation scores for each position (0-1 scale)
#'     \item metadata: List with calculation parameters and summary statistics
#'   }
#' @param plot_type Character string specifying the type of plot: "line" or "bar"
#'
#' @return A ggplot object showing the conservation scores.
#'
#' @examples
#' Dermcidin_Lacritin_Sequences <- readSequences(system.file("extdata", "IPR043557_Dermcidin_Lacritin.fasta",
#'                           package = "ProtConservR"))
#' Dermcidin_Lacritin_Alignment <- createAlignment(Dermcidin_Lacritin_Sequences, method = "ClustalW")
#'
#' # Calculate conservation scores
#' Dermcidin_Lacritin_Conservation <- calculateConservation(Dermcidin_Lacritin_Alignment)
#' conservationScorePlotLine <- plotConservationScores(Dermcidin_Lacritin_Conservation, plot_type = "line")
#' conservationScorePlotBar <- plotConservationScores(Dermcidin_Lacritin_Conservation, plot_type = "bar")
#'
#' @references
#' Wickham H (2016). \emph{ggplot2: Elegant Graphics for Data Analysis.}
#' Springer-Verlag New York. ISBN 978-3-319-24277-4, https://ggplot2.tidyverse.org.
#'
#' @export
#' @importFrom ggplot2 ggplot aes geom_line geom_point geom_bar labs theme_minimal
plotConservationScores <- function(conservation_data, plot_type = "line") {
  # Input validation
  if (!is.list(conservation_data) || !"scores" %in% names(conservation_data)) {
    stop("Input must be a list containing 'scores'")
  }

  # Extract conservation scores
  conservation_scores <- conservation_data$scores

  # Create a data frame for plotting
  plot_data <- data.frame(
    Position = 1:length(conservation_scores),
    Conservation = conservation_scores
  )

  # Generate the plot based on plot_type
  p <- ggplot(plot_data, aes(x = Position, y = Conservation)) +
    labs(title = "Position-Specific Conservation Scores",
         x = "Position in Alignment",
         y = "Conservation Score (0-1)") +
    theme_minimal()

  if (plot_type == "line") {
    p <- p + geom_line(color = "blue") + geom_point(color = "darkred")
  } else if (plot_type == "bar") {
    p <- p + geom_bar(stat = "identity", fill = "blue", color = "black")
  } else {
    stop("Invalid plot_type. Choose either 'line' or 'bar'.")
  }

  return(p)
}
