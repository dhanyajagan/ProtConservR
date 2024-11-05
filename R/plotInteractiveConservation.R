#' Visualize Protein Conservation Scores Interactively
#'
#' Creates an interactive visualization of protein conservation scores using plotly.
#' The function generates multiple visualization types including a heatmap view of
#' the conservation scores, a line plot showing conservation across residues, and
#' summary statistics.
#'
#' @param mapped_data List containing mapped conservation scores (output from mapConservation)
#' @param view_type Character, type of visualization: "all" (default), "heatmap", "line", or "stats"
#' @param color_scheme Character, color scheme for conservation: "blues" (default) or "rainbow"
#'
#' @return plotly object containing the interactive visualization
#'
#' @examples
#' # Read and process data
#' sequences <- readSequences(system.file("extdata", "IPR043557_Dermcidin_Lacritin.fasta",
#'                           package = "ProtConservR"))
#' alignment <- createAlignment(sequences, method = "ClustalW")
#' conservation <- calculateConservation(alignment)
#' structure <- readStructure(system.file("extdata", "Dermcidin.pdb",
#'                           package = "ProtConservR"))
#' mapped_data <- mapConservation(conservation, structure)
#'
#' # Create interactive visualization
#' plotInteractiveConservation(mapped_data)
#'
#' @export
#' @importFrom plotly plot_ly layout subplot
#' @importFrom stats density
plotInteractiveConservation <- function(mapped_data,
                                  view_type = c("all", "heatmap", "line", "stats"),
                                  color_scheme = c("blues", "rainbow")) {
  # Input validation
  view_type <- match.arg(view_type)
  color_scheme <- match.arg(color_scheme)

  if (!is.list(mapped_data) || !all(c("mapped_scores", "metadata") %in% names(mapped_data))) {
    stop("mapped_data must be output from mapConservation function")
  }

  # Extract data
  scores_df <- mapped_data$mapped_scores
  meta <- mapped_data$metadata

  # Define color scales
  if (color_scheme == "blues") {
    colors <- colorRamp(c("#F7FBFF", "#08519C"))
  } else {
    colors <- colorRamp(c("#FF0000", "#FFFF00", "#00FF00", "#00FFFF", "#0000FF"))
  }

  # Create heatmap plot
  heatmap <- plot_ly(
    x = scores_df$residue_number,
    y = 1,
    z = matrix(scores_df$conservation_score, nrow = 1),
    type = "heatmap",
    colorscale = "Viridis",
    showscale = TRUE,
    hoverongaps = FALSE
  ) %>%
    layout(
      title = "Conservation Score Heatmap",
      xaxis = list(title = "Residue Number"),
      yaxis = list(showticklabels = FALSE),
      margin = list(t = 50)
    )

  # Create line plot
  line <- plot_ly(
    x = scores_df$residue_number,
    y = scores_df$conservation_score,
    type = "scatter",
    mode = "lines+markers",
    marker = list(size = 6),
    line = list(width = 2),
    text = paste("Residue:", scores_df$residue_type,
                 "<br>Position:", scores_df$residue_number,
                 "<br>Conservation:", round(scores_df$conservation_score, 3)),
    hoverinfo = "text"
  ) %>%
    layout(
      title = "Conservation Score Profile",
      xaxis = list(title = "Residue Number"),
      yaxis = list(title = "Conservation Score",
                   range = c(0, 1)),
      margin = list(t = 50)
    )

  # Create statistics plot
  stats_data <- data.frame(
    stat = c("Mean", "Median", "Min", "Max"),
    value = c(
      mean(scores_df$conservation_score, na.rm = TRUE),
      median(scores_df$conservation_score, na.rm = TRUE),
      min(scores_df$conservation_score, na.rm = TRUE),
      max(scores_df$conservation_score, na.rm = TRUE)
    )
  )

  stats <- plot_ly(
    stats_data,
    x = ~stat,
    y = ~value,
    type = "bar",
    marker = list(color = "#1f77b4")
  ) %>%
    layout(
      title = "Conservation Statistics",
      xaxis = list(title = "Statistic"),
      yaxis = list(title = "Value",
                   range = c(0, 1)),
      margin = list(t = 50)
    )

  # Return appropriate plot based on view_type
  if (view_type == "heatmap") {
    return(heatmap)
  } else if (view_type == "line") {
    return(line)
  } else if (view_type == "stats") {
    return(stats)
  } else {
    # Combine all plots
    subplot(
      list(heatmap, line, stats),
      nrows = 3,
      heights = c(0.3, 0.4, 0.3),
      margin = 0.1
    ) %>%
      layout(
        title = "Protein Conservation Analysis",
        showlegend = FALSE,
        margin = list(t = 50)
      )
  }
}
