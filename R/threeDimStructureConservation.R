#' Visualize Protein Conservation in 3D
#'
#' Creates an interactive 3D visualization of protein structure colored by conservation
#' scores. Uses r3dmol to create an interactive viewer where conservation scores are
#' mapped to a color gradient from red (low conservation) to white (medium) to blue
#' (high conservation). NOTE:the conservation is not that apparent in the sphere
#' and stick visualization style. For best results, try the cartoon style
#'
#' @param input_pdb Path to input PDB file. Make sure to add a newline character to the PDB file if there is not one already
#' @param mapped_data List containing mapped conservation scores from mapConservation()
#' @param width Numeric, width of viewer in pixels (default: 800)
#' @param height Numeric, height of viewer in pixels (default: 600)
#' @param background Character, background color of viewer (default: "white")
#' @param style Character, visualization style ("cartoon", "sphere", "stick")
#'
#' @return r3dmol viewer object that can be displayed in R markdown or R console
#'
#' @details
#' The function creates an interactive 3D visualization where:
#' * Mouse left button rotates the structure
#' * Mouse right button translates
#' * Scroll wheel controls zoom
#' * Double click resets the view
#'
#' Conservation scores are mapped using color gradients where:
#' * Blue indicates high conservation (score close to 1)
#' * White indicates medium conservation (score around 0.5)
#' * Red indicates low conservation (score close to 0)
#'
#' @examples
#' # Read sequences and create alignment
#' sequences <- readSequences(system.file("extdata", "IPR043557_Dermcidin_Lacritin.fasta",
#'                           package = "ProtConservR"))
#' alignment <- createAlignment(sequences, method = "ClustalW")
#'
#' # Calculate conservation scores
#' conservation <- calculateConservation(alignment)
#'
#' # Read structure and map conservation
#' structure <- readStructure(system.file("extdata", "Dermcidin.pdb",
#'                          package = "ProtConservR"))
#' mapped_data <- mapConservation(conservation, structure)
#'
#' # Create interactive visualization
#' threeDimStructureConservation(
#'   system.file("extdata", "Dermcidin.pdb", package = "ProtConservR"),
#'   mapped_data = mapped_data
#' )
#'
#' @references
#' Nicholas Rego and David Koes, 3Dmol.js: molecular visualization with WebGL,
#' \emph{Bioinformatics}, Volume 31, Issue 8, April 2015, Pages 1322–1324,
#' \href{https://doi.org/10.1093/bioinformatics/btu829}{link}.
#'
#' @export
#' @importFrom r3dmol r3dmol m_add_model m_set_style m_style_cartoon m_zoom_to m_spin m_rotate m_sel m_add_label
#' @importFrom grDevices rgb colorRampPalette
threeDimStructureConservation <- function(input_pdb,
                                          mapped_data,
                                          width = 800,
                                          height = 600,
                                          background = "white",
                                          style = c("cartoon", "sphere", "stick")) {
  # Input validation
  style <- match.arg(style)

  if (!is.list(mapped_data) || !all(c("mapped_scores", "structure", "metadata") %in% names(mapped_data))) {
    stop("mapped_data must be the output from mapConservation function")
  }

  # Create color palette function (red -> white -> blue)
  color_palette <- colorRampPalette(c("#FF0000", "#FFFFFF", "#0000FF"))

  # Extract conservation scores
  scores_df <- mapped_data$mapped_scores

  # Initialize the viewer
  viewer <- r3dmol::r3dmol(
    id = "r3dmol_viewer",
    elementId = "r3dmol",
    width = width,
    height = height,
    backgroundColor = background
  ) %>%
    r3dmol::m_add_model(
      data = input_pdb,
      format = "pdb"
    ) %>%
    m_zoom_to()

  # Generate colors for each residue based on conservation scores
  n_colors <- 100
  color_map <- color_palette(n_colors)

  # Function to get color for a conservation score
  get_color <- function(score) {
    if (is.na(score)) return("#CCCCCC")  # Gray for NA values
    idx <- floor(score * (n_colors - 1)) + 1
    return(color_map[idx])
  }

  # Apply style and coloring based on conservation for each residue
  for (i in seq_len(nrow(scores_df))) {
    residue <- scores_df[i, ]
    color <- get_color(residue$conservation_score)

    # Create style based on selected visualization type
    style_params <- switch(style,
                           "cartoon" = m_style_cartoon(color = color),
                           "sphere" = m_style_sphere(color = color),
                           "stick" = m_style_stick(color = color)
    )

    # Apply style to specific residue
    viewer <- viewer %>%
      m_set_style(
        sel = m_sel(resi = residue$residue_number),
        style = style_params
      )
  }


  # Final viewer settings
  viewer <- viewer %>%
    m_rotate(angle = 90, axis = "y") %>%
    m_spin()

  return(viewer)
}
