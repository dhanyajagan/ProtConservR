#' Map Conservation Scores to Protein Structure
#'
#' This function maps sequence conservation scores onto protein structure by
#' aligning the structural residues with the corresponding positions in the
#' multiple sequence alignment. It handles potential gaps and mismatches between
#' the structure and alignment.
#'
#' @param conservation_data List containing conservation scores and metadata (output from calculateConservation)
#' @param structure_data List containing structural information (output from readStructure)
#' @param chain Character, chain ID to map conservation scores to (default: first chain)
#'
#' @return List containing:
#'   \itemize{
#'     \item mapped_scores: Data frame with residue numbers, residue types, and conservation scores
#'     \item structure: Updated structure object with conservation scores added
#'     \item metadata: List with mapping statistics and parameters
#'   }
#'
#' @examples
#' # Read sequences and create alignment
#' Dermcidin_Lacritin_Sequences <- readSequences(system.file("extdata", "IPR043557_Dermcidin_Lacritin.fasta",
#'                           package = "ProtConservR"))
#' Dermcidin_Lacritin_Alignment <- createAlignment(Dermcidin_Lacritin_Sequences, method = "ClustalW")
#'
#' # Calculate conservation scores
#' Dermcidin_Lacritin_Conservation <- calculateConservation(Dermcidin_Lacritin_Alignment)
#'
#' # Read structure
#' Dermcidin_Lacritin_structure <- readStructure(system.file("extdata", "Dermcidin.pdb", package = "ProtConservR"))
#'
#' # Map conservation scores to structure
#' Dermcidin_Lacritin_Mapped_Data <- mapConservation(Dermcidin_Lacritin_Conservation, Dermcidin_Lacritin_structure)
#'
#' @references
#' Grant, B.J. et al. (2006). Bio3d: an R package for the comparative analysis
#' of protein structures. \emph{Bioinformatics}, 22(21), 2695-2696.
#'
#' @export
#' @importFrom bio3d atom.select pdbseq
mapConservation <- function(conservation_data, structure_data, chain = NULL) {
  # Input validation
  if (!is.list(conservation_data) || !all(c("scores", "metadata") %in% names(conservation_data))) {
    stop("Conservation data must be a list containing 'scores' and 'metadata'")
  }
  if (!is.list(structure_data) || !all(c("atoms", "residues") %in% names(structure_data))) {
    stop("Structure data must be a list containing 'atoms' and 'residues'")
  }

  # If chain is not specified, use the first chain in the structure
  if (is.null(chain)) {
    chain <- structure_data$chains[1]
  }

  # Extract relevant data
  conservation_scores <- conservation_data$scores
  atom_data <- structure_data$atoms

  # Filter atoms for specified chain
  chain_atoms <- atom_data[atom_data$chain == chain, ]

  # Get unique residues for the chain
  chain_residues <- unique(chain_atoms[, c("resno", "resid")])
  chain_residues <- chain_residues[order(chain_residues$resno), ]

  # Create mapping between structure residues and conservation scores
  n_residues <- nrow(chain_residues)
  n_scores <- length(conservation_scores)

  # Handle case where number of residues doesn't match number of conservation scores
  if (n_residues != n_scores) {
    warning(sprintf("Number of residues (%d) doesn't match number of conservation scores (%d).
                   Scores will be mapped based on position.", n_residues, n_scores))
  }

  # Create mapping data frame
  mapped_data <- data.frame(
    residue_number = chain_residues$resno,
    residue_type = chain_residues$resid,
    conservation_score = NA_real_
  )

  # Map conservation scores to residues
  # Use minimum length to avoid index out of bounds
  max_idx <- min(n_residues, n_scores)
  mapped_data$conservation_score[1:max_idx] <- conservation_scores[1:max_idx]

  # Calculate mapping statistics
  mapping_stats <- list(
    chain_id = chain,
    n_residues_structure = n_residues,
    n_positions_alignment = n_scores,
    n_mapped_positions = max_idx,
    mean_conservation_mapped = mean(mapped_data$conservation_score, na.rm = TRUE),
    date = Sys.Date()
  )

  # Update atom data with conservation scores
  atom_data$conservation <- NA_real_
  for (i in 1:nrow(mapped_data)) {
    res_idx <- atom_data$resno == mapped_data$residue_number[i] &
      atom_data$chain == chain
    atom_data$conservation[res_idx] <- mapped_data$conservation_score[i]
  }

  # Update structure data
  structure_data$atoms <- atom_data

  return(list(
    mapped_scores = mapped_data,
    structure = structure_data,
    metadata = mapping_stats
  ))
}

#' Get Conservation Score for Specific Residue
#'
#' Helper function to retrieve the conservation score for a specific residue
#'
#' @param mapped_data Output from mapConservation function
#' @param residue_number Positive integer indicating the residue number to to retrieve the conservation score
#' @param chain Character chain ID (optional)
#'
#' @return Numeric conservation score or NA if residue not found
#'
#' @keywords internal
get_residue_conservation <- function(mapped_data, residue_number, chain = NULL) {
  if (is.null(chain)) {
    score <- mapped_data$mapped_scores$conservation_score[
      mapped_data$mapped_scores$residue_number == residue_number
    ]
  } else {
    score <- mapped_data$mapped_scores$conservation_score[
      mapped_data$mapped_scores$residue_number == residue_number &
        mapped_data$structure$atoms$chain == chain
    ]
  }
  return(if(length(score) == 0) NA_real_ else score[1])
}
