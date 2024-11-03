#' Calculate Position-Specific Conservation Scores
#'
#' This function calculates conservation scores for each position in a multiple sequence
#' alignment using Shannon's entropy. Lower entropy values indicate higher conservation
#' at that position. The scores are normalized to a 0-1 scale where 1 indicates
#' perfect conservation.
#'
#' @param alignment_data List containing alignment and metadata (output from createAlignment function)
#'
#' @return List containing:
#'   \itemize{
#'     \item scores: Numeric vector of conservation scores for each position (0-1 scale)
#'     \item metadata: List with calculation parameters and summary statistics
#'   }
#'
#' @examples
#' # Read sequences, create alignment, and calculate conservation scores.
#' Dermcidin_Lacritin_Sequences <- readSequences(system.file("extdata", "IPR043557_Dermcidin_Lacritin.fasta",
#'                           package = "ProtConservR"))
#' Dermcidin_Lacritin_Alignment <- createAlignment(Dermcidin_Lacritin_Sequences, method = "ClustalW")
#'
#' # Calculate conservation scores
#' Dermcidin_Lacritin_Conservation <- calculateConservation(Dermcidin_Lacritin_Alignment)
#'
#' @references
#' Valdar, W.S.J. (2002). Scoring residue conservation. \emph{Proteins: Structure,
#' Function, and Bioinformatics}, 48(2), 227-241.
#'
#' #' Barry J. Grant, Ana P. C. Rodrigues, Karim M. ElSawy, J. Andrew McCammon, Leo S.
#' D. Caves, Bio3d: an R package for the comparative analysis of protein structures,
#' \emph{Bioinformatics}, Volume 22, Issue 21, November 2006, Pages 2695–2696,
#' \href{https://doi.org/10.1093/bioinformatics/btl461}{link}.
#'
#' @export
#' @importFrom Biostrings consensusMatrix
calculateConservation <- function(alignment_data) {
  # Input validation
  if (!is.list(alignment_data) || !all(c("alignment", "metadata") %in% names(alignment_data))) {
    stop("Input must be a list containing 'alignment' and 'metadata' (output from createAlignment)")
  }

  # Extract alignment
  alignment <- alignment_data$alignment

  # Convert alignment to consensus matrix (probability matrix)
  cons_matrix <- Biostrings::consensusMatrix(alignment, as.prob = TRUE)

  # Remove gap row and renormalize probabilities
  if ("-" %in% rownames(cons_matrix)) {
    cons_matrix <- cons_matrix[rownames(cons_matrix) != "-", ]
    cons_matrix <- sweep(cons_matrix, 2, colSums(cons_matrix), "/")
  }

  # Calculate Shannon entropy for each position
  # Higher entropy = lower conservation
  entropy_scores <- apply(cons_matrix, 2, function(col) {
    # Remove zero probabilities to avoid log(0)
    col <- col[col > 0]
    -sum(col * log2(col))
  })

  # Convert entropy to conservation scores (0-1 scale)
  # Maximum possible entropy with 20 amino acids is -log2(1/20)
  # Invert so 1 = most conserved, 0 = least conserved
  max_entropy <- -log2(1/20)
  conservation_scores <- 1 - (entropy_scores / max_entropy)

  # Calculate summary statistics
  stats <- list(
    mean_conservation = mean(conservation_scores),
    median_conservation = median(conservation_scores),
    max_conservation = max(conservation_scores),
    min_conservation = min(conservation_scores),
    n_positions = length(conservation_scores)
  )

  return(list(
    scores = conservation_scores,
    metadata = stats
  ))
}
