#' Create Multiple Sequence Alignment
#'
#' Performs multiple sequence alignment on protein sequences using either MUSCLE or
#' ClustalW algorithm with default settings. Takes processed sequences from
#' readSequences() function as input.
#'
#' @param seq_data List containing sequences and metadata (output from readSequences function)
#' @param method Character, alignment method to use: "MUSCLE" (default) or "ClustalW"
#'
#' @return List containing:
#'   \itemize{
#'     \item alignment: The multiple sequence alignment object
#'     \item metadata: List with alignment statistics
#'   }
#'
#' @examples
#' # Read sequences: Dermcidin/Lacritin
#' Dermcidin_Lacritin_Sequences <- readSequences(system.file("extdata", "IPR043557_Dermcidin_Lacritin.fasta",
#'                                      package = "ProtConservR"))
#' # Create alignment using MUSCLE
#' alignment_muscle_DL <- createAlignment(Dermcidin_Lacritin_Sequences, method = "MUSCLE")
#' # Create alignment using ClustalW
#' alignment_clustal_DL <- createAlignment(Dermcidin_Lacritin_Sequences, method = "ClustalW")
#'
#' @references
#' Edgar, R.C. (2004) MUSCLE: multiple sequence alignment with high accuracy and
#' high throughput. \emph{Nucleic Acids Res.} 32(5):1792-1797.
#'
#' Thompson, J.D., Higgins, D.G., Gibson, T.J. (1994) CLUSTAL W: improving the
#' sensitivity of progressive multiple sequence alignment through sequence weighting,
#' position-specific gap penalties and weight matrix choice.
#' \emph{Nucleic Acids Res.} 22(22):4673-4680.
#'
#' Pagès H, Aboyoun P, Gentleman R, DebRoy S (2024). \emph{Biostrings: Efficient
#' manipulation of biological strings.} R package version 2.74.0,
#' \href{https://bioconductor.org/packages/Biostrings}{Link}.
#'
#' Barry J. Grant, Ana P. C. Rodrigues, Karim M. ElSawy, J. Andrew McCammon, Leo S.
#' D. Caves, Bio3d: an R package for the comparative analysis of protein structures,
#' \emph{Bioinformatics}, Volume 22, Issue 21, November 2006, Pages 2695–2696,
#' \href{https://doi.org/10.1093/bioinformatics/btl461}{link}.
#'
#' Bodenhofer U, Bonatesta E, Horejs-Kainrath C, Hochreiter S (2015).
#' “msa: an R package for multiple sequence alignment.” \emph{Bioinformatics},
#' 31(24), 3997–3999. doi:10.1093/bioinformatics/btv494
#'
#' @export
#' @importFrom Biostrings width
#' @importFrom muscle muscle
#' @importFrom msa msa
#' @import seqinr

createAlignment <- function(seq_data, method = c("MUSCLE", "ClustalW")) {

  # Input validation
  if (!is.list(seq_data) || !all(c("sequences", "metadata") %in% names(seq_data))) {
    stop("Input must be a list containing 'sequences' and 'metadata' (output from readSequences)")
  }

  method <- match.arg(method)

  # Extract sequences
  sequences <- seq_data$sequences

  # Perform alignment based on selected method
  if (method == "MUSCLE") {
    tryCatch({
      muscle_alignment <- msa::msa(sequences, method = "Muscle")
      # Convert to AAMultipleAlignment object and get consensus matrix
      alignment <- muscle_alignment@unmasked
    }, error = function(e) {
      stop("MUSCLE alignment failed: ", e$message)
    })
  } else if (method == "ClustalW") {
    tryCatch({
      # Perform ClustalW alignment
      clustal_alignment <- msa::msa(sequences, method = "ClustalW")
      # Convert to AAMultipleAlignment object and get consensus matrix
      alignment <- clustal_alignment@unmasked
    }, error = function(e) {
      stop("ClustalW alignment failed: ", e$message)
    })
  }

  # Get alignment dimensions
  if (method == "MUSCLE") {
    n_seq <- length(alignment)
    align_length <- width(alignment)[1]
  } else {
    n_seq <- nrow(as.matrix(alignment))
    align_length <- ncol(as.matrix(alignment))
  }

  # Calculate alignment statistics
  alignment_stats <- list(
    n_sequences = n_seq,
    alignment_length = align_length,
    gap_percentage = calculate_gap_percentage(alignment),
    alignment_method = method,
    date = Sys.Date()
  )

  return(list(
    alignment = alignment,
    metadata = alignment_stats
  ))
}



#' Calculate Gap Percentage in Alignment
#'
#' Helper function to calculate the percentage of gaps in the alignment
#'
#' @param alignment AAMultipleAlignment object or MsaAAMultipleAlignment object
#' @return Numeric value representing gap percentage
#'
#' @references Alignment scoring, Gaps & Similarities. Sequence Alignment:
#' Scores, Gaps and Gap Penalties. [accessed 2024 Nov 5].
#' https://proteinstructures.com/sequence/sequence-alignment/
#'
#' @keywords internal
calculate_gap_percentage <- function(alignment) {
  # Convert alignment to matrix, handling both types of alignment objects
  if (inherits(alignment, "MsaAAMultipleAlignment")) {
    align_matrix <- as.matrix(alignment@unmasked)
  } else {
    align_matrix <- as.matrix(alignment)
  }

  total_positions <- prod(dim(align_matrix))
  gap_count <- sum(align_matrix == "-")
  return((gap_count / total_positions) * 100)
}
