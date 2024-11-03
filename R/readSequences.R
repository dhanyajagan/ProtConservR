#' Read Protein Sequences from Fasta File
#'
#' Reads in unaligned protein sequences provided in a FASTA file. Protein
#' sequences should be of a particular protein family, like hemoglobin, insulin,
#' etc. User can use sample data located in inst/extdata/, or they can pull FASTA
#' files of protein families of their choice from InterPro website as input into
#' this function.
#'
#' @param input Character, Path to multifasta file containing protein family sequences
#' @param format Character, input format (default: "fasta")
#'
#' @return List containing processed sequences and metadata
#' @examples
#' # Reading from file: Insulin fasta
#' Insulin_Sequences <- read_sequences(system.file("extdata", "IPR004825_Insulin.fasta", package = "ProtConservR"))
#' head(Insulin_Sequences)
#'
#' # Reading from file: Dermcidin/Lacritin fasta
#' Dermcidin_Lacritin_Sequences <- read_sequences(system.file("extdata", "IPR043557_Dermcidin_Lacritin.fasta", package = "ProtConservR"))
#' head(Dermcidin_Lacritin_Sequences)
#'
#' @references
#' Pagès H, Aboyoun P, Gentleman R, DebRoy S (2024). \emph{Biostrings: Efficient manipulation of biological strings.} R package version 2.74.0,
#' \href{https://bioconductor.org/packages/Biostrings}{Link}.
#'
#' Barry J. Grant, Ana P. C. Rodrigues, Karim M. ElSawy, J. Andrew McCammon, Leo S.
#' D. Caves, Bio3d: an R package for the comparative analysis of protein structures,
#' \emph{Bioinformatics}, Volume 22, Issue 21, November 2006, Pages 2695–2696,
#' \href{https://doi.org/10.1093/bioinformatics/btl461}{link}.
#'
#' @export
#' @import Biostrings

read_sequences <- function(input, format = "fasta") {
  # Input validation
  if (!is.character(input)) {
    stop("Input must be a character string (file path)")
  }

  # Validate input format
  validate_input_format(input)

  # Read sequences from file
  if (!file.exists(input)) {
    stop("File does not exist: ", input)
  }
  sequences <- Biostrings::readAAStringSet(input, format = format)

  # Add metadata
  metadata <- list(
    input_id = input,
    date = Sys.Date(),
    n_sequences = length(sequences),
    avg_length = mean(width(sequences))
  )

  return(list(
    sequences = sequences,
    metadata = metadata
  ))
}


#' Validate Input Format
#'
#' Validate input format for above function read_sequences
#'
#' @param input Character string to validate
#'
#' @return Logical indicating if input is valid
#'
#' @export
validate_input_format <- function(input) {
  if (!grepl("\\.(fa|fasta|txt)$", input)) {
    stop("File must have .fa, .fasta, or .txt extension")
  }
}
