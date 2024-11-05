#' Read Protein Structure from PDB File
#'
#' This function reads and processes a protein structure from a specified PDB file.
#' Users must provide a file path to a PDB file.
#'
#' @param pdb_input Character. A file path to a PDB file.
#'
#' @return A list containing processed structural information, including atomic coordinates, chain IDs, and residue data.
#'
#' @examples
#' # Example with local PDB file: Insulin
#' Insulin_structure <- readStructure(system.file("extdata", "Insulin.pdb", package = "ProtConservR"))
#'
#' # Example with local PDB file: Dermcidin/Lacritin
#' Dermcidin_Lacritin_structure <- readStructure(system.file("extdata", "Dermcidin.pdb", package = "ProtConservR"))
#'
#' @references
#' #' Barry J. Grant, Ana P. C. Rodrigues, Karim M. ElSawy, J. Andrew McCammon, Leo S.
#' D. Caves, Bio3d: an R package for the comparative analysis of protein structures,
#' \emph{Bioinformatics}, Volume 22, Issue 21, November 2006, Pages 2695–2696,
#' \href{https://doi.org/10.1093/bioinformatics/btl461}{link}.
#'
#' @export
#' @importFrom bio3d read.pdb
readStructure <- function(pdb_input) {
  # Input validation
  if (!is.character(pdb_input)) {
    stop("PDB input must be a character string (file path to PDB file)")
  }

  # Check if the provided path is valid
  if (!file.exists(pdb_input)) {
    stop("File does not exist: ", pdb_input)
  }

  # Read the PDB structure from the provided file path
  structure <- bio3d::read.pdb(pdb_input)

  # Process and extract relevant data from structure
  structure_data <- list(
    atoms = structure$atom,
    residues = unique(structure$atom$resid),
    chains = unique(structure$atom$chain),
    n_atoms = nrow(structure$atom),
    n_residues = length(unique(structure$atom$resid)),
    date = Sys.Date()
  )

  return(structure_data)
}
