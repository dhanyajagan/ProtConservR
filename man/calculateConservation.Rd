% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/calculateConservation.R
\name{calculateConservation}
\alias{calculateConservation}
\title{Calculate Position-Specific Conservation Scores}
\usage{
calculateConservation(alignment_data)
}
\arguments{
\item{alignment_data}{List containing alignment and metadata (output from createAlignment() function)}
}
\value{
List containing:
  \itemize{
    \item scores: Numeric vector of conservation scores for each position (0-1 scale)
    \item metadata: List with calculation parameters and summary statistics
  }
}
\description{
This function calculates conservation scores for each position in a multiple sequence
alignment using Shannon's entropy. Lower entropy values indicate higher conservation
at that position. The scores are then normalized to a 0-1 scale where 1 indicates
perfect conservation.
}
\examples{
# Read sequences, create alignment, and calculate conservation scores.
Dermcidin_Lacritin_Sequences <- readSequences(system.file("extdata", "IPR043557_Dermcidin_Lacritin.fasta",
                          package = "ProtConservR"))
Dermcidin_Lacritin_Alignment <- createAlignment(Dermcidin_Lacritin_Sequences, method = "ClustalW")

# Calculate conservation scores
Dermcidin_Lacritin_Conservation <- calculateConservation(Dermcidin_Lacritin_Alignment)

}
\references{
Valdar, W.S.J. (2002). Scoring residue conservation. \emph{Proteins: Structure,
Function, and Bioinformatics}, 48(2), 227-241.

Barry J. Grant, Ana P. C. Rodrigues, Karim M. ElSawy, J. Andrew McCammon, Leo S.
D. Caves, Bio3d: an R package for the comparative analysis of protein structures,
\emph{Bioinformatics}, Volume 22, Issue 21, November 2006, Pages 2695–2696,
\href{https://doi.org/10.1093/bioinformatics/btl461}{link}.

Strait BJ, Dewey TG. The Shannon information entropy of protein sequences.
Biophys J. 1996 Jul;71(1):148-55. doi: 10.1016/S0006-3495(96)79210-X.
PMID: 8804598; PMCID: PMC1233466.
\href{https://pmc.ncbi.nlm.nih.gov/articles/PMC1233466/}{link}
}
