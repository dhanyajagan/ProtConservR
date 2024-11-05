# tests/testthat/test-calculateConservation.R
library(testthat)
library(ProtConservR)

test_that("calculateConservation function works as expected", {
  # Load sample data
  seq_data <- readSequences(system.file("extdata", "IPR043557_Dermcidin_Lacritin.fasta", package = "ProtConservR"))
  alignment_data <- createAlignment(seq_data, method = "ClustalW")

  # Test conservation calculation
  conservation_data <- calculateConservation(alignment_data)
  expect_true(is.list(conservation_data))
  expect_true("scores" %in% names(conservation_data))
  expect_true("metadata" %in% names(conservation_data))
  expect_equal(length(conservation_data$scores), alignment_data$metadata$alignment_length)
  expect_true(all(conservation_data$scores >= 0 & conservation_data$scores <= 1))

  # Test error handling
  expect_error(calculateConservation(NULL), "Input must be a list containing 'alignment' and 'metadata'")
})
