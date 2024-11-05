# tests/testthat/test-mapConservation.R
library(testthat)
library(ProtConservR)

test_that("mapConservation function works as expected", {
  # Load sample data
  seq_data <- readSequences(system.file("extdata", "IPR043557_Dermcidin_Lacritin.fasta", package = "ProtConservR"))
  alignment_data <- createAlignment(seq_data, method = "ClustalW")
  conservation_data <- calculateConservation(alignment_data)
  structure_data <- readStructure(system.file("extdata", "Dermcidin.pdb", package = "ProtConservR"))

  # Test mapping conservation to structure
  mapped_data <- mapConservation(conservation_data, structure_data)
  expect_true(is.list(mapped_data))
  expect_true("mapped_scores" %in% names(mapped_data))
  expect_true("structure" %in% names(mapped_data))
  expect_true("metadata" %in% names(mapped_data))

  # Test get_residue_conservation helper function
  expect_true(is.numeric(get_residue_conservation(mapped_data, 1)))
  expect_true(is.numeric(get_residue_conservation(mapped_data, 1, chain = structure_data$chains[1])))
  expect_true(is.na(get_residue_conservation(mapped_data, 999999)))

  # Test error handling
  expect_error(mapConservation(NULL, structure_data), "Conservation data must be a list")
  expect_error(mapConservation(conservation_data, NULL), "Structure data must be a list")
})
