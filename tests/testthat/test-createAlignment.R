# tests/testthat/test-createAlignment.R
library(testthat)
library(ProtConservR)

test_that("createAlignment function works as expected", {
  # Load sample data
  seq_data <- readSequences(system.file("extdata", "IPR043557_Dermcidin_Lacritin.fasta", package = "ProtConservR"))

  # Test MUSCLE alignment
  muscle_alignment <- createAlignment(seq_data, method = "MUSCLE")
  expect_true(is.list(muscle_alignment))
  expect_true("alignment" %in% names(muscle_alignment))
  expect_true("metadata" %in% names(muscle_alignment))
  expect_s4_class(muscle_alignment$alignment, "AAStringSet")

  # Test ClustalW alignment
  clustal_alignment <- createAlignment(seq_data, method = "ClustalW")
  expect_true(is.list(clustal_alignment))
  expect_true("alignment" %in% names(clustal_alignment))
  expect_true("metadata" %in% names(clustal_alignment))
  expect_s4_class(clustal_alignment$alignment, "AAStringSet")

})
