
<!-- README.md is generated from README.Rmd. Please edit that file -->

# ProtConservR

R package for Analyzing and Visualizing Protein Sequence Conservation.
<!-- badges: start --> <!-- badges: end -->

## Description

`ProtConservR` is an R package for analyzing and visualizing protein
sequence conservation. It provides functions to (1) Read and process
protein sequences from FASTA files, (2) Create multiple sequence
alignments from these protein sequence multifasta files using MUSCLE or
ClustalW. (3) Calculate position-specific conservation scores based on
shannon entropy. (4) Map sequence conservation scores to protein
structure data. (5) Create plots and interactive visualizations of
conservation data, including line plots of conservation scores across
residues, heatmaps of conservation scores on protein structure data
(residues), and 3D visualizations of protein structure colored by
conservation.

## Installation

You can install the latest version of ProtConservR like so: <br> <br>

``` r
install.packages("devtools")
library("devtools")
devtools::install_github("dhanyajagan/ProtConservR", build_vignettes = TRUE)
library("ProtConservR")
```

To run the Shiny app: `Under construction`

## Overview

<br> <br> <br>

``` r
ls("package:ProtConservR")
data(package = "ProtConservR") 
browseVignettes("ProtConservR")
```

`ProtConservR` contains 8 functions.

1.  ***readSequences*** Reads in unaligned protein sequences provided in
    a multifasta file. Protein sequences should be of a particular
    protein family, like hemoglobin, insulin,etc. User can use sample
    data located in inst/extdata/, or they can pull multifasta files of
    protein families of their choice from InterPro website as input into
    this function.

2.  ***createAlignment*** Performs multiple sequence alignment on
    protein sequences using either MUSCLE or ClustalW algorithm with
    default settings. Takes processed sequences from readSequences()
    function as input.

3.  ***calculateConservation*** This function calculates conservation
    scores for each position in a multiple sequence alignment using
    Shannon’s entropy. The scores are then normalized to a 0-1 scale
    where 1 indicates perfect conservation.

4.  ***plotConservationScores*** This function creates a plot of
    conservation scores for each position in a multiple sequence
    alignment (plotting sequence conservation scores).

5.  ***readStructure*** This function reads and processes a protein
    structure from a specified PDB file. Users must provide a file path
    to a PDB file.

6.  ***mapConservation*** This function maps sequence conservation
    scores onto protein structure by aligning the structural residues
    with the corresponding positions in the multiple sequence alignment.
    It handles potential gaps and mismatches between the structure and
    alignment in a very simple way.

7.  ***plotInteractiveConservation*** Creates an interactive
    visualization of protein conservation scores using plotly. This
    function takes as input the output of mapConservation() and
    generates multiple visualization types including a heatmap view of
    the conservation scores, a line plot showing conservation across
    residues, and summary statistics.

8.  ***threeDimStructureConservation*** Creates an interactive 3D
    visualization of the pdb protein structure colored by conservation
    scores. Takes as input the representative pdb file and output from
    mapConservation and uses r3dmol to create an interactive viewer
    where conservation scores are mapped to a color gradient from red
    (low conservation) to white (medium) to blue (high conservation).

## Contributions

## References

## Acknowledgements

This package was developed as part of an assessment for 2024 BCB410H:
Applied Bioinformatics course at the University of Toronto, Toronto,
CANADA. `ProtConservR` welcomes issues, enhancement requests, and other
contributions. To submit an issue, use the[GitHub
issues](https://github.com/dhanyajagan/ProtConservR/issues)
