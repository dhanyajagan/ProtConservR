library(shiny)
library(ProtConservR)
library(ggplot2)
library(plotly)
library(r3dmol)

# Define UI for random distribution app ----
ui <- fluidPage(
  titlePanel("ProtConservR: Protein Sequence Analysis"),

  sidebarLayout(

    sidebarPanel(
      tags$p("Description: This is a Shiny App that is part of the ProtConservR R
             package. Most of the functions available in the ProtConservR package
             are made available with Shiny App. ProtConservR is an R package for
             analyzing and visualizing protein sequence conservation. It provides
             functions to (1) Read and process protein sequences from FASTA files,
             (2) Create multiple sequence alignments from these protein sequence
             multifasta files using MUSCLE or ClustalW. (3) Calculate position-specific
             conservation scores based on shannon entropy. (4) Map sequence
             conservation scores to protein structure data. (5) Create plots and
             interactive visualizations of conservation data, including line plots
             of conservation scores across residues, heatmaps of conservation scores
             on protein structure data (residues), and 3D visualizations of protein
             structure colored by conservation. This app allows for the upload of
             a protein family multifasta file of choice, perform alignment of choice
             (ClustalW or MUSCLE), compute conservation scores and plot them as
             either a line plot or bar plot, and map these sequence based conservation
             scores onto protein structure (with a specified style). Note: for the
             3D protein structure visualization feature, selection cartoon results
             in the best visualization, as the other options include other colours
             for the atoms-- it can get confusing."),

      # br() element to introduce extra vertical spacing ----
      br(),
      br(),

      tags$p("Instructions: Below, select your fasta file with protein sequences
             of any protein family of your choice. You can find some on InterPro,
             or you can use the one provided: inst/extdata/IPR043557_Dermcidin_Lacritin.fasta"),

      br(),

      # File input for fasta file
      fileInput("fasta_file", "Choose Fasta File",
                accept = c(".fa", ".fasta", ".txt")),

      # Alignment method selection
      selectInput("alignment_method", "Alignment Method",
                  choices = c("MUSCLE", "ClustalW")),

      # File input for PDB file
      fileInput("pdb_file", "Choose PDB File (Optional)",
                accept = c(".pdb")),

      # Visualization type selection for conservation
      selectInput("conservation_plot_type", "Conservation Plot Type",
                  choices = c("Line Plot", "Bar Plot")),

      # 3D Structure visualization of conservation mapped onto proteins structure
      selectInput("structure_style", "3D Structure Style",
                  choices = c("Cartoon", "Sphere", "Stick")),

      # Action button to run analysis
      actionButton("run_analysis", "Run Analysis")
    ),

    mainPanel(
      tabsetPanel(
        tabPanel("Input Sequence Information",
                 verbatimTextOutput("sequence_info")),
        tabPanel("Alignment Details",
                 verbatimTextOutput("alignment_info")),
        tabPanel("Multiple Sequence Alignment (MSA) Conservation Scores",
                 uiOutput("conservation_visualization")),
        tabPanel("3D Structure Overlayed With Sequence Conservation Scores", r3dmolOutput("threeDim_conservation", width = "100%", height = "600px"))
      )
    )
  )
)

server <- function(input, output, session) {
  sequence_data <- reactiveVal(NULL)
  alignment_data <- reactiveVal(NULL)
  conservation_data <- reactiveVal(NULL)
  mapped_conservation_data <- reactiveVal(NULL)

  observeEvent(input$run_analysis, {
    req(input$fasta_file, input$alignment_method, input$pdb_file)

    tryCatch({
      # Load sequences
      seq_data <- readSequences(input$fasta_file$datapath)
      sequence_data(seq_data)

      # Show sequence info
      output$sequence_info <- renderPrint({
        cat("Number of Sequences:", seq_data$metadata$n_sequences, "\n")
        cat("Average Sequence Length:", round(seq_data$metadata$avg_length, 2), "\n")
      })

      # Align sequences
      align_data <- createAlignment(seq_data, method = input$alignment_method)
      alignment_data(align_data)

      # Show alignment info
      output$alignment_info <- renderPrint({
        str(align_data$metadata)
      })

      # Calculate conservation
      cons_data <- calculateConservation(align_data)
      conservation_data(cons_data)

      # Render conservation visualization based on selected type
      output$conservation_visualization <- renderUI({
        req(conservation_data())

        # Choose plot based on selected type
        plot <- switch(input$conservation_plot_type,
                       "Line Plot" = plotConservationScores(conservation_data(), plot_type = "line"),
                       "Bar Plot" = plotConservationScores(conservation_data(), plot_type = "bar")
        )

        # Render plot appropriately
        if (input$conservation_plot_type %in% c("Line Plot", "Bar Plot")) {
          plotOutput("static_conservation_plot", height = "500px")
        } else {
          plotlyOutput("interactive_conservation_plot", height = "500px")
        }
      })

      # Render static plot for line/bar options
      output$static_conservation_plot <- renderPlot({
        if (input$conservation_plot_type == "Line Plot") {
          plotConservationScores(conservation_data(), plot_type = "line")
        } else {
          plotConservationScores(conservation_data(), plot_type = "bar")
        }
      })


      # Map conservation to structure
      struct_data <- readStructure(input$pdb_file$datapath)
      mapped_data <- mapConservation(cons_data, struct_data)
      mapped_conservation_data(mapped_data)

      # Render 3D structure visualization
      output$threeDim_conservation <- renderR3dmol({
        threeDimStructureConservation(
          input_pdb = input$pdb_file$datapath,
          mapped_data = mapped_data,
          style = tolower(input$structure_style)
        )
      })

    }, error = function(e) {
      showNotification(paste("Error during analysis:", e$message), type = "error")
    })
  })
}

shinyApp(ui, server)


# [END]
