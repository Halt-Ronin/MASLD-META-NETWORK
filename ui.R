library(shiny)
library(readxl)
library(dplyr)
library(tidyr)
library(stringr)
library(visNetwork)
library(magick)
library(plotly)
library(shinyWidgets)
library(reactable)
library(shinyjs)


ui <- fluidPage(
  tags$head(
    tags$title("MASLD META NETWORK"),
    tags$style(HTML("
      body {
        background-color: #f9f9f9;
        color: #222222;
      }

      .sidebar {
        background-color: #ffffff;
        border-right: 2px solid #06d6a0;
        padding: 15px;
        border-radius: 0 10px 10px 0;
      }

      h2, h4 {
        color: #118ab2;
      }

      .btn-primary, .btn {
        background-color: #ef476f !important;
        color: white !important;
        border-color: #ef476f !important;
      }

      .btn-primary:hover, .btn:hover {
        background-color: #d43d66 !important;
        color: white !important;
      }

      .well {
        background-color: #f0f0f0;
        border: none;
        box-shadow: none;
      }

      .slider-animate-button {
        background-color: #ffd166 !important;
      }

      .selectize-input {
        background-color: #ffffff !important;
        color: #222 !important;
        border: 2px solid #06d6a0 !important;
        border-radius: 6px;
      }

      .selectize-dropdown {
        background-color: #ffffff !important;
        color: #222 !important;
      }

      .form-control {
        background-color: #ffffff !important;
        color: #222 !important;
        border-color: #06d6a0 !important;
      }

      table.dataTable {
        color: #222;
        background-color: #ffffff;
      }

      table.dataTable thead {
        background-color: #ffd166;
        color: #222;
      }

      table.dataTable tbody tr {
        background-color: #ffffff;
      }

      .dataTables_wrapper .dataTables_paginate .paginate_button {
        color: #ef476f !important;
      }

      #network select,
      #network input[type='text'],
      #network option {
        background-color: #ffffff !important;
        color: #222 !important;
        border: 1px solid #06d6a0 !important;
      }

      div.vis-tooltip {
        background-color: rgba(255, 255, 255, 0.95) !important;
        color: #222 !important;
        border: 1px solid #06d6a0 !important;
        font-size: 13px !important;
        padding: 8px !important;
        max-width: 300px;
        word-wrap: break-word;
        box-shadow: 0 0 10px #00000022;
      }
      
      table.dataTable {
      table-layout: fixed !important;
      width: 100% !important;
    }

    table.dataTable td {
      max-width: 300px;
      overflow: hidden;
      text-overflow: ellipsis;
      white-space: nowrap;
    }

    table.dataTable th {
      white-space: nowrap;
    }
    "))
  ),
  titlePanel(h2("MASLD META NETWORK", style = "color: #ef476f; font-weight: bold;")),
  
  tabsetPanel(
    tabPanel("Top-Score",
             sidebarLayout(
               sidebarPanel(
                 selectInput(
                   inputId = "metric_top_score",
                   label = "Histological Scoring System:",
                   choices = c("NAFLD Activity Score","Fibrosis Stage"),
                   selected = "NAFLD Activity Score"
                 ),
                 selectInput(
                   inputId = "direction_dropdown_top_score",
                   label = "Include Analysis With:",
                   choices = c("All Genes","Only Upregulated","Only Downregulated"),
                   selected = "All Genes"
                 ),
                 selectInput(
                   inputId = "category_dropdown_top_score",
                   label = "Select Category:",
                   choices = c("C2", "C3", "C4", "C5", "C6", "C7", "C8", "H"),
                   selected = "C5"
                 ),
                 sliderInput(
                   inputId = "qvalue_threshold_top_score",
                   label = "Minimum -log10(q-value):",
                   min = 1.3,   # ~0.05
                   max = 10,    # ~1e-10
                   value = 5,
                   step = 0.1
                 ),
                 selectInput(
                   inputId = "edge_filter_method_top_score",
                   label = "Cut network by:",
                   choices = c("Shared Genes","Jaccard Index"),
                   selected = "Jaccard Index"
                 ),
                 conditionalPanel(
                   condition = "input.edge_filter_method_top_score == 'Jaccard Index'",
                   sliderInput(
                     inputId = "jaccard_threshold_top_score",
                     label = "Minimum Jaccard index between pathways:",
                     min = 0,
                     max = 1,
                     step = 0.01,
                     value = 0.2
                   )
                 ),
                 conditionalPanel(
                   condition = "input.edge_filter_method_top_score == 'Shared Genes'",
                   sliderInput(
                     inputId = "shared_gene_threshold_top_score",
                     label = "Minimum number of shared genes between pathways:",
                     min = 1,
                     max = 20,
                     value = 3
                   )
                 ),
                 selectInput("clustering_method_top_score", "Clustering Method:",
                             choices = c("No Clustering", "Louvain", "Edge Betweenness", "Label Propagation"),
                             selected = "No Clustering"),
                 actionButton("open_cluster_rename_modal", "Rename Clusters"),
                 checkboxInput(
                   inputId = "include_unconnected_nodes_top_score",
                   label = "Keep only connected nodes",
                   value = FALSE  # or FALSE if you want it off by default
                 ), 
                 checkboxInput(
                   inputId = "enable_physics_top_score",
                   label = "Enable network motion",
                   value = TRUE
                 )
               ),
               
               mainPanel(
                 visNetworkOutput("top_score_network", height = "800px"),
                 div(
                   style = "width: 100%; overflow-x: auto;",
                   DT::dataTableOutput("network_summary_top_score")
                 ),
                 downloadButton("download_summary_top_score", "Download Summary Table (XLSX)"),
                 downloadButton("download_network_top_score", "Download Network as HTML")
               )
             )
    ),  
    
    tabPanel("Process-Gene",
      sidebarLayout(
        sidebarPanel(
            uiOutput("cluster_selector_process_gene"),
            sliderInput(
              inputId = "qvalue_threshold_process_gene",
              label = "Minimum -log10(q-value):",
              min = 1.3,   # ~0.05
              max = 10,    # ~1e-10
              value = 3,
              step = 0.1
            ),
            selectInput("centrality_method_process_gene", "Select Centrality Method:",
                        choices = c("Degree","Closeness","Betweenness","PageRank"),
                        selected = "Degree"),
            sliderInput("weight_total_process_gene", "Weight for total score:",
                        min = 0, max = 1, value = 0.5),
            sliderInput("threshold_process_gene", "Show genes with combined score above:",
                        min = 0, max = 1, value = 0),
            fileInput(
              inputId = "gene_logfc_input_process_gene",
              label = "Choose a CSV file",
              accept = c(".csv",".txt")  
            ),
            uiOutput("concordance_ui"),
            div(
              style = "text-align: center;",
              checkboxInput(
                inputId = "enable_physics_process_gene",
                label = "Enable network motion",
                value = TRUE
              ),
              actionButton("reset_network_process_gene", "Back to Full Network", class = "btn-primary")
             )
        ),
        mainPanel(
          visNetworkOutput("network_process_gene", height = "800px"),
          hr(),
          dataTableOutput("gene_table_process_gene")
        )
      )
    ),
    tabPanel(
      "Stage-Dependent Network",
      sidebarLayout(
        sidebarPanel(
          selectInput(
            inputId = "metric_temporal",
            label = "Histological Scoring System:",
            choices = c("NAFLD Activity Score","Fibrosis Stage"),
            selected = "NAFLD Activity Score"
          ),
          selectInput(
            inputId = "direction_dropdown_temporal",
            label = "Include Analysis With:",
            choices = c("All Genes","Only Upregulated","Only Downregulated"),
            selected = "All Genes"
          ),
          selectInput(
            inputId = "category_dropdown_temporal",
            label = "Select Category:",
            choices = c("C2","C3","C4","C5","C6","C7","C8","H"),
            selected = "C5"
          ),
          conditionalPanel(
            condition = "input.metric_temporal == 'NAFLD Activity Score'",
            sliderTextInput(
              inputId = "slider_temporal",
              label = "Choose NAS Group (Against NAS Group 0):",
              choices = c("1", "2", "3","4"),
              selected = "4"
            )
          ),
          conditionalPanel(
            condition = "input.metric_temporal == 'Fibrosis Stage'",
          sliderTextInput(
            inputId = "slider_temporal",
            label = "Choose Stage (Against Fibrosis Stages F0+F1):",
            choices = c("F2", "F3", "F4"),
            selected = "F2"
          )
          ),
          sliderInput(
            inputId = "qvalue_threshold_temporal",
            label = "Minimum -log10(q-value):",
            min = 1.3,   # ~0.05
            max = 10,    # ~1e-10
            value = 5,
            step = 0.1
          ),
          selectInput(
            inputId = "edge_filter_method_temporal",
            label = "Cut network by:",
            choices = c("Shared Genes","Jaccard Index"),
            selected = "Jaccard Index"
          ),
          conditionalPanel(
            condition = "input.edge_filter_method_temporal == 'Jaccard Index'",
            sliderInput(
              inputId = "jaccard_threshold_temporal",
              label = "Minimum Jaccard index between pathways:",
              min = 0,
              max = 1,
              step = 0.01,
              value = 0.2
            )
          ),
          conditionalPanel(
            condition = "input.edge_filter_method_temporal == 'Shared Genes'",
            sliderInput(
              inputId = "shared_gene_threshold_temporal",
              label = "Minimum number of shared genes between pathways:",
              min = 1,
              max = 20,
              value = 3
            )
          ),
          selectInput("clustering_method_temporal", "Clustering Method:",
                      choices = c("No Clustering", "Louvain", "Edge Betweenness", "Label Propagation"),
                      selected = "No Clustering"),
          #  actionButton("open_cluster_rename_modal_sex_aware", "Rename Clusters"),
          checkboxInput(
            inputId = "include_unconnected_nodes_temporal",
            label = "Keep only connected nodes",
            value = FALSE  # or FALSE if you want it off by default
          ), 
          checkboxInput(
            inputId = "enable_physics_temporal",
            label = "Enable network motion",
            value = TRUE
          ),
          checkboxInput(
            inputId = "previous_nodes",
            label = "Show lost nodes of this stage",
            value = FALSE
          )
        ),
        mainPanel(
          visNetworkOutput("temporal_net", height = "800px")
        )
    )
    ),
    tabPanel(
      "Sex-Aware Network",
      sidebarLayout(
        sidebarPanel(
          selectInput(
            inputId = "metric_sex_aware",
            label = "Histological Scoring System:",
            choices = c("NAFLD Activity Score","Fibrosis Stage"),
            selected = "NAFLD Activity Score"
          ),
          selectInput(
            inputId = "direction_dropdown_sex_aware",
            label = "Include Analysis With:",
            choices = c("All Genes","Only Upregulated","Only Downregulated"),
            selected = "All Genes"
          ),
          conditionalPanel(
            condition = "input.metric_sex_aware == 'NAFLD Activity Score'",
            selectInput(
              inputId = "category_dropdown_sex_aware",
              label = "Select Category:",
              choices = c("C2","C3","C4","C5","C6","C7","C8"),
              selected = "C5"
            )
          ),
          conditionalPanel(
            condition = "input.metric_sex_aware == 'Fibrosis Stage'",
            selectInput(
              inputId = "category_dropdown_sex_aware",
              label = "Select Category:",
              choices = c("C2","C3","C4","C5","C6","C7","C8","H"),
              selected = "C5"
            )
          ),
          sliderInput(
            inputId = "qvalue_threshold_sex_aware",
            label = "Minimum -log10(q-value):",
            min = 1.3,   # ~0.05
            max = 10,    # ~1e-10
            value = 5,
            step = 0.1
          ),
          selectInput(
            inputId = "edge_filter_method_sex_aware",
            label = "Cut network by:",
            choices = c("Shared Genes","Jaccard Index"),
            selected = "Jaccard Index"
          ),
          conditionalPanel(
            condition = "input.edge_filter_method_sex_aware == 'Jaccard Index'",
            sliderInput(
              inputId = "jaccard_threshold_sex_aware",
              label = "Minimum Jaccard index between pathways:",
              min = 0,
              max = 1,
              step = 0.01,
              value = 0.2
            )
          ),
          conditionalPanel(
            condition = "input.edge_filter_method_sex_aware == 'Shared Genes'",
            sliderInput(
              inputId = "shared_gene_threshold_sex_aware",
              label = "Minimum number of shared genes between pathways:",
              min = 1,
              max = 20,
              value = 3
            )
          ),
          selectInput("clustering_method_sex_aware", "Clustering Method:",
                      choices = c("No Clustering", "Louvain", "Edge Betweenness", "Label Propagation"),
                      selected = "No Clustering"),
        #  actionButton("open_cluster_rename_modal_sex_aware", "Rename Clusters"),
          checkboxInput(
            inputId = "include_unconnected_nodes_sex_aware",
            label = "Keep only connected nodes",
            value = FALSE  # or FALSE if you want it off by default
          ), 
          checkboxInput(
            inputId = "enable_physics_sex_aware",
            label = "Enable network motion",
            value = TRUE
          )
        ),
        mainPanel(
          visNetworkOutput("sexnet", height = "800px"),
          div(
            style = "width: 100%; overflow-x: auto;",
            DT::dataTableOutput("network_summary_sex_aware")
          ),
          downloadButton("download_summary_sex_aware", "Download Summary Table (XLSX)"),
          downloadButton("download_network_sex_aware", "Download Network as HTML")
        )
      )
    ),
    tabPanel("Transcriptome Browser",
             sidebarLayout(
               sidebarPanel(
                   selectInput(
                     inputId = "metric_browser",
                     label = "Histological Scoring System:",
                     choices = c("NAFLD Activity Score","Fibrosis Stage"),
                     selected = "NAFLD Activity Score"
                   ),
                   selectizeInput(
                     inputId = "selected_gene_browser",
                     label   = "Select Gene:",
                     choices = NULL,                # filled from server
                     multiple = FALSE,
                     options = list(placeholder = "Type a gene symbol…")
                   ),
                   fileInput("gene_txt", "Upload gene list to create dotplot (.txt)", accept = c(".txt")),
                   uiOutput("concordance_ui_browser")
               ),
               mainPanel(
                 fluidRow(
                   column(
                     width = 12,
                     h3("General Results"),
                     tags$hr(),
                     reactableOutput("gene_table_browser_pretty", height = "auto",width = "100%")
                   )
                 ),
                 br(),
                 div(
                   style = "background-color:#f8f9fa; padding:10px; border-radius:6px; 
               border:1px solid #e0e0e0; margin-bottom:10px;",
                   tags$p(
                     strong("Tip: "),
                     "Click on any section title below to expand or collapse its results."
                   )
                 ),
                 
                 tags$details(
                   tags$summary(tags$h3("Regression Results")),
                   tags$hr(),
                   uiOutput("regression_box_plot")
                 ),
                 
                 tags$details(
                   tags$summary(tags$h3("DGEA Results")),
                   tags$hr(),
                   plotlyOutput("logfc_heatmap", height = "520px")
                 ),
                 
                 tags$details(
                   tags$summary(tags$h3("Gene Disease Patterns Results")),
                   tags$hr(),
                   div(
                     style = "text-align: center; background-color: #f9f9f9;",
                     imageOutput("collage", height = "auto")
                   )
                 ),
                 
                 tags$details(
                   tags$summary(tags$h3("Gene List Dotplot")),
                   tags$hr(),
                   plotlyOutput("gene_logfc_dotplot", height = "600px")
                 )
               )
             )   
    )
  )
)
