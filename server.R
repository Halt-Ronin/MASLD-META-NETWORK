# Packages

library(shiny)
library(visNetwork)
library(dplyr)
library(stringr)
library(readxl)
library(tidyr)
library(igraph)
library(writexl)
library(htmlwidgets)
library(tidyverse)
library(scales)
library(plotly)
library(DT)
library(reactable)


server <- function(input, output, session) {
  
  # Data load 
  
  scores <- read_xlsx("overall_meta_analysis_overview_directions.xlsx")
  scores_fibrosis <- read_xlsx("overall_meta_analysis_overview_directions_fibrosis.xlsx") 
  orthologs_mouse_fibrosis <- read_xlsx("fibrosis_top_score_mouse.xlsx")
  orthologs_mouse_nas <- read_xlsx("nas_top_score_mouse.xlsx")
  orthologs_zebrafish_fibrosis <- read_xlsx("fibrosis_top_score_zebrafish.xlsx")
  orthologs_zebrafish_nas <- read_xlsx("nas_top_score_zebrafish.xlsx")
  orthologs_mouse_fibrosis <- orthologs_mouse_fibrosis %>%
    rename(`Mouse Symbol` = `Output Gene Symbol`)
  orthologs_mouse_nas <- orthologs_mouse_nas %>%
    rename(`Mouse Symbol` = `Output Gene Symbol`)
  orthologs_zebrafish_fibrosis <- orthologs_zebrafish_fibrosis %>%
    rename(`Zebrafish Symbol` = `Output Gene Symbol`)
  orthologs_zebrafish_nas <- orthologs_zebrafish_nas %>%
    rename(`Zebrafish Symbol` = `Output Gene Symbol`)
  string_nas <- read_xlsx("string_edges_detailed_nas.xlsx")
  string_fibrosis <-  read_xlsx("string_edges_detailed_fibrosis.xlsx")
  tcga_nas <-  read_xlsx("nas_tcga.xlsx")
  tcga_fibrosis <-  read_xlsx("fibrosis_tcga.xlsx")

  # words that mean "this first line is a header, not a gene"
  gene_header_words <- c("gene", "genes", "symbol", "symbols", "gene_symbol",
                         "gene symbol", "mapped_gene", "mapped gene", "id", "x")

  # Reads an uploaded gene list (.csv or .txt) and works out on its own whether the
  # first line is a header. Without this, .csv files without a header lose their first
  # gene (read.csv assumes header = TRUE) and .txt files with a header keep the header
  # as a gene and read the logFC column as text (read.table assumed header = FALSE).
  read_user_gene_file <- function(path, name) {
    ext <- tolower(sub(".*\\.", "", name))
    validate(need(ext %in% c("csv", "txt"), "Unsupported file type. Please upload .csv or .txt"))
    sep <- if (ext == "csv") "," else ""

    first <- read.table(path, sep = sep, nrows = 1, header = FALSE,
                        stringsAsFactors = FALSE, quote = "\"", comment.char = "")
    has_header <- if (ncol(first) >= 2) {
      # a header row has a non-numeric second column (e.g. "logFC")
      is.na(suppressWarnings(as.numeric(first[[2]][1])))
    } else {
      tolower(trimws(as.character(first[[1]][1]))) %in% gene_header_words
    }

    df <- read.table(path, sep = sep, header = has_header,
                     stringsAsFactors = FALSE, quote = "\"", comment.char = "")
    df[[1]] <- trimws(as.character(df[[1]]))
    if (ncol(df) >= 2) df[[2]] <- suppressWarnings(as.numeric(df[[2]]))
    df <- df[!is.na(df[[1]]) & nzchar(df[[1]]), , drop = FALSE]
    validate(need(nrow(df) > 0, "No gene found in the uploaded file."))
    df
  }

  # Guesses which species the uploaded symbols belong to and translates them into the
  # human symbols the app displays. `reference_hits` is how many uploaded symbols are
  # already present in the human data being shown: the switch to an ortholog table only
  # fires when that table matches *more* symbols than that. Without this guard a single
  # coincidence is enough to flip a human list (a few mouse symbols are spelled exactly
  # like human ones: C2/F11 for NAS, C3/C7/F3 for Fibrosis) and the whole list is then
  # replaced by that one gene.
  resolve_user_species <- function(user_genes, metric, reference_hits) {
    if (metric == "NAFLD Activity Score") {
      ortho_mouse <- orthologs_mouse_nas
      ortho_fish  <- orthologs_zebrafish_nas
    } else {
      ortho_mouse <- orthologs_mouse_fibrosis
      ortho_fish  <- orthologs_zebrafish_fibrosis
    }
    mouse_hits <- sum(user_genes %in% ortho_mouse$`Mouse Symbol`)
    fish_hits  <- sum(user_genes %in% ortho_fish$`Zebrafish Symbol`)

    if (mouse_hits > reference_hits & fish_hits == 0) {
      keep <- ortho_mouse$`Mouse Symbol` %in% user_genes
      list(species = "mouse", table = ortho_mouse, key = "Mouse Symbol",
           original = ortho_mouse$`Mouse Symbol`[keep],
           genes    = ortho_mouse$`Search Term`[keep])
    } else if (mouse_hits == 0 & fish_hits > reference_hits) {
      keep <- ortho_fish$`Zebrafish Symbol` %in% user_genes
      list(species = "zebrafish", table = ortho_fish, key = "Zebrafish Symbol",
           original = ortho_fish$`Zebrafish Symbol`[keep],
           genes    = ortho_fish$`Search Term`[keep])
    } else {
      # human list (or no upload at all): nothing to translate
      list(species = "human", table = NULL, key = NULL,
           original = user_genes, genes = user_genes)
    }
  }

  # Translates uploaded symbols into human symbols using an already resolved species.
  # Identity for human lists.
  translate_user_genes <- function(genes, species_info) {
    if (is.null(species_info$table)) return(genes)
    species_info$table$`Search Term`[species_info$table[[species_info$key]] %in% genes]
  }

  # Short recap shown under each upload widget so the user can see what was actually kept.
  upload_summary_msg <- function(n_uploaded, species, n_kept) {
    if (!n_uploaded) return("")
    sprintf("Uploaded: %d genes | detected species: %s | kept: %d | not found here: %d",
            n_uploaded, species, n_kept, n_uploaded - n_kept)
  }
  ################################################ Top-Score Over Representation Analysis ###########################
  
  summary_table_top_score <- reactiveVal(NULL) #reactive summary_table for the first tab
  vis_network_top_score <- reactiveVal(NULL) #reactive for downloading network as HTML
  upload_msg_process_gene <- reactiveVal("") #recap of the uploaded gene list, shown under the upload widget
  upload_msg_string <- reactiveVal("")

  output$upload_msg_process_gene <- renderText(upload_msg_process_gene())
  output$upload_msg_string <- renderText(upload_msg_string())
  cluster_rename_map_top_score <- reactiveVal(list()) #reactive list for renaming clusters
  
  observeEvent(input$use_paper_data_top_score, {
    
    # Only apply paper settings when checkbox is checked
    if (isTRUE(input$use_paper_data_top_score)) {
      
      updateSelectInput(
        session,
        "clustering_method_top_score",
        selected = "Louvain"
      )
      
      updateSelectInput(
        session,
        "direction_dropdown_top_score",
        selected = "All Genes"
      )
      
      updateSelectInput(
        session,
        "category_dropdown_top_score",
        selected = "C5"
      )
      
      updateSliderInput(
        session,
        "qvalue_threshold_top_score",
        value = 1.3
      )
      
      updateSelectInput(
        session,
        "edge_filter_method_top_score",
        selected = "Jaccard Index"
      )
      
      updateSliderInput(
        session,
        "jaccard_threshold_top_score",
        value = 0.2
      )
    }
  })
  
  selected_data_top_score <- reactive({
    req(input$category_dropdown_top_score, selected_data_top_score_xlsx()) ## require category_dropdown_top_score from the dropdown menu 
    if(input$metric_top_score == "NAFLD Activity Score"){
      if(input$direction_dropdown_top_score == "All Genes"){
        filename_rds <- paste0("precalculated_edges_overview_95_msigdb_", input$category_dropdown_top_score, ".rds") ## load the selected MSigDB result
        filepath_rds <- file.path("precalculated_edge_ora", filename_rds) ## get the file path for the file
      }else if(input$direction_dropdown_top_score == "Only Upregulated"){
        filename_rds <- paste0("precalculated_edges_overview_95_msigdb_", input$category_dropdown_top_score, "_upregulated.rds") ## load the selected MSigDB result
        filepath_rds <- file.path("precalculated_edge_ora", filename_rds) ## get the file path for the file
      }else{
        filename_rds <- paste0("precalculated_edges_overview_95_msigdb_", input$category_dropdown_top_score, "_downregulated.rds") ## load the selected MSigDB result
        filepath_rds <- file.path("precalculated_edge_ora", filename_rds) ## get the file path for the file
      } 
    }else{
      if(input$direction_dropdown_top_score == "All Genes"){
        filename_rds <- paste0("precalculated_edges_overview_95_msigdb_", input$category_dropdown_top_score, ".rds") ## load the selected MSigDB result
        filepath_rds <- file.path("precalculated_edge_ora/fibrosis", filename_rds) ## get the file path for the file
      }else if(input$direction_dropdown_top_score == "Only Upregulated"){
        filename_rds <- paste0("precalculated_edges_overview_95_msigdb_", "upregulated_", input$category_dropdown_top_score, ".rds") ## load the selected MSigDB result
        filepath_rds <- file.path("precalculated_edge_ora/fibrosis", filename_rds) ## get the file path for the file
      }else{
        filename_rds <- paste0("precalculated_edges_overview_95_msigdb_", "downregulated_", input$category_dropdown_top_score, ".rds") ## load the selected MSigDB result
        filepath_rds <- file.path("precalculated_edge_ora/fibrosis", filename_rds) ## get the file path for the file
      } 
    }
    
    if (!file.exists(filepath_rds)) {
      validate(need(FALSE, paste("No file found for", input$category_dropdown_top_score)))
      return(NULL)  # never actually reached because validate() stops
    }
    
    df_rds <- readRDS(filepath_rds) ## get the file
    
    df_xlsx <- selected_data_top_score_xlsx() #get msigdb output
    
    df_rds <- df_rds[df_rds$from %in% df_xlsx$ID, ] #filter for q value since some of msigdb will be eliminated
    df_rds <- df_rds[df_rds$to %in% df_xlsx$ID,] #filter for q value since some of msigdb will be eliminated
    
    return(df_rds)
  })
  
  selected_data_top_score_xlsx <- reactive({
    req(input$category_dropdown_top_score) ## require category_dropdown_top_score from the dropdown menu 
    if(input$metric_top_score == "NAFLD Activity Score"){
      if(input$direction_dropdown_top_score == "All Genes"){
        filename_xlsx <-  paste0("overview_95_msigdb_", input$category_dropdown_top_score, ".xlsx") ## load the selected MSigDB result
        filepath_xlsx <- file.path("top_score_ora", filename_xlsx) ## get the file path for the file
      }else if(input$direction_dropdown_top_score == "Only Upregulated"){
        filename_xlsx <- paste0("overview_95_msigdb_", input$category_dropdown_top_score, "_upregulated.xlsx") ## load the selected MSigDB result
        filepath_xlsx <- file.path("top_score_ora", filename_xlsx) ## get the file path for the file
      }else{
        filename_xlsx <- paste0("overview_95_msigdb_", input$category_dropdown_top_score, "_downregulated.xlsx") ## load the selected MSigDB result
        filepath_xlsx <- file.path("top_score_ora", filename_xlsx) ## get the file path for the file
      } 
    }else{
      if(input$direction_dropdown_top_score == "All Genes"){
        filename_xlsx <-  paste0("fibrosis_overview_95_msigdb_", input$category_dropdown_top_score, ".xlsx") ## load the selected MSigDB result
        filepath_xlsx <- file.path("top_score_ora/fibrosis", filename_xlsx) ## get the file path for the file
      }else if(input$direction_dropdown_top_score == "Only Upregulated"){
        filename_xlsx <- paste0("fibrosis_overview_95_msigdb_","upregulated_", input$category_dropdown_top_score, ".xlsx") ## load the selected MSigDB result
        filepath_xlsx <- file.path("top_score_ora/fibrosis", filename_xlsx) ## get the file path for the file
      }else{
        filename_xlsx <- paste0("fibrosis_overview_95_msigdb_", "downregulated_", input$category_dropdown_top_score, ".xlsx") ## load the selected MSigDB result
        filepath_xlsx <- file.path("top_score_ora/fibrosis", filename_xlsx) ## get the file path for the file
      } 
    }
    if (!file.exists(filepath_xlsx)) {
      validate(need(FALSE, paste("No file found for", input$category_dropdown_top_score)))
      return(NULL)  # never actually reached because validate() stops
    }
    
    df_xlsx <- read_xlsx(filepath_xlsx)
    req("qvalue" %in% colnames(df_xlsx))  # Check qvalue exists
    
    
    df_xlsx <- df_xlsx %>%
      filter(-log10(qvalue) >= input$qvalue_threshold_top_score) # Filter by threshold on -log10(qvalue)
    
    return(df_xlsx)
  })
  
  
  observeEvent(selected_data_top_score(), {
    df <- selected_data_top_score() #get the current rds
    if (is.null(df) || nrow(df) < 2) return() # return null if empty or 1 patheay
    max_weight <- max(df$weight, na.rm = TRUE) #get the maximum overlap
    updateSliderInput(
      session,
      inputId = "shared_gene_threshold_top_score",
      min = 1,
      max = max_weight, #update maximum of the shared gene slider to maximum current overlap.
      value = min(3, max_weight)
    )
  })
  
  
  output$top_score_network <- renderVisNetwork({
    req(selected_data_top_score(),selected_data_top_score_xlsx()) ## require selected_data_top_score reactive value
    
    edges <- selected_data_top_score() ## assign to a variable
    df_top_score <- selected_data_top_score_xlsx() ## assign to a variable
    df_top_score$gene_list <- strsplit(df_top_score$geneID, "/") ## seperate genes to get them as list
    
    if (nrow(edges) < 2) {
      showNotification("Not enough pathways to compute edges.", type = "warning")
      return(NULL)
    } #make sure there is enough data after filtering
    
    if (input$edge_filter_method_top_score == "Jaccard Index") {
      edges$title <- paste0(
        "Jaccard Index: ", round(edges$jaccard, 3),
        "\nUpregulated: ", edges$up,
        "\nDownregulated: ", edges$down
      ) #titles for jaccard option
    } else {
      edges$title <- paste0(
        edges$weight, " Shared genes",
        "\nUpregulated: ", edges$up,
        "\nDownregulated: ", edges$down
      )
    } #titles for shared gene option
    
    if (input$edge_filter_method_top_score == "Shared Genes") {
      edges <- edges[edges$weight >= input$shared_gene_threshold_top_score, ]
    } else if (input$edge_filter_method_top_score == "Jaccard Index") {
      edges <- edges[edges$jaccard >= input$jaccard_threshold_top_score, ]
    } #filterings forr edge cut methods
    
    df_top_score$qvalue[df_top_score$qvalue == 0] <- 1e-300 # Avoid log10(0)
    
    if (input$include_unconnected_nodes_top_score) {
      used_ids <- unique(c(edges$from, edges$to))
      df_top_score <- df_top_score[df_top_score$ID %in% used_ids, ]
    } #if checkbox is checked there is no unconnected nodes (I know the name indicates opposite I will fix the naming but it works)
    
    gene_counts <- lengths(df_top_score$gene_list)
    
    
    node_sizes <- -log10(df_top_score$qvalue) #node size = log transformed q value
    node_sizes_scaled <- scales::rescale(node_sizes, to = c(10, 40))  # adjust range if needed
    
    # ============================================================
    # CLUSTERING
    # ============================================================
    
    if (nrow(edges) > 0 && input$clustering_method_top_score != "No Clustering") {
      
      # ----------------------------------------------------------
      # PAPER CLUSTERING
      #
      # When the paper preset button has been clicked, the UI still
      # displays "Louvain", but the previously curated clustering
      # from the paper is used instead of recalculating Louvain.
      # ----------------------------------------------------------
      
      if (
        isTRUE(input$use_paper_data_top_score) &&
        input$clustering_method_top_score == "Louvain"
      ) {
        
        # Select clustering file according to histological score
        if (input$metric_top_score == "NAFLD Activity Score") {
          
          paper_cluster_file <- file.path(
            "paper_clustering",
            "NAS_top_score.xlsx"
          )
          
        } else {
          
          paper_cluster_file <- file.path(
            "paper_clustering",
            "fibrosis_top_score.xlsx"
          )
        }
        
        # Make sure file exists
        validate(
          need(
            file.exists(paper_cluster_file),
            paste("Paper clustering file not found:", paper_cluster_file)
          )
        )
        
        # Read clustering generated for the paper
        paper_clusters <- readxl::read_xlsx(
          paper_cluster_file
        )
        
        # Make sure required columns exist
        validate(
          need(
            all(c("Pathway ID", "Cluster Group") %in% colnames(paper_clusters)),
            "Paper clustering file must contain 'Pathway ID' and 'Cluster Group' columns."
          )
        )
        
        # Keep only pathway ID and its predefined paper cluster
        cluster_df <- paper_clusters %>%
          dplyr::transmute(
            id = as.character(`Pathway ID`),
            group = as.character(`Cluster Group`)
          ) %>%
          dplyr::distinct(id, .keep_all = TRUE)
        
        # Only retain pathways present in the current network
        cluster_df <- cluster_df %>%
          dplyr::filter(id %in% df_top_score$ID)
        
        # Add any pathway that is unexpectedly absent from the
        # paper clustering file as "Unclustered"
        missing_ids <- setdiff(
          df_top_score$ID,
          cluster_df$id
        )
        
        if (length(missing_ids) > 0) {
          
          cluster_df <- dplyr::bind_rows(
            cluster_df,
            data.frame(
              id = missing_ids,
              group = "Unclustered",
              stringsAsFactors = FALSE
            )
          )
        }
        
      } else {
        
        # --------------------------------------------------------
        # NORMAL APP CLUSTERING
        # --------------------------------------------------------
        
        edge_mat <- as.matrix(
          edges[, c("from", "to")]
        )
        
        g <- igraph::graph_from_edgelist(
          edge_mat,
          directed = FALSE
        )
        
        cluster_result <- switch(
          input$clustering_method_top_score,
          
          "Louvain" =
            igraph::cluster_louvain(g),
          
          "Edge Betweenness" =
            igraph::cluster_edge_betweenness(g),
          
          "Label Propagation" =
            igraph::cluster_label_prop(g)
        )
        
        membership <- igraph::membership(
          cluster_result
        )
        
        cluster_df <- data.frame(
          id = names(membership),
          group = paste0("Cluster ", membership),
          stringsAsFactors = FALSE
        )
      }
      
    } else {
      
      if (nrow(df_top_score) > 0) {
        
        cluster_df <- data.frame(
          id = df_top_score$ID,
          group = "Unclustered",
          stringsAsFactors = FALSE
        )
      }
    }
    
    if(nrow(df_top_score) == 0){
      return(NULL)
    }
    
    network_nodes <- unique(data.frame(
      id = df_top_score$ID,
      label = df_top_score$ID,
      size = node_sizes_scaled,
      title = paste(gene_counts, "genes in pathway"),
      stringsAsFactors = FALSE
    )) #create network_nodes
    
    network_nodes <- merge(network_nodes, cluster_df, by = "id", all.x = TRUE) #merge with clustering info previously obtained
    
    summary_table <- data.frame(
      `Pathway ID` = network_nodes$id,
      `Gene IDs` = df_top_score$geneID[match(network_nodes$id, df_top_score$ID)],
      `Cluster Group` = network_nodes$group,
      `q-value` = df_top_score$qvalue[match(network_nodes$id, df_top_score$ID)],
      `-log10(q)` = node_sizes[match(network_nodes$id, df_top_score$ID)],
      `Gene Count` = gene_counts[match(network_nodes$id, df_top_score$ID)],
      stringsAsFactors = FALSE
    ) #create the summary table
    
    # Save the original cluster labels before any renaming
    summary_table$`Original Cluster Group` <- network_nodes$group #put group information as well but in original vers.
    
    # Apply renaming to new 'Cluster Group' column
    rename_map <- cluster_rename_map_top_score()
    
    summary_table$`Cluster.Group` <- vapply(
      summary_table$`Original Cluster Group`,
      function(cl) {
        if (!is.null(rename_map[[cl]]) && nzchar(rename_map[[cl]])) {
          rename_map[[cl]]
        } else {
          cl
        }
      },
      character(1)
    ) #this part maps the new names to the cluster names assigned by default.
    
    network_nodes$group <- summary_table$`Cluster.Group` #change the group
    
    # Save table for reuse in UI/download
    summary_table_top_score(summary_table) #reactive
    updateCheckboxInput(session, "enable_physics_top_score", value = TRUE) #checkbox updater when there is a new network
    
    edges$color = ifelse( edges$color == 'green', "#D55E00", "#0072B2")
    network_object <- visNetwork(network_nodes, edges, height = "100%", width = "100%") %>%
      visNodes(shape = "dot", size = "size") %>%
      visEdges(smooth = FALSE) %>%
      visOptions(highlightNearest = TRUE, nodesIdSelection = TRUE, selectedBy = "group") %>%
      visOptions(highlightNearest = TRUE, nodesIdSelection = TRUE, manipulation = FALSE,selectedBy = "group") %>%
      visInteraction(dragView = TRUE, zoomView = TRUE, navigationButtons = TRUE) %>%
      visLayout(randomSeed = 123) %>%
      visPhysics(
        solver = "forceAtlas2Based",
        forceAtlas2Based = list(
          gravitationalConstant = -50,
          centralGravity = 0.001,
          springLength = 200,
          springConstant = 0.02,
          damping = 0.4,
          avoidOverlap = 1
        ),
        stabilization = list(
          enabled = TRUE,
          iterations = 200
        )
      ) #visualizes the network including its properties (physics etc.)
    
    # Save for download
    vis_network_top_score(network_object) #make network object reactive
    
    # Render to UI
    network_object #visualizes the network
  })
  
  observeEvent(input$enable_physics_top_score, {
    visNetworkProxy("top_score_network") %>%
      visPhysics(enabled = input$enable_physics_top_score)
  }) #if checked stop motion
  
  observeEvent({
    input$direction_dropdown_top_score
    input$category_dropdown_top_score
  }, {
    
    # Reset q-value threshold to default (e.g., 1.3 means q < 0.05)
    updateSliderInput(session, "qvalue_threshold_top_score", value = 5)
    
    # Reset edge filtering thresholds
    updateSliderInput(session, "shared_gene_threshold_top_score", value = 3)
    updateSliderInput(session, "jaccard_threshold_top_score", value = 0.2)
    updateSliderInput(session, "score_threshold_top_score", value = 0)
    
    # Reset edge filter method
    updateSelectInput(session, "edge_filter_method_top_score", selected = "Jaccard Index")
    
    # Reset clustering dropdown
    updateSelectInput(session, "clustering_method_top_score", selected = "No Clustering")
    
    # Optionally reset include_unconnected_nodes
    updateCheckboxInput(session, "include_unconnected_nodes_top_score", value = FALSE)
    
    updateCheckboxInput(session, "enable_physics_top_score", value = TRUE)
    updateCheckboxInput(
      session,
      "use_paper_data_top_score",
      value = FALSE
    )  }) #update everything to default when there is 
  
  observeEvent(input$metric_top_score, {
    
    updateCheckboxInput(
      session,
      "use_paper_data_top_score",
      value = FALSE
    )
    
  }, ignoreInit = TRUE)
  
  output$network_summary_top_score <- DT::renderDataTable({
    req(summary_table_top_score())
    
    df <- summary_table_top_score()
    df <- df[,-c(2,7)]
    colnames(df) <- c("Pathway ID", "Cluster Group", "q-value", "-log10(q)", "Gene Count")
    colnames(df) <- gsub("\\.", " ", colnames(df))  # Fix column names for display
    
    DT::datatable(
      df,
      options = list(
        pageLength = 10,
        autoWidth = TRUE,
        columnDefs = list(
          list(targets = 1, width = '200px', className = 'dt-wrap')
        )
      ),
      rownames = FALSE,
      escape = FALSE
    ) %>% DT::formatStyle(
      columns = 1,
      `white-space` = "normal",
      `word-wrap` = "break-word"
    )
  }) #visualize table output below network
  
  output$download_summary_top_score <- downloadHandler(
    filename = function() {
      paste0("Pathway_Summary_", Sys.Date(), ".xlsx")
    },
    content = function(file) {
      df <- summary_table_top_score()
      df <- df[,-c(2,7)]
      colnames(df) <- c("Pathway ID", "Cluster Group", "q-value", "-log10(q)", "Gene Count")
      colnames(df) <- gsub("\\.", " ", colnames(df))  # Fix column names for display
      writexl::write_xlsx(df, file)
    }
  ) #download the network clusterings etc. as a table.
  
  output$download_network_top_score <- downloadHandler(
    filename = function() {
      paste0("TopScore_Network_", Sys.Date(), ".html")
    },
    content = function(file) {
      temp_file <- tempfile(fileext = ".html")
      visSave(vis_network_top_score(), file = temp_file, selfcontained = TRUE)
      
      html <- readLines(temp_file)
      
      style_block <- '<style>
      html, body {
        margin: 0;
        padding: 0;
        width: 100%;
        height: 100%;
        overflow: hidden;
      }
      #htmlwidget_container {
        width: 100vw !important;
        height: 100vh !important;
        position: relative;
      }
      .vis-network {
        width: 100% !important;
        height: 100% !important;
      }
      .vis-manipulation,
      .vis-navigation {
        z-index: 9999 !important;
        position: absolute !important;
        top: 10px;
        right: 10px;
      }
    </style>'
      
      head_index <- grep("<head>", html, fixed = TRUE)
      if (length(head_index) > 0) {
        html <- append(html, style_block, after = head_index)
      }
      
      writeLines(html, file)
    }
  ) #download network as html
  
  observeEvent(input$open_cluster_rename_modal, {
    df <- summary_table_top_score()
    req(df)
    
    cluster_names <- unique(df$`Original Cluster Group`)  # Use stable original names
    
    showModal(modalDialog(
      title = "Rename Clusters",
      size = "l",
      easyClose = FALSE,
      footer = tagList(
        modalButton("Cancel"),
        actionButton("apply_cluster_renames_top_score", "Apply Renaming")
      ),
      fluidPage(
        lapply(cluster_names, function(cl) {
          inputId <- paste0("rename_", gsub(" ", "_", cl))
          current_input <- input[[inputId]]
          default_val <- cluster_rename_map_top_score()[[cl]] %||% cl
          textInput(
            inputId = inputId,
            label = paste("Rename", cl),
            value = if (!is.null(current_input) && nzchar(current_input)) current_input else default_val
          )
        })
      )
    ))
  }) #renaming option
  
  
  observeEvent(input$apply_cluster_renames_top_score, {
    df <- summary_table_top_score()
    req(df)
    cluster_names <- unique(df$`Original Cluster Group`)
    
    new_map <- setNames(
      lapply(cluster_names, function(cl) {
        input[[paste0("rename_", gsub(" ", "_", cl))]] %||% cl
      }),
      cluster_names
    )
    
    cluster_rename_map_top_score(new_map)
    removeModal()
  }) #apply renaming
  
  
  
  
  ################################################ Process-Gene #########################
  vis_network_process_gene <- reactiveVal(NULL) #reactive for downloading network as HTML
  
  
  observe({
    req(input$qvalue_threshold_top_score)
    
    updateSliderInput(
      session,
      inputId = "qvalue_threshold_process_gene",
      min = input$qvalue_threshold_top_score,
      max = 10,
      value = input$qvalue_threshold_top_score
    )
  }) #get the qvalue slider values from the previous tab
  
  user_data_process_gene <- reactive({
    req(input$gene_logfc_input_process_gene)  # wait until user uploads
    file <- input$gene_logfc_input_process_gene$datapath
    name <- input$gene_logfc_input_process_gene$name

    read_user_gene_file(file, name)
  }) #user gene list uploader
  
  output$concordance_ui <- renderUI({
    req(input$gene_logfc_input_process_gene)
    if(ncol(user_data_process_gene()) == 2){
      div(
        style = "text-align: center;",
        checkboxInput("concordance_process_gene", "Concordance filtering using logFC values", value = FALSE)
      )
    }
  }) # show concordance UI if it is not only a gene list
  
  df_process_gene <- reactive({
    if (#input$direction_dropdown_process_gene == "Use Clustering from Top-Score Tab" &&
        !is.null(summary_table_top_score())) {
      
      df_initial_gene <- summary_table_top_score()
      df_initial_gene <- df_initial_gene[,-7]
      colnames(df_initial_gene) <- c("Pathway","Gene","Cluster","q-value","-log10(q)","Gene Count")
      df_initial_gene <- df_initial_gene[df_initial_gene$'-log10(q)' >= input$qvalue_threshold_process_gene,]
      
      df_long_gene <- df_initial_gene %>%
        mutate(Gene = str_replace_all(Gene, "\\s*/\\s*", "/")) %>% # replace gaps in front or before spaces to "/"
        separate_rows(Gene, sep = "/") #convert to long format by seperating using /
      
      full_edges_gene <- df_long_gene %>%
        count(Gene, Cluster, name = "weight") # count the gene-category matches and call it as weight
      
      return(full_edges_gene)
      
    } else {
      NULL  # Or placeholder/alternative if desired
    }
  })
  
  all_processes_process_gene <- reactive({ 
    req(df_process_gene())
    
    df <- df_process_gene()
    df <- sort(unique(df$Cluster)) #get the category name list
    return(df)
    
  })
  
  output$cluster_selector_process_gene <- renderUI({
    req(all_processes_process_gene())
    selectInput("selected_cluster_process_gene",
                "Select Cluster(s):",
                choices = all_processes_process_gene(),
                selected = all_processes_process_gene(),
                multiple = TRUE)
  }) #This dynamically adjusts process_selector UI with 5 pre-selected categories and allowing multiple choicess
    
  show_full_network_process_gene <- reactiveVal(TRUE) #this creates a reactive value where the default is true
  
  observeEvent(input$network_process_gene_selected, {
    req(df_process_gene())
    df <- df_process_gene()
    selected <- input$network_process_gene_selected
    if (!is.null(selected) && selected %in% df$Cluster) {
      show_full_network_process_gene(FALSE)
      updateSelectInput(session, "selected_cluster_process_gene", selected = selected)
    }
  })
  
  observeEvent(input$reset_network_process_gene, {
    show_full_network_process_gene(TRUE)
  }) #when clicked on this button show_full_network_process_gene becomes TRUE
  
  output$network_process_gene <- renderVisNetwork({
    req(df_process_gene())
    req(input$selected_cluster_process_gene)
    
    if (show_full_network_process_gene()) {
      
      edges <- df_process_gene() 
      
      edges <- edges %>%
        filter(Cluster %in% input$selected_cluster_process_gene) #filters gene-category matches according to selection
      
      
      g <- graph_from_data_frame(edges, directed = FALSE) #creates an igraph object
      
      centrality_scores <- switch(
        input$centrality_method_process_gene,
        "Degree" = degree(g),
        "Closeness" = closeness(g, normalized = TRUE),
        "Betweenness" = betweenness(g, normalized = TRUE),
        "PageRank" = page_rank(g)$vector
      ) #switch-case for centrality_method calculation 
      #Degree measures how many edges each node has
      #closeness centrality is the inverse of the average shortest path length to all other nodes.
      #Betweenness is how often a node lies on the shortest path between other nodes.
      #PageRank Computes a version of importance based on links from important nodes
      
      centrality_df <- data.frame(
        id = names(centrality_scores),
        centrality = as.numeric(centrality_scores)
      ) %>%
        filter(id %in% edges$Gene) #creates a DF where id equals to node names and centrality equals to scores
      #then filters for only nodes in geneID (so only gene scores not process)
      
      current_gene_ids <- unique(edges$Gene) #list of genes on the network
      if(input$metric_top_score == "NAFLD Activity Score"){
        scores_filtered <- scores %>% filter(Gene %in% current_gene_ids) #score df for current genes
      }else{
        scores_filtered <- scores_fibrosis %>% filter(Gene %in% current_gene_ids) #score df for current genes
        
      }
      centrality_df_filtered <- centrality_df %>% filter(id %in% current_gene_ids) #centrality_df for current genes
      
      total_min <- min(scores_filtered$total, na.rm = TRUE) #minimum total score in existing genes
      total_max <- max(scores_filtered$total, na.rm = TRUE) #maximum total score in existing genes
      centrality_min <- min(centrality_df_filtered$centrality, na.rm = TRUE) #minimum centrality score in existing genes
      centrality_max <- max(centrality_df_filtered$centrality, na.rm = TRUE) #maximum centrality score in existing genes
      total_range <- total_max - total_min # range of total score
      centrality_range <- centrality_max - centrality_min # range of centrality score
      
      user_genes <- NULL
      if (!is.null(input$gene_logfc_input_process_gene)) {
        user_genes <- tryCatch({
          x <- user_data_process_gene()[[1]]         # first column
          x <- as.character(x)
          x <- trimws(x)
          unique(x[!is.na(x) & nzchar(x)])
        }, error = function(e) NULL)
      } ## get user genes
      n_uploaded <- length(user_genes)
      species_info <- resolve_user_species(user_genes, input$metric_top_score,
                                           sum(user_genes %in% current_gene_ids))
      user_genes_original_filtered <- species_info$original
      user_genes <- species_info$genes

      gene_nodes <- data.frame(
        id = current_gene_ids,
        label = current_gene_ids,
        group = "gene",
        stringsAsFactors = FALSE
      ) %>%
        { if (!is.null(user_genes)) dplyr::filter(., id %in% user_genes) else . } %>% ### user_gene filter
        left_join(scores_filtered, by = c("id" = "Gene")) %>%
        left_join(centrality_df_filtered, by = "id") %>%
        mutate(
          total = replace_na(total, 0),
          centrality = replace_na(centrality, 0),
          total_norm = if (total_range == 0) 0 else (total - total_min) / total_range,
          centrality_norm = if (centrality_range == 0) 0 else (centrality - centrality_min) / centrality_range,
          combined_score = input$weight_total_process_gene * total_norm + (1 - input$weight_total_process_gene) * centrality_norm
        ) %>%
        filter(combined_score >= input$threshold_process_gene) %>%
        mutate(
          color = case_when(
            total_upregulated > total_downregulated ~ "#D55E00",
            total_upregulated < total_downregulated ~ "#0072B2",
            TRUE ~ "blue"
          ),
          
          size = 10 + 20 * combined_score,
          title = paste0(
            "Gene: ", id, "<br>",
            "Total score: ", round(total, 2), "<br>",
            "Upregulated: ", total_upregulated, "<br>",
            "Downregulated: ", total_downregulated, "<br>",
            "Centrality: ", round(centrality, 3), "<br>",
            "Combined score: ", round(combined_score, 2)
          )
        )
      
      if(!is.null(input$concordance_process_gene) && isTRUE(input$concordance_process_gene)){
        user_data <- user_data_process_gene()
        user_data <- user_data[user_data[[1]] %in% user_genes_original_filtered, ]
        user_genes_up <- translate_user_genes(user_data[user_data[[2]] > 0,1], species_info)
        user_genes_down <- translate_user_genes(user_data[user_data[[2]] < 0,1], species_info)
        gene_nodes <- gene_nodes[(gene_nodes$color == "#D55E00" & gene_nodes$id %in% user_genes_up) | (gene_nodes$color == "#0072B2" & gene_nodes$id %in% user_genes_down),]
      }
      
      upload_msg_process_gene(upload_summary_msg(n_uploaded, species_info$species, nrow(gene_nodes)))
      
      output$gene_table_process_gene <- renderDataTable({
        gene_nodes %>%
          select(
            Gene = id,
            total,
            total_upregulated,
            total_downregulated,
            centrality,
            combined_score
          ) %>%
          arrange(desc(combined_score))
      })
      
      process_nodes <- data.frame(
        id = unique(edges$Cluster),
        label = unique(edges$Cluster),
        group = "process",
        color = "orange",
        size = 30,
        title = unique(edges$Cluster)
      )
      
      weight_min <- min(edges$weight, na.rm = TRUE)
      weight_max <- max(edges$weight, na.rm = TRUE)
      
      edges_formatted <- edges %>%
        filter(Gene %in% gene_nodes$id) %>%
        rename(from = Gene, to = Cluster) %>%
        mutate(
          width = if (weight_max == weight_min) 1 else 1 + 4 * (weight - weight_min) / (weight_max - weight_min),
          title = paste0("Gene and Process<br>Connections: ", weight)
        )
      
      process_nodes <- process_nodes %>%
        filter(id %in% edges_formatted$to)
      
      nodes <- bind_rows(gene_nodes, process_nodes)
    }else {
      selected_category <- input$selected_cluster_process_gene[1]
      
      df_initial_gene <- summary_table_top_score()
      df_initial_gene <- df_initial_gene[,-7]
      colnames(df_initial_gene) <- c("Pathway","Gene","Cluster","q-value","-log10(q)","Gene Count")
      df_initial_gene <- df_initial_gene[df_initial_gene$'-log10(q)' >= input$qvalue_threshold_process_gene,]
      df_initial_gene <- df_initial_gene[df_initial_gene$Cluster == selected_category, ]
      
      
      df_long_gene <- df_initial_gene %>%
        mutate(Gene = str_replace_all(Gene, "\\s*/\\s*", "/")) %>% # replace gaps in front or before spaces to "/"
        separate_rows(Gene, sep = "/") #convert to long format by seperating using /
      
      edges_sub <- df_long_gene %>%
        count(Gene, Pathway, name = "weight") # count the gene-category matches and call it as weight
      
      g <- graph_from_data_frame(edges_sub, directed = FALSE)
      centrality_scores <- switch(
        input$centrality_method_process_gene,
        "Degree" = degree(g),
        "Closeness" = closeness(g, normalized = TRUE),
        "Betweenness" = betweenness(g, normalized = TRUE),
        "PageRank" = page_rank(g)$vector
      )
      centrality_df <- data.frame(id = names(centrality_scores), centrality = as.numeric(centrality_scores))
      
      current_gene_ids <- unique(edges_sub$Gene)
      if(input$metric_top_score == "NAFLD Activity Score"){
        scores_filtered <- scores %>% filter(Gene %in% current_gene_ids) #score df for current genes
      }else{
        scores_filtered <- scores_fibrosis %>% filter(Gene %in% current_gene_ids) #score df for current genes
        
      }
      centrality_df_filtered <- centrality_df %>% filter(id %in% current_gene_ids)
      
      total_min <- min(scores_filtered$total, na.rm = TRUE)
      total_max <- max(scores_filtered$total, na.rm = TRUE)
      centrality_min <- min(centrality_df_filtered$centrality, na.rm = TRUE)
      centrality_max <- max(centrality_df_filtered$centrality, na.rm = TRUE)
      total_range <- total_max - total_min
      centrality_range <- centrality_max - centrality_min
      
      user_genes <- NULL
      if (!is.null(input$gene_logfc_input_process_gene)) {
        user_genes <- tryCatch({
          x <- user_data_process_gene()[[1]]         # first column
          x <- as.character(x)
          x <- trimws(x)
          unique(x[!is.na(x) & nzchar(x)])
        }, error = function(e) NULL)
      } ## get user genes
      n_uploaded <- length(user_genes)
      species_info <- resolve_user_species(user_genes, input$metric_top_score,
                                           sum(user_genes %in% current_gene_ids))
      user_genes_original_filtered <- species_info$original
      user_genes <- species_info$genes

      gene_nodes <- data.frame(
        id = current_gene_ids,
        label = current_gene_ids,
        group = "gene",
        stringsAsFactors = FALSE
      ) %>%
        { if (!is.null(user_genes)) dplyr::filter(., id %in% user_genes) else . } %>% ### user_gene filter
        left_join(scores_filtered, by = c("id" = "Gene")) %>%
        left_join(centrality_df_filtered, by = "id") %>%
        mutate(
          total = replace_na(total, 0),
          centrality = replace_na(centrality, 0),
          total_norm = if (total_range == 0) 0 else (total - total_min) / total_range,
          centrality_norm = if (centrality_range == 0) 0 else (centrality - centrality_min) / centrality_range,
          combined_score = input$weight_total_process_gene * total_norm + (1 - input$weight_total_process_gene) * centrality_norm
        ) %>%
        filter(combined_score >= input$threshold_process_gene) %>%
        mutate(
          color = case_when(
            total_upregulated > total_downregulated ~ "#D55E00",
            total_upregulated < total_downregulated ~ "#0072B2",
            TRUE ~ "blue"
          ),
          
          size = 10 + 20 * combined_score,
          title = paste0(
            "Gene: ", id, "<br>",
            "Total score: ", round(total, 2), "<br>",
            "Upregulated: ", total_upregulated, "<br>",
            "Downregulated: ", total_downregulated, "<br>",
            "Centrality: ", round(centrality, 3), "<br>",
            "Combined score: ", round(combined_score, 2)
          )
        )
      
      if(!is.null(input$concordance_process_gene) && isTRUE(input$concordance_process_gene)){
        user_data <- user_data_process_gene()
        user_data <- user_data[user_data[[1]] %in% user_genes_original_filtered, ]
        user_genes_up <- translate_user_genes(user_data[user_data[[2]] > 0,1], species_info)
        user_genes_down <- translate_user_genes(user_data[user_data[[2]] < 0,1], species_info)
        gene_nodes <- gene_nodes[(gene_nodes$color == "#D55E00" & gene_nodes$id %in% user_genes_up) | (gene_nodes$color == "#0072B2" & gene_nodes$id %in% user_genes_down),]
      }
      upload_msg_process_gene(upload_summary_msg(n_uploaded, species_info$species, nrow(gene_nodes)))
      
      output$gene_table <- renderDataTable({
        gene_nodes %>%
          select(
            Gene = id,
            total,
            total_upregulated,
            total_downregulated,
            centrality,
            combined_score
          ) %>%
          arrange(desc(combined_score))
      })
      
      process_nodes <- data.frame(
        id = unique(edges_sub$Pathway),
        label = unique(edges_sub$Pathway),
        group = "process",
        color = "orange",
        size = 30,
        title = unique(edges_sub$Pathway)
      )
      
      edges_formatted <- edges_sub %>%
        filter(Gene %in% gene_nodes$id) %>%
        rename(from = Gene, to = Pathway) %>%
        mutate(
          width = 1 + log1p(weight),
          title = paste0("Connections: ", weight)
        )
      
      process_nodes <- process_nodes %>%
        filter(id %in% edges_formatted$to)
      
      nodes <- bind_rows(gene_nodes, process_nodes)
    }
    
    gene_summary_table_process_gene <- reactive({
      gene_nodes %>%
        select(
          Gene = id,
          total,
          total_upregulated,
          total_downregulated,
          centrality,
          combined_score
        ) %>%
        arrange(desc(combined_score))
    })
    
    output$gene_table_process_gene <- renderDataTable({
      gene_summary_table_process_gene()
    })
    
    output$download_summary_process_gene <- downloadHandler(
      filename = function() {
        paste0("Process_Gene_Summary_", Sys.Date(), ".xlsx")
      },
      content = function(file) {
        df <- isolate(gene_summary_table_process_gene())
        writexl::write_xlsx(df, file)
      }
    )
    
    updateCheckboxInput(session, "enable_physics_process_gene", value = TRUE) #checkbox updater when there is a new network
    
    network_object <- visNetwork(nodes, edges_formatted) %>%
      visNodes(shape = "dot", size = "size") %>%
      visEdges(smooth = FALSE) %>%
      visOptions(highlightNearest = TRUE, nodesIdSelection = TRUE) %>%
      visInteraction(dragView = TRUE, zoomView = TRUE, navigationButtons = TRUE) %>%
      visLayout(randomSeed = 123) %>%
      visPhysics(
        solver = "forceAtlas2Based",
        forceAtlas2Based = list(
          gravitationalConstant = -50,
          centralGravity = 0.001,
          springLength = 200,
          springConstant = 0.02,
          damping = 0.4,
          avoidOverlap = 1
        ),
        stabilization = list(
          enabled = TRUE,
          iterations = 200
        )
      )
    
    vis_network_process_gene(network_object) #make network object reactive
    
    network_object
  })
  
  observeEvent(input$enable_physics_process_gene, {
    visNetworkProxy("network_process_gene") %>%
      visPhysics(enabled = input$enable_physics_process_gene)
  }) #if checked stop motion
  
  output$download_network_process_gene <- downloadHandler(
    filename = function() {
      paste0("Process_Gene_Network_", Sys.Date(), ".html")
    },
    content = function(file) {
      temp_file <- tempfile(fileext = ".html")
      visSave(vis_network_process_gene(), file = temp_file, selfcontained = TRUE)
      
      html <- readLines(temp_file)
      
      # Updated CSS block targeting standard htmlwidget classes
      style_block <- '<style>
      html, body {
        margin: 0;
        padding: 0;
        width: 100vw;
        height: 100vh;
        overflow: hidden;
      }
      /* Override the hardcoded inline styles from htmlwidgets */
      .html-widget, .visNetwork {
        width: 100vw !important;
        height: 100vh !important;
        position: absolute !important;
        top: 0;
        left: 0;
      }
      /* Ensure the inner vis.js canvas fills the new 100vw/vh container */
      .vis-network {
        width: 100% !important;
        height: 100% !important;
      }
      .vis-manipulation,
      .vis-navigation {
        z-index: 9999 !important;
        position: absolute !important;
        top: 10px;
        right: 10px;
      }
    </style>'
      
      head_index <- grep("<head>", html, fixed = TRUE)
      if (length(head_index) > 0) {
        html <- append(html, style_block, after = head_index)
      }
      
      writeLines(html, file)
    }
  )
  
  ################################################ STRING ###########################
  
  summary_table_string <- reactiveVal(NULL) #reactive summary_table for the first tab
  vis_network_string <- reactiveVal(NULL) #reactive for downloading network as HTML
  
  user_data_string<- reactive({
    req(input$gene_logfc_input_string)  # wait until user uploads
    file <- input$gene_logfc_input_string$datapath
    name <- input$gene_logfc_input_string$name

    read_user_gene_file(file, name)
  })
  
  output$concordance_ui_string <- renderUI({
    req(input$gene_logfc_input_string)
    if(ncol(user_data_string()) == 2){
      div(
        style = "text-align: center;",
        checkboxInput("concordance_string", "Concordance filtering using logFC values", value = FALSE)
      )
    }
  })
  
  
  output$string_network <- renderVisNetwork({
    
    if(input$string_top_score == "NAFLD Activity Score"){
      edges_string <- string_nas
    }else if(input$string_top_score == "Fibrosis Stage"){
      edges_string <- string_fibrosis
    }
    
    edges_string <- edges_string[
      edges_string$combined_score >= input$string_combined_score,
    ]
    
    if(input$toggle_advanced %% 2 == 1){
      edges_string <- edges_string[
        (edges_string$neighborhood >= input$neighborhood &
        edges_string$fusion >= input$fusion &
        edges_string$cooccurence >= input$cooccurence &
        edges_string$coexpression >= input$coexpression &
        edges_string$experimental >= input$experimental &
        edges_string$database >= input$database &
        edges_string$textmining >= input$textmining),
      ]
    }
    
    # if(input$pubmed_string == TRUE){
    #   if(input$string_top_score == "NAFLD Activity Score"){
    #     genes_remove <- tcga_nas[tcga_nas$pubmed_mentioned == TRUE,]$Gene
    #   }else if(input$string_top_score == "Fibrosis Stage"){
    #     genes_remove <- tcga_fibrosis[tcga_fibrosis$pubmed_mentioned == TRUE,]$Gene
    #   }
    #   edges_string <- edges_string[
    #     !(edges_string$from_gene %in% genes_remove |
    #         edges_string$to_gene %in% genes_remove),
    #   ]
    # }
    
    if(input$string_direction == "Only Upregulated"){
      if(input$string_top_score == "NAFLD Activity Score"){
        genes_keep <- tcga_nas[as.numeric(tcga_nas$total_upregulated) > as.numeric(tcga_nas$total_downregulated),]$Gene
        edges_string <- edges_string[
          (edges_string$from_gene %in% genes_keep &
              edges_string$to_gene %in% genes_keep),
        ]
      }else if(input$string_top_score == "Fibrosis Stage"){
        genes_keep <- tcga_fibrosis[as.numeric(tcga_fibrosis$total_upregulated) > as.numeric(tcga_fibrosis$total_downregulated),]$Gene
        edges_string <- edges_string[
          (edges_string$from_gene %in% genes_keep &
             edges_string$to_gene %in% genes_keep),
        ]
      }
    }else if(input$string_direction == "Only Downregulated"){
      if(input$string_top_score == "NAFLD Activity Score"){
        genes_keep <- tcga_nas[as.numeric(tcga_nas$total_downregulated) > as.numeric(tcga_nas$total_upregulated),]$Gene
        edges_string <- edges_string[
          (edges_string$from_gene %in% genes_keep &
             edges_string$to_gene %in% genes_keep),
        ]
      }else if(input$string_top_score == "Fibrosis Stage"){
        genes_keep <- tcga_fibrosis[as.numeric(tcga_fibrosis$total_downregulated) > as.numeric(tcga_fibrosis$total_upregulated),]$Gene
        edges_string <- edges_string[
          (edges_string$from_gene %in% genes_keep &
             edges_string$to_gene %in% genes_keep),
        ]
      }
    }
    
    if(input$string_top_score == "NAFLD Activity Score"){
      genes_remove <- tcga_nas[as.numeric(tcga_nas$total) < input$string_total_score,]$Gene
      edges_string <- edges_string[
        !(edges_string$from_gene %in% genes_remove |
            edges_string$to_gene %in% genes_remove),
      ]
    }else if(input$string_top_score == "Fibrosis Stage"){
      genes_remove <- tcga_fibrosis[as.numeric(tcga_fibrosis$total) < input$string_total_score,]$Gene
      edges_string <- edges_string[
        !(edges_string$from_gene %in% genes_remove |
            edges_string$to_gene %in% genes_remove),
      ]
    }
    
    
    if (nrow(edges_string) < 2) {
      showNotification("Not enough genes to compute edges.", type = "warning")
      return(NULL)
    } #make sure there is enough data after filtering
    
    if(input$string_top_score == "NAFLD Activity Score"){
      edges_meta <- scores[scores$Gene %in% c(edges_string$from_gene,edges_string$to_gene),]
    }else if(input$string_top_score == "Fibrosis Stage"){
      edges_meta <- scores_fibrosis[scores_fibrosis$Gene %in% c(edges_string$from_gene,edges_string$to_gene),]
    }
    
    edges <- edges_string[, c("from_gene", "to_gene","combined_score")]
    colnames(edges) <- c("from", "to","weight")
    
    edges$width <- scales::rescale(
      log(edges$weight),
      to = c(1,8)
    )
    
    edges$title <- paste0(
    "Combined Score: ", edges$weight
    ) 
    
    g_string <- igraph::graph_from_data_frame(as.matrix(edges), directed = FALSE)
    
    centrality_scores <- switch(
      input$centrality_method_string,
      "Degree" = degree(g_string),
      "Closeness" = closeness(g_string, normalized = TRUE),
      "Betweenness" = betweenness(g_string, normalized = TRUE),
      "PageRank" = page_rank(g_string)$vector
    )
    centrality_df_string <- data.frame(id = names(centrality_scores), centrality = as.numeric(centrality_scores))
    
    network_nodes <- data.frame(
      id = unique(c(edges$from, edges$to)),
      label = unique(c(edges$from, edges$to)),
      stringsAsFactors = FALSE
    )
    
    network_nodes <- merge(network_nodes, centrality_df_string, by = "id", all.x = TRUE)
    
    
    network_nodes <- merge(
      network_nodes,
      edges_meta[, c("Gene", "total", "total_upregulated", "total_downregulated")],
      by.x = "id",
      by.y = "Gene",
      all.x = TRUE
    )
    
    total_min <- min(network_nodes$total, na.rm = TRUE) #minimum total score in existing genes
    total_max <- max(network_nodes$total, na.rm = TRUE) #maximum total score in existing genes
    centrality_min <- min(network_nodes$centrality, na.rm = TRUE) #minimum centrality score in existing genes
    centrality_max <- max(network_nodes$centrality, na.rm = TRUE) #maximum centrality score in existing genes
    total_range <- total_max - total_min # range of total score
    centrality_range <- centrality_max - centrality_min # range of centrality score
    
    network_nodes %>% mutate(
      total = replace_na(total, 0),
      centrality = replace_na(centrality, 0))

    network_nodes$total_norm <- if (total_range == 0) 0 else (network_nodes$total - total_min) / total_range
    network_nodes$centrality_norm <- if (centrality_range == 0) 0 else (network_nodes$centrality - centrality_min) / centrality_range
    network_nodes$combined_score <- input$weight_total_string * network_nodes$total_norm + (1 - input$weight_total_string) * network_nodes$centrality_norm
      
    network_nodes <- network_nodes[network_nodes$combined_score >= input$threshold_string,]

    network_nodes$size <- scales::rescale(
      network_nodes$total,
      to = c(10, 40)
    )
    
    network_nodes$color <- ifelse(
      network_nodes$total_upregulated > network_nodes$total_downregulated,
      "#D55E00",   # #D55E00
      ifelse(
        network_nodes$total_downregulated > network_nodes$total_upregulated,
        "#0072B2", 
        "gray"  # gray
      )
    )
    
  
    
    # Convert edge list to igraph for clustering
    if (input$clustering_method_string != "No Clustering") {
      
      cluster_result <- switch(
        input$clustering_method_string,
        "Louvain" = igraph::cluster_louvain(g_string),
        "Edge Betweenness" = igraph::cluster_edge_betweenness(g_string),
        "Label Propagation" = igraph::cluster_label_prop(g_string)
      )
      
      membership <- igraph::membership(cluster_result)
      
      cluster_df <- data.frame(
        id = names(membership),
        group = paste0("Cluster ", membership),
        stringsAsFactors = FALSE
      )
      
      network_nodes <- merge(network_nodes, cluster_df, by = "id", all.x = TRUE)
      
    } else {
      network_nodes$group <- "No cluster"
    }
    
    network_nodes <- na.omit(network_nodes)
    
    user_genes <- NULL
    if (!is.null(input$gene_logfc_input_string)) {
      user_genes <- tryCatch({
        x <- user_data_string()[[1]]         # first column
        x <- as.character(x)
        x <- trimws(x)
        unique(x[!is.na(x) & nzchar(x)])
      }, error = function(e) NULL)
    } ## get user genes
    n_uploaded <- length(user_genes)
    species_info <- resolve_user_species(user_genes, input$string_top_score,
                                         sum(user_genes %in% network_nodes$id))
    user_genes_original_filtered <- species_info$original
    user_genes <- species_info$genes

    if (!is.null(user_genes)){
      network_nodes <- network_nodes[network_nodes$id %in% user_genes,]
    }


    if(!is.null(input$concordance_string) && isTRUE(input$concordance_string)){
      user_data <- user_data_string()
      user_data <- user_data[user_data[[1]] %in% user_genes_original_filtered, ]
      user_genes_up <- translate_user_genes(user_data[user_data[[2]] > 0,1], species_info)
      user_genes_down <- translate_user_genes(user_data[user_data[[2]] < 0,1], species_info)
      network_nodes <- network_nodes[(network_nodes$color ==  "#D55E00" & network_nodes$id %in% user_genes_up) | (network_nodes$color == "#0072B2" & network_nodes$id %in% user_genes_down),]
    }
    upload_msg_string(upload_summary_msg(n_uploaded, species_info$species, nrow(network_nodes)))
    
    summary_df_string <- data.frame(
      `Gene Name` = network_nodes$id,
      `Total Score` = network_nodes$total,
      `Degree` = network_nodes$centrality,
      `Combined Score` = network_nodes$combined_score,
      `Cluster Group` = network_nodes$group,
      stringsAsFactors = FALSE
    ) #create the summary table
    
    # Save table for reuse in UI/download
    summary_table_string(summary_df_string) #reactive
    
    
    updateCheckboxInput(session, "enable_physics_string", value = TRUE) #checkbox updater when there is a new network
    
    
    network_object <- visNetwork(network_nodes, edges, height = "100%", width = "100%") %>%
      visNodes(shape = "dot", size = "size") %>%
      visEdges(smooth = FALSE) %>%
      visOptions(highlightNearest = TRUE, nodesIdSelection = TRUE, selectedBy = "group") %>%
      visOptions(highlightNearest = TRUE, nodesIdSelection = TRUE, manipulation = FALSE,selectedBy = "group") %>%
      visInteraction(dragView = TRUE, zoomView = TRUE, navigationButtons = TRUE) %>%
      visLayout(randomSeed = 123) %>%
      visPhysics(
        solver = "forceAtlas2Based",
        forceAtlas2Based = list(
          gravitationalConstant = -50,
          centralGravity = 0.001,
          springLength = 200,
          springConstant = 0.02,
          damping = 0.4,
          avoidOverlap = 1
        ),
        stabilization = list(
          enabled = TRUE,
          iterations = 200
        )
      ) #visualizes the network including its properties (physics etc.)
    
    # Save for download
    vis_network_string(network_object) #make network object reactive
    
    # Render to UI
    network_object #visualizes the network
  })  
  
  observeEvent(input$enable_physics_string, {
    visNetworkProxy("string_network") %>%
      visPhysics(enabled = input$enable_physics_string)
  }) #if checked stop motion
  
  
  output$network_summary_string <- DT::renderDataTable({
    req(summary_table_string())
    
    df <- summary_table_string()
    colnames(df) <- c("Gene Name", "Total Score", "Degree", "Combined Score","Cluster Group")
    # colnames(df) <- gsub("\\.", " ", colnames(df))  # Fix column names for display
    
    DT::datatable(
      df,
      options = list(
        pageLength = 10,
        autoWidth = TRUE,
        columnDefs = list(
          list(targets = 1, width = '200px', className = 'dt-wrap')
        )
      ),
      rownames = FALSE,
      escape = FALSE
    ) %>% DT::formatStyle(
      columns = 1,
      `white-space` = "normal",
      `word-wrap` = "break-word"
    )
  }) #visualize table output below network
  
  output$download_summary_string <- downloadHandler(
    filename = function() {
      paste0("Pathway_Summary_STRING", Sys.Date(), ".xlsx")
    },
    content = function(file) {
      df <- summary_table_string()
      colnames(df) <- c("Gene Name", "Total Score", "Degree", "Combined Score","Cluster Group")
      colnames(df) <- gsub("\\.", " ", colnames(df))  # Fix column names for display
      writexl::write_xlsx(df, file)
    }
  ) #download the network clusterings etc. as a table.
  
  output$download_network_string <- downloadHandler(
    filename = function() {
      paste0("STRING_Network_", Sys.Date(), ".html")
    },
    content = function(file) {
      temp_file <- tempfile(fileext = ".html")
      visSave(vis_network_string(), file = temp_file, selfcontained = TRUE)
      
      html <- readLines(temp_file)
      
      style_block <- '<style>
      html, body {
        margin: 0;
        padding: 0;
        width: 100%;
        height: 100%;
        overflow: hidden;
      }
      #htmlwidget_container {
        width: 100vw !important;
        height: 100vh !important;
        position: relative;
      }
      .vis-network {
        width: 100% !important;
        height: 100% !important;
      }
      .vis-manipulation,
      .vis-navigation {
        z-index: 9999 !important;
        position: absolute !important;
        top: 10px;
        right: 10px;
      }
    </style>'
      
      head_index <- grep("<head>", html, fixed = TRUE)
      if (length(head_index) > 0) {
        html <- append(html, style_block, after = head_index)
      }
      
      writeLines(html, file)
    }
  ) #download network as html
  
  ############################### Stage-Dependent Networks
  
  
  summary_table_temporal <- reactiveVal(NULL) #reactive summary_table for the first tab
  vis_network_temporal <- reactiveVal(NULL) #reactive for downloading network as HTML

  selected_data_temporal <- reactive({
    req(input$category_dropdown_temporal, selected_data_temporal_xlsx()) ## require category_dropdown_top_score from the dropdown menu 
    if(input$metric_temporal == "NAFLD Activity Score"){
      if(input$direction_dropdown_temporal == "All Genes"){
        filename_rds <- paste0("NAS_", input$category_dropdown_temporal, "_total.rds") ## load the selected MSigDB result
        filepath_rds <- file.path("temporal/NAS/merged network", filename_rds) ## get the file path for the file
      }else if(input$direction_dropdown_temporal == "Only Upregulated"){
        filename_rds <- paste0("NAS_", input$category_dropdown_temporal, "_up.rds") ## load the selected MSigDB result
        filepath_rds <- file.path("temporal/NAS/merged network", filename_rds) ## get the file path for the file
      }else{
        filename_rds <- paste0("NAS_", input$category_dropdown_temporal, "_down.rds") ## load the selected MSigDB result
        filepath_rds <- file.path("temporal/NAS/merged network", filename_rds) ## get the file path for the file
      } 
    }else{
      if(input$direction_dropdown_temporal == "All Genes"){
        filename_rds <- paste0("fibrosis_", input$category_dropdown_temporal, "_total.rds") ## load the selected MSigDB result
        filepath_rds <- file.path("temporal/fibrosis/merged network", filename_rds) ## get the file path for the file
      }else if(input$direction_dropdown_temporal == "Only Upregulated"){
        filename_rds <- paste0("fibrosis_", input$category_dropdown_temporal, "_up.rds") ## load the selected MSigDB result
        filepath_rds <- file.path("temporal/fibrosis/merged network", filename_rds) ## get the file path for the file
      }else{
        filename_rds <- paste0("fibrosis_", input$category_dropdown_temporal, "_down.rds") ## load the selected MSigDB result
        filepath_rds <- file.path("temporal/fibrosis/merged network", filename_rds) ## get the file path for the file
      } 
    }
    
    if (!file.exists(filepath_rds)) {
      validate(need(FALSE, paste("No file found for", input$category_dropdown_temporal)))
      return(NULL)  # never actually reached because validate() stops
    }
    
    df_rds <- readRDS(filepath_rds) ## get the file
    
    df_xlsx <- selected_data_temporal_xlsx() #get msigdb output
    
    if (length(unique(df_xlsx$comparison)) > 1) {
      unique_comparisons <- unique(df_xlsx$comparison)
      df_rds_first <- df_rds[df_rds$comparison == unique_comparisons[1],]
      df_temporal_first <- df_xlsx[df_xlsx$comparison == unique_comparisons[1],]
      df_rds_first <- df_rds_first[(df_rds_first$from %in% df_temporal_first$ID), ] #filter for q value since some of msigdb will be eliminated
      df_rds_first <- df_rds_first[df_rds_first$to %in% df_temporal_first$ID,] #filter for q value since some of msigdb will be eliminated
      df_rds_second <- df_rds[df_rds$comparison == unique_comparisons[2],]
      df_temporal_second <- df_xlsx[df_xlsx$comparison == unique_comparisons[2],]
      df_rds_second <- df_rds_second[(df_rds_second$from %in% df_temporal_second$ID), ] #filter for q value since some of msigdb will be eliminated
      df_rds_second <- df_rds_second[df_rds_second$to %in% df_temporal_second$ID,] #filter for q value since some of msigdb will be eliminated
      df_rds <-  bind_rows(df_rds_first, df_rds_second)
    }else{
      df_rds <- df_rds[(df_rds$from %in% df_xlsx$ID), ] #filter for q value since some of msigdb will be eliminated
      df_rds <- df_rds[df_rds$to %in% df_xlsx$ID,] #filter for q value since some of msigdb will be eliminated
    }
    
    return(df_rds)
  })
  
  selected_data_temporal_xlsx <- reactive({
    req(input$category_dropdown_temporal) ## require category_dropdown_top_score from the dropdown menu 
    if(input$metric_temporal == "NAFLD Activity Score"){
      if(input$slider_temporal == "1"){
        current_time <- "1_VS_0"
        previous_time <- NULL
        current_folder <- "first"
        previous_folder <- NULL
      }else if(input$slider_temporal == "2"){
        current_time <- "2_VS_0"
        previous_time <- "1_VS_0"
        current_folder <- "second"
        previous_folder <- "first"
      }else if(input$slider_temporal == "3"){
        current_time <- "3_VS_0"
        previous_time <- "2_VS_0"
        current_folder <- "third"
        previous_folder <- "second"
      }else{
        current_time <- "4_VS_0"
        previous_time <- "3_VS_0"
        current_folder <- "fourth"
        previous_folder <- "third"
      }
      if(input$direction_dropdown_temporal == "All Genes"){
        if(is.null(previous_time)){
          current_filename_xlsx <- paste0(current_time,"_",input$category_dropdown_temporal,"_","total.xlsx")
          current_filepath_xlsx <- file.path("temporal/NAS",current_folder,current_filename_xlsx)
        }else{
          current_filename_xlsx <-  paste0(current_time,"_",input$category_dropdown_temporal,"_","total.xlsx")
          previous_filename_xlsx <-  paste0(previous_time,"_",input$category_dropdown_temporal,"_","total.xlsx")
          current_filepath_xlsx <- file.path("temporal/NAS",current_folder,current_filename_xlsx)
          previous_filepath_xlsx <- file.path("temporal/NAS",previous_folder,previous_filename_xlsx)
        }
      }else if(input$direction_dropdown_temporal == "Only Upregulated"){
        if(is.null(previous_time)){
          current_filename_xlsx <- paste0(current_time,"_",input$category_dropdown_temporal,"_","up.xlsx")
          current_filepath_xlsx <- file.path("temporal/NAS",current_folder,current_filename_xlsx)
        }else{
          current_filename_xlsx <-  paste0(current_time,"_",input$category_dropdown_temporal,"_","up.xlsx")
          previous_filename_xlsx <-  paste0(previous_time,"_",input$category_dropdown_temporal,"_","up.xlsx")
          current_filepath_xlsx <- file.path("temporal/NAS",current_folder,current_filename_xlsx)
          previous_filepath_xlsx <- file.path("temporal/NAS",previous_folder,previous_filename_xlsx)
        }
      }else{
        if(is.null(previous_time)){
          current_filename_xlsx <- paste0(current_time,"_",input$category_dropdown_temporal,"_","down.xlsx")
          current_filepath_xlsx <- file.path("temporal/NAS",current_folder,current_filename_xlsx)
        }else{
          current_filename_xlsx <-  paste0(current_time,"_",input$category_dropdown_temporal,"_","down.xlsx")
          previous_filename_xlsx <-  paste0(previous_time,"_",input$category_dropdown_temporal,"_","down.xlsx")
          current_filepath_xlsx <- file.path("temporal/NAS",current_folder,current_filename_xlsx)
          previous_filepath_xlsx <- file.path("temporal/NAS",previous_folder,previous_filename_xlsx)
        }
      } 
    }else{
      if(input$slider_temporal == "F4"){
        current_time <- "F4_VS_F0_F1"
        previous_time <- "F3_VS_F0_F1"
        current_folder <- "third"
        previous_folder <- "second"
      }else if(input$slider_temporal == "F3"){
        current_time <- "F3_VS_F0_F1"
        previous_time <- "F2_VS_F0_F1"
        current_folder <- "second"
        previous_folder <- "first"
      }else{
        current_time <- "F2_VS_F0_F1"
        previous_time <- NULL
        current_folder <- "first"
        previous_folder <- NULL
      }
      if(input$direction_dropdown_temporal == "All Genes"){
        if(is.null(previous_time)){
          current_filename_xlsx <- paste0(current_time,"_",input$category_dropdown_temporal,"_","total.xlsx")
          current_filepath_xlsx <- file.path("temporal/fibrosis",current_folder,current_filename_xlsx)
        }else{
          current_filename_xlsx <-  paste0(current_time,"_",input$category_dropdown_temporal,"_","total.xlsx")
          previous_filename_xlsx <-  paste0(previous_time,"_",input$category_dropdown_temporal,"_","total.xlsx")
          current_filepath_xlsx <- file.path("temporal/fibrosis",current_folder,current_filename_xlsx)
          previous_filepath_xlsx <- file.path("temporal/fibrosis",previous_folder,previous_filename_xlsx)
        }
      }else if(input$direction_dropdown_temporal == "Only Upregulated"){
        if(is.null(previous_time)){
          current_filename_xlsx <- paste0(current_time,"_",input$category_dropdown_temporal,"_","up.xlsx")
          current_filepath_xlsx <- file.path("temporal/fibrosis",current_folder,current_filename_xlsx)
        }else{
          current_filename_xlsx <-  paste0(current_time,"_",input$category_dropdown_temporal,"_","up.xlsx")
          previous_filename_xlsx <-  paste0(previous_time,"_",input$category_dropdown_temporal,"_","up.xlsx")
          current_filepath_xlsx <- file.path("temporal/fibrosis",current_folder,current_filename_xlsx)
          previous_filepath_xlsx <- file.path("temporal/fibrosis",previous_folder,previous_filename_xlsx)
        }
      }else{
        if(is.null(previous_time)){
          current_filename_xlsx <- paste0(current_time,"_",input$category_dropdown_temporal,"_","down.xlsx")
          current_filepath_xlsx <- file.path("temporal/fibrosis",current_folder,current_filename_xlsx)
        }else{
          current_filename_xlsx <-  paste0(current_time,"_",input$category_dropdown_temporal,"_","down.xlsx")
          previous_filename_xlsx <-  paste0(previous_time,"_",input$category_dropdown_temporal,"_","down.xlsx")
          current_filepath_xlsx <- file.path("temporal/fibrosis",current_folder,current_filename_xlsx)
          previous_filepath_xlsx <- file.path("temporal/fibrosis",previous_folder,previous_filename_xlsx)
        }
      } 
    }
    
    if (file.exists(current_filepath_xlsx)) {
      current_df_xlsx <- read_xlsx(current_filepath_xlsx)
      current_df_xlsx$comparison <- current_time
    }else{
      current_df_xlsx <- NULL
    }
    
    if(!is.null(previous_time)){
    if (file.exists(previous_filepath_xlsx)){
      previous_df_xlsx <- read_xlsx(previous_filepath_xlsx)
      previous_df_xlsx$comparison <- previous_time
    }else{
      previous_df_xlsx <- NULL
    }
    }else{
      previous_df_xlsx <- NULL
    }
    
    if(is.null(current_df_xlsx) & is.null(previous_df_xlsx)){
      validate(need(FALSE, paste("No file found for", input$category_dropdown_temporal)))
      return(NULL)  # never actually reached because validate() stops
    }
    
    df_list <- list(current_df_xlsx, previous_df_xlsx)
    df_list <- df_list[!sapply(df_list, is.null)]  # remove NULLs
    
    if (length(df_list) > 0) {
      df_xlsx <- bind_rows(df_list)
    }
    
    df_xlsx <- df_xlsx %>%
      filter(-log10(qvalue) >= input$qvalue_threshold_temporal) # Filter by threshold on -log10(qvalue)
    
    return(df_xlsx)
  })
  
  
  observeEvent(selected_data_temporal(), {
    df <- selected_data_temporal() #get the current rds
    if (is.null(df) || nrow(df) < 2) return() # return null if empty or 1 patheay
    max_weight <- max(df$weight, na.rm = TRUE) #get the maximum overlap
    updateSliderInput(
      session,
      inputId = "shared_gene_threshold_temporal",
      min = 1,
      max = max_weight, #update maximum of the shared gene slider to maximum current overlap.
      value = min(3, max_weight)
    )
  })
  
  
  output$temporal_net <- renderVisNetwork({
    req(selected_data_temporal(),selected_data_temporal_xlsx()) ## require selected_data_top_score reactive value
    edges <- selected_data_temporal() ## assign to a variable
    df_temporal <- selected_data_temporal_xlsx() ## assign to a variable
    df_temporal$gene_list <- strsplit(df_temporal$geneID, "/") ## seperate genes to get them as list
    
    if(input$metric_temporal == "NAFLD Activity Score"){
      if(input$slider_temporal == "1"){
        current_time <- "1_VS_0"
        previous_time <- NULL
        current_folder <- "first"
        previous_folder <- NULL
      }else if(input$slider_temporal == "2"){
        current_time <- "2_VS_0"
        previous_time <- "1_VS_0"
        current_folder <- "second"
        previous_folder <- "first"
      }else if(input$slider_temporal == "3"){
        current_time <- "3_VS_0"
        previous_time <- "2_VS_0"
        current_folder <- "third"
        previous_folder <- "second"
      }else{
        current_time <- "4_VS_0"
        previous_time <- "3_VS_0"
        current_folder <- "fourth"
        previous_folder <- "third"
      }
    }else{
      if(input$slider_temporal == "F4"){
        current_time <- "F4_VS_F0_F1"
        previous_time <- "F3_VS_F0_F1"
        current_folder <- "third"
        previous_folder <- "second"
      }else if(input$slider_temporal == "F3"){
        current_time <- "F3_VS_F0_F1"
        previous_time <- "F2_VS_F0_F1"
        current_folder <- "second"
        previous_folder <- "first"
      }else{
        current_time <- "F2_VS_F0_F1"
        previous_time <- NULL
        current_folder <- "first"
        previous_folder <- NULL
      }
    }
    
    if(is.null(previous_time)){
      edges <- edges[edges$comparison == current_time,]
    }else{
      edges <- edges[edges$comparison == current_time | edges$comparison == previous_time,]
    }
    
    if (nrow(edges) < 2) {
      showNotification("Not enough pathways to compute edges.", type = "warning")
      return(NULL)
    } #make sure there is enough data after filtering
    
    if (input$edge_filter_method_temporal == "Jaccard Index") {
      edges$title <- paste0(
        "Jaccard Index: ", round(edges$jaccard, 3),
        "\nUpregulated: ", edges$up,
        "\nDownregulated: ", edges$down
      ) #titles for jaccard option
    } else {
      edges$title <- paste0(
        edges$weight, " Shared genes",
        "\nUpregulated: ", edges$up,
        "\nDownregulated: ", edges$down
      )
    } #titles for shared gene option
    
    if (input$edge_filter_method_temporal == "Shared Genes") {
      edges <- edges[edges$weight >= input$shared_gene_threshold_temporal, ]
    } else if (input$edge_filter_method_temporal == "Jaccard Index") {
      edges <- edges[edges$jaccard >= input$jaccard_threshold_temporal, ]
    } #filterings forr edge cut methods
    
    df_temporal$qvalue[df_temporal$qvalue == 0] <- 1e-300 # Avoid log10(0)
    
    if (input$include_unconnected_nodes_temporal) {
      used_ids <- unique(c(edges$from, edges$to))
      df_temporal <- df_temporal[df_temporal$ID %in% used_ids, ]
    } #if checkbox is checked there is no unconnected nodes (I know the name indicates opposite I will fix the naming but it works)
    
    df_temporal <- df_temporal %>%
      arrange(desc(comparison == current_time)) %>%   # rows with "a" come first
      distinct(ID, .keep_all = TRUE)         # keep first occurrence of each ID
   
     gene_counts <- lengths(df_temporal$gene_list)
    node_sizes <- -log10(df_temporal$qvalue) #node size = log transformed q value
    node_sizes_scaled <- scales::rescale(node_sizes, to = c(10, 40))  # adjust range if needed
    
    # Convert edge list to igraph for clustering
    if (nrow(edges) > 0 && input$clustering_method_temporal != "No Clustering") {
      edge_mat <- as.matrix(edges[, c("from", "to")])
      g <- igraph::graph_from_edgelist(edge_mat, directed = FALSE)
      
      cluster_result <- switch(input$clustering_method_temporal,
                               "Louvain" = igraph::cluster_louvain(g),
                               "Edge Betweenness" = igraph::cluster_edge_betweenness(g),
                               "Label Propagation" = igraph::cluster_label_prop(g))
      
      membership <- igraph::membership(cluster_result)
      cluster_df <- data.frame(
        id = names(membership),
        group = paste0("Cluster ", membership),
        stringsAsFactors = FALSE
      )
    } else {
      if(nrow(df_temporal) > 0){
      cluster_df <- data.frame(
        id = df_temporal$ID,
        group = "Unclustered",
        stringsAsFactors = FALSE
      )
      }
    }
    
    if(nrow(df_temporal) == 0){
      return(NULL)
    }
    
    
    network_nodes <- unique(data.frame(
      id = df_temporal$ID,
      label = df_temporal$ID,
      size = node_sizes_scaled,
      title = paste(gene_counts, "genes in pathway"),
      color = NA,
      stringsAsFactors = FALSE
    )) #create network_nodes
    
    if(!is.null(previous_time)){
      edges$color <- ifelse(edges$comparison == previous_time, "gray", edges$color)
      edges <- edges %>%
        rowwise() %>%
        mutate(pair = paste(sort(c(from, to)), collapse = "_")) %>%
        ungroup() %>%
        group_by(pair) %>%
        filter(!(color == "gray" & any(color != "gray"))) %>%  # drop gray if a non-gray version exists
        ungroup() %>%
        select(-pair)
      ids_only_previous <- df_temporal %>%
        group_by(ID) %>%
        filter(all(comparison == previous_time)) %>%
        pull(ID) %>%
        unique()
       network_nodes$color <- ifelse(network_nodes$id %in% ids_only_previous, "gray", NA)
    }
    network_nodes <- merge(network_nodes, cluster_df, by = "id", all.x = TRUE) #merge with clustering info previously obtained
    network_nodes <- network_nodes %>%
      mutate(`Lost Node` = color == "gray")
    if(!is.null(previous_time)){
      network_nodes
      summary_table <- data.frame(
        `Pathway ID` = network_nodes$id,
        `Gene IDs` = df_temporal$geneID[match(network_nodes$id, df_temporal$ID)],
        `Cluster Group` = network_nodes$group,
        `q-value` = df_temporal$qvalue[match(network_nodes$id, df_temporal$ID)],
        `-log10(q)` = node_sizes[match(network_nodes$id, df_temporal$ID)],
        `Gene Count` = gene_counts[match(network_nodes$id, df_temporal$ID)],
        `Lost Node` = network_nodes$`Lost Node`,
        stringsAsFactors = FALSE
      ) #create the summary table
    }else{
      summary_table <- data.frame(
        `Pathway ID` = network_nodes$id,
        `Gene IDs` = df_temporal$geneID[match(network_nodes$id, df_temporal$ID)],
        `Cluster Group` = network_nodes$group,
        `q-value` = df_temporal$qvalue[match(network_nodes$id, df_temporal$ID)],
        `-log10(q)` = node_sizes[match(network_nodes$id, df_temporal$ID)],
        `Gene Count` = gene_counts[match(network_nodes$id, df_temporal$ID)],
        stringsAsFactors = FALSE
      ) #create the summary table
    }
    # Save table for reuse in UI/download
    summary_table_temporal(summary_table) #reactive
    updateCheckboxInput(session, "enable_physics_temporal", value = TRUE) #checkbox updater when there is a new network
    if(input$previous_nodes){
      edges <- edges[edges$color == "gray",]
      network_nodes <- network_nodes[network_nodes$color == "gray",]
      network_nodes <- na.omit(network_nodes)
      edges <- edges[edges$from %in% network_nodes$id & edges$to %in% network_nodes$id,]
      if(input$include_unconnected_nodes_temporal == TRUE){
        network_nodes <- network_nodes[network_nodes$id %in% c(edges$from, edges$to),]
      } #problem due to variable naming I know this seems counterintutitive will fix it!
    }
    edges$color <- ifelse(
      edges$color == "green",
      "#D55E00",
      ifelse(
        edges$color == "red",
        "#0072B2",
        edges$color
      )
    )    
    network_object <- visNetwork(network_nodes, edges, height = "100%", width = "100%") %>%
      visNodes(shape = "dot", size = "size") %>%
      visEdges(smooth = FALSE) %>%
      visOptions(highlightNearest = TRUE, nodesIdSelection = TRUE, selectedBy = "group") %>%
      visOptions(highlightNearest = TRUE, nodesIdSelection = TRUE, manipulation = FALSE,selectedBy = "group") %>%
      visInteraction(dragView = TRUE, zoomView = TRUE, navigationButtons = TRUE) %>%
      visLayout(randomSeed = 123) %>%
      visPhysics(
        solver = "forceAtlas2Based",
        forceAtlas2Based = list(
          gravitationalConstant = -50,
          centralGravity = 0.001,
          springLength = 200,
          springConstant = 0.02,
          damping = 0.4,
          avoidOverlap = 1
        ),
        stabilization = list(
          enabled = TRUE,
          iterations = 200
        )
      ) #visualizes the network including its properties (physics etc.)
    
    # Save for download
    vis_network_temporal(network_object) #make network object reactive
    
    # Render to UI
    network_object #visualizes the network
  })
  
  observeEvent(input$enable_physics_temporal, {
    visNetworkProxy("temporal_net") %>%
      visPhysics(enabled = input$enable_physics_temporal)
  }) #if checked stop motion
  
  observeEvent({
    input$direction_dropdown_temporal
    input$category_dropdown_temporal
  }, {
    # Reset q-value threshold to default (e.g., 1.3 means q < 0.05)
    updateSliderInput(session, "qvalue_threshold_temporal", value = 5)
    
    # Reset edge filtering thresholds
    updateSliderInput(session, "shared_gene_threshold_temporal", value = 3)
    updateSliderInput(session, "jaccard_threshold_temporal", value = 0.2)
    updateSliderInput(session, "score_threshold_temporal", value = 0)
    
    # Reset edge filter method
    updateSelectInput(session, "edge_filter_method_temporal", selected = "Jaccard Index")
    
    # Reset clustering dropdown
    updateSelectInput(session, "clustering_method_temporal", selected = "No Clustering")
    
    # Optionally reset include_unconnected_nodes
    updateCheckboxInput(session, "include_unconnected_nodes_temporal", value = FALSE)
    
    updateCheckboxInput(session, "enable_physics_temporal", value = TRUE)
    updateCheckboxInput(session, "previous_nodes", value = FALSE)
    
  }) #update everything to default when there is 
  
  output$network_summary_temporal <- DT::renderDataTable({
    req(summary_table_temporal())
    
    df <- summary_table_temporal()
    if(ncol(df) == 6){
      df <- df[,-c(2)]
      colnames(df) <- c("Pathway ID", "Cluster Group", "q-value", "-log10(q)", "Gene Count")
      colnames(df) <- gsub("\\.", " ", colnames(df))  # Fix column names for display
    }else{
      df <- df[,-c(2)]
      colnames(df) <- c("Pathway ID", "Cluster Group", "q-value", "-log10(q)", "Gene Count","Lost Node")
      colnames(df) <- gsub("\\.", " ", colnames(df))  # Fix column names for display
    }
    
    DT::datatable(
      df,
      options = list(
        pageLength = 10,
        autoWidth = TRUE,
        columnDefs = list(
          list(targets = 1, width = '200px', className = 'dt-wrap')
        )
      ),
      rownames = FALSE,
      escape = FALSE
    ) %>% DT::formatStyle(
      columns = 1,
      `white-space` = "normal",
      `word-wrap` = "break-word"
    )
  }) #visualize table output below network
  
  output$download_summary_temporal <- downloadHandler(
    filename = function() {
      paste0("Pathway_Summary_Temporal", Sys.Date(), ".xlsx")
    },
    content = function(file) {
      df <- summary_table_temporal()
      if(ncol(df) == 6){
        df <- df[,-c(2)]
        colnames(df) <- c("Pathway ID", "Cluster Group", "q-value", "-log10(q)", "Gene Count")
        colnames(df) <- gsub("\\.", " ", colnames(df))  # Fix column names for display
      }else{
        df <- df[,-c(2)]
        colnames(df) <- c("Pathway ID", "Cluster Group", "q-value", "-log10(q)", "Gene Count","Lost Node")
        colnames(df) <- gsub("\\.", " ", colnames(df))  # Fix column names for display
      }
      writexl::write_xlsx(df, file)
    }
  ) #download the network clusterings etc. as a table.
  
  output$download_network_temporal <- downloadHandler(
    filename = function() {
      paste0("Temporal_", Sys.Date(), ".html")
    },
    content = function(file) {
      temp_file <- tempfile(fileext = ".html")
      visSave(vis_network_temporal(), file = temp_file, selfcontained = TRUE)
      
      html <- readLines(temp_file)
      
      style_block <- '<style>
      html, body {
        margin: 0;
        padding: 0;
        width: 100%;
        height: 100%;
        overflow: hidden;
      }
      #htmlwidget_container {
        width: 100vw !important;
        height: 100vh !important;
        position: relative;
      }
      .vis-network {
        width: 100% !important;
        height: 100% !important;
      }
      .vis-manipulation,
      .vis-navigation {
        z-index: 9999 !important;
        position: absolute !important;
        top: 10px;
        right: 10px;
      }
    </style>'
      
      head_index <- grep("<head>", html, fixed = TRUE)
      if (length(head_index) > 0) {
        html <- append(html, style_block, after = head_index)
      }
      
      writeLines(html, file)
    }
  ) #download network as html
  
  
  
  ################################# Sex-Related Networks
  summary_table_sex_aware <- reactiveVal(NULL) #reactive summary_table for the first tab
  vis_network_sex_aware <- reactiveVal(NULL) #reactive for downloading network as HTML
  
  selected_data_sex_aware_xlsx <- reactive({
    req(input$category_dropdown_sex_aware) ## require category_dropdown_top_score from the dropdown menu 
    if(input$metric_sex_aware == "NAFLD Activity Score"){
      if(input$direction_dropdown_sex_aware == "All Genes"){
        female_filename_xlsx <- paste0("female_NAS_overview_msigdb_", input$category_dropdown_sex_aware, ".xlsx") ## load the selected MSigDB result
        female_filepath_xlsx <- file.path("gender/NAS", female_filename_xlsx) ## get the file path for the file
        male_filename_xlsx <- paste0("male_NAS_overview_msigdb_", input$category_dropdown_sex_aware, ".xlsx")
        male_filepath_xlsx <- file.path("gender/NAS", male_filename_xlsx)
      }else if(input$direction_dropdown_sex_aware == "Only Upregulated"){
        female_filename_xlsx <- paste0("female_NAS_overview_msigdb_", input$category_dropdown_sex_aware, "_upregulated.xlsx") 
        female_filepath_xlsx <- file.path("gender/NAS", female_filename_xlsx) 
        male_filename_xlsx <- paste0("male_NAS_overview_msigdb_", input$category_dropdown_sex_aware, "_upregulated.xlsx") 
        male_filepath_xlsx <- file.path("gender/NAS", male_filename_xlsx) 
      }else{
        female_filename_xlsx <- paste0("female_NAS_overview_msigdb_", input$category_dropdown_sex_aware, "_downregulated.xlsx") 
        female_filepath_xlsx <- file.path("gender/NAS", female_filename_xlsx) 
        male_filename_xlsx <- paste0("male_NAS_overview_msigdb_", input$category_dropdown_sex_aware, "_downregulated.xlsx") 
        male_filepath_xlsx <- file.path("gender/NAS", male_filename_xlsx) 
      } 
    }else{
      if(input$direction_dropdown_sex_aware == "All Genes"){
        female_filename_xlsx <- paste0("female_fibrosis_overview_msigdb_", input$category_dropdown_sex_aware, ".xlsx") ## load the selected MSigDB result
        female_filepath_xlsx <- file.path("gender/Fibrosis", female_filename_xlsx) ## get the file path for the file
        male_filename_xlsx <- paste0("male_fibrosis_overview_msigdb_", input$category_dropdown_sex_aware, ".xlsx")
        male_filepath_xlsx <- file.path("gender/Fibrosis", male_filename_xlsx)
      }else if(input$direction_dropdown_sex_aware == "Only Upregulated"){
        female_filename_xlsx <- paste0("female_fibrosis_overview_msigdb_", input$category_dropdown_sex_aware, "_upregulated.xlsx") 
        female_filepath_xlsx <- file.path("gender/Fibrosis", female_filename_xlsx) 
        male_filename_xlsx <- paste0("male_fibrosis_overview_msigdb_", input$category_dropdown_sex_aware, "_upregulated.xlsx") 
        male_filepath_xlsx <- file.path("gender/Fibrosis", male_filename_xlsx) 
      }else{
        female_filename_xlsx <- paste0("female_fibrosis_overview_msigdb_", input$category_dropdown_sex_aware, "_downregulated.xlsx") 
        female_filepath_xlsx <- file.path("gender/Fibrosis", female_filename_xlsx) 
        male_filename_xlsx <- paste0("male_fibrosis_overview_msigdb_", input$category_dropdown_sex_aware, "_downregulated.xlsx") 
        male_filepath_xlsx <- file.path("gender/Fibrosis", male_filename_xlsx) 
      } 
    }
    if (!file.exists(female_filepath_xlsx) | !file.exists(male_filepath_xlsx)) {
      validate(need(FALSE, paste("No file found for", input$category_dropdown_sex_aware)))
      return(NULL)  # never actually reached because validate() stops
    }

    female_df_xlsx <- read_xlsx(female_filepath_xlsx)
    male_df_xlsx <- read_xlsx(male_filepath_xlsx)
    req("qvalue" %in% colnames(female_df_xlsx))  # Check qvalue exists
    req("qvalue" %in% colnames(male_df_xlsx))  # Check qvalue exists
    
    female_df_xlsx <- female_df_xlsx %>%
      filter(-log10(qvalue) >= input$qvalue_threshold_sex_aware) # Filter by threshold on -log10(qvalue)
    male_df_xlsx <- male_df_xlsx %>%
      filter(-log10(qvalue) >= input$qvalue_threshold_sex_aware) # Filter by threshold on -log10(qvalue)
    female_df_xlsx$sex = "female"
    male_df_xlsx$sex = "male"
    df_xlsx <- bind_rows(female_df_xlsx, male_df_xlsx)
    return(df_xlsx)
  })
  
  selected_data_sex_aware <- reactive({
    req(input$category_dropdown_sex_aware, selected_data_sex_aware_xlsx()) ## require category_dropdown_top_score from the dropdown menu 
    if(input$metric_sex_aware == "NAFLD Activity Score"){
      if(input$direction_dropdown_sex_aware == "All Genes"){
        edge_filename_rds <- paste0("female_male_precalculated_edges_overview_msigdb_", input$category_dropdown_sex_aware, ".rds") ## load the selected MSigDB result
        edge_filepath_rds <- file.path("gender/NAS/network", edge_filename_rds) ## get the file path for the file
        node_filename_rds <- paste0("female_male_nodes_overview_msigdb_", input$category_dropdown_sex_aware, ".rds") ## load the selected MSigDB result
        node_filepath_rds <- file.path("gender/NAS/network", node_filename_rds) ## get the file path for the file
      }else if(input$direction_dropdown_sex_aware == "Only Upregulated"){
        edge_filename_rds <- paste0("female_male_precalculated_edges_overview_msigdb_", input$category_dropdown_sex_aware, "_upregulated.rds") ## load the selected MSigDB result
        edge_filepath_rds <- file.path("gender/NAS/network", edge_filename_rds) ## get the file path for the file
        node_filename_rds <- paste0("female_male_nodes_overview_msigdb_", input$category_dropdown_sex_aware, "_upregulated.rds") ## load the selected MSigDB result
        node_filepath_rds <- file.path("gender/NAS/network", node_filename_rds) ## get the file path for the file
      }else{
        edge_filename_rds <- paste0("female_male_precalculated_edges_overview_msigdb_", input$category_dropdown_sex_aware, "_downregulated.rds") ## load the selected MSigDB result
        edge_filepath_rds <- file.path("gender/NAS/network", edge_filename_rds) ## get the file path for the file
        node_filename_rds <- paste0("female_male_nodes_overview_msigdb_", input$category_dropdown_sex_aware, "_downregulated.rds") ## load the selected MSigDB result
        node_filepath_rds <- file.path("gender/NAS/network", node_filename_rds) ## get the file path for the file
      } 
    }else{
      if(input$direction_dropdown_sex_aware == "All Genes"){
        edge_filename_rds <- paste0("female_male_precalculated_edges_overview_msigdb_", input$category_dropdown_sex_aware, ".rds") ## load the selected MSigDB result
        edge_filepath_rds <- file.path("gender/Fibrosis/network", edge_filename_rds) ## get the file path for the file
        node_filename_rds <- paste0("female_male_nodes_overview_msigdb_", input$category_dropdown_sex_aware, ".rds") ## load the selected MSigDB result
        node_filepath_rds <- file.path("gender/Fibrosis/network", node_filename_rds) ## get the file path for the file
      }else if(input$direction_dropdown_sex_aware == "Only Upregulated"){
        edge_filename_rds <- paste0("female_male_precalculated_edges_overview_msigdb_", input$category_dropdown_sex_aware, "_upregulated.rds") ## load the selected MSigDB result
        edge_filepath_rds <- file.path("gender/Fibrosis/network", edge_filename_rds) ## get the file path for the file
        node_filename_rds <- paste0("female_male_nodes_overview_msigdb_", input$category_dropdown_sex_aware, "_upregulated.rds") ## load the selected MSigDB result
        node_filepath_rds <- file.path("gender/Fibrosis/network", node_filename_rds) ## get the file path for the file
      }else{
        edge_filename_rds <- paste0("female_male_precalculated_edges_overview_msigdb_", input$category_dropdown_sex_aware, "_downregulated.rds") ## load the selected MSigDB result
        edge_filepath_rds <- file.path("gender/Fibrosis/network", edge_filename_rds) ## get the file path for the file
        node_filename_rds <- paste0("female_male_nodes_overview_msigdb_", input$category_dropdown_sex_aware, "_downregulated.rds") ## load the selected MSigDB result
        node_filepath_rds <- file.path("gender/Fibrosis/network", node_filename_rds) ## get the file path for the file
      } 
    }
    
    if (!file.exists(edge_filepath_rds) | !file.exists(node_filepath_rds)) {
      validate(need(FALSE, paste("No file found for", input$category_dropdown_sex_aware)))
      return(NULL)  # never actually reached because validate() stops
    }
    
    edge_df_rds <- readRDS(edge_filepath_rds) ## get the file
    node_df_rds <- readRDS(node_filepath_rds) ## get the file
    
    df_xlsx <- selected_data_sex_aware_xlsx() #get msigdb output
    
    edge_df_rds <- edge_df_rds[edge_df_rds$from %in% df_xlsx$ID, ] #filter for q value since some of msigdb will be eliminated
    edge_df_rds <- edge_df_rds[edge_df_rds$to %in% df_xlsx$ID,] #filter for q value since some of msigdb will be eliminated
    node_df_rds <- node_df_rds[node_df_rds$id %in%  df_xlsx$ID,]
    #It is important to note that above code contains a filtering preference for pathways that exist in both female and male. For such pathways
    #Having one sex beeting the threshold is enough to keep the node
    
    list(
      edge_selected_data_sex_aware = edge_df_rds,
      node_selected_data_sex_aware = node_df_rds
    )
  })
  
  observeEvent(selected_data_sex_aware(), {
    df <- selected_data_sex_aware()$edge_selected_data_sex_aware #get the current rds
    if (is.null(df) || nrow(df) < 2) return() # return null if empty or 1 patheay
    max_weight <- max(df$weight, na.rm = TRUE) #get the maximum overlap
    updateSliderInput(
      session,
      inputId = "shared_gene_threshold_sex_aware",
      min = 1,
      max = max_weight, #update maximum of the shared gene slider to maximum current overlap.
      value = min(3, max_weight)
    )
  })
  
  half_svg_data_uri <- function(left_hex = "#ec4899", right_hex = "#3b82f6") {
    svg <- sprintf(
      '<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 100 100">
       <defs><clipPath id="c"><circle cx="50" cy="50" r="50"/></clipPath></defs>
       <g clip-path="url(#c)">
         <rect x="0" y="0" width="50" height="100" fill="%s"/>
         <rect x="50" y="0" width="50" height="100" fill="%s"/>
       </g>
       <circle cx="50" cy="50" r="49" fill="none" stroke="#333333" stroke-width="2"/>
     </svg>', left_hex, right_hex
    )
    paste0("data:image/svg+xml;utf8,", utils::URLencode(svg, reserved = TRUE))
  } #This function is used to create half blue half ping node coloring. It was made by ChatGPT
  
  output$sexnet <- renderVisNetwork({
    req(selected_data_sex_aware(), selected_data_sex_aware_xlsx()) ## require selected_data_top_score reactive value
    
    edges <- selected_data_sex_aware()$edge_selected_data_sex_aware  ## assign to a variable
    
    library(dplyr)
    
    edges <- edges %>%
      mutate(
        color = case_when(
          # male NAs → decide by female counts
          is.na(male_up) & is.na(male_down) & !is.na(female_up) & !is.na(female_down) & female_up > female_down ~ "#D55E00",
          is.na(male_up) & is.na(male_down) & !is.na(female_up) & !is.na(female_down) & female_down > female_up ~ "#0072B2",
          
          # female NAs → decide by male counts
          is.na(female_up) & is.na(female_down) & !is.na(male_up) & !is.na(male_down) & male_up > male_down ~ "#D55E00",
          is.na(female_up) & is.na(female_down) & !is.na(male_up) & !is.na(male_down) & male_down > male_up ~ "#0072B2",
          
          # otherwise keep whatever color was there
          TRUE ~ color
        )
      ) #new edge color logic on top of the preprocessed one!
  
    
    nodes <-selected_data_sex_aware()$node_selected_data_sex_aware
    df_sex_aware <- selected_data_sex_aware_xlsx() ## assign to a variable
    df_sex_aware$gene_list <- strsplit(df_sex_aware$geneID, "/") ## seperate genes to get them as list
    
    if (nrow(edges) < 2) {
      showNotification("Not enough pathways to compute edges.", type = "warning")
      return(NULL)
    } #make sure there is enough data after filtering

    if (input$edge_filter_method_sex_aware == "Jaccard Index") {
      edges$title <- paste0(
        "Jaccard Index (total): ", round(edges$jaccard, 3), "<br>",
        "Jaccard Index (female): ", round(edges$jaccard_female, 3), "<br>",
        "Jaccard Index (male): ",round(edges$jaccard_male, 3), "<br>",
         edges$female_up, " genes upregulated in females", "<br>",
         edges$female_down, " genes downregulated in females", "<br>",
         edges$male_up, " genes upregulated in males", "<br>",
         edges$male_down, " genes downregulated in males"
      ) #titles for jaccard option
    } else {
      edges$title <- paste0(
        edges$weight, " Shared genes", "<br>",
        edges$female_up, " genes upregulated in females", "<br>",
        edges$female_down, " genes downregulated in females", "<br>",
        edges$male_up, " genes upregulated in males", "<br>",
        edges$male_down, " genes downregulated in males"
      )
    } #titles for shared gene option
    
    if (input$edge_filter_method_sex_aware == "Shared Genes") {
      edges <- edges[edges$weight >= input$shared_gene_threshold_sex_aware, ]
    } else if (input$edge_filter_method_sex_aware == "Jaccard Index") {
      edges <- edges[edges$jaccard >= input$jaccard_threshold_sex_aware, ]
    } #the filtering is done over the unified jaccard value not to diff sex individually!
    
    df_sex_aware$qvalue[df_sex_aware$qvalue == 0] <- 1e-300 # Avoid log10(0)
    
    if (input$include_unconnected_nodes_sex_aware) {
      used_ids <- unique(c(edges$from, edges$to))
      df_sex_aware <- df_sex_aware[df_sex_aware$ID %in% used_ids, ]
    } #if checkbox is checked there is no unconnected nodes (I know the name indicates opposite I will fix the naming but it works)
    
    df_merged_sex_aware <- df_sex_aware %>%
      group_by(ID) %>%
      summarise(
        qvalue = if(n() > 1) max(qvalue) else qvalue,   # take max if duplicate
        sex = if(n_distinct(sex) > 1) "both" else first(sex),  # if both sexes present, label "both"
        gene_list = if(n() > 1)  list(unique(unlist(gene_list))) else list(first(gene_list))
      ) %>%
      ungroup() #CHAT GPT IS USED TO MERGE COMMON PATHWAYS KEEPING THE MAX QVALUE (SO IF EXIST IN BOTH CLOSER TO BEING NON-SIGNIFICANT IS USED)
    gene_counts <- lengths(df_merged_sex_aware$gene_list)
    
    node_sizes <- -log10(df_merged_sex_aware$qvalue) #node size = log transformed q value
    node_sizes_scaled <- scales::rescale(node_sizes, to = c(10, 40))  # adjust range if needed
    
    # Convert edge list to igraph for clustering
    if (nrow(edges) > 0 && input$clustering_method_sex_aware != "No Clustering") {
      edge_mat <- as.matrix(edges[, c("from", "to")])
      g <- igraph::graph_from_edgelist(edge_mat, directed = FALSE)
      
      cluster_result <- switch(input$clustering_method_sex_aware,
                               "Louvain" = igraph::cluster_louvain(g),
                               "Edge Betweenness" = igraph::cluster_edge_betweenness(g),
                               "Label Propagation" = igraph::cluster_label_prop(g))
      
      membership <- igraph::membership(cluster_result)
      cluster_df <- data.frame(
        id = names(membership),
        group = paste0("Cluster ", membership),
        stringsAsFactors = FALSE
      )
    } else {
      cluster_df <- data.frame(
        id = df_merged_sex_aware$ID,
        group = "Unclustered",
        stringsAsFactors = FALSE
      )
    }
    
    half_img <- half_svg_data_uri("#ec4899", "#3b82f6")
    network_nodes <- unique(data.frame(
      id = df_merged_sex_aware$ID,
      label = df_merged_sex_aware$ID,
      size = node_sizes_scaled,
      title = paste(gene_counts, "genes in pathway"),
      shape = ifelse(df_merged_sex_aware$sex == "both", "image", "dot"),
      image = ifelse(df_merged_sex_aware$sex == "both", half_img, NA),
      color.background = dplyr::case_when(
        df_merged_sex_aware$sex == "female" ~ "#ec4899",
        df_merged_sex_aware$sex == "male"   ~ "#3b82f6",
        TRUE                                 ~ NA_character_
      ),
      color.border = "#333333",
      stringsAsFactors = FALSE
    )) #create network_nodes
    
    # ensure half nodes not too tiny
    half_idx <- which(df_merged_sex_aware$sex == "both")
    if (length(half_idx)) {
      network_nodes$size[half_idx] <- pmax(network_nodes$size[half_idx], 22)
    }
    
    network_nodes <- merge(network_nodes, cluster_df, by = "id", all.x = TRUE) #merge with clustering info previously obtained
    
    # palette for cluster borders (no extra packages)
    clusters <- unique(network_nodes$group)
    clusters[is.na(clusters)] <- "Unclustered"
    network_nodes$group[is.na(network_nodes$group)] <- "Unclustered"
    
    nC <- length(unique(clusters))
    border_cols <- grDevices::hcl(
      h = seq(0, 360, length.out = nC + 1)[seq_len(nC)],
      c = 70, l = 45
    )
    border_map <- setNames(border_cols, sort(unique(clusters)))
    
    # assign border color & widths
    network_nodes$color.border <- unname(border_map[network_nodes$group])
    network_nodes$borderWidth <- 2
    network_nodes$borderWidthSelected <- 3
    
    summary_table_s <- data.frame(
      `Pathway ID` = network_nodes$id,
      `Gene IDs` = sapply(df_merged_sex_aware$gene_list[
        match(network_nodes$id, df_merged_sex_aware$ID)],
        function(v) paste(v, collapse = "/")),
      `Cluster Group` = network_nodes$group,
      `q-value` = df_merged_sex_aware$qvalue[match(network_nodes$id, df_merged_sex_aware$ID)],
      `-log10(q)` = node_sizes[match(network_nodes$id, df_merged_sex_aware$ID)],
      `Gene Count` = gene_counts[match(network_nodes$id, df_merged_sex_aware$ID)],
      `Sex` = df_merged_sex_aware$sex[match(network_nodes$id, df_merged_sex_aware$ID)],
      stringsAsFactors = FALSE
    ) #create the summary table
    
    # Save table for reuse in UI/download
    summary_table_sex_aware(summary_table_s) #reactive
    updateCheckboxInput(session, "enable_physics_sex_aware", value = TRUE) #checkbox updater when there is a new network
    
    
    network_object <- visNetwork(network_nodes, edges, height = "100%", width = "100%") %>%
      visNodes(size = "size") %>%
      visEdges(smooth = FALSE) %>%
      visOptions(highlightNearest = TRUE, nodesIdSelection = TRUE, selectedBy = "group") %>%
      visOptions(highlightNearest = TRUE, nodesIdSelection = TRUE, manipulation = FALSE,selectedBy = "group") %>%
      visInteraction(dragView = TRUE, zoomView = TRUE, navigationButtons = TRUE) %>%
      visLayout(randomSeed = 123) %>%
      visPhysics(
        solver = "forceAtlas2Based",
        forceAtlas2Based = list(
          gravitationalConstant = -50,
          centralGravity = 0.001,
          springLength = 200,
          springConstant = 0.02,
          damping = 0.4,
          avoidOverlap = 1
        ),
        stabilization = list(
          enabled = TRUE,
          iterations = 200
        )
      ) #visualizes the network including its properties (physics etc.)
    
    # Save for download
    vis_network_sex_aware(network_object) #make network object reactive
    
    # Render to UI
    network_object #visualizes the network
  })
  
  observeEvent(input$enable_physics_sex_aware, {
    visNetworkProxy("sexnet") %>%
      visPhysics(enabled = input$enable_physics_sex_aware)
  }) #if checked stop motion
  
  output$network_summary_sex_aware <- DT::renderDataTable({
    req(summary_table_sex_aware())
    
    df <- summary_table_sex_aware()
    df <- df[,-c(2)]
    colnames(df) <- c("Pathway ID", "Cluster Group", "q-value", "-log10(q)", "Gene Count","Sex")
    colnames(df) <- gsub("\\.", " ", colnames(df))  # Fix column names for display
    
    DT::datatable(
      df,
      options = list(
        pageLength = 10,
        autoWidth = TRUE,
        columnDefs = list(
          list(targets = 1, width = '200px', className = 'dt-wrap')
        )
      ),
      rownames = FALSE,
      escape = FALSE
    ) %>% DT::formatStyle(
      columns = 1,
      `white-space` = "normal",
      `word-wrap` = "break-word"
    )
  }) #visualize table output below network
  
  output$download_summary_sex_aware <- downloadHandler(
    filename = function() {
      paste0("Pathway_Summary_Sex_Aware", Sys.Date(), ".xlsx")
    },
    content = function(file) {
      df <- summary_table_sex_aware()
      df <- df[,-c(2)]
      colnames(df) <- c("Pathway ID", "Cluster Group", "q-value", "-log10(q)", "Gene Count","Sex")
      colnames(df) <- gsub("\\.", " ", colnames(df))  # Fix column names for display
      writexl::write_xlsx(df, file)
    }
  ) #download the network clusterings etc. as a table.
  
  output$download_network_sex_aware <- downloadHandler(
    filename = function() {
      paste0("SexAware_Network_", Sys.Date(), ".html")
    },
    content = function(file) {
      temp_file <- tempfile(fileext = ".html")
      visSave(vis_network_sex_aware(), file = temp_file, selfcontained = TRUE)
      
      html <- readLines(temp_file)
      
      style_block <- '<style>
      html, body {
        margin: 0;
        padding: 0;
        width: 100%;
        height: 100%;
        overflow: hidden;
      }
      #htmlwidget_container {
        width: 100vw !important;
        height: 100vh !important;
        position: relative;
      }
      .vis-network {
        width: 100% !important;
        height: 100% !important;
      }
      .vis-manipulation,
      .vis-navigation {
        z-index: 9999 !important;
        position: absolute !important;
        top: 10px;
        right: 10px;
      }
    </style>'
      
      head_index <- grep("<head>", html, fixed = TRUE)
      if (length(head_index) > 0) {
        html <- append(html, style_block, after = head_index)
      }
      
      writeLines(html, file)
    }
  ) #download network as html
  ############################################ Transcriptome Browser
  
  output$gene_table_browser_pretty <- renderReactable({
    gene_table_browser_df <- if (input$metric_browser == "NAFLD Activity Score") scores else scores_fibrosis
    gene_table_browser_df <- gene_table_browser_df[
      gene_table_browser_df$Gene == input$selected_gene_browser, 
      -ncol(gene_table_browser_df)
    ]
    
    reactable(
      gene_table_browser_df,
      fullWidth = TRUE,
      resizable  = TRUE,
      bordered   = TRUE,
      highlight  = TRUE,
      searchable = FALSE,
      pagination = TRUE,
      defaultPageSize = 5,
      # let columns expand; no maxWidth here
      defaultColDef = colDef(minWidth = 90, align = "right", style = list(whiteSpace = "normal")),
      columns = list(
        gene = colDef(align = "left", sticky = "left", minWidth = 160, style = list(whiteSpace = "normal")),
        global_upregulated             = colDef(name = "Up"),
        global_downregulated           = colDef(name = "Down"),
        individual_upregulated         = colDef(name = "Up"),
        individual_downregulated       = colDef(name = "Down"),
        regression_upregulated         = colDef(name = "Up"),
        regression_downregulated       = colDef(name = "Down"),
        GeneDiseasePattern_upregulated = colDef(name = "Up"),
        GeneDiseasePattern_downregulated=colDef(name = "Down"),
        total_upregulated              = colDef(name = "Up"),
        total_downregulated            = colDef(name = "Down"),
        total                          = colDef(name = "Total")
      ),
      columnGroups = list(
        colGroup(name = "Global",               columns = c("global_upregulated","global_downregulated")),
        colGroup(name = "Individual",           columns = c("individual_upregulated","individual_downregulated")),
        colGroup(name = "Regression",           columns = c("regression_upregulated","regression_downregulated")),
        colGroup(name = "Gene–Disease Pattern", columns = c("GeneDiseasePattern_upregulated","GeneDiseasePattern_downregulated")),
        colGroup(name = "Total",                columns = c("total_upregulated","total_downregulated"))
      )
    )
  })
  
  
  gse_ids <- c("GSE135251","GSE193066","GSE225740","GSE130970",
               "GSE192959","GSE207310","GSE162694","GSE174478","GSE185051")

  base_dir <- "."  # parent directory containing all the GSE folders

  NAS_datasets <- lapply(gse_ids, function(study){
    study_dir <- file.path(base_dir, study)

    f_vst  <- file.path(study_dir, paste0(study, "_vst_matrix_NAS.rds"))
    f_col  <- file.path(study_dir, paste0(study, "_colData_NAS.rds"))
    f_stat <- file.path(study_dir, paste0("NAS_", study, "_regression.csv"))

    if (!file.exists(f_vst) || !file.exists(f_col) || !file.exists(f_stat)) {
      warning("Missing files for ", study)
      return(NULL)
    }

    vst   <- readRDS(f_vst)
    cdat  <- readRDS(f_col)
    stats <- tryCatch(readr::read_csv(f_stat, show_col_types = FALSE), error = function(e) NULL)

    if (!is.null(stats) && !"adj_P_Value" %in% names(stats) && "P_Value" %in% names(stats)) {
      stats$adj_P_Value <- p.adjust(stats$P_Value, method = "BH")
    }

    list(vst = vst, colData = cdat, stats = stats)
  })
  names(NAS_datasets) <- gse_ids
  
  gse_ids <- c("GSE135251","GSE193066","GSE225740","GSE130970","GSE192959","GSE162694","GSE174478")
  
  base_dir <- "."  # parent directory containing all the GSE folders
  
  FIBROSIS_datasets <- lapply(gse_ids, function(study){
    study_dir <- file.path(base_dir, study)
    
    f_vst  <- file.path(study_dir, paste0(study, "_vst_matrix_FIBROSIS.rds"))
    f_col  <- file.path(study_dir, paste0(study, "_colData_FIBROSIS.rds"))
    f_stat <- file.path(study_dir, paste0("FIBROSIS_", study, "_regression.csv"))
    
    if (!file.exists(f_vst) || !file.exists(f_col) || !file.exists(f_stat)) {
      warning("Missing files for ", study)
      return(NULL)
    }
    
    vst   <- readRDS(f_vst)
    cdat  <- readRDS(f_col)
    stats <- tryCatch(readr::read_csv(f_stat, show_col_types = FALSE), error = function(e) NULL)
    
    if (!is.null(stats) && !"adj_P_Value" %in% names(stats) && "P_Value" %in% names(stats)) {
      stats$adj_P_Value <- p.adjust(stats$P_Value, method = "BH")
    }
    
    list(vst = vst, colData = cdat, stats = stats)
  })
  names(FIBROSIS_datasets) <- gse_ids
  
  NAS_gene_disease_patterns <- read_xlsx("NAS_genediseasepatterns_merged_results.xlsx")
  FIBROSIS_gene_disease_patterns <- read_xlsx("FIBROSIS_genediseasepatterns_merged_results.xlsx")
  
  # Precompute once (outside observers)
  genes_nas <- unique(scores$Gene)
  genes_fib <- unique(scores_fibrosis$Gene)
  
  default_nas <- scores$Gene[order(scores$total, decreasing = TRUE)][1]
  default_fib <- scores_fibrosis$Gene[order(scores_fibrosis$total, decreasing = TRUE)][1]
  
  # Initialize on first flush (so default metric fills the box fast)
  session$onFlushed(function() {
    updateSelectizeInput(
      session, "selected_gene_browser",
      choices = genes_nas,
      selected = default_nas,
      server = TRUE
    )
  }, once = TRUE)
  
  # Switch list when metric changes
  observeEvent(input$metric_browser, {
    if (input$metric_browser == "NAFLD Activity Score") {
      updateSelectizeInput(
        session, "selected_gene_browser",
        choices = genes_nas,
        selected = default_nas,
        server = TRUE
      )
    } else {
      updateSelectizeInput(
        session, "selected_gene_browser",
        choices = genes_fib,
        selected = default_fib,
        server = TRUE
      )
    }
  }, ignoreInit = TRUE)
  
  plot_df <- reactive({
    req(input$selected_gene_browser)
    gene <- input$selected_gene_browser

    if(input$metric_browser == "NAFLD Activity Score"){
      reg_dataset <- NAS_datasets
    }else{
      reg_dataset <- FIBROSIS_datasets
    }
    dplyr::bind_rows(lapply(names(reg_dataset), function(study){
      ds <- reg_dataset[[study]]
      if (is.null(ds) || !gene %in% rownames(ds$vst)) return(NULL)
      
      df <- data.frame(
        Sample = colnames(ds$vst),
        expr   = as.numeric(ds$vst[gene, ]),
        group  = ds$colData$score[match(colnames(ds$vst), rownames(ds$colData))],
        study  = study,
        stringsAsFactors = FALSE
      )
      df <- df[!is.na(df$group), ]
      if(input$metric_browser == "NAFLD Activity Score"){
        df$group <- factor(as.integer(as.character(df$group)), levels = 0:4)
      }else{
        df$group <- factor((as.character(df$group)), levels = c("F0_F1","F2","F3","F4"))
      }
      
      stat_row <- NULL
      if (!is.null(ds$stats) && "Gene" %in% names(ds$stats)) {
        stat_row <- ds$stats[ds$stats$Gene == gene, , drop = FALSE]
      }
      df$R2  <- if (!is.null(stat_row) && "R_squared" %in% names(stat_row)) stat_row$R_squared[1] else NA_real_
      df$adj <- if (!is.null(stat_row) && "adj_P_Value" %in% names(stat_row)) stat_row$adj_P_Value[1] else NA_real_
      
      df
    }))
  })

  output$regression_box_plot <- renderUI({
    plotOutput("boxplot_facets", height = 600)
  })

  output$boxplot_facets <- renderPlot({
    df <- plot_df(); req(nrow(df) > 0)

    lab_df <- df %>%
      dplyr::group_by(study) %>%
      dplyr::summarise(y = max(expr, na.rm = TRUE),
                       R2 = first(R2), adj = first(adj),
                       .groups = "drop") %>%
      dplyr::mutate(label = sprintf("R² = %s\nadj-p = %s",
                                    ifelse(is.na(R2), "NA", sprintf("%.3f", R2)),
                                    ifelse(is.na(adj), "NA", format(adj, digits = 2))),
                    x = 4.5)

    ggplot2::ggplot(df, ggplot2::aes(x = group, y = expr, fill = group)) +
      ggplot2::geom_boxplot(outlier.alpha = 0.35, color = "black") +
      ggplot2::facet_wrap(~ study, scales = "free_y") +
      ggplot2::scale_fill_brewer(palette = "Set2") +
      ggplot2::labs(
        title = paste0(input$selected_gene_browser, " VST expression across studies"),
        x = ifelse(input$metric_browser == "NAFLD Activity Score", "Group (0–4)", "Group (F0_F1 - F4)"), y = "VST"
      ) +
      ggplot2::theme_bw(base_size = 14) +
      ggplot2::theme(
        panel.background = element_rect(fill = "#f9f9f9", color = NA),
        plot.background  = element_rect(fill = "#f9f9f9", color = NA),
        strip.background = element_rect(fill = "#e6e6e6", color = NA),
        legend.position = "none"
      ) +
      ggplot2::geom_text(data = lab_df,
                         ggplot2::aes(x = x, y = y, label = label),
                         inherit.aes = FALSE, hjust = 1, vjust = 1, size = 3.5)
  })
  
  pattern_locations <- list(
    list(w = 130, h = 130, x = 105, y = 380), #Pattern 1
    list(w = 130, h = 130, x = 235, y = 380), #Pattern 2
    list(w = 130, h = 130, x = 365, y = 380), #Pattern 3
    list(w = 130, h = 130, x = 105, y = 250), #Pattern 4
    list(w = 130, h = 130, x = 235, y = 250), #Pattern 5
    list(w = 130, h = 130, x = 365, y = 250), #Pattern 6
    list(w = 130, h = 130, x = 105, y = 120), #Pattern 7
    list(w = 130, h = 130, x = 235, y = 120), #Pattern 8 
    list(w = 130, h = 130, x = 365, y = 120) #Pattern 9
  )

  
  # Tile & grid settings
  tile_w <- 256
  tile_h <- 256
  grid_cols <- 4
  grid_rows <- 2
  label_pointsize <- 18
  label_box <- "#00000080"
  label_color <- "white"
  max_cols <- 4   # up to 4 tiles per row
  
  # --------------- HELPERS ----------------
  
  is_blank_spec <- function(spec) {
    (is.numeric(spec) && length(spec) == 1 && spec == 0) ||
      (is.character(spec) && length(spec) == 1 && spec == "0")
  }
  
  # Convert a crop spec into a magick geometry object
  as_geometry <- function(spec) {
    if (is.character(spec)) {
      return(spec)  # already "WxH+X+Y"
    }
    if (is.list(spec) || (is.data.frame(spec) && nrow(spec) == 1)) {
      w <- as.numeric(spec$w); h <- as.numeric(spec$h)
      x <- as.numeric(spec$x); y <- as.numeric(spec$y)
      return(geometry_area(width = w, height = h, x_off = x, y_off = y))
    }
    stop("Unsupported crop spec format. Use \"WxH+X+Y\", list/data.frame {w,h,x,y}, or 0.")
  }
  
  # Make one tile:
  # - if crop_spec == 0 -> blank tile
  # - else crop, resize to tile size, and annotate with label
  make_tile <- function(path, crop_spec, label, tw, th) {
    if (is_blank_spec(crop_spec)) {
      # empty tile, no label
      return(image_blank(tw, th, color = "#f9f9f9"))
    }
    
    img <- image_read(path)
    geom <- as_geometry(crop_spec)
    
    cropped <- image_crop(img, geom)
    cropped <- image_resize(
      cropped,
      geometry_size_pixels(width = tw, height = th, preserve_aspect = FALSE)
    )
    
    image_annotate(
      cropped,
      text = label,
      gravity = "south",
      size = label_pointsize,
      color = label_color,
      boxcolor = label_box
    )
  }
  
  # Append a vector of images horizontally (one row)
  append_row <- function(imgs) image_append(image_join(imgs), stack = FALSE)
  # Append rows vertically to make the final grid
  append_col <- function(rows) image_append(image_join(rows), stack = TRUE)
  
  # --------------- BUILD COLLAGE ----------------
  
  output$collage <- renderImage({
    req(input$selected_gene_browser)
    if (input$metric_browser == "NAFLD Activity Score"){
      gene_disease_selected <- NAS_gene_disease_patterns[NAS_gene_disease_patterns$total_vector == input$selected_gene_browser,3:10]
      img_names <- colnames(gene_disease_selected)
      img_paths <- paste0("gene_disease_patterns_images/", img_names, "_SOM_results.png")
      crop_specs <- ifelse(is.na(gene_disease_selected[1, ]), 0, pattern_locations[as.integer(sub("Patern ", "", gene_disease_selected[1, ]))])
    }else{
      gene_disease_selected <- FIBROSIS_gene_disease_patterns[FIBROSIS_gene_disease_patterns$total_vector == input$selected_gene_browser,3:9]
      img_names <- colnames(gene_disease_selected)
      img_paths <- paste0("gene_disease_patterns_images/fibrosis/", img_names, "_SOM_results.png")
      crop_specs <- ifelse(is.na(gene_disease_selected[1, ]), 0, pattern_locations[as.integer(sub("Patern ", "", gene_disease_selected[1, ]))])
    }
    # --- keep only real tiles (crop_spec != 0) ---
is_blank_spec <- function(spec) {
  (is.numeric(spec) && length(spec) == 1 && spec == 0) ||
    (is.character(spec) && length(spec) == 1 && spec == "0")
}
keep <- !vapply(crop_specs, is_blank_spec, logical(1))
img_paths_kept <- img_paths[keep]
img_names_kept <- img_names[keep]
crop_specs_kept <- crop_specs[keep]

m <- length(img_paths_kept)
if (m == 0) {
  final_img <- image_annotate(
    image_blank(800, 300, color = "#f9f9f9"),
    text = "No images to display",
    gravity = "center", size = 24, color = "gray20"
  )
} else {
  # --- auto-fit grid: up to max_cols per row, rows computed automatically ---
  cols <- min(max_cols, m)
  rows <- ceiling(m / cols)

  # build tiles
  tiles <- lapply(seq_len(m), function(i)
    make_tile(img_paths_kept[i], crop_specs_kept[[i]], img_names_kept[i], tile_w, tile_h))

  # split into rows (no padding—last row may be shorter)
  split_idx <- rep(seq_len(rows), each = cols)[seq_len(m)]
  row_list  <- split(tiles, split_idx)

  # join rows then stack
  row_imgs  <- lapply(row_list, append_row)  # append_row(image_join(imgs), stack=FALSE)
  final_img <- append_col(row_imgs)          # append_col(image_join(rows),  stack=TRUE)
}
    
    tmp <- tempfile(fileext = ".png")
    image_write(final_img, path = tmp, format = "png", density = "600x600")
    list(src = tmp, contentType = "image/png")
  }, deleteFile = TRUE)

# DGEA
x <- readRDS("logfc_upper_matrices.rds")
levels_ord <- x$levels
Mats <- x$matrices

x <- readRDS("fibrosis_logfc_upper_matrices.rds")
levels_ord_fibrosis <- x$levels
Mats_fibrosis <- x$matrices


make_heatmap_upper <- function(mat, clip = 2, show_vals = FALSE) {
  stopifnot(is.matrix(mat))
  
  # --- force numeric matrix safely ---
  rn <- rownames(mat); cn <- colnames(mat)
  mat_num <- suppressWarnings(matrix(
    as.numeric(mat),
    nrow = nrow(mat), ncol = ncol(mat),
    dimnames = list(rn, cn)
  ))
  mat <- mat_num
  mat[!is.finite(mat)] <- NA  # handle Inf, -Inf, NaN
  
  # --- clip guard ---
  if (is.null(clip) || !is.numeric(clip) || length(clip) != 1L) clip <- 2
  
  # Keep only upper triangle if square
  if (nrow(mat) == ncol(mat)) {
    mat[lower.tri(mat, diag = FALSE)] <- NA
  }
  
  # Symmetric clipping around 0
  z <- pmax(pmin(mat, clip), -clip)
  
  # Build hover text
  rv <- rownames(z); cv <- colnames(z)
  hover <- matrix(
    paste0("Case: ", rep(rv, times = ncol(z)),
           "<br>Control: ", rep(cv, each = nrow(z)),
           "<br>logFC: ", sprintf("%.3f", as.vector(z))),
    nrow = nrow(z)
  )
  
  p <- plot_ly(
    x = colnames(z), y = rownames(z),
    z = z, type = "heatmap",
    colorscale = "RdBu", reversescale = TRUE,
    zmid = 0, zmin = -clip, zmax = clip,
    text = hover, hoverinfo = "text",
    showscale = TRUE
  ) %>%
    layout(
      xaxis = list(title = "Controls", tickangle = 45, side = "top"),
      yaxis = list(title = "Cases", autorange = "reversed"),
      plot_bgcolor = "#f9f9f9",
      paper_bgcolor = "#f9f9f9",
      margin = list(l = 80, r = 20, b = 80, t = 40)
    )
  
  if (isTRUE(show_vals)) {
    lab <- ifelse(is.na(z), "", sprintf("%.2f", z))
    p <- p %>% add_annotations(
      x = rep(colnames(z), each = nrow(z)),
      y = rep(rownames(z), times = ncol(z)),
      text = as.vector(lab),
      showarrow = FALSE, font = list(size = 10)
    )
  }
  p
}

# ---- SERVER example ----
output$logfc_heatmap <- renderPlotly({
  req(input$selected_gene_browser)
  if (input$metric_browser == "NAFLD Activity Score"){
    mat <- Mats[[input$selected_gene_browser]]
  }else{
    mat <- Mats_fibrosis[[input$selected_gene_browser]]
  }
  
  
  if (is.null(mat) || !is.matrix(mat)) {
    return(plotly_empty(type = "heatmap") %>%
             layout(title = paste("No valid matrix for", input$selected_gene_browser)))
  }
  
  
  rownames(mat) <- as.character(rownames(mat))
  colnames(mat) <- as.character(colnames(mat))
  make_heatmap_upper(mat, clip = input$hm_clip, show_vals = input$hm_show_vals)
})

## gene list dotplot
indiv_10 <- read_xlsx("nas_individual/Individual_Analysis_1_VS_0.xlsx")
indiv_10 <- indiv_10[indiv_10$DE.invnorm == 1 & indiv_10$DE.fishercomb == 1,]
indiv_10 <- indiv_10[,c("Gene","average_logFC")]
indiv_20 <- read_xlsx("nas_individual/Individual_Analysis_2_VS_0.xlsx")
indiv_20 <- indiv_20[indiv_20$DE.invnorm == 1 & indiv_20$DE.fishercomb == 1,]
indiv_20 <- indiv_20[,c("Gene","average_logFC")]
indiv_30 <- read_xlsx("nas_individual/Individual_Analysis_3_VS_0.xlsx")
indiv_30 <- indiv_30[indiv_30$DE.invnorm == 1 & indiv_30$DE.fishercomb == 1,]
indiv_30 <- indiv_30[,c("Gene","average_logFC")]
indiv_40 <- read_xlsx("nas_individual/Individual_Analysis_4_VS_0.xlsx")
indiv_40 <- indiv_40[indiv_40$DE.invnorm == 1 & indiv_40$DE.fishercomb == 1,]
indiv_40 <- indiv_40[,c("Gene","average_logFC")]

indiv_F2_F0 <- read_xlsx("fibrosis_individual/Individual_Analysis_F2_VS_F0_F1.xlsx")
indiv_F2_F0 <- indiv_F2_F0[indiv_F2_F0$DE.invnorm == 1 & indiv_F2_F0$DE.fishercomb == 1,]

indiv_F3_F0 <- read_xlsx("fibrosis_individual/Individual_Analysis_F3_VS_F0_F1.xlsx")
indiv_F3_F0 <- indiv_F3_F0[indiv_F3_F0$DE.invnorm == 1 & indiv_F3_F0$DE.fishercomb == 1,]

indiv_F4_F0 <- read_xlsx("fibrosis_individual/Individual_Analysis_F4_VS_F0_F1.xlsx")
indiv_F4_F0 <- indiv_F4_F0[indiv_F4_F0$DE.invnorm == 1 & indiv_F4_F0$DE.fishercomb == 1,]


# 1) Read genes from .txt (keep order)
gene_upload_info_browser <- reactive({
  req(input$gene_txt)
  x <- readLines(input$gene_txt$datapath, warn = FALSE, encoding = "UTF-8")
  x <- trimws(x)
  x <- x[nzchar(x)]
  if(length(x) && tolower(x[1]) %in% gene_header_words) x <- x[-1] # drop a header line such as MAPPED_GENE
  x <- unique(x)
  n_uploaded <- length(x)

  human_reference <- if(input$metric_browser == "NAFLD Activity Score") scores$Gene else scores_fibrosis$Gene
  species_info <- resolve_user_species(x, input$metric_browser, sum(x %in% human_reference))

  list(genes = unique(species_info$genes), species = species_info$species, n_uploaded = n_uploaded)
})

genes_uploaded <- reactive({
  gene_upload_info_browser()$genes
})

output$concordance_ui_browser <- renderUI({
  req(input$gene_txt)
    div(
      style = "text-align: center;",
      checkboxInput("concordance_gene_browser", "Include only top scoring genes", value = TRUE)
      # selectInput(
      #   inputId = "heatmap_distance_method",
      #   label = "Distance metric:",
      #   choices = c(
      #     "Pearson correlation" = "correlation",
      #     "Euclidean"           = "euclidean",
      #     "Manhattan"           = "manhattan"
      #   ),
      #   selected = "correlation",
      #   width = "100%"
      # ),
      # selectInput(
      #   inputId = "heatmap_clustering_method",
      #   label = "Clustering method:",
      #   choices = c(
      #     "Average linkage"  = "average",
      #     "Complete linkage" = "complete",
      #     "Single linkage"   = "single"
      #   ),
      #   selected = "average",
      #   width = "100%"
      # )
    )
})

# 2) Gather available contrasts (skip missing objects)
long_logfc <- reactive({
  req(genes_uploaded())

  if (input$metric_browser == "NAFLD Activity Score"){
    dfs <- list(
      "4 vs 0" = get0("indiv_40", ifnotfound = NULL),
      "3 vs 0" = get0("indiv_30", ifnotfound = NULL),
      "2 vs 0" = get0("indiv_20", ifnotfound = NULL),
      "1 vs 0" = get0("indiv_10", ifnotfound = NULL)
    )
  }else{
    dfs <- list(
      "F4 vs 0" = get0("indiv_F4_F0", ifnotfound = NULL),
      "F3 vs 0" = get0("indiv_F3_F0", ifnotfound = NULL),
      "F2 vs 0" = get0("indiv_F2_F0", ifnotfound = NULL)
    )
  }
  
  dfs <- dfs[!vapply(dfs, is.null, logical(1))]
  validate(need(length(dfs) > 0, "No result found."))
  
  bind_rows(lapply(names(dfs), function(nm){
    df <- dfs[[nm]]
    # Expect columns: Gene, average_logFC
    if (!all(c("Gene","average_logFC") %in% names(df))) return(NULL)
    transmute(df,
              Gene     = as.character(.data$Gene),
              contrast = nm,
              logFC    = as.numeric(.data$average_logFC)
    )
  })) |>
    distinct(Gene, contrast, .keep_all = TRUE)
})

# 3) Prepare dotplot data (only rows with values; preserve gene order)
dot_df <- reactive({
  req(long_logfc(), genes_uploaded())
  gs <- genes_uploaded()
  if(!is.null(input$concordance_gene_browser) && isTRUE(input$concordance_gene_browser)){
    if(input$metric_browser == "NAFLD Activity Score"){
      gene_scores <- scores
      cutoff <- quantile(gene_scores$total, 0.95)
      gene_scores <- gene_scores[gene_scores$total >= cutoff,]$Gene
      gs <- gs[gs %in% gene_scores]
    }else if(input$metric_browser == "Fibrosis Stage"){
      gene_scores <- scores_fibrosis
      cutoff <- quantile(gene_scores$total, 0.95)
      gene_scores <- gene_scores[gene_scores$total >= cutoff,]$Gene
      gs <- gs[gs %in% gene_scores]
    }
  }
  print(gs)
  if(length(gs) > 20){
    showNotification("Recommended amount of genes is 20 for this plot!", type = "message")
  }
  # Keep only requested genes; set factor to enforce the y order
  long_logfc() |>
    filter(Gene %in% gs) |>
    mutate(
      Gene = factor(Gene, levels =  unique(gs)),
      abs_logFC = abs(logFC)
    )
})

# 4) Missing gene report
output$missing_genes_msg <- renderText({
  info <- gene_upload_info_browser()
  present <- unique(as.character(dot_df()$Gene))
  missing <- setdiff(info$genes, present)
  msg <- upload_summary_msg(info$n_uploaded, info$species, length(present))
  if (length(missing)) {
    paste0(msg, "\nGenes not found in our datasets: ", paste(missing, collapse = ", "))
  } else msg
})

# 5) Dotplot (bubble plot): color = logFC, size = |logFC|
# helper: sort "a vs b" by numeric a, then b
order_contrasts_numeric <- function(x) {
  # extract first two integers from each string
  nums <- lapply(x, function(s) {
    m <- regmatches(s, gregexpr("-?\\d+", s))[[1]]
    if (length(m) >= 2) c(as.integer(m[1]), as.integer(m[2])) else c(Inf, Inf)
  })
  mat <- do.call(rbind, nums)
  x[order(mat[,1], mat[,2], na.last = TRUE)]
}

output$gene_logfc_dotplot_ui <- renderUI({
  
  req(dot_df())
  
  df <- dot_df()
  
  n_genes <- length(
    unique(as.character(df$Gene))
  )
  
  plot_height <- max(
    450,
    n_genes * 38 + 180
  )
  
  plotlyOutput(
    "gene_logfc_dotplot",
    height = paste0(plot_height, "px")
  )
})


output$gene_logfc_dotplot <- renderPlotly({
  
  req(dot_df())
  
  df <- dot_df()
  
  validate(
    need(nrow(df) > 0, "No matching genes to plot.")
  )
  
  x_levels <- order_contrasts_numeric(
    unique(df$contrast)
  )
  
  lim <- max(
    abs(df$logFC),
    na.rm = TRUE
  )
  
  # Smaller maximum bubble so neighboring genes do not overlap
  size_range <- c(4, 12)
  
  p <- ggplot(
    df,
    aes(
      x = factor(contrast, levels = x_levels),
      y = Gene
    )
  ) +
    
    geom_point(
      aes(
        size = abs_logFC,
        color = logFC
      ),
      alpha = 0.9
    ) +
    
    scale_size(
      range = size_range,
      name = "|logFC|",
      breaks = scales::breaks_pretty(n = 4)
    ) +
    
    scale_color_gradient2(
      name = "logFC",
      low = "#4575b4",
      mid = "#f7f7f7",
      high = "#d73027",
      midpoint = 0,
      limits = c(-lim, lim),
      oob = scales::squish
    ) +
    
    labs(
      x = NULL,
      y = NULL
    ) +
    
    theme_minimal(
      base_size = 13
    ) +
    
    theme(
      panel.background = element_rect(
        fill = NA,
        color = NA
      ),
      plot.background = element_rect(
        fill = NA,
        color = NA
      ),
      axis.text.x = element_text(
        angle = 45,
        hjust = 1,
        vjust = 1
      ),
      panel.grid.major.y = element_blank(),
      panel.grid.minor = element_blank(),
      legend.title = element_text(size = 11),
      legend.text = element_text(size = 10)
    )
  
  ggplotly(
    p,
    tooltip = c(
      "y",
      "x",
      "colour",
      "size"
    )
  ) %>%
    layout(
      autosize = TRUE,
      legend = list(
        orientation = "v"
      ),
      paper_bgcolor = "rgba(0,0,0,0)",
      plot_bgcolor = "rgba(0,0,0,0)"
    )
})
  
# ============================================================
# FULL LOGFC TABLES FOR CLUSTERED GENE HEATMAP
# Do NOT filter for significance here.
# ============================================================

heat_nas_10 <- read_xlsx("nas_individual/Individual_Analysis_1_VS_0.xlsx") %>%
  dplyr::select(Gene, average_logFC)

heat_nas_20 <- read_xlsx("nas_individual/Individual_Analysis_2_VS_0.xlsx") %>%
  dplyr::select(Gene, average_logFC)

heat_nas_30 <- read_xlsx("nas_individual/Individual_Analysis_3_VS_0.xlsx") %>%
  dplyr::select(Gene, average_logFC)

heat_nas_40 <- read_xlsx("nas_individual/Individual_Analysis_4_VS_0.xlsx") %>%
  dplyr::select(Gene, average_logFC)


heat_fib_F2 <- read_xlsx(
  "fibrosis_individual/Individual_Analysis_F2_VS_F0_F1.xlsx"
) %>%
  dplyr::select(Gene, average_logFC)

heat_fib_F3 <- read_xlsx(
  "fibrosis_individual/Individual_Analysis_F3_VS_F0_F1.xlsx"
) %>%
  dplyr::select(Gene, average_logFC)

heat_fib_F4 <- read_xlsx(
  "fibrosis_individual/Individual_Analysis_F4_VS_F0_F1.xlsx"
) %>%
  dplyr::select(Gene, average_logFC)

# ============================================================
# CLUSTERED GENE LIST HEATMAP DATA
# ============================================================

heatmap_long_logfc <- reactive({
  
  req(genes_uploaded())
  
  if (input$metric_browser == "NAFLD Activity Score") {
    
    dfs <- list(
      "1 vs 0" = heat_nas_10,
      "2 vs 0" = heat_nas_20,
      "3 vs 0" = heat_nas_30,
      "4 vs 0" = heat_nas_40
    )
    
  } else {
    
    dfs <- list(
      "F2 vs F0/F1" = heat_fib_F2,
      "F3 vs F0/F1" = heat_fib_F3,
      "F4 vs F0/F1" = heat_fib_F4
    )
    
  }
  
  dplyr::bind_rows(
    lapply(names(dfs), function(nm) {
      
      df <- dfs[[nm]]
      
      df %>%
        dplyr::transmute(
          Gene = as.character(Gene),
          stage = nm,
          logFC = as.numeric(average_logFC)
        )
    })
  ) %>%
    dplyr::distinct(Gene, stage, .keep_all = TRUE)
})

heatmap_df <- reactive({
  
  req(heatmap_long_logfc(), genes_uploaded())
  
  gs <- genes_uploaded()
  
  # Same optional top scoring gene filter as dotplot
  if (
    !is.null(input$concordance_gene_browser) &&
    isTRUE(input$concordance_gene_browser)
  ) {
    
    if (input$metric_browser == "NAFLD Activity Score") {
      
      gene_scores <- scores
      cutoff <- quantile(
        gene_scores$total,
        0.95,
        na.rm = TRUE
      )
      
      top_genes <- gene_scores$Gene[
        gene_scores$total >= cutoff
      ]
      
    } else {
      
      gene_scores <- scores_fibrosis
      cutoff <- quantile(
        gene_scores$total,
        0.95,
        na.rm = TRUE
      )
      
      top_genes <- gene_scores$Gene[
        gene_scores$total >= cutoff
      ]
    }
    
    gs <- gs[gs %in% top_genes]
  }
  
  heatmap_long_logfc() %>%
    dplyr::filter(Gene %in% gs)
})

# ============================================================
# CLUSTERED GENE LIST HEATMAP
# ============================================================
output$gene_logfc_heatmap_ui <- renderUI({
  
  req(heatmap_df())
  
  df <- heatmap_df()
  
  # Count only genes that have at least one real logFC value
  n_genes <- df %>%
    dplyr::group_by(Gene) %>%
    dplyr::summarise(
      has_value = any(is.finite(logFC)),
      .groups = "drop"
    ) %>%
    dplyr::filter(has_value) %>%
    nrow()
  
  # Enough vertical space for gene labels + dendrogram
  plot_height <- max(
    500,
    n_genes * 28 + 220
  )
  
  plotlyOutput(
    "gene_logfc_heatmap",
    height = paste0(plot_height, "px")
  )
})

output$gene_logfc_heatmap <- renderPlotly({
  
  df <- heatmap_df()
  
  validate(
    need(nrow(df) > 0, "No matching genes to plot.")
  )
  
  # ==========================================================
  # STAGE ORDER
  # ==========================================================
  
  if (input$metric_browser == "NAFLD Activity Score") {
    
    stage_order <- c(
      "1 vs 0",
      "2 vs 0",
      "3 vs 0",
      "4 vs 0"
    )
    
  } else {
    
    stage_order <- c(
      "F2 vs F0/F1",
      "F3 vs F0/F1",
      "F4 vs F0/F1"
    )
  }
  
  
  # ==========================================================
  # LONG -> WIDE
  # ==========================================================
  
  wide <- df %>%
    dplyr::select(Gene, stage, logFC) %>%
    tidyr::pivot_wider(
      names_from = stage,
      values_from = logFC
    )
  
  available_stages <- stage_order[
    stage_order %in% colnames(wide)
  ]
  
  validate(
    need(
      length(available_stages) > 1,
      "Not enough stages available for clustering."
    )
  )
  
  wide <- wide %>%
    dplyr::select(
      Gene,
      dplyr::all_of(available_stages)
    )
  
  
  # ==========================================================
  # RAW logFC MATRIX
  # ==========================================================
  
  mat_raw <- as.matrix(
    wide[, available_stages, drop = FALSE]
  )
  
  storage.mode(mat_raw) <- "numeric"
  
  rownames(mat_raw) <- wide$Gene
  
  
  # ==========================================================
  # DROP GENES THAT ARE NA AT EVERY STAGE
  # ==========================================================
  
  keep <- rowSums(is.finite(mat_raw)) > 0
  
  mat_raw <- mat_raw[
    keep,
    ,
    drop = FALSE
  ]
  
  validate(
    need(
      nrow(mat_raw) > 0,
      "None of the selected genes have stage-specific logFC values."
    )
  )
  
  
  # ==========================================================
  # PREPARE MATRIX ONLY FOR CLUSTERING
  #
  # Missing individual stage values are replaced by the
  # gene-specific mean ONLY for calculation of distances.
  #
  # They remain NA in the displayed heatmap.
  # ==========================================================
  
  mat_cluster <- mat_raw
  
  for (i in seq_len(nrow(mat_cluster))) {
    
    vals <- mat_cluster[i, ]
    
    missing <- !is.finite(vals)
    
    if (any(missing)) {
      
      replacement <- mean(
        vals[is.finite(vals)],
        na.rm = TRUE
      )
      
      vals[missing] <- replacement
    }
    
    mat_cluster[i, ] <- vals
  }
  
  
  # ==========================================================
  # Z SCORE EACH GENE ACROSS STAGES
  #
  # This is used ONLY to decide which genes have similar
  # expression patterns.
  # ==========================================================
  
  mat_z <- t(
    scale(
      t(mat_cluster)
    )
  )
  
  mat_z[!is.finite(mat_z)] <- 0
  
  
  # ==========================================================
  # HIERARCHICAL CLUSTERING
  # ==========================================================
  
  if (nrow(mat_z) > 1) {
    
    # --------------------------------------------------------
    # DISTANCE METRIC
    # --------------------------------------------------------
    
    if (input$heatmap_distance_method == "correlation") {
      
      # Pearson correlation distance:
      # identical expression patterns -> distance 0
      # opposite expression patterns  -> distance 2
      
      cor_mat <- cor(
        t(mat_z),
        method = "pearson",
        use = "pairwise.complete.obs"
      )
      
      # Protect against numerical issues
      cor_mat[!is.finite(cor_mat)] <- 0
      
      gene_dist <- as.dist(
        1 - cor_mat
      )
      
    } else {
      
      gene_dist <- dist(
        mat_z,
        method = input$heatmap_distance_method
      )
    }
    
    
    # --------------------------------------------------------
    # LINKAGE METHOD
    # --------------------------------------------------------
    
    gene_hclust <- hclust(
      gene_dist,
      method = input$heatmap_clustering_method
    )
    
    gene_dendrogram <- as.dendrogram(
      gene_hclust
    )
    
  } else {
    
    gene_dendrogram <- FALSE
  }
  
  
  # ==========================================================
  # SYMMETRIC COLOR LIMIT
  #
  # BLUE = downregulated
  # WHITE = unchanged
  # RED = upregulated
  # ==========================================================
  
  lim <- max(
    abs(mat_raw),
    na.rm = TRUE
  )
  
  if (!is.finite(lim) || lim == 0) {
    lim <- 1
  }
  
  
  # ==========================================================
  # HEATMAP + DENDROGRAM
  # ==========================================================
  p <- heatmaply::heatmaply(
    
    mat_raw,
    
    Rowv = gene_dendrogram,
    Colv = FALSE,
    
    colors = colorRampPalette(
      c(
        "#4575b4",
        "#f7f7f7",
        "#d73027"
      )
    )(256),
    
    limits = c(
      -lim,
      lim
    ),
    
    scale = "none",
    
    na.value = "#d9d9d9",
    
    dendrogram = "row",
    
    branches_lwd = 1.2,
    
    xlab = "",
    ylab = "",
    
    main = "Clustered Gene Expression Patterns Across Stages",
    
    label_names = c(
      "Gene",
      "Stage",
      "logFC"
    ),
    
    plot_method = "plotly"
  )
  
  
  # Transparent background
  p <- p %>%
    plotly::layout(
      paper_bgcolor = "rgba(0,0,0,0)",
      plot_bgcolor  = "rgba(0,0,0,0)"
    )
  
  
  # Rename heatmap colorbar to logFC
  p <- plotly::plotly_build(p)
  
  for (i in seq_along(p$x$data)) {
    
    if (
      !is.null(p$x$data[[i]]$type) &&
      p$x$data[[i]]$type == "heatmap"
    ) {
      
      p$x$data[[i]]$colorbar$title <- list(
        text = "logFC"
      )
    }
  }
  
  p
})


  
}

