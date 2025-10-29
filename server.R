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


server <- function(input, output, session) {
  
  scores <- read_xlsx("overall_meta_analysis_overview_directions.xlsx") #load scores for each gene
  scores_fibrosis <- read_xlsx("overall_meta_analysis_overview_directions_fibrosis.xlsx") #load scores for each gene - fibrosis
  
  ################################################ Top-Score Over Representation Analysis ###########################
  
  summary_table_top_score <- reactiveVal(NULL) #reactive summary_table for the first tab
  vis_network_top_score <- reactiveVal(NULL) #reactive for downloading network as HTML
  cluster_rename_map_top_score <- reactiveVal(list()) #reactive list for renaming clusters
  
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
    
    # Convert edge list to igraph for clustering
    if (nrow(edges) > 0 && input$clustering_method_top_score != "No Clustering") {
      edge_mat <- as.matrix(edges[, c("from", "to")])
      g <- igraph::graph_from_edgelist(edge_mat, directed = FALSE)
      
      cluster_result <- switch(input$clustering_method_top_score,
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
        id = df_top_score$ID,
        group = "Unclustered",
        stringsAsFactors = FALSE
      )
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
  }) #update everything to default when there is 
  
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
  
  observe({
    req(input$qvalue_threshold_top_score)
    
    updateSliderInput(
      session,
      inputId = "qvalue_threshold_process_gene",
      min = input$qvalue_threshold_top_score,
      max = 10,
      value = input$qvalue_threshold_top_score
    )
  })
  
  df_process_gene <- reactive({
    if (input$direction_dropdown_process_gene == "Use Clustering from Top-Score Tab" &&
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
      
      gene_nodes <- data.frame(
        id = current_gene_ids,
        label = current_gene_ids,
        group = "gene",
        stringsAsFactors = FALSE
      ) %>%
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
            total_upregulated > total_downregulated ~ "green",
            total_upregulated < total_downregulated ~ "red",
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
      
      gene_nodes <- data.frame(
        id = current_gene_ids,
        label = current_gene_ids,
        group = "gene",
        stringsAsFactors = FALSE
      ) %>%
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
            total_upregulated > total_downregulated ~ "green",
            total_upregulated < total_downregulated ~ "red",
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
    
    gene_summary_table_process_gene <- gene_nodes %>%
      select(
        Gene = id,
        total,
        total_upregulated,
        total_downregulated,
        centrality,
        combined_score
      ) %>%
      arrange(desc(combined_score))
    
    output$gene_table_process_gene <- renderDataTable({
      gene_summary_table_process_gene
    })
    
    updateCheckboxInput(session, "enable_physics_process_gene", value = TRUE) #checkbox updater when there is a new network
    
    visNetwork(nodes, edges_formatted) %>%
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
    
  })
  
  observeEvent(input$enable_physics_process_gene, {
    visNetworkProxy("network_process_gene") %>%
      visPhysics(enabled = input$enable_physics_process_gene)
  }) #if checked stop motion
  
  
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
          is.na(male_up) & is.na(male_down) & !is.na(female_up) & !is.na(female_down) & female_up > female_down ~ "green",
          is.na(male_up) & is.na(male_down) & !is.na(female_up) & !is.na(female_down) & female_down > female_up ~ "red",
          
          # female NAs → decide by male counts
          is.na(female_up) & is.na(female_down) & !is.na(male_up) & !is.na(male_down) & male_up > male_down ~ "green",
          is.na(female_up) & is.na(female_down) & !is.na(male_up) & !is.na(male_down) & male_down > male_up ~ "red",
          
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
genes_uploaded <- reactive({
  req(input$gene_txt)
  x <- readLines(input$gene_txt$datapath, warn = FALSE, encoding = "UTF-8")
  x <- trimws(x)
  x <- x[nzchar(x)]
  # unique() preserves first appearance -> keeps user order
  unique(x)
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
  
  # Keep only requested genes; set factor to enforce the y order
  long_logfc() |>
    filter(Gene %in% gs) |>
    mutate(
      Gene = factor(Gene, levels = gs),
      abs_logFC = abs(logFC)
    )
})

# 4) Missing gene report
output$missing_genes_msg <- renderText({
  req(genes_uploaded())
  present <- unique(as.character(dot_df()$Gene))
  missing <- setdiff(genes_uploaded(), present)
  if (length(missing)) {
    paste0("Genes not found in our datasets: ",
           paste(missing, collapse = ", "))
  } else ""
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

output$gene_logfc_dotplot <- renderPlotly({
  req(dot_df())
  df <- dot_df()
  validate(need(nrow(df) > 0, "No matching genes to plot."))
  
  # x-axis order like: 1 vs 0, 2 vs 0, 3 vs 0, 4 vs 0, 4 vs 1, 4 vs 2, 4 vs 3, ...
  x_levels <- order_contrasts_numeric(unique(df$contrast))
  
  # symmetric color limits around 0
  lim <- max(abs(df$logFC), na.rm = TRUE)
  
  # size range for bubbles
  size_range <- c(4, 18)
  
  p <- ggplot(df, aes(
    x = factor(contrast, levels = x_levels),
    y = Gene
  )) +
    geom_point(aes(size = abs_logFC, color = logFC), alpha = 0.9) +
    scale_size(
      range = size_range, name = "|logFC|",
      breaks = scales::breaks_pretty(n = 4)
    ) +
    scale_color_gradient2(
      name = "logFC",
      low = "#4575b4", mid = "#f7f7f7", high = "#d73027",
      midpoint = 0, limits = c(-lim, lim),
      oob = scales::squish
    ) +
    labs(x = NULL, y = NULL) +
    theme_minimal(base_size = 13) +
    theme(
      panel.background = element_rect(fill = "#f9f9f9", color = NA),
      plot.background  = element_rect(fill = "#f9f9f9", color = NA),
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
      panel.grid.major.y = element_blank(),
      panel.grid.minor = element_blank(),
      legend.title = element_text(size = 11),
      legend.text  = element_text(size = 10)
    )
  
  # make height depend on #genes (about 28px per gene, with padding)
  n_genes <- length(levels(df$Gene))
  h <- max(350, 28 * n_genes + 140)
  
  ggplotly(p, tooltip = c("y","x","colour","size")) %>%
    layout(legend = list(orientation = "v"), height = h)
})
  
  
}

