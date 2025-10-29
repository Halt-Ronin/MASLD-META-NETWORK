library(readxl)
library(writexl)
library(scales)

scores <- read_xlsx("overall_meta_analysis_overview_directions.xlsx") #load scores for each gene
#scores$total <- scores$total_downregulated + scores$total_upregulated #calculate the total score
#write_xlsx(scores,"overall_meta_analysis_overview_directions.xlsx")

#C6_downregulated is only 1 pathway

for (y in c("C2","C2_downregulated","C2_upregulated","C3","C3_upregulated","C3_downregulated","C4","C4_downregulated","C4_upregulated","C5","C5_downregulated","C5_upregulated","C6","C6_upregulated","C7","C7_downregulated","C7_upregulated","C8","C8_upregulated","C8_downregulated","H","H_upregulated","H_downregulated")){
  df_top_score <- read_xlsx(paste0("top_score_ora/overview_95_msigdb_", y, ".xlsx"))
  df_top_score$gene_list <- strsplit(df_top_score$geneID, "/") ## seperate genes to get them as list
  
  edge_list <- list()
  k <- 1
  
  for (i in 1:(nrow(df_top_score) - 1)) {
    genes_i <- df_top_score$gene_list[[i]]
    
    for (j in (i + 1):nrow(df_top_score)) {
      genes_j <- df_top_score$gene_list[[j]]
      
      shared_genes <- intersect(genes_i, genes_j)
      union_genes <- union(genes_i, genes_j)
      overlap <- length(shared_genes)
      
      if (overlap > 0) {
        
        # Calculate Jaccard Index only if needed
        jaccard <- overlap / length(union_genes)
        
        
        # Get regulation info for shared genes
        gene_scores <- scores[scores$Gene %in% shared_genes, ]
        
        # Calculate average score only if Score method is selected
        avg_score <- mean(gene_scores$total, na.rm = TRUE)
        
        n_up <- sum(gene_scores$total_upregulated > gene_scores$total_downregulated, na.rm = TRUE)
        n_down <- sum(gene_scores$total_downregulated > gene_scores$total_upregulated, na.rm = TRUE)
        
        # Determine edge color
        edge_color <- if (n_up > n_down) {
          "green"
        } else if (n_down > n_up) {
          "red"
        } else {
          "gray"  # tie or no data
        }
        
        # Build edge
        edge_list[[k]] <- data.frame(
          from = df_top_score$ID[i],
          to = df_top_score$ID[j],
          weight = overlap,
          jaccard = jaccard,
          avg_score = avg_score,
          color = edge_color,
          up = n_up,
          down = n_down,
          stringsAsFactors = FALSE
        )
        
        k <- k + 1
      }
    }
  }
  
  edge_df <- do.call(rbind, edge_list)
  edge_df$width <- scales::rescale(edge_df$weight, to = c(1, 10))
  saveRDS(edge_df, file = paste0("precalculated_edge_ora/precalculated_edges_overview_95_msigdb_", y, ".rds"))
}

############################# Same but with fibrosis ###############

scores <- read_xlsx("overall_meta_analysis_overview_directions_fibrosis.xlsx") #load scores for each gene

#C1_downregulated is only 1 pathway

for (y in c("C2","C3","C4","C5","C6","C7","C8","downregulated_C2","downregulated_C3","downregulated_C4","downregulated_C5","downregulated_C8","downregulated_H","H","upregulated_C2","upregulated_C3","upregulated_C4","upregulated_C5","upregulated_C6","upregulated_C7","upregulated_C8","upregulated_H")){
  df_top_score <- read_xlsx(paste0("top_score_ora/fibrosis/fibrosis_overview_95_msigdb_", y, ".xlsx"))
  df_top_score$gene_list <- strsplit(df_top_score$geneID, "/") ## seperate genes to get them as list
  
  edge_list <- list()
  k <- 1
  
  for (i in 1:(nrow(df_top_score) - 1)) {
    genes_i <- df_top_score$gene_list[[i]]
    
    for (j in (i + 1):nrow(df_top_score)) {
      genes_j <- df_top_score$gene_list[[j]]
      
      shared_genes <- intersect(genes_i, genes_j)
      union_genes <- union(genes_i, genes_j)
      overlap <- length(shared_genes)
      
      if (overlap > 0) {
        
        # Calculate Jaccard Index only if needed
        jaccard <- overlap / length(union_genes)
        
        
        # Get regulation info for shared genes
        gene_scores <- scores[scores$Gene %in% shared_genes, ]
        
        # Calculate average score only if Score method is selected
        avg_score <- mean(gene_scores$total, na.rm = TRUE)
        
        n_up <- sum(gene_scores$total_upregulated > gene_scores$total_downregulated, na.rm = TRUE)
        n_down <- sum(gene_scores$total_downregulated > gene_scores$total_upregulated, na.rm = TRUE)
        
        # Determine edge color
        edge_color <- if (n_up > n_down) {
          "green"
        } else if (n_down > n_up) {
          "red"
        } else {
          "gray"  # tie or no data
        }
        
        # Build edge
        edge_list[[k]] <- data.frame(
          from = df_top_score$ID[i],
          to = df_top_score$ID[j],
          weight = overlap,
          jaccard = jaccard,
          avg_score = avg_score,
          color = edge_color,
          up = n_up,
          down = n_down,
          stringsAsFactors = FALSE
        )
        
        k <- k + 1
      }
    }
  }
  
  edge_df <- do.call(rbind, edge_list)
  edge_df$width <- scales::rescale(edge_df$weight, to = c(1, 10))
  saveRDS(edge_df, file = paste0("precalculated_edge_ora/fibrosis/precalculated_edges_overview_95_msigdb_", y, ".rds"))
}


################## GENDER
common_lists_nas <- c("C2","C3","C4","C5","C7","C8")
common_lists_nas_upregulated <- c("C2","C3","C4","C5","C7","C8")
common_lists_nas_downregulated <- c("C2","C3","C4","C5","C8")
#### NAS
scores <- read_xlsx("overall_meta_analysis_overview_directions.xlsx") 

# small helper: split "A/B/C" -> unique trimmed character vector
split_genes <- function(x) {
  if (is.na(x) || length(x)==0) return(character(0))
  u <- unique(strsplit(x, "/")[[1]])
  trimws(u[nzchar(u)])
}

for (y in common_lists_nas_downregulated) {
  # --- Load female & male tables
  df_f <- read_xlsx(paste0("gender/NAS/female_NAS_overview_msigdb_", y, "_downregulated.xlsx"))
  df_m <- read_xlsx(paste0("gender/NAS/male_NAS_overview_msigdb_",   y, "_downregulated.xlsx"))
  
  # Prepare gene lists
  df_f$gene_list <- lapply(df_f$geneID, split_genes)
  df_m$gene_list <- lapply(df_m$geneID, split_genes)
  
  # Index by pathway ID for quick lookup
  gf <- setNames(df_f$gene_list, df_f$ID)
  gm <- setNames(df_m$gene_list, df_m$ID)
  
  # --- Node colors
  ids_f <- names(gf)
  ids_m <- names(gm)
  all_ids <- sort(unique(c(ids_f, ids_m)))
  
  node_color <- sapply(all_ids, function(id) {
    in_f <- id %in% ids_f
    in_m <- id %in% ids_m
    if (in_f && in_m) "half" else if (in_f) "pink" else "blue"
  }, USE.NAMES = TRUE)
  
  node_df <- data.frame(
    id = all_ids,
    color = unname(node_color),
    stringsAsFactors = FALSE
  )
  
  # --- Build edges over all unique ID pairs
  edge_list <- list()
  k <- 1L
  
  # convenience to fetch genes (empty if missing)
  getf <- function(id) if (id %in% ids_f) gf[[id]] else character(0)
  getm <- function(id) if (id %in% ids_m) gm[[id]] else character(0)
  
  n_all <- length(all_ids)
  if (n_all >= 2) {
    for (ii in 1:(n_all - 1)) {
      for (jj in (ii + 1):n_all) {
        id_i <- all_ids[ii]
        id_j <- all_ids[jj]
        
        gi_f <- getf(id_i); gj_f <- getf(id_j)
        gi_m <- getm(id_i); gj_m <- getm(id_j)
        
        # sex-specific overlaps
        overlap_f <- intersect(gi_f, gj_f)
        overlap_m <- intersect(gi_m, gj_m)
        
        # total overlap = union of sex-specific overlaps (prevents double counting)
        overlap_total <- union(overlap_f, overlap_m)
        weight <- length(overlap_total)
        
        if (weight == 0) next  # skip edges with no shared genes in any sex
        
        # unions for Jaccard denominators
        union_total <- union(union(gi_f, gi_m), union(gj_f, gj_m))
        jaccard_total <- if (length(union_total) > 0) weight / length(union_total) else NA_real_
        
        jaccard_female <- if ((length(gi_f) > 0) && (length(gj_f) > 0)) {
          length(overlap_f) / length(union(gi_f, gj_f))
        } else NA_real_
        
        jaccard_male <- if ((length(gi_m) > 0) && (length(gj_m) > 0)) {
          length(overlap_m) / length(union(gi_m, gj_m))
        } else NA_real_
        
        # scores over total overlap
        gene_scores_total <- scores[scores$Gene %in% overlap_total, , drop = FALSE]
        avg_score <- if (nrow(gene_scores_total) > 0) mean(gene_scores_total$total, na.rm = TRUE) else NA_real_
        
        # female up/down over female overlap
        if (length(overlap_f) > 0) {
          gs_f <- scores[scores$Gene %in% overlap_f, , drop = FALSE]
          female_up   <- sum(gs_f$total_upregulated > gs_f$total_downregulated, na.rm = TRUE)
          female_down <- sum(gs_f$total_downregulated > gs_f$total_upregulated, na.rm = TRUE)
        } else {
          female_up <- NA_integer_; female_down <- NA_integer_
        }
        
        # male up/down over male overlap
        if (length(overlap_m) > 0) {
          gs_m <- scores[scores$Gene %in% overlap_m, , drop = FALSE]
          male_up   <- sum(gs_m$total_upregulated > gs_m$total_downregulated, na.rm = TRUE)
          male_down <- sum(gs_m$total_downregulated > gs_m$total_upregulated, na.rm = TRUE)
        } else {
          male_up <- NA_integer_; male_down <- NA_integer_
        }
        
        # edge color rule
        both_up   <- (!is.na(female_up) && !is.na(female_down) && female_up > female_down) &&
          (!is.na(male_up)   && !is.na(male_down)   && male_up   > male_down)
        both_down <- (!is.na(female_up) && !is.na(female_down) && female_down > female_up) &&
          (!is.na(male_up)   && !is.na(male_down)   && male_down   > male_up)
        
        edge_color <- if (both_up) "green" else if (both_down) "red" else "gray"
        
        edge_list[[k]] <- data.frame(
          from = id_i,
          to   = id_j,
          weight = weight,
          jaccard = jaccard_total,
          jaccard_female = jaccard_female,
          jaccard_male   = jaccard_male,
          avg_score = avg_score,
          female_up = female_up,
          female_down = female_down,
          male_up = male_up,
          male_down = male_down,
          color = edge_color,
          stringsAsFactors = FALSE
        )
        k <- k + 1L
      }
    }
  }
  
  edge_df <- if (length(edge_list)) do.call(rbind, edge_list) else {
    data.frame(from=character(), to=character(), weight=integer(),
               jaccard=numeric(), jaccard_female=numeric(), jaccard_male=numeric(),
               avg_score=numeric(), female_up=integer(), female_down=integer(),
               male_up=integer(), male_down=integer(), color=character(),
               stringsAsFactors = FALSE)
  }
  
  # scale widths from total overlap
  if (nrow(edge_df) > 0) {
    edge_df$width <- rescale(edge_df$weight, to = c(1, 10))
  }
  
  # Save outputs
  dir.create("gender/NAS/network", showWarnings = FALSE, recursive = TRUE)
  saveRDS(edge_df, file = paste0("gender/NAS/network/female_male_precalculated_edges_overview_msigdb_", y, "_downregulated.rds"))
  saveRDS(node_df, file = paste0("gender/NAS/network/female_male_nodes_overview_msigdb_", y, "_downregulated.rds"))
}

### Fibrosis

common_lists_fibrosis <- c("C2","C3","C4","C5","C7","C8","H")
common_lists_fibrosis_upregulated <- c("C2","C3","C4","C5","C6","C7","C8","H")
common_lists_fibrosis_downregulated <- c("C2","C3","C5","C8")


scores <- read_xlsx("overall_meta_analysis_overview_directions_fibrosis.xlsx") #load scores for each gene

# small helper: split "A/B/C" -> unique trimmed character vector
split_genes <- function(x) {
  if (is.na(x) || length(x)==0) return(character(0))
  u <- unique(strsplit(x, "/")[[1]])
  trimws(u[nzchar(u)])
}

for (y in common_lists_fibrosis_downregulated) {
  # --- Load female & male tables
  df_f <- read_xlsx(paste0("gender/Fibrosis/female_fibrosis_overview_msigdb_", y, "_downregulated.xlsx"))
  df_m <- read_xlsx(paste0("gender/Fibrosis/male_fibrosis_overview_msigdb_",   y, "_downregulated.xlsx"))
  
  # Prepare gene lists
  df_f$gene_list <- lapply(df_f$geneID, split_genes)
  df_m$gene_list <- lapply(df_m$geneID, split_genes)
  
  # Index by pathway ID for quick lookup
  gf <- setNames(df_f$gene_list, df_f$ID)
  gm <- setNames(df_m$gene_list, df_m$ID)
  
  # --- Node colors
  ids_f <- names(gf)
  ids_m <- names(gm)
  all_ids <- sort(unique(c(ids_f, ids_m)))
  
  node_color <- sapply(all_ids, function(id) {
    in_f <- id %in% ids_f
    in_m <- id %in% ids_m
    if (in_f && in_m) "half" else if (in_f) "pink" else "blue"
  }, USE.NAMES = TRUE)
  
  node_df <- data.frame(
    id = all_ids,
    color = unname(node_color),
    stringsAsFactors = FALSE
  )
  
  # --- Build edges over all unique ID pairs
  edge_list <- list()
  k <- 1L
  
  # convenience to fetch genes (empty if missing)
  getf <- function(id) if (id %in% ids_f) gf[[id]] else character(0)
  getm <- function(id) if (id %in% ids_m) gm[[id]] else character(0)
  
  n_all <- length(all_ids)
  if (n_all >= 2) {
    for (ii in 1:(n_all - 1)) {
      for (jj in (ii + 1):n_all) {
        id_i <- all_ids[ii]
        id_j <- all_ids[jj]
        
        gi_f <- getf(id_i); gj_f <- getf(id_j)
        gi_m <- getm(id_i); gj_m <- getm(id_j)
        
        # sex-specific overlaps
        overlap_f <- intersect(gi_f, gj_f)
        overlap_m <- intersect(gi_m, gj_m)
        
        # total overlap = union of sex-specific overlaps (prevents double counting)
        overlap_total <- union(overlap_f, overlap_m)
        weight <- length(overlap_total)
        
        if (weight == 0) next  # skip edges with no shared genes in any sex
        
        # unions for Jaccard denominators
        union_total <- union(union(gi_f, gi_m), union(gj_f, gj_m))
        jaccard_total <- if (length(union_total) > 0) weight / length(union_total) else NA_real_
        
        jaccard_female <- if ((length(gi_f) > 0) && (length(gj_f) > 0)) {
          length(overlap_f) / length(union(gi_f, gj_f))
        } else NA_real_
        
        jaccard_male <- if ((length(gi_m) > 0) && (length(gj_m) > 0)) {
          length(overlap_m) / length(union(gi_m, gj_m))
        } else NA_real_
        
        # scores over total overlap
        gene_scores_total <- scores[scores$Gene %in% overlap_total, , drop = FALSE]
        avg_score <- if (nrow(gene_scores_total) > 0) mean(gene_scores_total$total, na.rm = TRUE) else NA_real_
        
        # female up/down over female overlap
        if (length(overlap_f) > 0) {
          gs_f <- scores[scores$Gene %in% overlap_f, , drop = FALSE]
          female_up   <- sum(gs_f$total_upregulated > gs_f$total_downregulated, na.rm = TRUE)
          female_down <- sum(gs_f$total_downregulated > gs_f$total_upregulated, na.rm = TRUE)
        } else {
          female_up <- NA_integer_; female_down <- NA_integer_
        }
        
        # male up/down over male overlap
        if (length(overlap_m) > 0) {
          gs_m <- scores[scores$Gene %in% overlap_m, , drop = FALSE]
          male_up   <- sum(gs_m$total_upregulated > gs_m$total_downregulated, na.rm = TRUE)
          male_down <- sum(gs_m$total_downregulated > gs_m$total_upregulated, na.rm = TRUE)
        } else {
          male_up <- NA_integer_; male_down <- NA_integer_
        }
        
        # edge color rule
        both_up   <- (!is.na(female_up) && !is.na(female_down) && female_up > female_down) &&
          (!is.na(male_up)   && !is.na(male_down)   && male_up   > male_down)
        both_down <- (!is.na(female_up) && !is.na(female_down) && female_down > female_up) &&
          (!is.na(male_up)   && !is.na(male_down)   && male_down   > male_up)
        
        edge_color <- if (both_up) "green" else if (both_down) "red" else "gray"
        
        edge_list[[k]] <- data.frame(
          from = id_i,
          to   = id_j,
          weight = weight,
          jaccard = jaccard_total,
          jaccard_female = jaccard_female,
          jaccard_male   = jaccard_male,
          avg_score = avg_score,
          female_up = female_up,
          female_down = female_down,
          male_up = male_up,
          male_down = male_down,
          color = edge_color,
          stringsAsFactors = FALSE
        )
        k <- k + 1L
      }
    }
  }
  
  edge_df <- if (length(edge_list)) do.call(rbind, edge_list) else {
    data.frame(from=character(), to=character(), weight=integer(),
               jaccard=numeric(), jaccard_female=numeric(), jaccard_male=numeric(),
               avg_score=numeric(), female_up=integer(), female_down=integer(),
               male_up=integer(), male_down=integer(), color=character(),
               stringsAsFactors = FALSE)
  }
  
  # scale widths from total overlap
  if (nrow(edge_df) > 0) {
    edge_df$width <- rescale(edge_df$weight, to = c(1, 10))
  }
  
  # Save outputs
  dir.create("gender/Fibrosis/network", showWarnings = FALSE, recursive = TRUE)
  saveRDS(edge_df, file = paste0("gender/Fibrosis/network/female_male_precalculated_edges_overview_msigdb_", y, "_downregulated.rds"))
  saveRDS(node_df, file = paste0("gender/Fibrosis/network/female_male_nodes_overview_msigdb_", y, "_downregulated.rds"))
}
