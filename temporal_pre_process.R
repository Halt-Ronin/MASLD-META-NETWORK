library(readxl)
library(writexl)
library(msigdbr)
library(clusterProfiler)
library(dplyr)
modes <- c("total","up","down")
##### NAS

folders <- c("first","second","third","fourth")
comparsions <- c("1_VS_0","2_VS_0","3_VS_0","4_VS_0")
genesets <- c("H","C1","C2","C3","C4","C5","C6","C7","C8")

for (q in 1:4){
  dge_table <- read_xlsx(paste("temporal/NAS/","Individual_Analysis_",comparsions[q], ".xlsx",sep = ""))
  dge_table <- dge_table[dge_table$DE.fishercomb == 1 & dge_table$DE.invnorm == 1,]
  for (j in modes){
    if(j =="up"){
      dge_table_last <- dge_table[dge_table$signFC == 1,]
    }else if(j == "down"){
      dge_table_last <- dge_table[dge_table$signFC == -1,]
    }else{
      dge_table_last <- dge_table
    }
for (i in genesets){
  msig_hallmark <- msigdbr(species = "Homo sapiens", category = i)
  hallmark_list <- msig_hallmark[, c("gs_name", "gene_symbol")]
  ora_hallmark <- enricher(gene = dge_table_last$Gene,
                           TERM2GENE = hallmark_list,
                           pAdjustMethod = "BH",
                           qvalueCutoff = 0.05)
  data_framed_ora <- as.data.frame(ora_hallmark)
  if (nrow(data_framed_ora) > 1){
    write_xlsx(data_framed_ora, paste("temporal/NAS/",folders[q],"/",comparsions[q],"_",i,"_",j,".xlsx", sep=""))
  }
}
  }
}

##### Fibrosis 


folders <- c("first","second","third")
comparsions <- c("F2_VS_F0_F1","F3_VS_F0_F1","F4_VS_F0_F1")
genesets <- c("H","C1","C2","C3","C4","C5","C6","C7","C8")

for (q in 1:3){
  dge_table <- read_xlsx(paste("temporal/fibrosis/","Individual_Analysis_",comparsions[q], ".xlsx",sep = ""))
  dge_table <- dge_table[dge_table$DE.fishercomb == 1 & dge_table$DE.invnorm == 1,]
  for (j in modes){
    if(j =="up"){
      dge_table_last <- dge_table[dge_table$signFC == 1,]
    }else if(j == "down"){
      dge_table_last <- dge_table[dge_table$signFC == -1,]
    }else{
      dge_table_last <- dge_table
    }
    for (i in genesets){
      msig_hallmark <- msigdbr(species = "Homo sapiens", category = i)
      hallmark_list <- msig_hallmark[, c("gs_name", "gene_symbol")]
      ora_hallmark <- enricher(gene = dge_table_last$Gene,
                               TERM2GENE = hallmark_list,
                               pAdjustMethod = "BH",
                               qvalueCutoff = 0.05)
      data_framed_ora <- as.data.frame(ora_hallmark)
      if (nrow(data_framed_ora) > 1){
        write_xlsx(data_framed_ora, paste("temporal/fibrosis/",folders[q],"/",comparsions[q],"_",i,"_",j,".xlsx", sep=""))
      }
    }
  }
}

############ NETWORK

modes <- c("total","up","down")

###### NAS
folders <- c("first","second","third","fourth")
comparisons <- c("1_VS_0","2_VS_0","3_VS_0","4_VS_0")
genesets <- c("H","C1","C2","C3","C4","C5","C6","C7","C8")


for (y in genesets){
  for (x in modes){
    for (z in 1:4){
      file_path <- paste0("temporal/NAS/", folders[z], "/", comparisons[z], "_", y, "_", x, ".xlsx")
      out_rds   <- paste0("temporal/NAS/network/precalculated_", comparisons[z], "_", y, "_", x, ".rds")
      
      if (file.exists(file_path) && !file.exists(out_rds)) {
        
        print(comparisons[z]); print(x); print(y)
        
        df_enrichment <- read_xlsx(file_path)
        dge_table     <- read_xlsx(paste0("temporal/NAS/","Individual_Analysis_", comparisons[z], ".xlsx"))
        
        # expect columns: ID, geneID (slash-separated)
        df_enrichment$gene_list <- strsplit(df_enrichment$geneID, "/", fixed = TRUE)
        
        edge_list <- list()
        k <- 1
        
        # pairwise over rows
        if (nrow(df_enrichment) >= 2) {
          for (i in 1:(nrow(df_enrichment) - 1)) {
            genes_i <- df_enrichment$gene_list[[i]]
            
            for (j in (i + 1):nrow(df_enrichment)) {
              genes_j <- df_enrichment$gene_list[[j]]
              
              shared_genes <- intersect(genes_i, genes_j)
              if (length(shared_genes) == 0) next
              
              union_genes <- union(genes_i, genes_j)
              overlap     <- length(shared_genes)
              jaccard     <- overlap / length(union_genes)
              
              # regulation of shared genes (expect columns: Gene, signFC)
              gene_scores <- dge_table[dge_table$Gene %in% shared_genes, ]
              
              n_up   <- sum(gene_scores$signFC == 1,  na.rm = TRUE)
              n_down <- sum(gene_scores$signFC == -1, na.rm = TRUE)
              
              edge_color <- if (n_up > n_down) "green" else if (n_down > n_up) "red" else "gray"
              
              edge_list[[k]] <- data.frame(
                from    = df_enrichment$ID[i],
                to      = df_enrichment$ID[j],
                weight  = overlap,
                jaccard = jaccard,
                color   = edge_color,
                up      = n_up,
                down    = n_down,
                stringsAsFactors = FALSE
              )
              k <- k + 1
            }
          }
        }
        
        # Save if any edges found
        if (length(edge_list)) {
          edge_df <- do.call(rbind, edge_list)
          edge_df$width <- scales::rescale(edge_df$weight, to = c(1, 10))
          saveRDS(edge_df, file = out_rds)
        }
      }
    }
  }
}

### Fibrosis

folders <- c("first","second","third")
comparisons <- c("F2_VS_F0_F1","F3_VS_F0_F1","F4_VS_F0_F1")
genesets <- c("H","C1","C2","C3","C4","C5","C6","C7","C8")

for (y in genesets){
  for (x in modes){
    for (z in 1:3){
      file_path <- paste0("temporal/fibrosis/", folders[z], "/", comparisons[z], "_", y, "_", x, ".xlsx")
      out_rds   <- paste0("temporal/fibrosis/network/precalculated_", comparisons[z], "_", y, "_", x, ".rds")
      
      if (file.exists(file_path) && !file.exists(out_rds)) {
        
        print(comparisons[z]); print(x); print(y)
        
        df_enrichment <- read_xlsx(file_path)
        dge_table     <- read_xlsx(paste0("temporal/fibrosis/","Individual_Analysis_", comparisons[z], ".xlsx"))
        
        # expect columns: ID, geneID (slash-separated)
        df_enrichment$gene_list <- strsplit(df_enrichment$geneID, "/", fixed = TRUE)
        
        edge_list <- list()
        k <- 1
        
        # pairwise over rows
        if (nrow(df_enrichment) >= 2) {
          for (i in 1:(nrow(df_enrichment) - 1)) {
            genes_i <- df_enrichment$gene_list[[i]]
            
            for (j in (i + 1):nrow(df_enrichment)) {
              genes_j <- df_enrichment$gene_list[[j]]
              
              shared_genes <- intersect(genes_i, genes_j)
              if (length(shared_genes) == 0) next
              
              union_genes <- union(genes_i, genes_j)
              overlap     <- length(shared_genes)
              jaccard     <- overlap / length(union_genes)
              
              # regulation of shared genes (expect columns: Gene, signFC)
              gene_scores <- dge_table[dge_table$Gene %in% shared_genes, ]
              
              n_up   <- sum(gene_scores$signFC == 1,  na.rm = TRUE)
              n_down <- sum(gene_scores$signFC == -1, na.rm = TRUE)
              
              edge_color <- if (n_up > n_down) "green" else if (n_down > n_up) "red" else "gray"
              
              edge_list[[k]] <- data.frame(
                from    = df_enrichment$ID[i],
                to      = df_enrichment$ID[j],
                weight  = overlap,
                jaccard = jaccard,
                color   = edge_color,
                up      = n_up,
                down    = n_down,
                stringsAsFactors = FALSE
              )
              k <- k + 1
            }
          }
        }
        
        # Save if any edges found
        if (length(edge_list)) {
          edge_df <- do.call(rbind, edge_list)
          edge_df$width <- scales::rescale(edge_df$weight, to = c(1, 10))
          saveRDS(edge_df, file = out_rds)
        }
      }
    }
  }
}



######## RDS MERGE 

modes <- c("total","up","down")

###### NAS
folders <- c("first","second","third","fourth")
comparisons <- c("1_VS_0","2_VS_0","3_VS_0","4_VS_0")
genesets <- c("H","C1","C2","C3","C4","C5","C6","C7","C8")

for (y in genesets){
  for (x in modes){
    temp_list <- list()  # store data for all comparisons
    for (z in 1:4){
      file_path   <-  paste0("temporal/NAS/network/precalculated_", comparisons[z], "_", y, "_", x, ".rds")
      if (file.exists(file_path)){
        a <- readRDS(file_path)
        a$comparison <- comparisons[z]
        temp_list[[length(temp_list) + 1]] <- a
      }
    }
    if (length(temp_list) > 0) {
      merged_df <- bind_rows(temp_list)
      saveRDS(merged_df, paste0("temporal/NAS/merged network/NAS_", y, "_", x, ".rds"))
    }
  }
}
       

### Fibrosis
folders <- c("first","second","third")
comparisons <- c("F2_VS_F0_F1","F3_VS_F0_F1","F4_VS_F0_F1")
genesets <- c("H","C1","C2","C3","C4","C5","C6","C7","C8")


for (y in genesets){
  for (x in modes){
    temp_list <- list()  # store data for all comparisons
    for (z in 1:3){
      file_path   <-  paste0("temporal/fibrosis/network/precalculated_", comparisons[z], "_", y, "_", x, ".rds")
      if (file.exists(file_path)){
        a <- readRDS(file_path)
        a$comparison <- comparisons[z]
        temp_list[[length(temp_list) + 1]] <- a
      }
    }
    if (length(temp_list) > 0) {
      merged_df <- bind_rows(temp_list)
      saveRDS(merged_df, paste0("temporal/fibrosis/merged network/fibrosis_", y, "_", x, ".rds"))
    }
  }
}