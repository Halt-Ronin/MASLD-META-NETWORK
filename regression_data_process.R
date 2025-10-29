library(data.table)
library(org.Hs.eg.db)
library(GEOquery)
library(Rtsne)
library(ggfortify)
library(dplyr)
library(gridExtra)
library(sva)
library(limma)
library(DESeq2)
library(edgeR)
library(plotly)
library(rmarkdown)
library(vegan)
library(tidyr)
library(pheatmap)
set.seed(42)

load_dataset <- function(id){
  a <- fread(paste("raw_counts/", id, ".tsv", sep="")) #load the count optained from NCBI-generated RNA-seq
  gene_names <- mapIds(
    org.Hs.eg.db,
    keys = as.character(a$GeneID),
    column = "SYMBOL", 
    keytype = "ENTREZID",
    multiVals = "first"
  ) # Get the gene symbols from ENTREZIDs
  a$GeneID <- gene_names #assign gene names instead of ENTREZID
  colnames(a)[1] <- "Gene" #change the column name
  a <- a[!duplicated(a$Gene),] #remove duplicated gene names
  a <- a[!is.na(a$Gene),] #remove any NAs in gene name
  a <- as.data.frame(a) #turn to df 
  rownames(a) <- a$Gene #assign gene names as row names.
  a$Gene <- NULL #remove column
  return(a) #return fully modified df
}

load_meta <- function(id){
  a  <- getGEO(id, GSEMatrix = TRUE) #get the meta data via accession.
  a  <- pData(a[[1]]) #get the df portion
  return(a) #return it
}

regression <- function(data, meta, col_index, name){
  meta <- meta[rownames(meta) %in% colnames(data),] # only samples where count data is available
  colData <- as.data.frame(meta[,col_index])        # create the design df
  rownames(colData) <- rownames(meta)               # assign sample names as rownames
  print(sum(rownames(colData) != colnames(data)))   # should be 0
  colnames(colData) <- "score"                      # name column
  
  # drop NAs
  colData <- colData[!colData[, 1] == "NA", , drop = FALSE]
  
  # ---- CHANGE: recode scores to 0..4 ----
  colData <- colData %>%
    mutate(
      score = suppressWarnings(as.integer(score)),
      score = case_when(
        score %in% c(0,1,2) ~ 0L,   # collapsed to 0
        score == 3          ~ 1L,
        score == 4          ~ 2L,
        score == 5          ~ 3L,
        score %in% c(6,7,8) ~ 4L,   # collapsed to 4
        TRUE                ~ NA_integer_
      )
    )
  colData <- colData[!is.na(colData$score), , drop = FALSE]
  colData$score <- factor(colData$score, levels = 0:4)
  # ---------------------------------------
  
  # align order
  data <- data[, colnames(data) %in% rownames(colData)] 
  colData <- colData[colnames(data), , drop = FALSE]
  
  
  dds <- DESeqDataSetFromMatrix(countData = data, colData = colData, design = ~ score)
  dds <- dds[rowSums(counts(dds)) > 10, ]
  dds <- DESeq(dds) 
  
  vst_counts <- vst(dds, blind = FALSE)
  vst_matrix <- assay(vst_counts)
  vst_matrix_df <- as.data.frame(vst_matrix)
  
  # --- Save for Shiny ---
  if (!dir.exists(name)) dir.create(name, recursive = TRUE)
  saveRDS(vst_matrix, file = file.path(name, paste0(name, "_vst_matrix_NAS.rds")), compress = "xz")
  saveRDS(colData,    file = file.path(name, paste0(name, "_colData_NAS.rds")),    compress = "xz")
}

GSE <- c("GSE135251","GSE193066","GSE225740","GSE130970","GSE192959","GSE207310","GSE162694","GSE174478","GSE185051")
index <- c(44,45,38,46,41,46,44,42,46)
for(i in 8:9){
  data <- load_dataset(GSE[i])
  meta <- load_meta(GSE[i])
  col_index <- index[i]
  print(colnames(meta)[col_index])
  regression(data, meta, col_index, GSE[i])
}

regression <- function(data, meta, col_index, name){
  meta <- meta[rownames(meta) %in% colnames(data),] # only samples where count data is available
  colData <- as.data.frame(meta[,col_index])        # create the design df
  rownames(colData) <- rownames(meta)               # assign sample names as rownames
  print(sum(rownames(colData) != colnames(data)))   # should be 0
  colnames(colData) <- "score"                      # name column
  
  # drop NAs
  colData <- colData[!colData[, 1] == "NA", , drop = FALSE]
  
  # ---- CHANGE: recode scores to 0..4 ----
colData <- colData %>% mutate(
  score = case_when(
      score %in% c(0,1,"0","1") ~ "F0_F1",
      score %in% c(2,"2") ~ "F2",
      score %in% c(3,"3") ~ "F3",
      score %in% c(4,"4") ~ "F4",
      TRUE ~ "NA"
    ))
  colData <- colData[!colData[, 1] == "NA", , drop = FALSE]
  colData <- colData[!is.na(colData$score), , drop = FALSE]
  colData$score <- factor(colData$score, levels = c("F0_F1","F2","F3","F4"))
  # ---------------------------------------
  
  # align order
  data <- data[, colnames(data) %in% rownames(colData)] 
  colData <- colData[colnames(data), , drop = FALSE]

  dds <- DESeqDataSetFromMatrix(countData = data, colData = colData, design = ~ score)
  dds <- dds[rowSums(counts(dds)) > 10, ]
  dds <- DESeq(dds) 
  
  vst_counts <- vst(dds, blind = FALSE)
  vst_matrix <- assay(vst_counts)
  vst_matrix_df <- as.data.frame(vst_matrix)
  
  # --- Save for Shiny ---
  if (!dir.exists(name)) dir.create(name, recursive = TRUE)
  saveRDS(vst_matrix, file = file.path(name, paste0(name, "_vst_matrix_FIBROSIS.rds")), compress = "xz")
  saveRDS(colData,    file = file.path(name, paste0(name, "_colData_FIBROSIS.rds")),    compress = "xz")
}

GSE <- c("GSE135251","GSE193066","GSE225740","GSE130970","GSE192959","GSE162694","GSE174478")
index <- c(42,44,37,44,40,43,41)
for(i in 7:7){
  data <- load_dataset(GSE[i])
  meta <- load_meta(GSE[i])
  col_index <- index[i]
  print(colnames(meta)[col_index])
  regression(data, meta, col_index, GSE[i])
}
