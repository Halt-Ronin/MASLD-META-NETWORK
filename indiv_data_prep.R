library(readxl)
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
indiv_31 <- read_xlsx("nas_individual/Individual_Analysis_3_VS_1.xlsx")
indiv_31 <- indiv_31[indiv_31$DE.invnorm == 1 & indiv_31$DE.fishercomb == 1,]
indiv_31 <- indiv_31[,c("Gene","average_logFC")]
indiv_32 <- read_xlsx("nas_individual/Individual_Analysis_3_VS_2.xlsx")
indiv_32 <- indiv_32[indiv_32$DE.invnorm == 1 & indiv_32$DE.fishercomb == 1,]
indiv_32 <- indiv_32[,c("Gene","average_logFC")]
indiv_41 <- read_xlsx("nas_individual/Individual_Analysis_4_VS_1.xlsx")
indiv_41 <- indiv_41[indiv_41$DE.invnorm == 1 & indiv_41$DE.fishercomb == 1,]
indiv_41 <- indiv_41[,c("Gene","average_logFC")]
indiv_42 <- read_xlsx("nas_individual/Individual_Analysis_4_VS_2.xlsx")
indiv_42 <- indiv_42[indiv_42$DE.invnorm == 1 & indiv_42$DE.fishercomb == 1,]
indiv_42 <- indiv_42[,c("Gene","average_logFC")]
## 2) Bind into one long table with explicit case/control
pairs <- list(
  `1_0` = indiv_10,
  `2_0` = indiv_20,
  `3_0` = indiv_30,
  `4_0` = indiv_40,
  `3_1` = indiv_31,
  `3_2` = indiv_32,
  `4_1` = indiv_41,
  `4_2` = indiv_42
)

long <- bind_rows(lapply(names(pairs), function(k) {
  df <- pairs[[k]]
  if (!all(c("Gene", "average_logFC") %in% names(df))) return(NULL)
  parts <- strsplit(k, "_", fixed = TRUE)[[1]]
  case <- parts[1]; control <- parts[2]
  df %>%
    transmute(
      Gene = as.character(Gene),
      contrast = k,
      case = case,
      control = control,
      logFC = as.numeric(average_logFC)
    )
}))

# In case of duplicates (same Gene-case-control), keep first or average:
long <- long %>% 
  group_by(Gene, case, control) %>% 
  summarise(logFC = mean(logFC, na.rm = TRUE), .groups = "drop")

## 3) Builder: square matrix with ONLY upper triangle filled
build_matrix_for_gene <- function(gene, long_tbl,
                                  levels = c("4","3","2","1","0"),
                                  diag_value = 0,
                                  triangle = c("upper","lower","full")) {
  triangle <- match.arg(triangle)
  mat <- matrix(NA_real_, nrow = length(levels), ncol = length(levels),
                dimnames = list(levels, levels))
  diag(mat) <- diag_value
  
  sub <- long_tbl[long_tbl$Gene == gene, ]
  if (nrow(sub) > 0) {
    for (k in seq_len(nrow(sub))) {
      i <- match(sub$case[k], levels)
      j <- match(sub$control[k], levels)
      if (!is.na(i) && !is.na(j)) mat[i, j] <- sub$logFC[k]
    }
  }
  
  if (triangle == "upper")  mat[lower.tri(mat, diag = FALSE)] <- NA
  if (triangle == "lower")  mat[upper.tri(mat, diag = FALSE)] <- NA
  mat
}

## 4) Build for all genes and save
levels_ord <- c("4","3","2","1","0")        # rows = case, cols = control
genes <- sort(unique(long$Gene))

matrices <- setNames(
  lapply(genes, function(g)
    build_matrix_for_gene(g, long, levels = levels_ord, diag_value = 0, triangle = "upper")),
  genes
)


saveRDS(
  list(levels = levels_ord, matrices = matrices),
  file = "logfc_upper_matrices.rds"
)


#same for fibrosis


indiv_F2_F0 <- read_xlsx("fibrosis_individual/Individual_Analysis_F2_VS_F0_F1.xlsx")
indiv_F2_F0 <- indiv_F2_F0[indiv_F2_F0$DE.invnorm == 1 & indiv_F2_F0$DE.fishercomb == 1,]

indiv_F3_F0 <- read_xlsx("fibrosis_individual/Individual_Analysis_F3_VS_F0_F1.xlsx")
indiv_F3_F0 <- indiv_F3_F0[indiv_F3_F0$DE.invnorm == 1 & indiv_F3_F0$DE.fishercomb == 1,]

indiv_F3_F2 <- read_xlsx("fibrosis_individual/Individual_Analysis_F3_VS_F2.xlsx")
indiv_F3_F2 <- indiv_F3_F2[indiv_F3_F2$DE.invnorm == 1 & indiv_F3_F2$DE.fishercomb == 1,]

indiv_F4_F0 <- read_xlsx("fibrosis_individual/Individual_Analysis_F4_VS_F0_F1.xlsx")
indiv_F4_F0 <- indiv_F4_F0[indiv_F4_F0$DE.invnorm == 1 & indiv_F4_F0$DE.fishercomb == 1,]

indiv_F4_F2 <- read_xlsx("fibrosis_individual/Individual_Analysis_F4_VS_F2.xlsx")
indiv_F4_F2 <- indiv_F4_F2[indiv_F4_F2$DE.invnorm == 1 & indiv_F4_F2$DE.fishercomb == 1,]

indiv_F4_F3 <- read_xlsx("fibrosis_individual/Individual_Analysis_F4_VS_F3.xlsx")
indiv_F4_F3 <- indiv_F4_F3[indiv_F4_F3$DE.invnorm == 1 & indiv_F4_F3$DE.fishercomb == 1,]

# 2) Bind into one long table with explicit case/control
pairs <- list(
  list(case = "F2", control = "F0_F1", df = indiv_F2_F0),
  list(case = "F3", control = "F0_F1", df = indiv_F3_F0),
  list(case = "F3", control = "F2",    df = indiv_F3_F2),
  list(case = "F4", control = "F0_F1", df = indiv_F4_F0),
  list(case = "F4", control = "F2",    df = indiv_F4_F2),
  list(case = "F4", control = "F3",    df = indiv_F4_F3)
)

long <- bind_rows(lapply(pairs, function(x) {
  df <- x$df
  if (!all(c("Gene", "average_logFC") %in% names(df))) return(NULL)
  transmute(df,
            Gene    = as.character(Gene),
            case    = x$case,
            control = x$control,
            logFC   = as.numeric(average_logFC))
}))

# In case of duplicates (same Gene-case-control), average them
long <- long %>%
  group_by(Gene, case, control) %>%
  summarise(logFC = mean(logFC, na.rm = TRUE), .groups = "drop")

# 3) Builder: square matrix with ONLY upper triangle filled
build_matrix_for_gene <- function(gene, long_tbl,
                                  levels = c("F4","F3","F2","F0_F1"),
                                  diag_value = 0,
                                  triangle = c("upper","lower","full")) {
  triangle <- match.arg(triangle)
  mat <- matrix(NA_real_, nrow = length(levels), ncol = length(levels),
                dimnames = list(levels, levels))
  diag(mat) <- diag_value
  
  sub <- long_tbl[long_tbl$Gene == gene, ]
  if (nrow(sub) > 0) {
    for (k in seq_len(nrow(sub))) {
      i <- match(sub$case[k],    levels)
      j <- match(sub$control[k], levels)
      if (!is.na(i) && !is.na(j)) mat[i, j] <- sub$logFC[k]
    }
  }
  
  if (triangle == "upper")  mat[lower.tri(mat, diag = FALSE)] <- NA
  if (triangle == "lower")  mat[upper.tri(mat, diag = FALSE)] <- NA
  mat
}

# 4) Build for all genes and save
levels_ord <- c("F4","F3","F2","F0_F1")  # rows = case, cols = control (yields upper-tri for your contrasts)
genes <- sort(unique(long$Gene))

matrices <- setNames(
  lapply(genes, function(g)
    build_matrix_for_gene(g, long, levels = levels_ord, diag_value = 0, triangle = "upper")),
  genes
)

saveRDS(
  list(levels = levels_ord, matrices = matrices),
  file = "fibrosis_logfc_upper_matrices.rds"
)