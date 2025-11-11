library(homologene)
library(readxl)
library(writexl)

########## Fibrosis

fibrosis_top_score <- read_xlsx("overall_meta_analysis_overview_directions_fibrosis.xlsx")
cutoff <- quantile(fibrosis_top_score$total, 0.95)
print(cutoff)
fibrosis_top_score <- fibrosis_top_score[fibrosis_top_score$total >= cutoff,]$Gene
orthologs_zebrafish <- diopt(fibrosis_top_score, inTax = 9606, outTax = 7955, delay = 1)
orthologs_zebrafish <- orthologs_zebrafish[orthologs_zebrafish$`Best Score` == "Yes" | orthologs_zebrafish$`Best Score` =="Yes_Adjusted",]
orthologs_mouse <- diopt(fibrosis_top_score, inTax = 9606, outTax = 10090, delay = 1)
orthologs_mouse <- orthologs_mouse[orthologs_mouse$`Best Score` == "Yes" | orthologs_mouse$`Best Score` =="Yes_Adjusted",]

write_xlsx(orthologs_zebrafish, "fibrosis_top_score_zebrafish.xlsx")
write_xlsx(orthologs_mouse, "fibrosis_top_score_mouse.xlsx")

############ NAS

nas_top_score <- read_xlsx("overall_meta_analysis_overview_directions.xlsx")
cutoff <- quantile(nas_top_score$total, 0.95)
print(cutoff)
nas_top_score <- nas_top_score[nas_top_score$total >= cutoff,]$Gene
orthologs_zebrafish <- diopt(nas_top_score, inTax = 9606, outTax = 7955, delay = 1)
orthologs_zebrafish <- orthologs_zebrafish[orthologs_zebrafish$`Best Score` == "Yes" | orthologs_zebrafish$`Best Score` =="Yes_Adjusted",]
orthologs_mouse <- diopt(nas_top_score, inTax = 9606, outTax = 10090, delay = 1)
orthologs_mouse <- orthologs_mouse[orthologs_mouse$`Best Score` == "Yes" | orthologs_mouse$`Best Score` =="Yes_Adjusted",]

write_xlsx(orthologs_zebrafish, "nas_top_score_zebrafish.xlsx")
write_xlsx(orthologs_mouse, "nas_top_score_mouse.xlsx")

