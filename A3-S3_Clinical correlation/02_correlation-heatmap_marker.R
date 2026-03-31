rm(list=ls())

getwd()

setwd("D:/Analysis/Machine Learning/output/Correlation_log2_and_normalization_data/")

#install.packages("tidyr")

library(pheatmap)
library(RColorBrewer)
library(dplyr)
library(tibble)   
library(tidyr)

Sys.glob("*.csv")

### Data loading
## 1. Severe 
# To analyze data from a specific time point, users can modify the variable name (df <-> df2)
df <- data.frame(read.csv("correlation_sev_markers_and_clinical_first.csv"))
df2 <- data.frame(read.csv("correlation_sev_markers_and_clinical_last.csv"))

r_matrix <- df %>% select(Marker, Clinical, r) %>% pivot_wider(names_from= Clinical, values_from = r) %>%
  column_to_rownames("Marker") %>% as.matrix()

p_matrix <- df %>% select(Marker, Clinical, p) %>% pivot_wider(names_from= Clinical, values_from = p) %>%
  column_to_rownames("Marker") %>% round(., digits=3) %>% as.matrix()

group_df <- data.frame(read.csv("marker_group_information.csv"))

annotation_row <- group_df %>% filter(Marker %in% rownames(r_matrix)) %>%
                                        column_to_rownames("Marker")
annot_col <- list(Analysis=c(NMF = "lightgreen", MOFA="skyblue", ML="salmon"))

r_value <- formatC(r_matrix, format="f", digits=2)
p_value <- formatC(p_matrix, format="f", digits=3)

combined_value <- matrix(paste0(r_value, "\n(", p_value, ")"),
                         nrow = nrow(r_matrix), ncol = ncol(r_matrix), 
                         dimnames = dimnames(r_matrix))

p <- pheatmap(r_matrix, display_numbers = combined_value, 
         cluster_rows=T, cluster_cols=T, 
         color = colorRampPalette(rev(brewer.pal(n=9, name="RdBu")))(100),
         breaks = seq(-1, 1, length.out=101), 
         annotation_row = annotation_row, annotation_colors = annot_col,
         annotation_legend = F,
         #main = "Correlation between markers and clinical variables (Severe)", 
         fontsize_row=10, fontsize_col=10, fontsize_number=8,
         filename="test2_heatmap.png", dpi=600, width=11.8, height=8.8)

#, number_color="black")

p

## 2. Target 
df <- data.frame(read.csv("correlation_tar_markers_and_clinical_first.csv"))
df2 <- data.frame(read.csv("correlation_tar_markers_and_clinical_last.csv"))

r_matrix <- df %>% select(Marker, Clinical, r) %>% pivot_wider(names_from= Clinical, values_from = r) %>%
  column_to_rownames("Marker") %>% as.matrix()

p_matrix <- df %>% select(Marker, Clinical, p) %>% pivot_wider(names_from= Clinical, values_from = p) %>%
  column_to_rownames("Marker") %>% round(., digits=3) %>% as.matrix()

group_df <- data.frame(read.csv("marker_group_information.csv"))

annotation_row <- group_df %>% filter(Marker %in% rownames(r_matrix)) %>%
  column_to_rownames("Marker")
annot_col <- list(Analysis=c(NMF = "lightgreen", MOFA="skyblue", ML="salmon"))

r_value <- formatC(r_matrix, format="f", digits=2)
p_value <- formatC(p_matrix, format="f", digits=3)

combined_value <- matrix(paste0(r_value, "\n(", p_value, ")"),
                         nrow = nrow(r_matrix), ncol = ncol(r_matrix), 
                         dimnames = dimnames(r_matrix))

p2 <- pheatmap(r_matrix, display_numbers = combined_value, 
         cluster_rows=T, cluster_cols=T, 
         color = colorRampPalette(rev(brewer.pal(n=9, name="RdBu")))(100),
         breaks = seq(-0.5, 0.5, length.out=101), 
         annotation_row = annotation_row, annotation_colors = annot_col,
         annotation_legend = F,
         #main = "Correlation between markers and clinical variables (Target)", 
         fontsize_row=10, fontsize_col=10, fontsize_number=8,
         filename="test2_heatmap.png", dpi=600, width=11.8, height=8.8)
#, number_color="black")

p2

## 3. Normal 
df <- data.frame(read.csv("correlation_nor_markers_and_clinical_first.csv"))
df2 <- data.frame(read.csv("correlation_nor_markers_and_clinical_last.csv"))

r_matrix <- df %>% select(Marker, Clinical, r) %>% pivot_wider(names_from= Clinical, values_from = r) %>%
  column_to_rownames("Marker") %>% as.matrix()

p_matrix <- df %>% select(Marker, Clinical, p) %>% pivot_wider(names_from= Clinical, values_from = p) %>%
  column_to_rownames("Marker") %>% round(., digits=3) %>% as.matrix()

group_df <- data.frame(read.csv("marker_group_information.csv"))

annotation_row <- group_df %>% filter(Marker %in% rownames(r_matrix)) %>%
  column_to_rownames("Marker")
annot_col <- list(Analysis=c(NMF = "lightgreen", MOFA="skyblue", ML="salmon"))

## what if the correlation value and p-value is combinated?
r_value <- formatC(r_matrix, format="f", digits=2)
p_value <- formatC(p_matrix, format="f", digits=3)

combined_value <- matrix(paste0(r_value, "\n(", p_value, ")"),
                         nrow = nrow(r_matrix), ncol = ncol(r_matrix), 
                         dimnames = dimnames(r_matrix))

p3 <- pheatmap(r_matrix, display_numbers = combined_value, 
         cluster_rows=T, cluster_cols=T, 
         color = colorRampPalette(rev(brewer.pal(n=9, name="RdBu")))(100),
         breaks = seq(-0.5, 0.5, length.out=101), 
         annotation_row = annotation_row, annotation_colors = annot_col,
         annotation_legend = F,
         #main = "Correlation between markers and clinical variables (Normal)", 
         fontsize_row=10, fontsize_col=10, fontsize_number=8,
         filename="test2_heatmap.png", dpi=600, width=11.8, height=8.8)
#, number_color="black")

p3