rm(list=ls())

## The library for dep analysis (2 installation methods)
devtools::install_github("mildpiggy/DEP2", dependencies = TRUE)  # https://github.com/mildpiggy/DEP2
BiocManager::install("DEP") #https://bioconductor.org/packages/devel/bioc/vignettes/DEP/inst/doc/DEP.html#interactive-analysis-using-the-dep-shiny-apps

library(readxl)
library(DEP)
library(dplyr)
library(ggplot2)

getwd()
setwd("D:/NMF_Clustering/input/dep")
Sys.glob("*")

data <- read.csv("proteome_last_update_dep_data.csv")
data <- data.frame(data, stringsAsFactors = FALSE)
#data$X <- NULL

dim(data)
colnames(data)

data$Gene_name %>% duplicated() %>% any()  # TRUE -> overlapped gene
data %>% group_by(Gene_name) %>% summarize(frequency = n()) %>%
  arrange(desc(frequency)) %>% filter(frequency > 1)

# Unique identifiers for each protein are required to further analysis
data_unique <- make_unique(data, "Gene_name", "Protein_ID", delim=";")
## added name and protein ID

data$name %>% duplicated() %>% any()   # FALSE -> no overlapped gene

## data preparation ##
# in Excel, first. change label and group to integer values for analysis
# in Excel, second. change column name of cluster to condition 
clinical <- read.csv("clinical_multi-omics_L-exposed_cluster_5_data.csv", check.names=FALSE)
clinical <- clinical %>% select(label, condition, Time, Age, Sex, group)

LFQ_columns <- grep("LFQ_", colnames(data_unique))
colnames(data_unique) <- sub("^LFQ_", "\\", colnames(data_unique))

clinical <- clinical[order(clinical$group), ]

clinical <- clinical %>% group_by(condition) %>% 
  mutate(replicate = row_number()) %>% ungroup()

clinical <- clinical %>% mutate(condition = recode(
  condition, "Control"="control", "High-risk"="high_risk"))

data_se <- make_se(data_unique, LFQ_columns, clinical)

# missing value check
plot_frequency(data_se) 
plot_numbers(data_se)   
plot_coverage(data_se)   

## normalization
# background corrected & VSN (Variance Stabilizing Transformation) Normalization 
data_norm <- normalize_vsn(data_se)

plot_normalization(data_se, data_norm)

## Differential enrichment analysis - linear model and empirical bayes statistics
# test_diff : using limma 
# all, control(every vs control), manual 
#data_diff <- test_diff(data_norm, type="control", control="Control")
data_diff <- test_diff(data_norm, type="manual", test=c("high_risk_vs_control"))

# Test all possible comparisons
data_diff_all <- test_diff(data_norm, type='all')

# The definition using user-drived cutoff values
dep <- add_rejections(data_diff, alpha=0.05, lfc=log2(1.5))
dep_all <- add_rejections(data_diff_all, alpha=0.05, lfc=log2(1.5))

## Visualization
help(plot_pca)
# PCA plot
plot_pca(dep, x=1, y=2, n=500, point_size=4)
plot_pca(dep, x=1, y=2,indicate=c("condition", "group"), n=500, point_size=4)
plot_pca(dep_all, x=1, y=2, n=500, point_size=4)  

# correlation matrix
plot_cor(dep, significant = TRUE, lower=0, upper=1, pal="Reds")

help(plot_heatmap)
# Heatmap
plot_heatmap(dep, type="centered", kmeans=TRUE, k=6,
             col_limit=4, show_row_names=FALSE, 
             indicate="condition")

plot_heatmap(dep, type="contrast", kmeans=TRUE, 
             k=6, col_limit=10, show_row_names=TRUE)

# added group information
plot_heatmap(dep, type="centered", kmeans=TRUE, k=6,
             col_limit=10, show_row_names=FALSE, 
             indicate="group")

# volcano plot 
plot_volcano(dep, contrast="high_risk_vs_control", label_size=2, add_names=TRUE)  

result_df <- get_results(dep)

DEP_info <- result_df %>% filter(high_risk_vs_control_significant == TRUE)

DEP_info <- DEP_info %>% select(name, ID, high_risk_vs_control_p.val, high_risk_vs_control_p.adj,
                                high_risk_vs_control_significant, high_risk_vs_control_ratio, high_risk_centered)

colnames(DEP_info) <- c("name", "ID", "p.value", "p.value_adj", "significance", "log2FC", "centered")

result_df <- result_df %>% select(name, ID, high_risk_vs_control_p.val, high_risk_vs_control_p.adj,
                                  high_risk_vs_control_significant, high_risk_vs_control_ratio, high_risk_centered)

colnames(result_df) <- c("name", "ID", "p.value", "p.value_adj", "sig", "log2FC", "centered")

lfc_cutoff <- log2(1.5)
alpha_cutoff <- 0.05

result_df <- result_df %>%
  mutate(
    significance = case_when(
      p.value_adj < alpha_cutoff & log2FC > lfc_cutoff ~ "High-risk",
      p.value_adj < alpha_cutoff & log2FC < -lfc_cutoff ~ "Control",
      TRUE ~ "insig"
    )
  )

p <- ggplot(result_df, aes(x = log2FC, y = -log10(p.value_adj), color = significance)) +
  geom_point(alpha = 0.8, size = 2) +
  scale_color_manual(values = c("High-risk" = "red", "Control" = "blue", "insig" = "gray")) +
  geom_vline(xintercept = c(-lfc_cutoff, lfc_cutoff), linetype = "dashed", color = "black") +
  geom_hline(yintercept = c(-log10(alpha_cutoff)), linetype = "dashed", color = "black") +
  geom_text_repel(data = DEP_info,
                  aes(label = name),
                  size = 4,
                  box.padding = 0.4,
                  max.overlaps = 100,
                  segment.color = "black") +
  theme_minimal()+
  labs(
    #title = "High-risk vs Control (Last)",
    x = expression(log[2]~Fold~change),
    y = expression(-log[10]~Adjusted~P~value)
  ) +
  theme(
    plot.title = element_text(
      hjust = 0.5,          
      face = "bold",       
      family = "sans",     
      size = 14           
    ),
    axis.title.x = element_text(size = 14, face = "bold"),
    axis.title.y = element_text(size = 14, face = "bold"),
    legend.title = element_text(size = 13),
    legend.text = element_text(size = 12)
  )

getwd()
setwd("../../output/")

p

ggsave("dep_volcano_plot_high-res_last_cluster_5.png", plot=p,
       width=10.2, height=6,
       dpi=300, units="in", bg="white")

# Adjusted for alignment with the format in Figure 2B. 
dep$condition <- recode(dep$condition, "control" = "Cluster1", "high_risk" = "Cluster3")
dep_all$condition <- recode(dep_all$condition, "control" = "Cluster1", "high_risk" = "Cluster3")

# barplot for interested proteins
plot_single(dep, proteins=c("TPM4", "TLN1", "ITGA2B", "ACTB", "ACTA2", "FL"), type="centered")
plot_single(dep, proteins=c("ATP5F1A", "ATP2A3", "GPX1", "ARPC4", "CTTN", "CALD1"), type="centered")

p2 <- plot_single(dep_all, proteins=c("ATP5F1A", "ATP2A3", "GPX1", "TPM4", "TLN1", "ITGA2B"), type="centered")

p2

ggsave("top6_dep_last_rank_05.png", plot=p2,
       width=9.5, height=6.4,
       dpi=300, units="in", bg="white")

p3 <- plot_single(dep_all, proteins=c("PPBP", "PF4", "FERMT3", "RTN4", "CD36"), type="centered")

p3

ggsave("nmf-marker_dep_last_rank_05.png", plot=p3,
       width=9.5, height=6.0,
       dpi=300, units="in", bg="white")

plot_single(dep, proteins=c("PPBP", "PF4", "FERMT3", "RTN4", "CD36"), type="centered", plot=FALSE)

# Frequency plot
plot_cond(dep_all)

# result table
data_results <- get_results(dep)

data_results %>% filter(significant) %>% nrow() 
colnames(data_results)

write.csv(data_results, "../../output/dep_result/dep_analysis_result_rank_05.csv")
save.image("../rdata/publication/dep_analysis-last_rank_05.RData")
load("D:/NMF_Clustering/rdata/dep_analysis-last_rank_05.RData")

# transform result to dataframe
df_wide <- get_df_wide(dep)
df_long <- get_df_long(dep)

getwd()
setwd("D:/NMF_Clustering/output/dep_result/")

write.csv(df_wide, "./dep_analysis_result_prot-cluster.csv")
write.csv(df_long, "./dep_analysis_result_information.csv")

save(data_se, data_norm, data_diff, dep, file="./dep_analysis_cluster_05.RData")
load("example.RData")