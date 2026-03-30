# install packages
install.packages("ggplot2")
install.packages("gridExtra")

library(BiocManager)
BiocManager::install("NMF")
#BiocManager::install("Biobase")

suppressMessages({
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(gridExtra)
})

library(NMF)

getwd()
setwd("D:/Analysis/NMF_Clustering/output")

Sys.glob("*")  

# load data
omics <- data.frame(fread("multiomics_F-exposed_log2_norm.csv"), row.names=1)
dim(omics)
omics[1:5, 1:5]

clinical <- read.csv("first_all_information.csv", row.names=1)
dim(clinical)
clinical[1:5, 1:ncol(clinical)]

clinical <- clinical %>% filter(Group != "Group5") # excluded non-exposed group
clinical <- clinical[clinical$Label %in% colnames(omics), ]

# save sample list to apply same samples
sample_list <- intersect(colnames(omics), clinical$Label)
length(sample_list) 

typeof(omics)                         
matrix_omics <- as.matrix(omics)      
typeof(matrix_omics)                

## 110,180 element (pro + metabo)

test_nmf <- nmf(matrix_omics, rank=10, seed=1218220)  # test run : 40min (local)
# Seed for reproducible analysis, 1218 is birthday and 220 is lab office number
test_nmf

w_matrics <- test_nmf@fit@W # Omics x 3 matrics => fitting values to 3 rank by genes
head(w_matrics)

h_matrics <- test_nmf@fit@H # 3 matrics x Samples => find cluster numbers by samples
head(h_matrics)

# check cluster
h_matrics[1:3, 1:5] 
apply(h_matrics, 2, which.max)[1:10]  
# mapping to cluster (max value)

getwd()

# actual NMF analysis (run by server) 
nmf_res <- nmf(matrix_omics, rank=2:10, method="brunet", nrun=50, seed=1218220) 

nmf_res 

## Proteome + Metabolome (First)
load("D:/NMF_Clustering/rdata/nmf_first_data.RData")

getwd()

pdf("./check_nmf_optimal_k_value.pdf")
plot(nmf_res)
dev.off()  

pdf("./nmf_cluster_consensusmap_best.pdf")
#pdf("./nmf_cluster_consensusmap_extend.pdf", width = 16, height = 16)
consensusmap(nmf_res)
dev.off()

# save figure (option)
png("NMF_consensusmap_high-res_first.png", width=3000, height=3000, res=300)
consensusmap(nmf_res)
dev.off()

### first-data ### 
## select to optimal cluster number
# 1. Optimal nmf cluster = 3
h_matrics <- nmf_res$fit$'3'@fit@H
nmf_cluster <- apply(h_matrics, 2, which.max)
nmf_cluster 

# ordered data (after 181, the row mixed up)
clinical$Label
head(clinical)
clinical_f <- clinical[order(match(clinical$Label, names(nmf_cluster))), ]
#clinical_f: combine the samples & results (ordered data)
head(clinical_f)
clinical_f$Cluster <- nmf_cluster
head(clinical_f)
write.csv(clinical_f, "D:/NMF_Clustering/input/clinical_info_F-exposed_cluster_3.csv")

nmf_rank3 <- clinical_f

# Second Optimal nmf cluster = 4
h_matrics <- nmf_res$fit$'4'@fit@H
nmf_cluster <- apply(h_matrics, 2, which.max)
nmf_cluster # match sample to cluster

head(clinical)
clinical_f <- clinical[order(match(clinical$Label, names(nmf_cluster))), ]
#clinical_f: combine the samples & results (ordered data)
clinical_f$Cluster <- nmf_cluster
head(clinical_f)

write.csv(clinical_f, "D:/NMF_Clustering/input/clinical_info_F-exposed_cluster_4.csv")

nmf_rank4 <- clinical_f

write.csv(nmf_rank4, "./NMF_First_rank4_data.csv")

# NMF result analysis
group_distribution <- nmf_rank4 %>%
  group_by(Group, Cluster) %>%
  summarize(count = n(), .groups='drop')

print(group_distribution)

cluster_counts <- nmf_rank4 %>%
  group_by(Cluster) %>%
  summarize(total_samples = n(), .groups='drop')

print(cluster_counts)

## visualization for publication is graphic by Graphpad (Prism) ##
# transform to factor data type for boxplot of data distribution
nmf_rank4$Cluster <- factor(nmf_rank4$Cluster)

cluster_colors <- c("#B3CDE3", "#FBB4AE", "#CCEBC5", "#DECBE4")

# FEV1
pdf("./FEV1_boxplot_cluster_03.pdf")
ggplot(nmf_rank4, aes(x = Cluster, y = FEV1_z_F, fill=Cluster)) +
  geom_boxplot(outlier.color="red") +
  scale_fill_manual(values = cluster_colors) +        
  labs(title = "Boxplot of FEV1 score by Cluster",
       x = "Cluster",
       y = "FEV1 (z-score)") +
  theme_minimal()

dev.off()

# FEV1/FVC
pdf("./FEV1FVC_boxplot_cluster_03.pdf")
ggplot(nmf_rank4, aes(x = Cluster, y = FEV1FVC_z_F, fill=Cluster)) +
  geom_boxplot(outlier.color="red") +
  scale_fill_manual(values = cluster_colors) +
  labs(title = "Boxplot of FEV1/FVC score by Cluster",
       x = "Cluster",
       y = "FEV1/FVC (z-score)") +
  theme_minimal()

dev.off()

# FVC
pdf("./FVC_boxplot_cluster_04.pdf")
ggplot(nmf_rank4, aes(x = Cluster, y = FVC_z_F, fill=Cluster)) +
  geom_boxplot(outlier.color="red") +
  scale_fill_manual(values = cluster_colors) +           
  labs(title = "Boxplot of FVC score by Cluster",
       x = "Cluster",
       y = "FVC (z-score)") +
  theme_minimal()

dev.off()

# FEF
pdf("./FEF_boxplot_cluster_04.pdf")
ggplot(nmf_rank4, aes(x = Cluster, y = FEF_z_F, fill=Cluster)) +
  geom_boxplot(outlier.color="red") +
  scale_fill_manual(values = cluster_colors) +           
  labs(title = "Boxplot of FEF score by Cluster",
       x = "Cluster",
       y = "FEF (z-score)") +
  theme_minimal()

dev.off()

## Outlier ##
# FEV1FVC: -2.247 -1.657 (Cluster 1) -2.540 2.506 (Cluster 3)
# FEF: -2.95 (Cluster 1)
# FVC, FEV1: none
