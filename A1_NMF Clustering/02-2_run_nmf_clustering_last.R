# isntall library
install.packages("ggplot2")
install.packages("gridExtra")

library(BiocManager)
BiocManager::install("NMF")

suppressMessages({
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(gridExtra)
})

library(NMF)

getwd()
setwd("D:/NMF_Clustering/output")

Sys.glob("*")  

# Raw data
omics <- data.frame(fread("multiomics_L-exposed_tri_log2_norm.csv"), row.names=1)
dim(omics)
omics[1:5, 1:5]

clinical <- read.csv("Last_all_information.csv", row.names=1)
dim(clinical)
clinical[1:5, 1:ncol(clinical)]

clinical <- clinical %>% filter(Group != "Group5") # excluded non-exposed group
clinical <- clinical[clinical$Label %in% colnames(omics), ]
clinical$Label == "GA249"
# GA249 sample has not transcriptome data

sample_list <- intersect(colnames(omics), clinical$Label)
length(sample_list)  

typeof(omics)                        
matrix_omics <- as.matrix(omics)      
typeof(matrix_omics)            

## 2,256,921 element (trans + pro + metabo)

test_nmf <- nmf(matrix_omics, rank=10, seed=1218220) 
test_nmf

w_matrics <- test_nmf@fit@W # Omics x 3 matrics 
head(w_matrics)

h_matrics <- test_nmf@fit@H # 3 matrics x Samples 
head(h_matrics)

# check cluster
h_matrics[1:3, 1:5] 
apply(h_matrics, 2, which.max)[1:10]  

getwd()
# if you want to backup data, save present environment to rdata before NMF analysis.
save.image("../nmf_data.RData")
load("D:/NMF_Clustering/nmf_data.RData")

# actual NMF analysis -> server
nmf_res <- nmf(matrix_omics, rank=2:10, method="brunet", nrun=50, seed=1218220) 
nmf_res

## Proteome + Transcriptome + Metabolome (2nd analysis)
save.image("../nmf_last-data.RData")
load("D:/NMF_Clustering/rdata/nmf_last-data.RData")

getwd()
#dir.create("../../NMF_Clustering/", showWarnings = FALSE)

pdf("./check_nmf_tri-omics_optimal_k_value.pdf")
plot(nmf_res)
dev.off()  

pdf("./nmf_cluster_tri-omics_consensusmap.pdf")
consensusmap(nmf_res)   
dev.off()

# save figure
png("NMF_consensusmap_high-res_last.png", width=3000, height=3000, res=300)
consensusmap(nmf_res)
dev.off()

# downstream analysis is same to 02-1_run_nmf_clustering_first.R
# However, in each cases, optimal cluster is different 
# Optimal Cluster = 5
h_matrics <- nmf_res$fit$'5'@fit@H
nmf_cluster <- apply(h_matrics, 2, which.max)   
nmf_cluster 

clinical$Label
head(clinical)
clinical_f <- clinical[order(match(clinical$Label, names(nmf_cluster))), ]
head(clinical_f)

clinical_f$Cluster <- nmf_cluster
head(clinical_f)
write.csv(clinical_f, "D:/NMF_Clustering/input/clinical_info_tri_L-exposed_cluster_5.csv")

nmf_rank5 <- clinical_f

# Second Optimal nmf cluster = 4
h_matrics <- nmf_res$fit$'4'@fit@H
nmf_cluster <- apply(h_matrics, 2, which.max)
nmf_cluster 

head(clinical)
clinical_f <- clinical[order(match(clinical$Label, names(nmf_cluster))), ]
clinical_f$Cluster <- nmf_cluster
head(clinical_f)
write.csv(clinical_f, "D:/NMF_Clustering/input/clinical_info_tri_L-exposed_cluster_4.csv")

nmf_rank4 <- clinical_f

### Cluster=5에 대한 결과 분석
group_distribution <- nmf_rank5 %>%
  group_by(Group, Cluster) %>%
  summarize(count = n(), .groups='drop') %>%
  pivot_wider(names_from = Cluster, values_from = count) %>%
  rename(Cluster = Group)

cluster_counts <- nmf_rank5 %>%
  group_by(Cluster) %>%
  summarize(total_samples = n(), .groups='drop')

print(group_distribution)
print(cluster_counts)

## visualization for publication is graphic by Graphpad (Prism) ##
# transform to factor data type for boxplot of data distribution
nmf_rank5$Cluster <- factor(nmf_rank5$Cluster)

cluster_colors <- c("#B3CDE3", "#FBB5AE", "#CCEBC5", "#DECBE5", "#FED9A6")

pdf("./FEV1_tri-omics_boxplot_cluster_05.pdf")
ggplot(nmf_rank5, aes(x = Cluster, y = FEV1_z_L, fill=Cluster)) +
  geom_boxplot(outlier.color="red") +
  scale_fill_manual(values = cluster_colors) +           
  labs(title = "Boxplot of FEV1 score by Cluster",
       x = "Cluster",
       y = "FEV1 (z-score)") +
  theme_minimal()

dev.off()  

# FEV1/FVC
pdf("./FEV1FVC_tri-omics_boxplot_cluster_05.pdf")
ggplot(nmf_rank5, aes(x = Cluster, y = FEV1FVC_z_L, fill=Cluster)) +
  geom_boxplot(outlier.color="red") +
  scale_fill_manual(values = cluster_colors) +
  labs(title = "Boxplot of FEV1/FVC score by Cluster",
       x = "Cluster",
       y = "FEV1/FVC (z-score)") +
  theme_minimal()

dev.off()

# FVC
pdf("./FVC_tri-omics_boxplot_cluster_05.pdf")
ggplot(nmf_rank5, aes(x = Cluster, y = FVC_z_L, fill=Cluster)) +
  geom_boxplot(outlier.color="red") +
  scale_fill_manual(values = cluster_colors) +           
  labs(title = "Boxplot of FVC score by Cluster",
       x = "Cluster",
       y = "FVC (z-score)") +
  theme_minimal()

dev.off()

# FEF
pdf("./FEF_tri-omics_boxplot_cluster_05.pdf")
ggplot(nmf_rank5, aes(x = Cluster, y = FEF_z_L, fill=Cluster)) +
  geom_boxplot(outlier.color="red") +
  scale_fill_manual(values = cluster_colors) +           
  labs(title = "Boxplot of FEF score by Cluster",
       x = "Cluster",
       y = "FEF (z-score)") +
  theme_minimal()

dev.off()

## Outlier ##
# FEV1: -4.415 -3.952 (Cluster 2) -4.489 (Cluster 3)
# FEv1FVC: 4.748 (Cluster 4)
# FVC: -4.771 3.142 (Cluster 2) 3.436 (Cluster 3)
# FEF: none