#Packages
#install.packages("factoextra")

library(factoextra)
library(cluster)
library(stats)
library(data.table)
library(dplyr)
library(ggplot2)

setwd("D:/")    # Change the directory to your data file storage
getwd()

Sys.glob("*")

# Change the file name to your downloading files' name
omics <- data.frame(fread("./multiomics_F-exposed_log2_norm.csv"), stringsAsFactors = FALSE)
rownames(omics) <- omics$V1
omics[1:5, 1:5]

omics$V1 <- NULL
colnames(omics)

summary(omics)

### Find Optimal K
## Elbow Method
Elbow_withinss <- c()

for(i in 1:12) {
  set.seed(1218220)
  kmeans_cluster <- kmeans(omics, centers=i, iter.max=1000)
  Elbow_withinss[i] <- kmeans_cluster$tot.withinss
}

pdf("../output/optimal_k_elbow_f-exposed_rank4.pdf")
plot(c(1:12), Elbow_withinss, type="b", main="Optimal number of clusters", xlab="Number of clusters", 
     ylab="Total within-cluster sum of squares")
dev.off()

## Silhouette Method
set.seed(1218220)
res <- kmeans(omics, centers=4)

# visualization
sil <- silhouette(res$cluster, dist(omics))

p2 <- fviz_silhouette(sil)

ggsave("../output/optimal_k_silhouette_f-exposed_rank4_resize.png", plot=p2, dpi=300, width=5.5, height=5.1)

# Last data
omics <- data.frame(fread("./multiomics_L-exposed_log2_norm.csv"), stringsAsFactors = FALSE)
rownames(omics) <- omics$V1
omics[1:5, 1:5]

omics$V1 <- NULL
colnames(omics)

summary(omics)

Elbow_withinss <- c()

for(i in 1:10) {
  set.seed(1218220)
  kmeans_cluster <- kmeans(omics, centers=i, iter.max=1000)
  Elbow_withinss[i] <- kmeans_cluster$tot.withinss
}

pdf("../output/optimal_k_elbow_l-exposed_rank5.pdf")
plot(c(1:10), Elbow_withinss, type="b", main="Optimal number of clusters", xlab="Number of clusters", 
     ylab="Total within-cluster sum of squares")
dev.off()

# silhouette score for evaluate k-means clusters
set.seed(1218220)
res <- kmeans(omics, centers=4)

# visualization
sil <- silhouette(res$cluster, dist(omics))

p3 <- fviz_silhouette(sil)

ggsave("../output/optimal_k_silhouette_l-exposed_rank4_resize.png", plot=p3, dpi=300, width=5.5, height=5.1)
