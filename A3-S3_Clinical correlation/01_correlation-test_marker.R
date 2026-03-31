rm(list=ls())

#######################
# All marker analysis #
#######################

# first
library(ggplot2)
library(dplyr)
library(psych)
library(grid)

setwd("check_clinical_correlation/")
getwd()
Sys.glob("*")

fdata <- read.csv("multiome_clinical_norm_comp4_first.csv", header=T)
fdata <- as.data.frame(fdata)

rownames(fdata) <- fdata$Label

clinical_list <- c("FEV1_z", "FVC_z", "FEV1FVC_z", "FEF_z")
marker_list <- c("P02775", "P02776", "P16671", "P48059.2", "Q15942", "Q86UX7.Q86UX7.2","Q9NQC3.3",  # NMF
                 "P01023", "P12814", "O95477", "Q15848", "P01019", "P02647", "P02652",     # MOFA
                 "P01814", "met294", "met342")                                             # ML

clinical <- subset(fdata, select=clinical_list)

nor_data <- fdata[fdata$Group == 1, ]   # 20
tar_data <- fdata[fdata$Group == 2, ]   # 40
sev_data <- fdata[fdata$Group == 3, ]   # 10

all_marker <- subset(fdata, select=marker_list)

normal_marker <- subset(nor_data, select=marker_list)
target_marker <- subset(tar_data, select=marker_list)
severe_marker <- subset(sev_data, select=marker_list)

normal_clinical <- subset(nor_data, select=clinical_list)
target_clinical <- subset(tar_data, select=clinical_list)
severe_clinical <- subset(sev_data, select=clinical_list)

corr.test(all_marker, clinical, use='complete', method='spearman', adjust='BH')
corr.test(normal_marker, normal_clinical, use='complete', method='spearman', adjust='BH')
corr.test(target_marker, target_clinical, use='complete', method='spearman', adjust='BH')
corr.test(severe_marker, severe_clinical, use='complete', method='spearman', adjust='BH')

res_all = corr.test(all_marker, clinical, use='complete', method='spearman', adjust='BH')
res_nor = corr.test(normal_marker, normal_clinical, use='complete', method='spearman', adjust='BH')
res_tar = corr.test(target_marker, target_clinical, use='complete', method='spearman', adjust='BH')
res_sev = corr.test(severe_marker, severe_clinical, use='complete', method='spearman', adjust='BH')

res_all$stars

# save result as dataframe type
result_all <- data.frame(
  Marker = rep(rownames(res_all$r), ncol(res_all$r)),
  Clinical = rep(colnames(res_all$r), each=nrow(res_all$r)),
  r = as.vector(res_all$r),
  t = as.vector(res_all$t),
  p = as.vector(res_all$p),
  p.adj = as.vector(res_all$p.adj),
  se = as.vector(res_all$se),
  stars = as.vector(res_all$stars)
)

result_nor <- data.frame(
  Marker = rep(rownames(res_nor$r), ncol(res_nor$r)),
  Clinical = rep(colnames(res_nor$r), each=nrow(res_nor$r)),
  r = as.vector(res_nor$r),
  t = as.vector(res_nor$t),
  p = as.vector(res_nor$p),
  p.adj = as.vector(res_nor$p.adj),
  se = as.vector(res_nor$se),
  stars = as.vector(res_nor$stars)
)

result_tar <- data.frame(
  Marker = rep(rownames(res_tar$r), ncol(res_tar$r)),
  Clinical = rep(colnames(res_tar$r), each=nrow(res_tar$r)),
  r = as.vector(res_tar$r),
  t = as.vector(res_tar$t),
  p = as.vector(res_tar$p),
  p.adj = as.vector(res_tar$p.adj),
  se = as.vector(res_tar$se),
  stars = as.vector(res_tar$stars)
)

result_sev <- data.frame(
  Marker = rep(rownames(res_sev$r), ncol(res_sev$r)),
  Clinical = rep(colnames(res_sev$r), each=nrow(res_sev$r)),
  r = as.vector(res_sev$r),
  t = as.vector(res_sev$t),
  p = as.vector(res_sev$p),
  p.adj = as.vector(res_sev$p.adj),
  se = as.vector(res_sev$se),
  stars = as.vector(res_sev$stars)
)

write.csv(result_all, file="../output/correlation_all_markers_and_clinical_first.csv")
write.csv(result_nor, file="../output/correlation_nor_markers_and_clinical_first.csv")
write.csv(result_tar, file="../output/correlation_tar_markers_and_clinical_first.csv")
write.csv(result_sev, file="../output/correlation_sev_markers_and_clinical_first.csv")

# last data
ldata <- read.csv("multiome_clinical_norm_comp4_last.csv", header=T)
ldata <- as.data.frame(ldata)

rownames(ldata) <- ldata$Label

clinical_list <- c("FEV1_z", "FVC_z", "FEV1FVC_z", "FEF_z")
marker_list <- c("P02775", "P02776", "P16671", "P48059.2", "Q15942", "Q86UX7.Q86UX7.2","Q9NQC3.3",  # NMF
                 "P01023", "P12814", "O95477", "Q15848", "P01019", "P02647", "P02652",     # MOFA
                 "P01814", "met294", "met342")                                             # ML

clinical <- subset(ldata, select=clinical_list)

nor_data <- ldata[ldata$Group == 1, ]   # 20
tar_data <- ldata[ldata$Group == 2, ]   # 40
sev_data <- ldata[ldata$Group == 3, ]   # 10

all_marker <- subset(ldata, select=marker_list)

normal_marker <- subset(nor_data, select=marker_list)
target_marker <- subset(tar_data, select=marker_list)
severe_marker <- subset(sev_data, select=marker_list)

normal_clinical <- subset(nor_data, select=clinical_list)
target_clinical <- subset(tar_data, select=clinical_list)
severe_clinical <- subset(sev_data, select=clinical_list)

corr.test(all_marker, clinical, use='complete', method='spearman', adjust='BH')
corr.test(normal_marker, normal_clinical, use='complete', method='spearman', adjust='BH')
corr.test(target_marker, target_clinical, use='complete', method='spearman', adjust='BH')
corr.test(severe_marker, severe_clinical, use='complete', method='spearman', adjust='BH')

res_all = corr.test(all_marker, clinical, use='complete', method='spearman', adjust='BH')
res_nor = corr.test(normal_marker, normal_clinical, use='complete', method='spearman', adjust='BH')
res_tar = corr.test(target_marker, target_clinical, use='complete', method='spearman', adjust='BH')
res_sev = corr.test(severe_marker, severe_clinical, use='complete', method='spearman', adjust='BH')

res_all$stars

# save result as dataframe type
result_all <- data.frame(
  Marker = rep(rownames(res_all$r), ncol(res_all$r)),
  Clinical = rep(colnames(res_all$r), each=nrow(res_all$r)),
  r = as.vector(res_all$r),
  t = as.vector(res_all$t),
  p = as.vector(res_all$p),
  p.adj = as.vector(res_all$p.adj),
  se = as.vector(res_all$se),
  stars = as.vector(res_all$stars)
)

result_nor <- data.frame(
  Marker = rep(rownames(res_nor$r), ncol(res_nor$r)),
  Clinical = rep(colnames(res_nor$r), each=nrow(res_nor$r)),
  r = as.vector(res_nor$r),
  t = as.vector(res_nor$t),
  p = as.vector(res_nor$p),
  p.adj = as.vector(res_nor$p.adj),
  se = as.vector(res_nor$se),
  stars = as.vector(res_nor$stars)
)

result_tar <- data.frame(
  Marker = rep(rownames(res_tar$r), ncol(res_tar$r)),
  Clinical = rep(colnames(res_tar$r), each=nrow(res_tar$r)),
  r = as.vector(res_tar$r),
  t = as.vector(res_tar$t),
  p = as.vector(res_tar$p),
  p.adj = as.vector(res_tar$p.adj),
  se = as.vector(res_tar$se),
  stars = as.vector(res_tar$stars)
)

result_sev <- data.frame(
  Marker = rep(rownames(res_sev$r), ncol(res_sev$r)),
  Clinical = rep(colnames(res_sev$r), each=nrow(res_sev$r)),
  r = as.vector(res_sev$r),
  t = as.vector(res_sev$t),
  p = as.vector(res_sev$p),
  p.adj = as.vector(res_sev$p.adj),
  se = as.vector(res_sev$se),
  stars = as.vector(res_sev$stars)
)

write.csv(result_all, file="../output/correlation_all_markers_and_clinical_last.csv")
write.csv(result_nor, file="../output/correlation_nor_markers_and_clinical_last.csv")
write.csv(result_tar, file="../output/correlation_tar_markers_and_clinical_last.csv")
write.csv(result_sev, file="../output/correlation_sev_markers_and_clinical_last.csv") 