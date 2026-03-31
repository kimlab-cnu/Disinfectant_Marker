# Convert clinical data for use as MOFA metadata
rm(list=ls())

if (!requireNamespace("preprocessCore", quietly = TRUE)) {
  install.packages("preprocessCore")
}

# loading Library
suppressMessages({
  library(data.table)
  library(readxl)
  library(dplyr)
  library(impute)
  library(stringr)
})

setwd("D:/Desktop/Analysis/MOFA/input/")
getwd()
Sys.glob("*")

sampleinfo <- read.csv("SampleInfo.csv")
dim(sampleinfo)
sampleinfo[1:5,]

sampleinfo <- data.frame(sampleinfo, stringsAsFactors = FALSE)  # needed group factoring
sampleinfo$Group <- as.factor(sampleinfo$Group)
sampleinfo$Sex <- ifelse(
  sampleinfo$Sex == 1,
  sampleinfo$Sex <- "M",
  sampleinfo$Sex <- "F"
)
sampleinfo$Sex <- as.factor(sampleinfo$Sex)
sampleinfo$Time <- as.factor(sampleinfo$Time)

str(sampleinfo)

sample_label <- sampleinfo %>% select(Label, Prot_ID, Metabo_ID, Trans_ID, Methyl_ID)
sample_label[1:5, ]

F_info <- sampleinfo %>% filter(Time == 'F') %>% 
  select(Label, Group, Time, Age, Sex, FEV1_z_F, FVC_z_F, FEV1FVC_z_F, FEF_z_F)

L_info <- sampleinfo %>% filter(Time == 'L') %>% 
  select(Label, Group, Time, Age, Sex, FEV1_z_L, FVC_z_L, FEV1FVC_z_L, FEF_z_L)

F_cont_exp <- F_info %>% filter(Group == 'Group5') 
F_exposed <- F_info %>% filter(Group != 'Group5')    # not including group5 (control)

L_cont_exp <- L_info %>% filter(Group == 'Group5') 
L_exposed <- L_info %>% filter(Group != 'Group5')    # not including group5 (control)

write.csv(F_info, "clinical_first_data.csv")
write.csv(L_info, "clinical_last_data.csv")
write.csv(F_exposed, "clinical_first_exposed_data.csv")
write.csv(L_exposed, "clinical_last_exposed_data.csv")

## We need to data alignment (for ordering and downstream analysis)
# Issue: misaligned groups between sample and clinical data
library(dplyr)
library(preprocessCore)

metabo <- data.frame(fread("metabolome_last_norm_data.csv"), stringsAsFactors = FALSE)
methyl <- data.frame(fread("methylome_last_norm_data_t.csv"), stringsAsFactors = FALSE)
prot <- data.frame(fread("proteome_last_impute-norm_data.csv"), stringsAsFactors = FALSE)
trans <- data.frame(fread("transcriptome_last_norm_data_t.csv"), stringsAsFactors = FALSE)

L_clinical <- read.csv("clinical_last_ordered_data.csv")
F_clinical <- read.csv("clinical_first_ordered_data.csv")

#methyl <- t(methyl)
#colnames(methyl) <- methyl[1,]
#write.csv(methyl, "./methylome_last_norm_data_t.csv")

#colnames(trans) <- sub("^X(\\d)", "\\1", colnames(trans))
#trans <- t(trans)
#colnames(trans) <- trans[1,]
#write.csv(trans, "./transcriptome_last_norm_data_t.csv")

rownames(prot) <- prot$Label
rownames(methyl) <- methyl$Label
rownames(trans) <- trans$Label
rownames(metabo) <- metabo$Label

prot$Label <- NULL
methyl$Label <- NULL
trans$Label <- NULL
metabo$Label <- NULL

## sample_label의 Label에 맞게 정렬해주는 코드
reorder_data <- function(data, sample_labels) {
  valid_labels <- sample_labels[sample_labels %in% colnames(data)]
  ordered_data <- data[, match(valid_labels, colnames(data))]
  return(ordered_data)
}

L_proteome <- reorder_data(prot, L_clinical$Label)
L_methylome <- reorder_data(methyl, L_clinical$Label)
L_transcriptome <- reorder_data(trans, L_clinical$Label)
L_metabolome <- reorder_data(metabo, L_clinical$Label)

write.csv(L_proteome, "./proteome_last_impute-norm_ordered_data.csv")
write.csv(L_methylome, "./methylome_last_norm_ordered_data.csv")
write.csv(L_transcriptome, "./transcriptome_last_norm_ordered_data.csv")
write.csv(L_metabolome, "./metabolome_last_norm_ordered_data.csv")

selected_cols_by_standard <- function(df, standard) {
  labels <- colnames(standard)
  
  selected_data <- df[, colnames(df) %in% labels]
  
  return(selected_data)
}

# Equalize the number of samples across all layers based on the minimum count (n)
# GA749 (Group4, in Last data, not including transcriptome data)
L_pro <- selected_cols_by_standard(L_proteome, L_transcriptome)
L_methyl <- selected_cols_by_standard(L_methylome, L_transcriptome)
L_metabo <- selected_cols_by_standard(L_metabolome, L_transcriptome)

write.csv(L_pro, "./proteome_last_impute-norm_final_data.csv")
write.csv(L_methyl, "./methylome_last_norm_final_data.csv")
write.csv(L_transcriptome, "./transcriptome_last_norm_final_data.csv")   
write.csv(L_metabo, "./metabolome_last_norm_final_data.csv")

## First data
metabo <- data.frame(fread("metabolome_first_norm_data.csv"), stringsAsFactors = FALSE)
methyl <- data.frame(fread("methylome_first_norm_data_t.csv"), stringsAsFactors = FALSE)
prot <- data.frame(fread("proteome_first_impute-norm_data.csv"), stringsAsFactors = FALSE)

#methyl <- t(methyl)
#colnames(methyl) <- methyl[1,]
#write.csv(methyl, "./methylome_first_norm_data_t.csv")

rownames(prot) <- prot$Label
rownames(methyl) <- methyl$Label
rownames(metabo) <- metabo$Label

prot$Label <- NULL
methyl$Label <- NULL
metabo$Label <- NULL

F_proteome <- reorder_data(prot, F_clinical$Label)
F_methylome <- reorder_data(methyl, F_clinical$Label)
F_metabolome <- reorder_data(metabo, F_clinical$Label)

write.csv(F_proteome, "./proteome_first_impute-norm_ordered_data.csv")
write.csv(F_methylome, "./methylome_first_norm_ordered_data.csv")
write.csv(F_metabolome, "./metabolome_first_norm_ordered_data.csv")

selected_cols_by_standard <- function(df, standard) {
  labels <- colnames(standard)
  
  selected_data <- df[, colnames(df) %in% labels]
  
  return(selected_data)
}

# Equalize the number of samples across all layers based on the minimum count (n)
F_pro <- selected_cols_by_standard(F_proteome, F_methylome)
F_metabo <- selected_cols_by_standard(F_metabolome, F_methylome)
F_methyl <- F_methylome

write.csv(F_pro, "./proteome_first_impute-norm_final_data.csv")
write.csv(F_methyl, "./methylome_first_norm_final_data.csv")             
write.csv(F_metabo, "./metabolome_first_norm_final_data.csv")
