# Data Preprocessing
rm(list=ls())

# Package download 
install.packages("data.table")
install.packages("readxl")
install.packages("dplyr")
install.packages("stringr")
BiocManager::install("impute")   

# loading Library
suppressMessages({
  library(data.table)
  library(readxl)
  library(dplyr)
  library(impute)
  library(stringr)
})

getwd() 
setwd("D:/NMF_Clustering/input") # Change the directory to your data file storage
Sys.glob("*") 

## load multiomics data
# methylome data
methylome <- data.frame(fread("Methylome_data.csv"), stringsAsFactors = FALSE)
methylome

# delete X infront of column name (start numbers) by R Name Rule
colnames(methylome) <- sub("^X(\\d)", "\\1", colnames(methylome))
colnames(methylome)

methylome_info <- methylome$TargetID
methylome <- methylome[,-1]

dim(methylome)
colnames(methylome)  # Sample ID
rownames(methylome) <- methylome_info
rownames(methylome)
summary(methylome)

# proteome data
# 999 proteins -> not 7 annotated proteins -> 992 proteins
proteome <- read.csv("Proteome_annotated_data.csv")
proteome[1:5, 1:5]
dim(proteome)

proteome <- data.frame(proteome, stringsAsFactors = FALSE)
proteome[1:5,]
#summary(proteome)

# delete X infront of column name (start numbers) by R Name Rule
colnames(proteome) <- sub("^X(\\d)", "\\1", colnames(proteome))
colnames(proteome)

# Can't use varaible name that starts number by R Name Rule
proteome <- proteome %>% rename('First.Protein' = '1st.Protein')
#proteome <- proteome %>% relocate(c(Gene,First.Protein), .after=Accession)
colnames(proteome)

proteome_info <- proteome %>% select(Accession, Gene, First.Protein, Description)
proteome <- proteome %>% select(-c(Accession, Gene, First.Protein, Description))

rownames(proteome) <- proteome_info$Accession
head(proteome_info, 3)
head(proteome, 3)

summary(proteome)

transcriptome <- data.frame(fread("Transcriptome_data.csv"), stringsAsFactors = FALSE)
transcriptome[1:5,]

# delete X infront of column name (start numbers) by R Name Rule
colnames(transcriptome) <- sub("^X(\\d)", "\\1", colnames(transcriptome))

transcriptome$Gene_ID <- as.character(transcriptome$Gene_ID)  # id is reflected to string.
str(transcriptome)

transcriptome_info <- transcriptome %>% select(c(1, 81:118))
transcriptome <- transcriptome %>% select(-c(1, 81:118))

transcriptome_info$Gene_Symbol <- as.character(transcriptome_info$Gene_Symbol)

rownames(transcriptome) <- transcriptome_info$probeset_id
#rownames(transcriptome) <- transcriptome_info$Gene_ID
colnames(transcriptome_info)
rownames(transcriptome)

metabolome <- data.frame(fread("Metabolomics_data.csv"), stringsAsFactors = FALSE)
metabolome[1:5,]

colnames(metabolome)
metabolome_info <- metabolome %>% select(Features, MW_RT_mode, Name, Formula, KEGG.ID)
metabolome <- metabolome %>% select(-c(Features, MW_RT_mode, Name, Formula, KEGG.ID))

colnames(metabolome_info)
rownames(metabolome) <- metabolome_info$Features
rownames(metabolome)

sampleinfo <- read.csv("SampleInfo.csv")
dim(sampleinfo)
sampleinfo[1:5,]

sampleinfo <- data.frame(sampleinfo, stringsAsFactors = FALSE)  # Group factoring 
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

# check the negative value - none
summary(methylome)              # 0.006 ~ 0.99
summary(proteome)               # min 0.7 ~ max 59000
summary(transcriptome)          # min 1 ~ max 14
summary(metabolome)             # min 299 ~ 90000

sum(is.na(methylome))
sum(is.na(proteome))           # Proteome have many NA value -> preprocessing
sum(is.na(transcriptome))
sum(is.na(metabolome))

### integrate column name to total sample label
methyl_name <- list()
prot_name <- list()
trans_name <- list()
metabo_name <- list()

for (i in colnames(methylome)) {
  for (j in 1:nrow(sample_label)) {
    if (i %in% sample_label[j,]) {
      k <- sample_label[j,1]
      methyl_name <- append(k, methyl_name)
    }
  }
}

for (i in colnames(proteome)) {
  for (j in 1:nrow(sample_label)) {
    if (i %in% sample_label[j,]) {
      k <- sample_label[j,1]
      prot_name <- append(k, prot_name)
    }
  }
}

for (i in colnames(transcriptome)) {
  for (j in 1:nrow(sample_label)) {
    if (i %in% sample_label[j,]) {
      k <- sample_label[j,1]
      trans_name <- append(k, trans_name)
    }
  }
}

for (i in colnames(metabolome)) {
  for (j in 1:nrow(sample_label)) {
    if (i %in% sample_label[j,]) {
      k <- sample_label[j,1]
      metabo_name <- append(k, metabo_name)
    }
  }
}

methyl_name <- rev(methyl_name)
prot_name <- rev(prot_name)
trans_name <- rev(trans_name)
metabo_name <- rev(metabo_name)

colnames(methylome)
colnames(proteome)
colnames(transcriptome)
colnames(metabolome)

methyl_name
prot_name
trans_name
metabo_name

# Store original omics IDs separately, then create an integration label to enable analysis.
methyl_id <- colnames(methylome)
prot_id <- colnames(proteome)
trans_id <- colnames(transcriptome)
metabo_id <- colnames(metabolome)

colnames(methylome) <- methyl_name
colnames(proteome) <- prot_name
colnames(transcriptome) <- trans_name
colnames(metabolome) <- metabo_name

## align data to sample label 
reorder_data <- function(data, sample_labels) {
  valid_labels <- sample_labels[sample_labels %in% colnames(data)]
  ordered_data <- data[, match(valid_labels, colnames(data))]
  return(ordered_data)
}

proteome <- reorder_data(proteome, sample_label$Label)
methylome <- reorder_data(methylome, sample_label$Label)
transcriptome <- reorder_data(transcriptome, sample_label$Label)
metabolome <- reorder_data(metabolome, sample_label$Label)

### input data for RapidMiner
#t_proteome <- t(proteome)
#t_methylome <- t(methylome)
#t_transcriptome <- t(transcriptome)
#t_metabolome <- t(metabolome)

#write.csv(t_methylome, "../methylome_data.csv")
#write.csv(t_proteome, "../proteome_data.csv")
#write.csv(t_transcriptome, "../transriptome_data.csv")
#write.csv(t_metabolome, "../metabolome_data.csv")

sampleinfo$Group
sampleinfo$Time

## Data Preparation for analysis
# Case 1. Group5 vs Group1-4 in First (Control vs HD-exposed)
# Case 2. Group4 vs Group1-3 in First (Normal vs impaired lung function)
# Case 3. Group5 vs Group1-4 in Last (Control vs HD-exposed)
# Case 4. Group4 vs Group1-3 in Last (Normal vs impaired lung function)

F_info <- sampleinfo %>% filter(Time == 'F') %>% 
  select(Label, Group, Time, Age, Sex, FEV1_z_F, FVC_z_F, FEV1FVC_z_F, FEF_z_F)

L_info <- sampleinfo %>% filter(Time == 'L') %>% 
  select(Label, Group, Time, Age, Sex, FEV1_z_L, FVC_z_L, FEV1FVC_z_L, FEF_z_L)

F_cont_exp <- F_info %>% filter(Group == 'Group5') 
F_exposed <- F_info %>% filter(Group != 'Group5')    # not included control

F_cont_sym <- F_info %>% filter(Group == 'Group4')
F_symptom <- F_exposed %>% filter(Group != 'Group4')  # not included normal

L_cont_exp <- L_info %>% filter(Group == 'Group5') 
L_exposed <- L_info %>% filter(Group != 'Group5')    # not included control

L_cont_sym <- L_info %>% filter(Group == 'Group4')
L_symptom <- L_exposed %>% filter(Group != 'Group4')  # not included normal

tmp1 <- methylome %>% log2()
summary(tmp1)
summary(methylome)
### log2 normalization - fitting distribution across multiple datasets
## Make log2 matrix for NMF
missing_cut = 0.6
neighbors = 7     # Default: k=10
option_subset = 'All'

# Proteome (have NA values)
mean(is.na(proteome))   # NA ratio: 21.8% (high, but it's common in real world data)
sum(is.na(proteome))    # 3

ncol(proteome)   # sample count
rowSums(is.na(proteome))  # NA distribution by proteins
ncol(proteome)*0.6        # NA >= 84 (60%)

tmp1 <- proteome[rowSums(is.na(proteome)) < (ncol(proteome)*missing_cut),] 
tmp1 <- impute.knn(tmp1 %>% as.matrix(), k=neighbors, rowmax=0, rng.seed=1218220) 
tmp1 <- tmp1[["data"]]
sum(is.na(tmp1))     
# after filtering, The sum of NA values is 0

tmp1 <- tmp1 %>% log2 %>% data.frame()   
head(tmp1)
summary(tmp1)  # remaining negative value

### negative data processing
tmp2 <- tmp1
row.names(tmp1) <- paste0(rownames(tmp1), '1', sep='_')
row.names(tmp2) <- paste0(rownames(tmp2), '2', sep='_')

# Step 1. change all negative value to zero
tmp1[tmp1 <0] <- 0

# Step 2. change all positive value to zero (tmp1, tmp2 - parallel)
tmp2[tmp2 >0] <- 0 
tmp2 <- tmp2 * -1

# Step 3. Concatenation
res <- rbind(tmp1, tmp2)
dim(res)

# Step 4. Remove Null value
res = res[rowSums(res)!=0, ]
dim(res)
res[1:5, 1:5]

rowSums(res)

# Step 5. Remove zero value (using index)
res = res[-c(823:849), ]

View(res)
summary(res)

## confirmed data -> all negative values are disappeared
new_names <- sub('1_$', "", rownames(res))
rownames(res) <- new_names

log2_pro <- res

## methylome data
summary(methylome)  # raw data

tmp1 <- methylome %>% log2()
tmp2 <- tmp1
row.names(tmp1) <- paste0(rownames(tmp1), '1', sep='_')
row.names(tmp2) <- paste0(rownames(tmp2), '2', sep='_')

# Step 1
tmp1[tmp1 <0] <- 0

# Step 2
tmp2[tmp2 >0] <- 0 
tmp2 <- tmp2 * -1

# Step 3
res <- rbind(tmp1, tmp2)
dim(res)

# Step 4
res = res[rowSums(res)!=0, ]
dim(res)
res[1:5, 1:5]
View(res)

summary(res)
rownames(res)
rownames(methylome)
rownames(res)

new_names <- sub('2_$', "", rownames(res))
rownames(res) <- new_names

log2_methyl <- res

## Transcriptome data
summary(transcriptome)

tmp1 <- transcriptome %>% log2()
tmp2 <- tmp1
row.names(tmp1) <- paste0(rownames(tmp1), '1', sep='_')
row.names(tmp2) <- paste0(rownames(tmp2), '2', sep='_')

# Step 1
tmp1[tmp1 <0] <- 0

# Step 2
tmp2[tmp2 >0] <- 0 
tmp2 <- tmp2 * -1

# Step 3
res <- rbind(tmp1, tmp2)
dim(res)

# Step 4
res = res[rowSums(res)!=0, ]
dim(res)
res[1:5, 1:5]
View(res)

# Step 5
res[31160,]
res = res[-c(31136:31160), ]

View(res)
summary(res)

rownames(res)

new_names <- sub('1_$', "", rownames(res))
rownames(res) <- new_names

log2_trans <- res

## metabolome -> not negative value, so skip the negative value processing step.
tmp1 <- metabolome %>% log2()

rownames(tmp1)
log2_metabo <- tmp1

# data -> log2_metabo, log2_pro, log2_trans, log2_methyl
select_cols_by_label <- function(data, ref) {
  labels <- ref$Label
  selected_data <- data[, colnames(data) %in% labels]
  return(selected_data)
}

save.image("D:/NMF_Clustering/rdata/preprocessing_final.RData")
load("D:/NMF_Clustering/rdata/preprocessing_final.RData")

# Data preparation
# F_cont_exp vs F_exposed
# L_cont_exp vs L_exposed 

F_pro <- select_cols_by_label(log2_pro, F_info)
F_trans <- select_cols_by_label(log2_trans, F_info)     # only Transcriptome in Last 
F_metabo <- select_cols_by_label(log2_metabo, F_info)
F_methyl <- select_cols_by_label(log2_methyl, F_info)

L_pro <- select_cols_by_label(log2_pro, L_info)
L_trans <- select_cols_by_label(log2_trans, L_info)     
L_metabo <- select_cols_by_label(log2_metabo, L_info)
L_methyl <- select_cols_by_label(log2_methyl, L_info)

####### For ML analysis, Rapid Miner ######
F_pro <- select_cols_by_label(impute_pro, F_info)
F_trans <- select_cols_by_label(transcriptome, F_info)    
F_metabo <- select_cols_by_label(metabolome, F_info)
F_methyl <- select_cols_by_label(methylome, F_info)

L_pro <- select_cols_by_label(impute_pro, L_info)
L_trans <- select_cols_by_label(transcriptome, L_info)     
L_metabo <- select_cols_by_label(metabolome, L_info)
L_methyl <- select_cols_by_label(methylome, L_info)

write.csv(F_pro, "../proteome_first_impute-norm_data.csv")
write.csv(F_trans, "../transcriptome_first_norm_data.csv")
write.csv(F_metabo, "../metabolome_first_norm_data.csv")
write.csv(F_methyl, "../methylome_first_norm_data.csv")

write.csv(L_pro, "../proteome_last_impute-norm_data.csv")
write.csv(L_trans, "../transcriptome_last_norm_data.csv")
write.csv(L_metabo, "../metabolome_last_norm_data.csv")
write.csv(L_methyl, "../methylome_last_norm_data.csv")
##########

selected_cols_by_standard <- function(df, standard) {
  labels <- colnames(standard)
  selected_data <- df[, colnames(df) %in% labels]
  return(selected_data)
}

## First
F_pro <- selected_cols_by_standard(F_pro, F_trans)
F_methyl <- selected_cols_by_standard(F_methyl, F_trans)
F_metabo <- selected_cols_by_standard(F_metabo, F_trans)

# Final stage
F_concat_exp <- rbind(F_trans, F_pro)
F_concat_exp <- rbind(F_concat_exp, F_methyl)
F_concat_exp <- rbind(F_concat_exp, F_metabo)

print( paste('After log2 transform: Total Data - the number of rows:', nrow(F_concat_exp)))

setwd("D:/NMF_Clustering/")

write.csv(F_concat_exp, "./output/multiomics_F-exposed_log2_norm.csv")

## Last
# GA749 (Group2) don't include transcriptome data (exceptional 1 case)
L_pro <- selected_cols_by_standard(L_pro, L_trans)
L_methyl <- selected_cols_by_standard(L_methyl, L_trans)
L_metabo <- selected_cols_by_standard(L_metabo, L_trans)

# Final stage
L_concat_exp <- rbind(L_trans, L_pro)
L_concat_exp <- rbind(L_concat_exp, L_methyl)
L_concat_exp <- rbind(L_concat_exp, L_metabo)

print( paste('After log2 transform: Total Data - the number of rows:', nrow(L_concat_exp)))

## Save to file
write.csv(L_concat_exp, "./output/multiomics_L-exposed_tri_log2_norm.csv")

# data and clinical information 
write.csv(F_info, "D:/NMF_Clustering/output/First_information.csv")
write.csv(L_info, "D:/NMF_Clustering/output/Last_information.csv")

dim(L_concat_exp)  

L_concat_exp[1:5, 1:5] 
# RNA(1~31135 rows), Protein(31136~31958)
# methylome(31958, 732367), metabolome(732368, 733119)
View(L_concat_exp)
