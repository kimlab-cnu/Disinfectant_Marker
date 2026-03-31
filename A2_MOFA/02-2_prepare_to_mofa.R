rm(list = ls())

setwd("D:/Analysis/MOFA/input/")
getwd() 

### library
library(BiocManager)
library(R.utils)
library(data.table)
library(devtools)
library(ggplot2)
library(ggpubr)
library(reticulate)
library(psych)
library(tidyverse)
library(survival)
library(survminer)
library(MOFA2)
library(dplyr)

#### first_data & last_data 
multiome <- data.frame(fread("multiome_first_comp1_data.csv"), stringsAsFactors = FALSE)

multiome$V1 <- NULL

## Update group labels
multiome$group <- as.character(multiome$group)
multiome$group <- factor(multiome$group,
                         levels = c("group_1", "group_2", "group_3"),
                         labels = c("Normal", "Target", "Severe"))

multiome$sample <- sub("(_group_\\d+)$", "", multiome$sample)

# multi-group framework => the number of groups must be specified 
#groups = c(rep("group_1", 20), rep("group_2", 40), rep("group_3",10))

# update group labels..
groups = c(rep("Normal", 20), rep("Target", 40), rep("Severe",10))

head(multiome)
MOFAobject = create_mofa(multiome, groups=groups)   # ,groups=groups  # 15:18-

print(MOFAobject)                 # Number of View : number of omics layer, Views name : omics type
plot_data_overview(MOFAobject)    # Feature number: by omics layer, group number & name & sample by group 

pdf("../output/mofa_first_comp1_distribution.pdf")
plot_data_overview(MOFAobject)
dev.off()

## Adjust based on user preference
# Define options-1 (data)
data_opts = get_default_data_options(MOFAobject)
head(data_opts) 
#data_opts$scale_groups = TRUE

# Define options-1 (model)
model_opts = get_default_model_options(MOFAobject)
head(model_opts)  
model_opts$num_factors = 10

# Define options-1 (train)
train_opts = get_default_training_options(MOFAobject)
head(train_opts)
# convergence_mode => slow or 
# test =  train_opts$convergence_mode
# train_opts$convergence_mode <- "slow"         # fast -> yields satisfactory results

# Build and train the MOFA object
MOFAobject = prepare_mofa(object = MOFAobject, data_options = data_opts, 
                          model_options = model_opts, training_options = train_opts)

outfile = file.path(getwd(),"multiome_first_comp1.hdf5")

# Model training
MOFAobject.trained = run_mofa(MOFAobject, outfile, use_basilisk = TRUE)

test = MOFAobject.trained

test@expectations$Z$Normal[1:5,]

test@expectations$Z$group_1[1:5,]
MOFAobject.trained@expectations$Z$group_1[1:5,]

# Check in group_2
test@expectations$Z$group_2[1:5,]
MOFAobject.trained@expectations$Z$group_2[1:5,]

# Check in group_3
test@expectations$Z$group_3[1:5,]
MOFAobject.trained@expectations$Z$group_3[1:5,]

#### last_data ####

multiome <- data.frame(fread("multiome_last_comp1_data.csv"), stringsAsFactors = FALSE)
multiome$V1 <- NULL

## Update group labels
multiome$group <- as.character(multiome$group)
multiome$group <- factor(multiome$group,
                         levels = c("group_1", "group_2", "group_3"),
                         labels = c("Normal", "Target", "Severe"))

multiome$sample <- sub("(_group_\\d+)$", "", multiome$sample)

# multi-group framework => the number of groups must be specified 
#groups = c(rep("group_1", 20), rep("group_2", 40), rep("group_3",10))

# update group labels..
groups = c(rep("Normal", 20), rep("Target", 40), rep("Severe",10))

head(multiome)
MOFAobject = create_mofa(multiome, groups=groups)   # ,groups=groups  # 15:18-

print(MOFAobject)                 # Number of View : number of omics layer, Views name : omics type
plot_data_overview(MOFAobject)    # Feature number: by omics layer, group number & name & sample by group 

pdf("../output/mofa_last_comp1_distribution.pdf")
plot_data_overview(MOFAobject)
dev.off()

## Adjust based on user preference
# Define options-1 (data)
data_opts = get_default_data_options(MOFAobject)
head(data_opts) 
#data_opts$scale_groups = TRUE

# Define options-1 (model)
model_opts = get_default_model_options(MOFAobject)
head(model_opts)  
model_opts$num_factors = 15

# Define options-1 (train)
train_opts = get_default_training_options(MOFAobject)
head(train_opts)
# convergence_mode => slow or 
# test =  train_opts$convergence_mode
# train_opts$convergence_mode <- "slow"         # fast -> yields satisfactory results

# Build and train the MOFA object
MOFAobject = prepare_mofa(object = MOFAobject, data_options = data_opts, 
                          model_options = model_opts, training_options = train_opts)

outfile = file.path(getwd(),"multiome_last_comp1.hdf5")

# Model training
MOFAobject.trained = run_mofa(MOFAobject, outfile, use_basilisk = TRUE)

test = MOFAobject.trained

test@expectations$Z$group_1[1:5,]
MOFAobject.trained@expectations$Z$group_1[1:5,]