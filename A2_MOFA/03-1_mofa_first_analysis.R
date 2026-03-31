rm(list = ls()) 
getwd() 
setwd("D:/Analysis/MOFA/input/")

library(ggplot2)
library(data.table)
library(MOFA2)

#### first_data (15 factors -> 9 factors; other factors explain no variance, so they're removed)
### comp1 (normal vs target vs critical, three groups)
model = load_model("multiome_first_comp1.hdf5") 

# Add metadata to loaded MOFA object
# samples_metadata -> only group and sample information
# Link group and sample information using metadata

setwd("../output/")
metadata = data.frame(fread("clinical_first_exposed_comp1_data.csv"), stringsAsFactors = FALSE)
metadata$V1 = NULL

samples_metadata(model) <- metadata    
# sample, condition, age, sex + PFT variables / change sample ID
# multi-group mode: need to group column
model@samples_metadata$sample 

## Model overview
# Note the contents of trained object, focusing on 'data' and 'expectations'
slotNames(model)           
plot_data_overview(model)  
names(model@data)        

# save figure
p1 <- plot_data_overview(model)
ggsave("MOFA_high-res_first_distribution_comp1.png", plot=p1, 
       width=7.6, height=7.3, dpi=300)

model@data$metabolome[[1]][1:5,1:5]
model@data$methylome[[1]][1:5,1:5]
model@data$proteome[[1]][1:5,1:5]

model@samples_metadata[1:5, 1:9]

# expectations
dim(model@expectations$Z$Normal)    # Check the values for the 20 samples x 9 factors
model@expectations$Z$Normal[1:5,]

dim(model@expectations$Z$Target)    # Check the values for the 40 samples x 9 factors
model@expectations$Z$Target[1:5,]

dim(model@expectations$Z$Severe)    # Check the values for the 10 samples x 9 factors
model@expectations$Z$Severe[1:5,]

plot_factor_cor(model)              # Identify the correlation between factors.

pdf("../output/mofa_correlation_factor_first_comp1_data.pdf")
plot_factor_cor(model)
dev.off()

dim(model@expectations$W$metabolome) 
model@expectations$W$metabolome[1:5,]
model@expectations$W$methylome[1:5,]
model@expectations$W$proteome[1:5,]

## variance explained (R square) 
# How well the inferred factors explain the variance in the input views (omics) and groups
# Total variance explained per view and group
head(model@cache$variance_explained$r2_total[[1]]) 
head(model@cache$variance_explained$r2_per_factor[[1]]) 

model@cache$variance_explained$r2_per_factor

plot_variance_explained(model, x="view", y="factor")  # variance explained by each factor
plot_variance_explained(model, x="group", y="factor", plot_total = T)[[2]]  # total variance explained 

pdf("../output/mofa_variance_explain_factor_first_comp1.pdf")
plot_variance_explained(model, x="view", y="factor")  
dev.off()

pdf("../output/mofa_variance_explain_total_first_comp1.pdf")
plot_variance_explained(model, x="group", y="factor", plot_total = T)[[2]]  
dev.off()

# save figure
p2 <- plot_variance_explained(model, x="view", y="factor") 
ggsave("MOFA_high-res_variance_explain_factor_first_comp1.png", plot=p2, 
       width=9, height=5.8, dpi=300)

help("correlate_factors_with_covariates")
correlate_factors_with_covariates(model,
                                  covariates = c("sex","condition", "age"), 
                                  plot="log_pval")

# factor 1, 8 associated with 'age'

pdf("../output/correlation_factors_with_PFT_by-group_first_comp1.pdf")

correlate_factors_with_covariates(model,
                                  covariates = c("FEV1FVC_z_F", "FEV1_z_F", "FVC_z_F", "FEF_z_F"),
                                  groups="Normal",
                                  plot="log_pval")

correlate_factors_with_covariates(model,
                                  covariates = c("FEV1FVC_z_F", "FEV1_z_F", "FVC_z_F", "FEF_z_F"),
                                  groups="Target",
                                  plot="log_pval")

correlate_factors_with_covariates(model,
                                  covariates = c("FEV1FVC_z_F", "FEV1_z_F", "FVC_z_F", "FEF_z_F"),
                                  groups="Severe",
                                  plot="log_pval")

dev.off()

save.image("../rdata/mofa_first_comp1.RData")
load("../rdata/mofa_first_comp1.RData")

p <- plot_factor(model, 
                 factors = c(1, 3, 6, 7, 8),
                 color_by = "condition",
                 dot_size = 3,       
                 dodge = T,        
                 legend = T,          
                 add_violin = T,     
                 violin_alpha = 0.25  
)

# The output of plot_factor is a ggplot2 object that we can edit
p <- p + 
  scale_color_manual(values=c("normal"="black", "obstructive"="blue", "restrictive"="purple", "critical"="red")) +
  scale_fill_manual(values=c("normal"="black", "obstructive"="blue", "restrictive"="purple", "critical"="red"))

print(p)

#### Change feature name (if you want)
backup <- features_names(model)[["proteome"]]
new_name <- sapply(strsplit(features_names(model)[["proteome"]], "_"), `[`, 1)
features_names(model)[["proteome"]] <- new_name

backup2 <- features_names(model)[["metabolome"]]
new_name <- sapply(strsplit(features_names(model)[["metabolome"]], "_"), `[`, 1)
features_names(model)[["metabolome"]] <- new_name

backup3 <- features_names(model)[["methylome"]]
new_name <- sapply(strsplit(features_names(model)[["methylome"]], "_"), `[`, 1)
features_names(model)[["methylome"]] <- new_name

#### plotting by factor
## factor 1
plot_weights(model,view = "metabolome", factor = 1, nfeatures = 10, scale = T, text_size=3.0)
plot_top_weights(model, view = "metabolome", factor = 1, nfeatures = 10, scale = T) 

plot_weights(model,view = "methylome", factor = 1, nfeatures = 10, scale = T, text_size=3.0)
plot_top_weights(model, view = "methylome", factor = 1, nfeatures = 10, scale = T) 

plot_weights(model,view = "proteome", factor = 1, nfeatures = 10, scale = T, text_size=2.5)
plot_top_weights(model, view = "proteome", factor = 1, nfeatures = 10, scale = T) 

pdf("../output/mofa_factor1_weight-metabo_name-clear_first_comp1.pdf")
plot_weights(model,view = "metabolome", factor = 1, nfeatures = 10, scale = T, text_size=3.0)
plot_top_weights(model, view = "metabolome", factor = 1, nfeatures = 10, scale = T) 
dev.off()

pdf("../output/mofa_factor1_weight-methyl_first_comp1.pdf")
plot_weights(model,view = "methylome", factor = 1, nfeatures = 10, scale = T, text_size=3.0)
plot_top_weights(model, view = "methylome", factor = 1, nfeatures = 10, scale = T) 
dev.off()

pdf("../output/mofa_factor1_weight-prot_first_comp1.pdf")
plot_weights(model,view = "proteome", factor = 1, nfeatures = 10, scale = T, text_size=2.5)
plot_top_weights(model, view = "proteome", factor = 1, nfeatures = 10, scale = T) 
dev.off()

# heatmap - association with specific factors by each factor
pdf("../output/mofa_Factor1_metabolome_corr-heatmap_first_comp1.pdf")
plot_data_heatmap(model, view = "metabolome", factor = 1, features = 20, denoise = TRUE, cluster_rows = TRUE, 
                  cluster_cols = FALSE, show_rownames = TRUE, show_colnames = FALSE, scale = "row")
dev.off()

help(plot_data_scatter)
# scatter plot
pdf("../output/mofa_Factor1_group1-metabo_scatter_first_comp1.pdf")
plot_data_scatter(model, view = "metabolome", factor = 1, features = 6, groups = "Normal", sign = "positive")
plot_data_scatter(model, view = "metabolome", factor = 1, features = 6, groups = "Normal", sign = "negative")
dev.off()

pdf("../output/mofa_Factor1_group2-metabo_scatter_first_comp1.pdf")
plot_data_scatter(model, view = "metabolome", factor = 1, features = 6, groups = "Target", sign = "positive")
plot_data_scatter(model, view = "metabolome", factor = 1, features = 6, groups = "Target", sign = "negative")
dev.off()

pdf("../output/mofa_Factor1_group3-metabo_scatter_first_comp1.pdf")
plot_data_scatter(model, view = "metabolome", factor = 1, features = 4, groups = "Severe", sign = "positive")
plot_data_scatter(model, view = "metabolome", factor = 1, features = 4, groups = "Severe", sign = "negative")
dev.off()

pdf("../output/first_data/comp1/mofa_Factor1_top9-weight_metabo-by-group_scatter_first_comp1.pdf")
plot_data_scatter(model,
                  view = "metabolome",
                  factor = 1,      
                  features = 9,         
                  groups = "Normal",
                  add_lm = TRUE          # add linear regression
)

plot_data_scatter(model,
                  view = "metabolome",        
                  factor = 1,           
                  features = 9,           
                  groups = "Target",
                  add_lm = TRUE,          
                  color_by = "condition"
)

plot_data_scatter(model,
                  view = "metabolome",      
                  factor = 1,            
                  features = 9,          
                  groups = "Severe",
                  add_lm = TRUE         
)

dev.off()

## factor 3
plot_weights(model,view = "metabolome", factor = 3, nfeatures = 10, scale = T, text_size=3.0)
plot_top_weights(model, view = "metabolome", factor = 3, nfeatures = 10, scale = T) 

plot_weights(model,view = "methylome", factor = 3, nfeatures = 10, scale = T, text_size=3.0)
plot_top_weights(model, view = "methylome", factor = 3, nfeatures = 10, scale = T) 

plot_weights(model,view = "proteome", factor = 3, nfeatures = 10, scale = T, text_size=2.5)
plot_top_weights(model, view = "proteome", factor = 3, nfeatures = 10, scale = T) 

pdf("../output/mofa_factor3_weight-metabo_name-clear_first_comp1.pdf")
plot_weights(model,view = "metabolome", factor = 3, nfeatures = 10, scale = T, text_size=3.0)
plot_top_weights(model, view = "metabolome", factor = 3, nfeatures = 10, scale = T) 
dev.off()

pdf("../output/mofa_factor3_weight-methyl_first_comp1.pdf")
plot_weights(model,view = "methylome", factor = 3, nfeatures = 10, scale = T, text_size=3.0)
plot_top_weights(model, view = "methylome", factor = 3, nfeatures = 10, scale = T) 
dev.off()

pdf("../output/mofa_factor3_weight-prot_first_comp1.pdf")
plot_weights(model,view = "proteome", factor = 3, nfeatures = 10, scale = T, text_size=2.5)
plot_top_weights(model, view = "proteome", factor = 3, nfeatures = 10, scale = T) 
dev.off()

pdf("../output/mofa_factor3_metabolome_corr-heatmap_first_comp1.pdf")
plot_data_heatmap(model, view = "metabolome", factor = 3, features = 20, denoise = TRUE, cluster_rows = TRUE, 
                  cluster_cols = FALSE, show_rownames = TRUE, show_colnames = FALSE, scale = "row")
dev.off()

#help(plot_data_scatter)
# scatter plot
pdf("../output/mofa_factor3_group1-metabo_scatter_first_comp1.pdf")
plot_data_scatter(model, view = "metabolome", factor = 3, features = 4, groups = "Normal", sign = "positive")
plot_data_scatter(model, view = "metabolome", factor = 3, features = 4, groups = "Normal", sign = "negative")
dev.off()

pdf("../output/mofa_factor3_group2-metabo_scatter_first_comp1.pdf")
plot_data_scatter(model, view = "metabolome", factor = 3, features = c("met589_metabolome", "met239_metabolome", "met974_metabolome", "met1036_metabolome", "met1020_metabolome"), 
                  groups = "Target", sign = "positive")
plot_data_scatter(model, view = "metabolome", factor = 3, features = 4, groups = "Target", sign = "negative")
dev.off()

pdf("../output/mofa_factor3_group3-metabo_scatter_first_comp1.pdf")
plot_data_scatter(model, view = "metabolome", factor = 3, features = 4, groups = "Severe", sign = "positive")
plot_data_scatter(model, view = "metabolome", factor = 3, features = 4, groups = "Severe", sign = "negative")
dev.off()  

pdf("../output/first_data/comp1/mofa_Factor3_top9-weight_metabo-by-group_scatter_first_comp1.pdf")
plot_data_scatter(model,
                  view = "metabolome",         
                  factor = 3,          
                  features = 9,          
                  groups = "Normal",
                  add_lm = TRUE    
)

plot_data_scatter(model,
                  view = "metabolome",   
                  factor = 3,           
                  features = 9,          
                  groups = "Target",
                  add_lm = TRUE,         
                  color_by = "condition"
)

plot_data_scatter(model,
                  view = "metabolome",    
                  factor = 3,            
                  features = 9,           
                  groups = "Severe",
                  add_lm = TRUE         
)

dev.off()

## factor 5
plot_weights(model,view = "metabolome", factor = 5, nfeatures = 10, scale = T, text_size=3.0)
plot_top_weights(model, view = "metabolome", factor = 5, nfeatures = 10, scale = T) 

plot_weights(model,view = "methylome", factor = 5, nfeatures = 10, scale = T, text_size=3.0)
plot_top_weights(model, view = "methylome", factor = 5, nfeatures = 10, scale = T) 

plot_weights(model,view = "proteome", factor = 5, nfeatures = 10, scale = T, text_size=2.5)
plot_top_weights(model, view = "proteome", factor = 5, nfeatures = 10, scale = T) 

pdf("../output/mofa_factor5_weight-metabo_first_comp1.pdf")
plot_weights(model,view = "metabolome", factor = 5, nfeatures = 10, scale = T, text_size=3.0)
plot_top_weights(model, view = "metabolome", factor = 5, nfeatures = 10, scale = T) 
dev.off()

pdf("../output/mofa_factor5_weight-methyl_first_comp1.pdf")
plot_weights(model,view = "methylome", factor = 5, nfeatures = 10, scale = T, text_size=3.0)
plot_top_weights(model, view = "methylome", factor = 5, nfeatures = 10, scale = T) 
dev.off()

pdf("../output/mofa_factor5_weight-prot_first_comp1.pdf")
plot_weights(model,view = "proteome", factor = 5, nfeatures = 10, scale = T, text_size=2.5)
plot_top_weights(model, view = "proteome", factor = 5, nfeatures = 10, scale = T) 
dev.off()

pdf("../output/mofa_factor5_methylome_corr-heatmap_first_comp1.pdf")
plot_data_heatmap(model, view = "methylome", factor = 5, features = 20, denoise = TRUE, cluster_rows = TRUE, 
                  cluster_cols = FALSE, show_rownames = TRUE, show_colnames = FALSE, scale = "row")
dev.off()

#help(plot_data_scatter)
# scatter plot
pdf("../output/mofa_factor5_group1-methyl_scatter_first_comp1.pdf")
plot_data_scatter(model, view = "methylome", factor = 5, features = 9, groups = "Normal", sign = "positive")
plot_data_scatter(model, view = "methylome", factor = 5, features = c("cg20483857_methylome", "cg10387956_methylome"),
                  groups = "Normal", sign = "negative")
dev.off()

pdf("../output/mofa_factor5_group2-methyl_scatter_first_comp1.pdf")
plot_data_scatter(model, view = "methylome", factor = 5, features = 9, groups = "Target", sign = "positive")
plot_data_scatter(model, view = "methylome", factor = 5, features = c("cg25052156_methylome", "cg20483857_methylome",
                                                                      "cg08857906_methylome"), groups = "Target", sign = "negative")
dev.off()

pdf("../output/mofa_factor5_group3-methyl_scatter_first_comp1.pdf")
plot_data_scatter(model, view = "methylome", factor = 5, features = 9, groups = "Severe", sign = "positive")
plot_data_scatter(model, view = "methylome", factor = 5, features = c("cg02002523_methylome",
                                                                      "cg20483857_methylome"), groups = "Severe", sign = "negative")
dev.off()

save.image("../rdata/mofa_first_comp1_half.RData")
load("../rdata/mofa_first_comp1_half.RData")

## factor 6
plot_weights(model,view = "metabolome", factor = 6, nfeatures = 10, scale = T, text_size=3.0)
plot_top_weights(model, view = "metabolome", factor = 6, nfeatures = 10, scale = T) 

plot_weights(model,view = "methylome", factor = 6, nfeatures = 10, scale = T, text_size=3.0)
plot_top_weights(model, view = "methylome", factor = 6, nfeatures = 10, scale = T) 

plot_weights(model,view = "proteome", factor = 6, nfeatures = 10, scale = T, text_size=2.5)
plot_top_weights(model, view = "proteome", factor = 6, nfeatures = 10, scale = T) 

pdf("../output/mofa_factor6_weight-metabo_first_comp1.pdf")
plot_weights(model,view = "metabolome", factor = 6, nfeatures = 10, scale = T, text_size=3.0)
plot_top_weights(model, view = "metabolome", factor = 6, nfeatures = 10, scale = T) 
dev.off()

pdf("../output/mofa_factor6_weight-methyl_first_comp1.pdf")
plot_weights(model,view = "methylome", factor = 6, nfeatures = 10, scale = T, text_size=3.0)
plot_top_weights(model, view = "methylome", factor = 6, nfeatures = 10, scale = T) 
dev.off()

pdf("../output/first_data/comp1/mofa_factor6_weight-prot_name-clear_first_comp1.pdf")
plot_weights(model,view = "proteome", factor = 6, nfeatures = 10, scale = T, text_size=2.5)
plot_top_weights(model, view = "proteome", factor = 6, nfeatures = 10, scale = T) 
dev.off()

pdf("../output/mofa_factor6_proteome_corr-heatmap_first_comp1.pdf")
plot_data_heatmap(model, view = "proteome", factor = 6, features = 20, groups = "Severe", denoise = TRUE, cluster_rows = TRUE, 
                  cluster_cols = FALSE, show_rownames = TRUE, show_colnames = FALSE, scale = "row")
dev.off()

#help(plot_data_scatter)
# scatter plot
pdf("../output/mofa_factor6_group1-prot_scatter_first_comp1.pdf")
plot_data_scatter(model, view = "proteome", factor = 6, features = c("P07339_proteome","A0A075B6J1_proteome", "P40967;P40967-2;P40967-4;P40967-5_proteome"
                                                                     ), groups = "Normal", sign = "positive")
plot_data_scatter(model, view = "proteome", factor = 6, features = 6, groups = "Normal", sign = "negative")
dev.off()

pdf("../output/mofa_factor6_group2-prot_scatter_first_comp1.pdf")
plot_data_scatter(model, view = "proteome", factor = 6, features = c("Q86XI6_proteome", "P0DJI8_proteome", "Q8TAW3_proteome", "Q96KN2_proteome", "P40967;P40967-2;P40967-4;P40967-5_proteome"
                                                                     ), groups = "Target", sign = "positive")
plot_data_scatter(model, view = "proteome", factor = 6, features = 6, groups = "Target", sign = "negative")
dev.off()

pdf("../output/first_data/comp1/mofa_factor6_group3-prot_scatter_first_comp1.pdf")
plot_data_scatter(model, view = "proteome", factor = 6, features = 9, groups = "Severe", sign = "positive")
plot_data_scatter(model, view = "proteome", factor = 6, features = c("A0A075B6P5;P01615_proteome", "A0A075B6S9;P0DSN7_proteome",
                                                                     "A0A0A0MT89_proteome", "A0A0C4DH29_proteome"), groups = "Severe", sign = "negative")
dev.off()

## factor 7
plot_weights(model,view = "metabolome", factor = 7, nfeatures = 10, scale = T, text_size=3.0)
plot_top_weights(model, view = "metabolome", factor = 7, nfeatures = 10, scale = T) 

plot_weights(model,view = "methylome", factor = 7, nfeatures = 10, scale = T, text_size=3.0)
plot_top_weights(model, view = "methylome", factor = 7, nfeatures = 10, scale = T) 

plot_weights(model,view = "proteome", factor = 7, nfeatures = 10, scale = T, text_size=2.5)
plot_top_weights(model, view = "proteome", factor = 7, nfeatures = 10, scale = T) 

pdf("../output/mofa_factor7_weight-metabo_first_comp1.pdf")
plot_weights(model,view = "metabolome", factor = 7, nfeatures = 10, scale = T, text_size=3.0)
plot_top_weights(model, view = "metabolome", factor = 7, nfeatures = 10, scale = T) 
dev.off()

pdf("../output/mofa_factor7_weight-methyl_first_comp1.pdf")
plot_weights(model,view = "methylome", factor = 7, nfeatures = 10, scale = T, text_size=3.0)
plot_top_weights(model, view = "methylome", factor = 7, nfeatures = 10, scale = T) 
dev.off()

pdf("../output/first_data/comp1/mofa_factor7_weight-prot_name-clear_first_comp1.pdf")
plot_weights(model,view = "proteome", factor = 7, nfeatures = 10, scale = T, text_size=2.5)
plot_top_weights(model, view = "proteome", factor = 7, nfeatures = 10, scale = T) 
dev.off()

pdf("../output/mofa_factor7_proteome_corr-heatmap_first_comp1.pdf")
plot_data_heatmap(model, view = "proteome", factor = 7, features = 20, groups="Severe", denoise = TRUE, cluster_rows = TRUE, 
                  cluster_cols = FALSE, show_rownames = TRUE, show_colnames = FALSE, scale = "row")
dev.off()

#help(plot_data_scatter)
# scatter plot
pdf("../output/mofa_factor7_group1-prot_scatter_first_comp1.pdf")
plot_data_scatter(model, view = "proteome", factor = 7, features = c("Q96Q89;Q96Q89-2;Q96Q89-3;Q96Q89-4_proteome", "P12109_proteome", "Q16515-2_proteome",
                                                                     "P17483_proteome"), groups = "Normal", sign = "positive")
plot_data_scatter(model, view = "proteome", factor = 7, features = 6, groups = "Normal", sign = "negative")
dev.off()

pdf("../output/mofa_factor7_group2-prot_scatter_first_comp1.pdf")
plot_data_scatter(model, view = "proteome", factor = 7, features = c("Q14376;Q14376-2_proteome", "P12109_proteome", "A0A1B0GVH6_proteome", "P17483_proteome"), 
                  groups = "Target", sign = "positive")
plot_data_scatter(model, view = "proteome", factor = 7, features = 5, groups = "Target", sign = "negative")
dev.off()

pdf("../output/mofa_factor7_group3-prot_scatter_first_comp1.pdf")
plot_data_scatter(model, view = "proteome", factor = 7, features = c("Q96Q89;Q96Q89-2;Q96Q89-3;Q96Q89-4_proteome", "A0A1B0GVH6_proteome", "Q16515-2_proteome", "P17483_proteome"),
                  groups = "Severe", sign = "positive")
plot_data_scatter(model, view = "proteome", factor = 7, features = 6, groups = "Severe", sign = "negative")
dev.off()

## factor 8
plot_weights(model,view = "metabolome", factor = 8, nfeatures = 10, scale = T, text_size=3.0)
plot_top_weights(model, view = "metabolome", factor = 8, nfeatures = 10, scale = T) 

plot_weights(model,view = "methylome", factor = 8, nfeatures = 10, scale = T, text_size=3.0)
plot_top_weights(model, view = "methylome", factor = 8, nfeatures = 10, scale = T) 

plot_weights(model,view = "proteome", factor = 8, nfeatures = 10, scale = T, text_size=2.5)
plot_top_weights(model, view = "proteome", factor = 8, nfeatures = 10, scale = T) 

pdf("../output/mofa_factor8_weight-metabo_first_comp1.pdf")
plot_weights(model,view = "metabolome", factor = 8, nfeatures = 10, scale = T, text_size=3.0)
plot_top_weights(model, view = "metabolome", factor = 8, nfeatures = 10, scale = T) 
dev.off()

pdf("../output/mofa_factor8_weight-methyl_first_comp1.pdf")
plot_weights(model,view = "methylome", factor = 8, nfeatures = 10, scale = T, text_size=3.0)
plot_top_weights(model, view = "methylome", factor = 8, nfeatures = 10, scale = T) 
dev.off()

pdf("../output/mofa_factor8_weight-prot_first_comp1.pdf")
plot_weights(model,view = "proteome", factor = 8, nfeatures = 10, scale = T, text_size=2.5)
plot_top_weights(model, view = "proteome", factor = 8, nfeatures = 10, scale = T) 
dev.off()

pdf("../output/mofa_factor8_metabolome_corr-heatmap_first_comp1.pdf")
plot_data_heatmap(model, view = "metabolome", factor = 8, features = 20, groups="Target", denoise = TRUE, cluster_rows = TRUE, 
                  cluster_cols = FALSE, show_rownames = TRUE, show_colnames = FALSE, scale = "row")

plot_data_heatmap(model, view = "proteome", factor = 8, features = 20, groups="Target", denoise = TRUE, cluster_rows = TRUE, 
                  cluster_cols = FALSE, show_rownames = TRUE, show_colnames = FALSE, scale = "row")
dev.off()

#help(plot_data_scatter)
# scatter plot
pdf("../output/mofa_factor8_group1-metabo_scatter_first_comp1.pdf")
plot_data_scatter(model, view = "metabolome", factor = 8, features = 9, groups = "Normal", sign = "positive")
plot_data_scatter(model, view = "metabolome", factor = 8, features = c("met719_metabolome", "met736_metabolome"), groups = "Normal", sign = "negative")
dev.off()

pdf("../output/mofa_factor8_group2-metabo_scatter_first_comp1.pdf")
plot_data_scatter(model, view = "metabolome", factor = 8, features = "met548_metabolome", groups = "Target", sign = "positive")
plot_data_scatter(model, view = "metabolome", factor = 8, features = c("met363_metabolome", "met364_metabolome"), groups = "Target", sign = "negative")
dev.off()

pdf("../output/mofa_factor8_group3-metabo_scatter_first_comp1.pdf")
plot_data_scatter(model, view = "metabolome", factor = 8, features = c("met115_metabolome", "met548_metabolome"),
                  groups = "Severe", sign = "positive")
plot_data_scatter(model, view = "metabolome", factor = 8, features = c("met364_metabolome", "met348_metabolome"),
                  groups = "Severe", sign = "negative")
dev.off()

pdf("../output/mofa_factor8_group1-prot_scatter_first_comp1.pdf")
plot_data_scatter(model, view = "proteome", factor = 8, features = 4, groups = "Normal", sign = "positive")
plot_data_scatter(model, view = "proteome", factor = 8, features = c("O95452_proteome", "Q9H257;Q9H257-2;Q9H257-3_proteome",
                                                                     "O95477_proteome", "Q9NZ43-3_proteome", "Q96QB1;Q96QB1-3;Q96QB1-5_proteome", "P78563;P78563-2;P78563-4;P78563-5_proteome", "P05019;P05019-2;P05019-3;P05019-4_proteome"), groups = "Normal", sign = "negative")
dev.off()

pdf("../output/mofa_factor8_group2-prot_scatter_first_comp1.pdf")
plot_data_scatter(model, view = "proteome", factor = 8, features = 6, groups = "Target", sign = "positive")
plot_data_scatter(model, view = "proteome", factor = 8, features = 9, groups = "Target", sign = "negative")
dev.off()

pdf("../output/mofa_factor8_group3-prot_scatter_first_comp1.pdf")
plot_data_scatter(model, view = "proteome", factor = 8, features = c("P52701;P52701-3;P52701-4_proteome", "P02763_proteome",
                                                                     "P00738_proteome", "Q8NI27_proteome"), groups = "Severe", sign = "positive")
plot_data_scatter(model, view = "proteome", factor = 8, features = "Q9H257;Q9H257-2;Q9H257-3_proteome", groups = "Severe", sign = "negative")
dev.off()

## factor 9
plot_weights(model,view = "metabolome", factor = 9, nfeatures = 10, scale = T, text_size=3.0)
plot_top_weights(model, view = "metabolome", factor = 9, nfeatures = 10, scale = T) 

plot_weights(model,view = "methylome", factor = 9, nfeatures = 10, scale = T, text_size=3.0)
plot_top_weights(model, view = "methylome", factor = 9, nfeatures = 10, scale = T) 

plot_weights(model,view = "proteome", factor = 9, nfeatures = 10, scale = T, text_size=2.5)
plot_top_weights(model, view = "proteome", factor = 9, nfeatures = 10, scale = T) 

pdf("../output/mofa_factor9_weight-metabo_first_comp1.pdf")
plot_weights(model,view = "metabolome", factor = 9, nfeatures = 10, scale = T, text_size=3.0)
plot_top_weights(model, view = "metabolome", factor = 9, nfeatures = 10, scale = T) 
dev.off()

pdf("../output/mofa_factor9_weight-methyl_first_comp1.pdf")
plot_weights(model,view = "methylome", factor = 9, nfeatures = 10, scale = T, text_size=3.0)
plot_top_weights(model, view = "methylome", factor = 9, nfeatures = 10, scale = T) 
dev.off()

pdf("../output/mofa_factor9_weight-prot_first_comp1.pdf")
plot_weights(model,view = "proteome", factor = 9, nfeatures = 10, scale = T, text_size=2.5)
plot_top_weights(model, view = "proteome", factor = 9, nfeatures = 10, scale = T) 
dev.off()

pdf("../output/mofa_factor9_metabolome_corr-heatmap_first_comp1.pdf")
plot_data_heatmap(model, view = "metabolome", factor = 9, features = 20, groups = "Target",denoise = TRUE, cluster_rows = TRUE, 
                  cluster_cols = FALSE, show_rownames = TRUE, show_colnames = FALSE, scale = "row")

plot_data_heatmap(model, view = "proteome", factor = 9, features = 20, groups="Severe", denoise = TRUE, cluster_rows = TRUE, 
                  cluster_cols = FALSE, show_rownames = TRUE, show_colnames = FALSE, scale = "row")
dev.off()

#help(plot_data_scatter)
# scatter plot
pdf("../output/mofa_factor9_group1-metabo_scatter_first_comp1.pdf")
plot_data_scatter(model, view = "metabolome", factor = 9, features="met1074_metabolome", groups = "Normal", sign = "positive")
plot_data_scatter(model, view = "metabolome", factor = 9, features=c("met73_metabolome", "met303_metabolome",
                  "met943_metabolome", "met467_metabolome", "met441_metabolome"), groups = "Normal", sign = "negative")
dev.off()

pdf("../output/mofa_factor9_group2-metabo_scatter_first_comp1.pdf")
plot_data_scatter(model, view = "metabolome", factor = 9, features = 4, groups = "Target", sign = "positive")
plot_data_scatter(model, view = "metabolome", factor = 9, features = c("met347_metabolome", "met254_metabolome", "met996_metabolome"), groups = "Target", sign = "negative")
dev.off()

pdf("../output/mofa_factor9_group1-prot_scatter_first_comp1.pdf")
plot_data_scatter(model, view = "proteome", factor = 9, features = 9, groups = "Normal", sign = "positive")
plot_data_scatter(model, view = "proteome", factor = 9, features = 9, groups = "Normal", sign = "negative")
dev.off()

pdf("../output/mofa_factor9_group2-prot_scatter_first_comp1.pdf")
plot_data_scatter(model, view = "proteome", factor = 9, features = 6, groups = "Target", sign = "positive")
plot_data_scatter(model, view = "proteome", factor = 9, features = 6, groups = "Target", sign = "negative")
dev.off()

pdf("../output/mofa_factor9_group3-prot_scatter_first_comp1.pdf")
plot_data_scatter(model, view = "proteome", factor = 9, features = "Q86UX2;Q86UX2-3_proteome", groups = "Severe", sign = "positive")
plot_data_scatter(model, view = "proteome", factor = 9, features = , groups = "Severe", sign = "negative")
dev.off()

## best scatter plots
pdf("../output/mofa_factor9_group2-methyl_scatter_by-condition_first_comp1.pdf")
plot_data_scatter(model, view = "methylome", factor = 9, features = 12, groups = "Target", 
                  add_lm = TRUE, color_by = "condition")
dev.off()

pdf("../output/mofa_factor1_group2-metabo_scatter_by-condition_first_comp1.pdf")
plot_data_scatter(model, view = "metabolome", factor = 1, features = 12, groups = "Target", 
                  add_lm = TRUE, color_by = "condition")
dev.off()

pdf("../output/mofa_factor3_group2-metabo_scatter_by-condition_first_comp1.pdf")
plot_data_scatter(model, view = "metabolome", factor = 3, features = 3, groups = "Target", 
                  add_lm = TRUE, color_by = "condition", text_size = 3.2)
dev.off()

pdf("../output/mofa_factor5_group2-methyl_scatter_by-condition_first_comp1.pdf")
plot_data_scatter(model, view = "methylome", factor = 5, features = 12, groups = "Target", 
                  add_lm = TRUE, color_by = "condition")
dev.off()

pdf("../output/mofa_factor6_group2-prot_scatter_by-condition_first_comp1.pdf")
plot_data_scatter(model, view = "proteome", factor = 6, features = 12, groups = "Target", 
                  add_lm = TRUE, color_by = "condition")
dev.off()

pdf("../output/mofa_factor7_group2-prot_scatter_by-condition_first_comp1.pdf")
plot_data_scatter(model, view = "proteome", factor = 7, features = 12, groups = "Target", 
                  add_lm = TRUE, color_by = "condition")
dev.off()

pdf("../output/mofa_factor8_group2-methyl_scatter_by-condition_first_comp1.pdf")
plot_data_scatter(model, view = "methylome", factor = 8, features = 6, groups = "Target", 
                  add_lm = TRUE, color_by = "condition")
dev.off()

## extracting data
factors <- get_factors(model, factors="all")
factors_target <- get_factors(model, groups = "Target", factors="all")
lapply(factors, dim)
lapply(factors_target, dim)

weights <- get_weights(model, views = "all", factors= "all")
lapply(weights, dim)

data <- get_data(model)
data_target <- get_data(model, groups="Target")

features_names(model)[["proteome"]]

save(factors, weights, data, model, metadata, 
     file="../rdata/mofa_first_comp1_add_study.RData")