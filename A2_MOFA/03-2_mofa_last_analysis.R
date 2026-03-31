rm(list = ls()) 
getwd() 
setwd("D:/Analysis/MOFA/input/")

library(ggplot2)
library(MOFA2)

# visualization
library(corrplot)
library(RColorBrewer)
library(psych)

### comp1 (normal vs target vs critical, three groups)
model = load_model("multiome_last_comp1.hdf5")  # warning : factor 2

# Add metadata 
setwd("../output/")
metadata = data.frame(fread("clinical_last_exposed_comp1_data.csv"), stringsAsFactors = FALSE)
metadata$V1 = NULL

samples_metadata(model) <- metadata    
model@samples_metadata$sample 

## Model overview
slotNames(model)
plot_data_overview(model)  
names(model@data) 

p1 <- plot_data_overview(model)
ggsave("MOFA_high-res_last_distribution_comp1.png", plot=p1, 
       width=7.6, height=7.3, dpi=300)

model@data$metabolome[[1]][1:5,1:5]
model@data$methylome[[1]][1:5,1:5]
model@data$proteome[[1]][1:5,1:5]
model@data$transcriptome[[1]][1:5,1:5]

model@samples_metadata

# expectation
dim(model@expectations$Z$Normal)    # Check the values for the 20 samples x 11 factors
model@expectations$Z$Normal[1:5,]

dim(model@expectations$Z$Target)    # Check the values for the 40 samples x 11 factors
model@expectations$Z$Target[1:5,]

dim(model@expectations$Z$Severe)    # Check the values for the 10 samples x 11 factors
model@expectations$Z$Severe[1:5,]

plot_factor_cor(model)               # Identify the correlation between factors.

pdf("../output/mofa_correlation_factor_last_comp1_data.pdf")
plot_factor_cor(model)
dev.off()

# Examine the feature loadings (weight) for the latent facotrs. 
dim(model@expectations$W$metabolome) 
model@expectations$W$metabolome[1:5,]
model@expectations$W$methylome[1:5,]
model@expectations$W$proteome[1:5,]
model@expectations$W$transcriptome[1:5,]

## variance explained (R square) 
# How well the inferred factors explain the variance in the input views (omics) and groups
# Total variance explained per view and group
head(model@cache$variance_explained$r2_total[[1]]) 
head(model@cache$variance_explained$r2_per_factor[[1]])

model@cache$variance_explained$r2_per_factor

plot_variance_explained(model, x="view", y="factor")  
plot_variance_explained(model, x="group", y="factor", plot_total = T)[[2]] 

pdf("../output/mofa_variance_explain_factor_last_comp1.pdf")
plot_variance_explained(model, x="view", y="factor")  
dev.off()

pdf("../output/mofa_variance_factor-omics_last_comp1.pdf", height=7, width=13)
plot_variance_explained(model, x="factor", y="view")  
dev.off()

pdf("../output/mofa_variance_explain_total_last_comp1.pdf")
plot_variance_explained(model, x="group", y="factor", plot_total = T)[[2]]  
dev.off()

# save figure
p2 <- plot_variance_explained(model, x="view", y="factor") 
ggsave("MOFA_high-res_variance_explain_factor_last_comp1.png", plot=p2, 
       width=10.5, height=7.3, dpi=300)

correlate_factors_with_covariates(model,
                                  covariates = c("sex", "age"), 
                                  plot="log_pval")

correlate_factors_with_covariates(model,
                                  covariates = c("FEV1FVC_z_L", "FEV1_z_L", "FVC_z_L", "FEF_z_L"), 
                                  plot="log_pval")

correlate_factors_with_covariates(model,
                                  covariates = c("FEV1FVC_z_L", "FEV1_z_L", "FVC_z_L", "FEF_z_L"),
                                  groups="Normal",
                                  plot="log_pval")

correlate_factors_with_covariates(model,
                                  covariates = c("FEV1FVC_z_L", "FEV1_z_L", "FVC_z_L", "FEF_z_L"),
                                  groups="Target",
                                  plot="log_pval")

correlate_factors_with_covariates(model,
                                  covariates = c("FEV1FVC_z_L", "FEV1_z_L", "FVC_z_L", "FEF_z_L"),
                                  groups="Severe",
                                  plot="log_pval")

save.image("../rdata/mofa_last_comp1.RData")
load("../rdata/mofa_last_comp1.RData")

# correlation analysis and visualization 
groups <- model@samples_metadata$group
target_index <- which(groups=="Target") 

normal_index <- which(groups=="Normal")
severe_index <- which(groups=="Severe")

factors <- model@expectations$Z
target_factors <- factors[["Target"]]  

normal_factors <- factors[["Normal"]] 
severe_factors <- factors[["Severe"]] 

covariates <- model@samples_metadata

filter_cols <- c("age", "FEV1_z_L", "FVC_z_L", "FEV1FVC_z_L", "FEF_z_L")
target_covariates <- covariates[target_index, filter_cols]

normal_covariates <- covariates[normal_index, filter_cols]
severe_covariates <- covariates[severe_index, filter_cols]

cor_matrix <- cor(target_factors, target_covariates, use="pairwise.complete.obs") # pearson
cor_target <- corr.test(target_factors, target_covariates, use='pairwise', method='pearson', adjust='BH')

cor_matrix2 <- cor(normal_factors, normal_covariates, use="pairwise.complete.obs") # pearson
cor_normal <- corr.test(normal_factors, normal_covariates, use='pairwise', method='pearson', adjust='BH')

cor_matrix3 <- cor(severe_factors, severe_covariates, use="pairwise.complete.obs") # pearson
cor_severe <- corr.test(severe_factors, severe_covariates, use='pairwise', method='pearson', adjust='BH')

corr.test(target_factors, target_covariates, use='pairwise', method='pearson', adjust='BH')
p_matrix <- cor_target$p

p_matrix2 <- cor_normal$p
p_matrix3 <- cor_severe$p

cor_matrix <- t(cor_matrix)
p_matrix <- t(p_matrix)

cor_matrix2 <- t(cor_matrix2)
p_matrix2 <- t(p_matrix2)

cor_matrix3 <- t(cor_matrix3)
p_matrix3 <- t(p_matrix3)

## cor_matrix = target
png("correlation_target_last_comp1.png", width=2535, height=1260, res=300)
par(cex=0.9)
corrplot(cor_matrix, method='circle', col=rev(brewer.pal(n=7, name="RdYlBu")),
         p.mat=p_matrix, insig='label_sig',
         tl.col="black", tl.srt=45, tl.cex=1.1, tl.pos="lt", col.lim=c(-0.5, 0.5), cl.pos='b', cl.cex=1.1, cl.ratio=0.2)
dev.off()

## cor_matrix2 = normal
corrplot(cor_matrix2, method='circle', col=rev(brewer.pal(n=7, name="RdYlBu")),
         p.mat=p_matrix2, insig='label_sig',
         tl.col="black", tl.srt=45, col.lim=c(-1.0, 1.0), cl.length=5, cl.pos='b', cl.cex=1.0, cl.ratio=0.2)

## cor_matrix3 = severe
corrplot(cor_matrix3, method='circle', col=rev(brewer.pal(n=7, name="RdYlBu")),
         p.mat=p_matrix3, insig='label_sig',
         tl.col="black", tl.srt=45, col.lim=c(-1.0, 1.0), cl.length=5, cl.pos='b', cl.cex=1.0, cl.ratio=0.2)

write.csv(cor_matrix, "../output/target_correlation_factor_and_clinical_comp1_last.csv")
write.csv(t(p_matrix), "../output/target_p-value_factor_and_clinical_comp1_last.csv")

## change feature name
backup <- features_names(model)[["proteome"]]
new_name <- sapply(strsplit(features_names(model)[["proteome"]], "_"), `[`, 1)
features_names(model)[["proteome"]] <- new_name

backup2 <- features_names(model)[["metabolome"]]
new_name <- sapply(strsplit(features_names(model)[["metabolome"]], "_"), `[`, 1)
features_names(model)[["metabolome"]] <- new_name

backup3 <- features_names(model)[["methylome"]]
new_name <- sapply(strsplit(features_names(model)[["methylome"]], "_"), `[`, 1)
features_names(model)[["methylome"]] <- new_name

backup4 <- features_names(model)[["transcriptome"]]
new_name <- sapply(strsplit(features_names(model)[["transcriptome"]], "_"), `[`, 1)
features_names(model)[["transcriptome"]] <- new_name

# Factor's distribution by condition
p <- plot_factor(model, 
                 factors = c(2, 6),
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

## factor 2
pdf("../output/last_data/mofa_factor2_weight-metabo_last_comp1.pdf")
plot_weights(model,view = "metabolome", factor = 2, nfeatures = 10, scale = T, text_size=3.0)
plot_top_weights(model, view = "metabolome", factor = 2, nfeatures = 10, scale = T) 
dev.off()

pdf("../output/last_data/mofa_factor2_weight-methyl_last_comp1.pdf")
plot_weights(model,view = "methylome", factor = 2, nfeatures = 10, scale = T, text_size=3.0)
plot_top_weights(model, view = "methylome", factor = 2, nfeatures = 10, scale = T) 
dev.off()

pdf("../output/last_data/mofa_factor2_weight-prot_last_comp1.pdf")
plot_weights(model,view = "proteome", factor = 2, nfeatures = 10, scale = T, text_size=2.5)
plot_top_weights(model, view = "proteome", factor = 2, nfeatures = 10, scale = T) 
dev.off()

pdf("../output/last_data/mofa_factor2_weight-trans_last_comp1.pdf")
plot_weights(model,view = "transcriptome", factor = 2, nfeatures = 10, scale = T, text_size=2.5)
plot_top_weights(model, view = "transcriptome", factor = 2, nfeatures = 10, scale = T) 
dev.off()

plot_data_heatmap(model, view = "transcriptome", factor = 2, groups = "Target", features = 20, denoise = TRUE, cluster_rows = TRUE, 
                  cluster_cols = FALSE, show_rownames = TRUE, show_colnames = FALSE, scale = "row")

pdf("../output/last_data/mofa_factor2_group1-trans_scatter_last_comp1.pdf")
plot_data_scatter(model, view = "transcriptome", factor = 2, features = 12, groups = "Normal", sign = "positive")
plot_data_scatter(model, view = "transcriptome", factor = 2, features = 12, groups = "Normal", sign = "negative")
dev.off()

pdf("../output/last_data/mofa_factor2_group2-trans_scatter_last_comp1.pdf")
plot_data_scatter(model, view = "transcriptome", factor = 2, features = 12, groups = "Target", sign = "positive")
plot_data_scatter(model, view = "transcriptome", factor = 2, features = 12, groups = "Target", sign = "negative")
dev.off()

pdf("../output/last_data/mofa_factor2_group3-trans_scatter_last_comp1.pdf")
plot_data_scatter(model, view = "transcriptome", factor = 2, features = 9, groups = "Severe", sign = "positive")
plot_data_scatter(model, view = "transcriptome", factor = 2, features = 12, groups = "Severe", sign = "negative")
dev.off()

## factor 5
plot_weights(model,view = "metabolome", factor = 5, nfeatures = 10, scale = T, text_size=3.0)
plot_top_weights(model, view = "metabolome", factor = 5, nfeatures = 10, scale = T) 

plot_weights(model,view = "methylome", factor = 5, nfeatures = 10, scale = T, text_size=3.0)
plot_top_weights(model, view = "methylome", factor = 5, nfeatures = 10, scale = T) 

plot_weights(model,view = "proteome", factor = 5, nfeatures = 10, scale = T, text_size=2.5)
plot_top_weights(model, view = "proteome", factor = 5, nfeatures = 10, scale = T) 

plot_weights(model,view = "transcriptome", factor = 5, nfeatures = 10, scale = T, text_size=2.5)
plot_top_weights(model, view = "transcriptome", factor = 5, nfeatures = 10, scale = T) 

pdf("../output/last_data/comp1/mofa_factor5_weight-metabo_last_comp1.pdf")
plot_weights(model,view = "metabolome", factor = 5, nfeatures = 10, scale = T, text_size=3.0)
plot_top_weights(model, view = "metabolome", factor = 5, nfeatures = 10, scale = T) 
dev.off()

pdf("../output/last_data/comp1/mofa_factor5_weight-methyl_last_comp1.pdf")
plot_weights(model,view = "methylome", factor = 5, nfeatures = 10, scale = T, text_size=3.0)
plot_top_weights(model, view = "methylome", factor = 5, nfeatures = 10, scale = T) 
dev.off()

pdf("../output/last_data/comp1/mofa_factor5_weight-prot_last_comp1.pdf")
plot_weights(model,view = "proteome", factor = 5, nfeatures = 10, scale = T, text_size=2.5)
plot_top_weights(model, view = "proteome", factor = 5, nfeatures = 10, scale = T) 
dev.off()

pdf("../output/last_data/comp1/mofa_factor5_weight-trans_last_comp1.pdf")
plot_weights(model,view = "transcriptome", factor = 5, nfeatures = 10, scale = T, text_size=2.5)
plot_top_weights(model, view = "transcriptome", factor = 5, nfeatures = 10, scale = T) 
dev.off()

help(plot_data_heatmap)
# heatmap - association with specific factors by each factor
pdf("../output/last_data/comp1/mofa_factor5_proteome_corr-heatmap_last_comp1.pdf")
plot_data_heatmap(model, view = "proteome", factor = 5, features = 20, denoise = TRUE, cluster_rows = TRUE, 
                  cluster_cols = FALSE, show_rownames = TRUE, show_colnames = FALSE, scale = "row")
dev.off()

#help(plot_data_scatter)
# scatter plot
pdf("../output/last_data/comp1/mofa_factor5_group1-prot_scatter_last_comp1.pdf")
plot_data_scatter(model, view = "proteome", factor = 5, features = 12, groups = "Normal", sign = "positive")
plot_data_scatter(model, view = "proteome", factor = 5, features = 12, groups = "Normal", sign = "negative")
dev.off()

pdf("../output/last_data/comp1/mofa_factor5_group2-prot_scatter_last_comp1.pdf")
plot_data_scatter(model, view = "proteome", factor = 5, features = 12, groups = "Target", sign = "positive")
plot_data_scatter(model, view = "proteome", factor = 5, features = 12, groups = "Target", sign = "negative")
dev.off()

pdf("../output/last_data/comp1/mofa_factor5_group3-prot_scatter_last_comp1.pdf")
plot_data_scatter(model, view = "proteome", factor = 5, features = 9, groups = "Severe", sign = "positive")
plot_data_scatter(model, view = "proteome", factor = 5, features = 12, groups = "Severe", sign = "negative")
dev.off()

## factor 6
plot_weights(model,view = "metabolome", factor = 6, nfeatures = 10, scale = T, text_size=3.0)
plot_top_weights(model, view = "metabolome", factor = 6, nfeatures = 10, scale = T) 

plot_weights(model,view = "methylome", factor = 6, nfeatures = 10, scale = T, text_size=3.0)
plot_top_weights(model, view = "methylome", factor = 6, nfeatures = 10, scale = T) 

plot_weights(model,view = "proteome", factor = 6, nfeatures = 10, scale = T, text_size=2.5)
plot_top_weights(model, view = "proteome", factor = 6, nfeatures = 10, scale = T) 

plot_weights(model,view = "transcriptome", factor = 6, nfeatures = 10, scale = T, text_size=2.5)
plot_top_weights(model, view = "transcriptome", factor = 6, nfeatures = 10, scale = T) 

pdf("../output/last_data/mofa_factor6_weight-metabo_last_comp1.pdf")
plot_weights(model,view = "metabolome", factor = 6, nfeatures = 10, scale = T, text_size=3.0)
plot_top_weights(model, view = "metabolome", factor = 6, nfeatures = 10, scale = T) 
dev.off()

pdf("../output/last_data/mofa_factor6_weight-methyl_last_supp.pdf")
plot_weights(model,view = "methylome", factor = 6, nfeatures = 10, scale = T, text_size=3.0)
plot_top_weights(model, view = "methylome", factor = 6, nfeatures = 10, scale = T) 
dev.off()

pdf("../output/last_data/mofa_factor6_weight-prot_name-clear_last_supp.pdf")
plot_weights(model,view = "proteome", factor = 6, nfeatures = 10, scale = T, text_size=2.5)
plot_top_weights(model, view = "proteome", factor = 6, nfeatures = 10, scale = T) 
dev.off()

pdf("../output/last_data/mofa_factor6_weight-trans_last_supp.pdf")
plot_weights(model,view = "transcriptome", factor = 6, nfeatures = 10, scale = T, text_size=2.5)
plot_top_weights(model, view = "transcriptome", factor = 6, nfeatures = 10, scale = T) 
dev.off()

## save figure (figure 4D)
label_map <- c(
  "P39060;P39060-1"="COL18A1 (isoform 1)",
  "P02750"="LRG1",
  "P51884"="LUM",
  "P22891;P22891-2"="PROZ (isoform 1-2)",
  "P80108"="GPLD1",
  "P02652"="APOA2",
  "P00742"="F10",
  "P29622"="SERPI4",
  "P02675"="FGB",
  "P02679;P02679-2"="FGG (isoform 1-2)"
)

p3 <- plot_top_weights(model, view = "proteome", factor = 6, nfeatures = 10, scale = T) 

p3 <- p3 + scale_x_discrete(labels = function(x) ifelse(x %in% names(label_map), 
                                                  label_map[x], x))
p3 <- p3 + theme(
  axis.text.y=element_text(size=11)
)

p3

ggsave("mofa_factor6_weight-prot_high-res2_last_comp1.png", plot=p3,
       width=7.7, height=7.7, 
       dpi=300)

help(plot_data_heatmap)

pdf("../output/last_data/comp1/mofa_factor6_metabolome_corr-heatmap_last_comp1.pdf")
plot_data_heatmap(model, view = "metabolome", groups = 'Target', factor = 6, features = 20, denoise = TRUE, cluster_rows = TRUE, 
                  cluster_cols = FALSE, show_rownames = TRUE, show_colnames = FALSE, scale = "row")
dev.off()

#help(plot_data_scatter)
# scatter plot
pdf("../output/last_data/mofa_factor6_group1-metabo_scatter_last_comp1.pdf")
plot_data_scatter(model, view = "metabolome", factor = 6, features = 12, groups = "Normal", sign = "positive")
plot_data_scatter(model, view = "metabolome", factor = 6, features = 12, groups = "Normal", sign = "negative")
dev.off()

pdf("../output/last_data/mofa_factor6_group2-metabo_scatter_last_comp1.pdf")
plot_data_scatter(model, view = "metabolome", factor = 6, features = c("met1079", "met800", "met527", "met1053", "met1037", "met188", "met1047", "met1001", "met589"),
                  groups = "Target", sign = "positive")
plot_data_scatter(model, view = "metabolome", factor = 6, features = 12, groups = "Target", sign = "negative")
dev.off()

pdf("../output/last_data/mofa_factor6_group3-metabo_scatter_last_comp1.pdf")
plot_data_scatter(model, view = "metabolome", factor = 6, features = c("met1092", "met1070", "met1053", "met1037", "met1047", "met1001", "met589"),
                  groups = "Severe", sign = "positive")
plot_data_scatter(model, view = "metabolome", factor = 6, features = 12, groups = "Severe", sign = "negative")
dev.off()

## factor 9
plot_weights(model,view = "metabolome", factor = 9, nfeatures = 10, scale = T, text_size=3.0)
plot_top_weights(model, view = "metabolome", factor = 9, nfeatures = 10, scale = T) 

plot_weights(model,view = "methylome", factor = 9, nfeatures = 10, scale = T, text_size=3.0)
plot_top_weights(model, view = "methylome", factor = 9, nfeatures = 10, scale = T) 

plot_weights(model,view = "proteome", factor = 9, nfeatures = 10, scale = T, text_size=2.5)
plot_top_weights(model, view = "proteome", factor = 9, nfeatures = 10, scale = T) 

plot_weights(model,view = "transcriptome", factor = 9, nfeatures = 10, scale = T, text_size=2.5)
plot_top_weights(model, view = "transcriptome", factor = 9, nfeatures = 10, scale = T) 

pdf("../output/last_data/mofa_factor9_weight-metabo_last_comp1.pdf")
plot_weights(model,view = "metabolome", factor = 9, nfeatures = 10, scale = T, text_size=3.0)
plot_top_weights(model, view = "metabolome", factor = 9, nfeatures = 10, scale = T) 
dev.off()

pdf("../output/last_data/mofa_factor9_weight-methyl_last_comp1.pdf")
plot_weights(model,view = "methylome", factor = 9, nfeatures = 10, scale = T, text_size=3.0)
plot_top_weights(model, view = "methylome", factor = 9, nfeatures = 10, scale = T) 
dev.off()

pdf("../output/last_data/mofa_factor9_weight-prot_clear-name_last_comp1.pdf")
plot_weights(model,view = "proteome", factor = 9, nfeatures = 10, scale = T, text_size=2.5)
plot_top_weights(model, view = "proteome", factor = 9, nfeatures = 10, scale = T) 
dev.off()

pdf("../output/last_data/mofa_factor9_weight-trans_last_comp1.pdf")
plot_weights(model,view = "transcriptome", factor = 9, nfeatures = 10, scale = T, text_size=2.5)
plot_top_weights(model, view = "transcriptome", factor = 9, nfeatures = 10, scale = T) 
dev.off()

pdf("../output/mofa_factor9_proteome_corr-heatmap_last_comp1.pdf")
plot_data_heatmap(model, view = "proteome", factor = 9, features = 20, groups="Severe", denoise = TRUE, cluster_rows = TRUE, 
                  cluster_cols = FALSE, show_rownames = TRUE, show_colnames = FALSE, scale = "row")
dev.off()

# scatter plot
pdf("../output/last_data/mofa_factor9_group1-prot_scatter_last_comp1.pdf")
#plot_data_scatter(model, view = "proteome", factor = 9, features = 12, groups = "Normal", sign = "positive")
plot_data_scatter(model, view = "proteome", factor = 9, features = c("P05156", "P00751", "P02743", "P18428"),
                  groups = "Normal", sign = "positive")
plot_data_scatter(model, view = "proteome", factor = 9, features = 12, groups = "Normal", sign = "negative")
#plot_data_scatter(model, view = "proteome", factor = 9, features = c("Q9BUN1_proteome", "Q9NZJ4;Q9NZJ4-2_proteome", "Q5TCS8_proteome", "P49747_proteome", "Q96Q40-3;Q96Q40-4;Q96Q40-5_proteome", "P13591-1;P13591-3;P13591-4_proteome", "Q9NZ43-3_proteome", "Q96QB1;Q96QB1-3;Q96QB1-5_proteome",
                                                                     #"Q6UY14;Q6UY14-2;Q6UY14-3_proteome", "Q6UXB8_proteome", "P51884_proteome"), groups = "Normal", sign = "negative")
dev.off()

pdf("../output/last_data/mofa_factor9_group2-prot_scatter_last_comp1.pdf")
plot_data_scatter(model, view = "proteome", factor = 9, features = c("P02763", "Q92954-3", "P05156", "Q96S38;Q96S38-2", "P00751", "P02743", "P18428"), groups = "Target", sign = "positive")
plot_data_scatter(model, view = "proteome", factor = 9, features = c("Q9BUN1","P49747", "P13591;P13591-1;P13591-3;P13591-4", "Q96QB1;Q96QB1-3;Q96QB1-5", "Q6UY14;Q6UY14-2;Q6UY14-3", "Q6UXB8", "P04278", "P51884"), groups = "Target", sign = "negative")
dev.off()

pdf("../output/last_data/mofa_factor9_group2-prot_scatter_last_comp1.pdf")
plot_data_scatter(model, view = "proteome", factor = 9, features = c("P02763", "Q92954-3", "P05156", "Q96S38;Q96S38-2", "P00751", "P02743", "P18428"), groups = "Target", sign = "positive")
plot_data_scatter(model, view = "proteome", factor = 9, features = c("Q9BUN1","P49747", "P13591;P13591-1;P13591-3;P13591-4", "Q96QB1;Q96QB1-3;Q96QB1-5", "Q6UY14;Q6UY14-2;Q6UY14-3", "Q6UXB8", "P04278", "P51884"), groups = "Target", sign = "negative")
dev.off() 

pdf("../output/last_data/mofa_factor9_group3-prot_scatter_last_comp1.pdf")
plot_data_scatter(model, view = "proteome", factor = 9, features = c("P05156", "P0DJI8", "Q03591", "P16930", "Q92954-3", "P18428", "P02743", "P00738", "Q8NI27"),
                  groups = "Severe", sign = "positive")
plot_data_scatter(model, view = "proteome", factor = 9, features = 12, groups = "Severe", sign = "negative")
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
     file="../rdata/mofa_last_comp1_add_study.RData")