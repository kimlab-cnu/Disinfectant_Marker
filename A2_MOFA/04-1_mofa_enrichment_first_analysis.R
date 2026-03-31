rm(list=ls())

library(ggplot2)
library(data.table)
library(cowplot)
library(purrr)
library(MOFA2)
library(dplyr)

setwd("D:/Analysis/MOFA/input/")
getwd()

load("../rdata/mofa_first_comp1_add_study.RData")

model

#write.csv(factors$group_1, "../output/factor_value_normal_first_comp1.csv")
#write.csv(factors$group_2, "../output/factor_value_target_first_comp1.csv")
#write.csv(factors$group_3, "../output/factor_value_critical_first_comp1.csv")

progene_map <- data.frame(fread("protein-gene_mapping.csv"), stringsAsFactors = FALSE)

features_names(model)$proteome
prot_name <- features_names(model)$proteome

## Process redundant gene identifiers
duplicated_gene <- duplicated(progene_map$Gene)
progene_map$Gene[duplicated_gene]

progene_map$new_Gene <- ave(progene_map$Gene, progene_map$Gene, FUN = function(x) {
  if (length(x) > 1) {
    paste0(x, "_", seq_along(x))
  } else {
    x
  }
})         # If they're overlapped, received the number by row orders (ex. APOM_1, APOM_2 ....)

duplicated_gene <- duplicated(progene_map$new_Gene)
progene_map$new_Gene[duplicated_gene]  # Identify no overlap (unique)

features_names(model)$proteome <- progene_map$new_Gene
features_names(model)$proteome

## enrichment analysis
library(msigdbr)  
library(BiocParallel)
library(tibble)

#MsigDB HALLMARK gene sets
hall_gene = msigdbr(species="human", category="H")
hall_list = split(x = hall_gene$gene_symbol, f = hall_gene$gs_name)
hall_list %>% head() %>% lapply(head)

# Update subcollection and names due to recent updates. 
#MsigDB KEGG pathway
kegg_gene <- msigdbr(species = "human", category = "C2", subcategory = "CP:KEGG")
kegg_list = split(x = kegg_gene$gene_symbol, f = kegg_gene$gs_name)

#MsigDB reactome pathway
reactome_gene <- msigdbr(species = "human", category = "C2", subcategory = "CP:REACTOME")
reactome_list = split(x = reactome_gene$gene_symbol, f = reactome_gene$gs_name)

#MsigDB GO pathway
GO_gene <- msigdbr(species = "human", category = "C5", subcategory = "GO:BP")
GO_list = split(x = GO_gene$gene_symbol, f = GO_gene$gs_name)

all_gene_list <- c(hall_list, kegg_list, reactome_list, GO_list)
length(all_gene_list) 

## prot_genes <- features_names(model)$proteome
all_gene_names <- names(all_gene_list)

binary_matrix <- sapply(all_gene_list, function(finder_genes) {
  as.integer(prot_genes %in% finder_genes)
})

rownames(binary_matrix) <- prot_genes
colnames(binary_matrix) <- all_gene_names
binary_matrix
binary_matrix2 <- t(binary_matrix)

##

### negative
enrichment.parametric <- run_enrichment(model,
                                        view = "proteome", factors = 1:9,
                                        feature.sets = binary_matrix2,
                                        sign = "negative",
                                        statistical.test = "parametric"
)

enrichment.parametric.adj <- run_enrichment(model,
                                            view = "proteome", factors = 1:9,
                                            feature.sets = binary_matrix2,
                                            sign = "negative",
                                            statistical.test = "cor.adj.parametric"
)

### positive
enrichment.parametric <- run_enrichment(model,
                                        view = "proteome", factors = 1:9,
                                        feature.sets = binary_matrix2,
                                        sign = "positive",
                                        statistical.test = "parametric"
)

enrichment.parametric.adj <- run_enrichment(model,
                                            view = "proteome", factors = 1:9,
                                            feature.sets = binary_matrix2,
                                            sign = "positive",
                                            statistical.test = "cor.adj.parametric"
)

##

names(enrichment.parametric)
enrichment.parametric$set.statistics[1:5, 1]
enrichment.parametric$pval.adj[1:5, 1]

pdf("../output/enrichment_pos_heatmap_factor_first_comp1.pdf")
plot_enrichment_heatmap(enrichment.parametric)
dev.off()

## comp1 (first, positive); factor 1-3,6,8

#enrichment_top10-pathways_factor1_first_comp1

plot_enrichment(enrichment.parametric, 
                factor=1, max.pathways=10)

plot_enrichment(enrichment.parametric, 
                factor=3, max.pathways=10)

plot_enrichment(enrichment.parametric, 
                factor=4, max.pathways=10)  # comp1(first+pos); little

plot_enrichment(enrichment.parametric, 
                factor=5, max.pathways=10)

plot_enrichment(enrichment.parametric, 
                factor=6, max.pathways=10)

plot_enrichment(enrichment.parametric, 
                factor=7, max.pathways=10)

plot_enrichment(enrichment.parametric, 
                factor=8, max.pathways=10)

plot_enrichment(enrichment.parametric, 
                factor=9, max.pathways=10)  # comp1(first+pos); none

#enrichment_weight_factor1_first_comp1
plot_enrichment_detailed(enrichment.parametric,
                         factor=1, max.genes=5, max.pathways=7)

plot_enrichment_detailed(enrichment.parametric,
                         factor=3, max.genes=5, max.pathways=7)

plot_enrichment_detailed(enrichment.parametric,
                         factor=4, max.genes=5, max.pathways=7)

plot_enrichment_detailed(enrichment.parametric,
                         factor=5, max.genes=5, max.pathways=7)

plot_enrichment_detailed(enrichment.parametric,
                         factor=6, max.genes=5, max.pathways=7)

plot_enrichment_detailed(enrichment.parametric,
                         factor=7, max.genes=5, max.pathways=7)

plot_enrichment_detailed(enrichment.parametric,
                         factor=8, max.genes=5, max.pathways=7)

plot_enrichment_detailed(enrichment.parametric,
                         factor=9, max.genes=5, max.pathways=7)

h <- as.data.frame(enrichment.parametric$feature.sets)
h <- t(h)
write.csv(h, "../output/first_data/comp1/enrichment_pathway-result_pos-first_comp1.csv")

### Comp1 - negative
h_weight <- as.data.frame(enrichment.parametric$feature.statistics)
h_weight$GOBP_ADAPTIVE_IMMUNE_RESPONSE <- h[,"GOBP_ADAPTIVE_IMMUNE_RESPONSE"]
h_weight$REACTOME_FCGR_ACTIVATION <- h[,"REACTOME_FCGR_ACTIVATION"]
h_weight$REACTOME_ROLE_OF_PHOSPHOLIPIDS_IN_PHAGOCYTOSIS <- h[,"REACTOME_ROLE_OF_PHOSPHOLIPIDS_IN_PHAGOCYTOSIS"]
h_weight$GOBP_B_CELL_RECEPTOR_SIGNALING_PATHWAY <- h[,"GOBP_B_CELL_RECEPTOR_SIGNALING_PATHWAY"]
h_weight$REACTOME_FCGR3A_MEDIATED_IL10_SYNTHESIS <- h[,"REACTOME_FCGR3A_MEDIATED_IL10_SYNTHESIS"]
h_weight$REACTOME_PARASITE_INFECTION <- h[,"REACTOME_PARASITE_INFECTION"]
h_weight$REACTOME_FCERI_MEDIATED_MAPK_ACTIVATION <- h[,"REACTOME_FCERI_MEDIATED_MAPK_ACTIVATION"]

### comp1 - positive
h_weight <- as.data.frame(enrichment.parametric$feature.statistics)

colnames(h_weight)
write.csv(h_weight, "../output/enrichment_weight-result_first_comp1.csv")

##### enriched gene analysis to specific factor or pathway
genes <- list("IGKV1D-3", "IGHV1-3")  # adaptive_immune_response; top2 weight genes. 

genes %>% map(~ plot_factors(model, 
                             factors = 6, 
                             color_by = 'condition',       #'group', 
                             #shape_by = 'condition',
                             scale = T,
                             legend = T
)) %>% cowplot::plot_grid(plotlist=., nrow=1)

help(plot_factors)
#####

dt <- rbind(
  enrichment.parametric$pval[,c(5, 7, 8, 9)] %>% as.data.table %>% .[,c("test","pathway"):=list("parametric",1:.N)],
  enrichment.parametric.adj$pval[,c(5, 7, 8, 9)] %>% as.data.table %>% .[,c("test","pathway"):=list("parametric.adj",1:.N)]
) %>% melt(id.vars=c("test","pathway"), variable.name="factor")

ggplot(dt, aes(x=value, fill=test)) +
  facet_wrap(~factor, scales="free_y", nrow=1) +
  geom_histogram() +
  theme_bw() +
  theme(
    legend.position = "top",
    legend.title = element_blank()
  )

dt2 <- dt %>% dcast(factor+pathway~test)

ggplot(dt2, aes(x=parametric, y=parametric.adj)) +
  geom_point(size=0.5) +
  geom_abline(slope=1, intercept=0, color="orange") +
  facet_wrap(~factor, scales="free_y", nrow=1) +
  labs(x="Parametric p-value", y="Adjusted parametric p-value") +
  theme_bw() +
  theme(
    legend.position = "top"
  )

save.image("../rdata/mofa_enrichment_first_comp1_neg.RData")
load("../rdata/mofa_enrichment_first_comp1_neg.RData")

save.image("../rdata/mofa_enrichment_first_comp1_pos.RData")  
load("../rdata/mofa_enrichment_first_comp1_pos.RData")
