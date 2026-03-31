rm(list=ls())

library(ggplot2)
library(data.table)
library(cowplot)
library(purrr)
library(MOFA2)
library(dplyr)

setwd("D:/Analysis/MOFA/input/")
getwd()

load("../rdata/mofa_last_comp1_add_study.RData")

model

#write.csv(factors$group_1, "../output/factor_value_normal_last_comp1.csv")
#write.csv(factors$group_2, "../output/factor_value_target_last_comp1.csv")
#write.csv(factors$group_3, "../output/factor_value_critical_last_comp1.csv")

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
})    

duplicated_gene <- duplicated(progene_map$new_Gene)
progene_map$new_Gene[duplicated_gene]  

features_names(model)$proteome <- progene_map$new_Gene
features_names(model)$proteome

## enrichment analysis
library(msigdbr) 
library(BiocParallel)
library(tibble)

# Update subcollection and names due to recent updates. 
#MsigDB HALLMARK gene sets
hall_gene = msigdbr(species="human", category="H")
hall_list = split(x = hall_gene$gene_symbol, f = hall_gene$gs_name)
hall_list %>% head() %>% lapply(head)

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

prot_genes <- features_names(model)$proteome
all_gene_names <- names(all_gene_list)

binary_matrix <- sapply(all_gene_list, function(finder_genes) {
  as.integer(prot_genes %in% finder_genes)
})

rownames(binary_matrix) <- prot_genes
colnames(binary_matrix) <- all_gene_names
binary_matrix
binary_matrix2 <- t(binary_matrix)

##

## 

enrichment.parametric <- run_enrichment(model,
                                        view = "proteome", factors = 1:11,
                                        feature.sets = binary_matrix2,
                                        sign = "negative",
                                        statistical.test = "parametric"
)

enrichment.parametric.adj <- run_enrichment(model,
                                            view = "proteome", factors = 1:11,
                                            feature.sets = binary_matrix2,
                                            sign = "negative",
                                            statistical.test = "cor.adj.parametric"
)

# positive
enrichment.parametric <- run_enrichment(model,
                                        view = "proteome", factors = 1:11,
                                        feature.sets = binary_matrix2,
                                        sign = "positive",
                                        statistical.test = "parametric"
)

names(enrichment.parametric)
enrichment.parametric$set.statistics[1:5, 1]
enrichment.parametric$pval.adj[1:5, 1]

pdf("../output/enrichment_heatmap_factor_last_comp1.pdf")   
plot_enrichment_heatmap(enrichment.parametric)
dev.off()

#enrichment_top10-pathways_factor1_last_comp1

plot_enrichment(enrichment.parametric, 
                factor=1, max.pathways=10)

plot_enrichment(enrichment.parametric, 
                factor=2, max.pathways=10)

plot_enrichment(enrichment.parametric, 
                factor=3, max.pathways=10)

plot_enrichment(enrichment.parametric, 
                factor=4, max.pathways=10)

plot_enrichment(enrichment.parametric, 
                factor=5, max.pathways=10)

plot_enrichment(enrichment.parametric, 
                factor=6, max.pathways=10)

plot_enrichment(enrichment.parametric, 
                factor=7, max.pathways=10)

plot_enrichment(enrichment.parametric, 
                factor=8, max.pathways=10)

plot_enrichment(enrichment.parametric, 
                factor=9, max.pathways=10)

plot_enrichment(enrichment.parametric, 
                factor=10, max.pathways=10)

plot_enrichment(enrichment.parametric, 
                factor=11, max.pathways=10)

# save figure
setwd("../output/")

p1 <- plot_enrichment(enrichment.parametric, factor=6, max.pathways=10)

ggsave("enrichment_top10_pathway_factor6_neg_high-res_last_comp1.png",plot=p1,
       width=8.4, height=5.3, dpi=300)

#ggsave("enrichment_top10_pathway_factor6_pos_high-res_last_comp1.png",plot=p2,
#       width=8.4, height=5.3, dpi=300)

#enrichment_weight_factor1_last_comp2
plot_enrichment_detailed(enrichment.parametric,
                         factor=1, max.genes=5, max.pathways=7)

plot_enrichment_detailed(enrichment.parametric,
                         factor=2, max.genes=5, max.pathways=7)

plot_enrichment_detailed(enrichment.parametric,
                         factor=3, max.genes=5, max.pathways=7)

plot_enrichment_detailed(enrichment.parametric,
                         factor=4, max.genes=5, max.pathways=7)

plot_enrichment_detailed(enrichment.parametric,
                         factor=5, max.genes=5, max.pathways=7)

plot_enrichment_detailed(enrichment.parametric,
                         factor=6, max.genes=5, max.pathways=10)

plot_enrichment_detailed(enrichment.parametric,
                         factor=8, max.genes=5, max.pathways=7)

plot_enrichment_detailed(enrichment.parametric,
                         factor=9, max.genes=5, max.pathways=7)

plot_enrichment_detailed(enrichment.parametric,
                         factor=10, max.genes=5, max.pathways=7)

plot_enrichment_detailed(enrichment.parametric,
                         factor=11, max.genes=5, max.pathways=7)


h <- as.data.frame(enrichment.parametric$feature.sets)
h <- t(h)
write.csv(h, "../output/enrichment_pathway-positive-result_last_comp1.csv")

### comp1 - negative
## factor 9
h_weight <- as.data.frame(enrichment.parametric$feature.statistics)
h_weight$REACTOME_ECM_PROTEOGLYCANS <- h[,"REACTOME_ECM_PROTEOGLYCANS"]

write.csv(h_weight, "../output/enrichment_factor9_weight-result_last_comp1.csv")

## factor 6
h_weight2 <- as.data.frame(enrichment.parametric$feature.statistics)
h_weight2$GOBP_POSITIVE_REGULATION_OF_LIPID_METABOLIC_PROCESS <- h[,"GOBP_POSITIVE_REGULATION_OF_LIPID_METABOLIC_PROCESS"]
h_weight2$GOBP_REGULATION_OF_LIPID_METABOLIC_PROCESS <- h[,"GOBP_REGULATION_OF_LIPID_METABOLIC_PROCESS"]
h_weight2$GOBP_PROTEIN_CONTAINING_COMPLEX_REMODELING <- h[,"GOBP_PROTEIN_CONTAINING_COMPLEX_REMODELING"]
h_weight2$GOBP_REGULATION_OF_PLASMA_LIPOPROTEIN_PARTICLE_LEVELS <- h[,"GOBP_REGULATION_OF_PLASMA_LIPOPROTEIN_PARTICLE_LEVELS"]
h_weight2$GOBP_PROTEIN_LIPID_COMPLEX_SUBUNIT_ORGANIZATION <- h[,"GOBP_PROTEIN_LIPID_COMPLEX_SUBUNIT_ORGANIZATION"]
h_weight2$GOBP_LIPID_MODIFICATION <- h[,"GOBP_LIPID_MODIFICATION"]
h_weight2$GOBP_PHOSPHATIDYLCHOLINE_METABOLIC_PROCESS <- h[,"GOBP_PHOSPHATIDYLCHOLINE_METABOLIC_PROCESS"]
h_weight2$GOBP_LIPID_LOCALIZATION <- h[,"GOBP_LIPID_LOCALIZATION"]
h_weight2$GOBP_CHOLESTEROL_EFFLUX <- h[,"GOBP_CHOLESTEROL_EFFLUX"]
h_weight2$GOBP_ORGANIC_HYDROXY_COMPUND_TRANSPORT <- h[,"GOBP_ORGANIC_HYDROXY_COMPOUND_TRANSPORT"]

## factor 2
h_weight3 <- as.data.frame(enrichment.parametric$feature.statistics)
h_weight3$GOBP_COMPLEMENT_ACTIVATION_ALTERNATIVE_PATHWAY <- h[,"GOBP_COMPLEMENT_ACTIVATION_ALTERNATIVE_PATHWAY"]
h_weight3$KEGG_SYSTEMIC_LUPUS_ERYTHEMATOSUS <- h[,"KEGG_SYSTEMIC_LUPUS_ERYTHEMATOSUS"]
h_weight3$KEGG_PRION_DISEASES <- h[,"KEGG_PRION_DISEASES"]
h_weight3$GOBP_POSITIVE_REGULATION_OF_IMMUNE_RESPONSE <- h[,"GOBP_POSITIVE_REGULATION_OF_IMMUNE_RESPONSE"]
h_weight3$GOBP_MYELOID_CELL_DIFFERENTIATION <- h[,"GOBP_MYELOID_CELL_DIFFERENTIATION"]
h_weight3$GOBP_IMMUNE_EFFECTOR_PROCESS <- h[,"GOBP_IMMUNE_EFFECTOR_PROCESS"]
h_weight3$GOBP_LEUKOCYTE_MEDIATED_IMMUNITY <- h[,"GOBP_LEUKOCYTE_MEDIATED_IMMUNITY"]

write.csv(h_weight2, "../output/enrichment_factor6_negative_weight-result_last_comp1.csv")
write.csv(h_weight3, "../output/enrichment_factor2_weight-result_last_comp1.csv")

## comp1 - positive
# factor 6
h_weight <- as.data.frame(enrichment.parametric$feature.statistics)
h_weight$REACTOME_INTEGRIN_CELL_SURFACE_INTERACTIONS <- h[,"REACTOME_INTEGRIN_CELL_SURFACE_INTERACTIONS"]
h_weight$REACTOME_EXTRACELLULAR_MATRIX_ORGANIZATION <- h[,"REACTOME_EXTRACELLULAR_MATRIX_ORGANIZATION"]
h_weight$GOBP_HETEROTYPIC_CELL_CELL_ADHESION <- h[,"GOBP_HETEROTYPIC_CELL_CELL_ADHESION"]
h_weight$GOBP_EPITHELIAL_CELL_APOPTOTIC_PROCESS <- h[,"GOBP_EPITHELIAL_CELL_APOPTOTIC_PROCESS"]
h_weight$REACTOME_MAP2K_AND_MAPK_ACTIVATION <- h[,"REACTOME_MAP2K_AND_MAPK_ACTIVATION"]
h_weight$REACTOME_SIGNALING_BY_MODERATE_KINASE_ACTIVITY_BRAF_MUTANTS <- h[,"REACTOME_SIGNALING_BY_MODERATE_KINASE_ACTIVITY_BRAF_MUTANTS"]
h_weight$REACTOME_TOLL_LIKE_RECEPTOR_TLR1_TLR2_CASCADE <- h[,"REACTOME_TOLL_LIKE_RECEPTOR_TLR1_TLR2_CASCADE"]
h_weight$GOBP_PROTEIN_ACTIVATION_CASCADE <- h[,"GOBP_PROTEIN_ACTIVATION_CASCADE"]
h_weight$REACTOME_ONCOGENIC_MAPK_SIGNALING <- h[,"REACTOME_ONCOGENIC_MAPK_SIGNALING"]
h_weight$REACTOME_INNATE_IMMUNE_SYSTEM <- h[,"REACTOME_INNATE_IMMUNE_SYSTEM"]

write.csv(h_weight, "../output/enrichment_factor6_positive_weight-result_last_comp1.csv")

## last-only
genes <- list("LUM", "NCAM1", "COMP", "TNXB", "HSPG2")  # platelet_activation_signaling_and_aggregation; TOP-4 Weights
genes %>% map(~ plot_factors(model, 
                             factors = 9, 
                             color_by = 'sex',       #'group', 
                             shape_by = 'condition',
                             scale = T,
                             legend = T
)) %>% cowplot::plot_grid(plotlist=., nrow=2)

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

save.image("../rdata/mofa_enrichment_last_comp1_neg.RData")
load("../rdata/mofa_enrichment_last_comp1_neg.RData")

save.image("../rdata/mofa_enrichment_last_comp1_pos.RData")
load("../rdata/mofa_enrichment_last_comp1_pos.RData")