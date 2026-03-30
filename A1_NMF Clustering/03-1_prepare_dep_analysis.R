library(dplyr)

# prepare to dep analysis 
prot_first <- read.csv("./dep/proteome_first_impute_data.csv")
prot_last <- read.csv("./dep/Proteome_last_impute_data.csv")

proteome_info

## in 01_data_preprocessing.R ## 
# Rows with a high frequency of NA values are excluded via imputation.
# match to proteome data (*822 proteins) after imputation
prot_list <- intersect(prot_first$Protein_ID, proteome_info$Accession)
length(prot_list)

rownames(proteome_info) <- proteome_info$Accession
rownames(proteome_info)

rownames(prot_first) <- prot_first$Protein_ID
rownames(prot_last) <- prot_last$Protein_ID

# fixed 822 protein, other 7 protein is non-annotaiton 
prot_info <- proteome_info[prot_list, ]
prot_first <- prot_first[prot_list, ]
prot_last <- prot_last[prot_list, ]

# Incorporated 3 additional protein attributes (ref to. DEP2 library tutorial)
prot_first$First_Protein_ID <- prot_info$First.Protein
prot_first$Gene_name <- prot_info$Gene
prot_first$Description <- prot_info$Description

prot_first <- prot_first %>% relocate(c(First_Protein_ID, Gene_name, Description), .after=Protein_ID)

write.csv(prot_first, "./dep/proteome_first_dep_data.csv")

prot_last$First_Protein_ID <- prot_info$First.Protein
prot_last$Gene_name <- prot_info$Gene
prot_last$Description <- prot_info$Description

prot_last <- prot_last %>% relocate(c(First_Protein_ID, Gene_name, Description), .after=Protein_ID)

write.csv(prot_last, "./dep/proteome_last_dep_data.csv")
