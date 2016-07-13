# 03/23/2016
#
# Download annotation for KEGG pathways using KEGGREST

# source("https://bioconductor.org/biocLite.R")
# biocLite("KEGGREST")

library(KEGGREST)
# ===== KEGG Human Pathways =====
# get all pathway ids
pathways_list = keggList("pathway", "hsa")
# make pathway ids into KEGG-style human pathway identifiers
pathway_codes <- sub("path:", "", names(pathways_list))
# get gene ids for each pathway
# subsetting by c(TRUE, FALSE) -- which repeats
# as many times as needed, sorts through some
# unexpected packaging of geneIDs in the GENE element
# of each pw[[n]]
genes_by_pathway = sapply(pathway_codes,
                          function(pwid){
                            pw = keggGet(pwid)
                            pw[[1]]$GENE[c(TRUE, FALSE)]
                          })

# save constituent genes of each pathway, their KEGG ids and names
names(pathways_list) = pathway_codes
KEGG_hsa = list(gs=genes_by_pathway,pathway_names=pathways_list)
saveRDS(KEGG_hsa,"Hide Lab/PCxN/data/Gene Sets/KEGG_hsa.RDS")

