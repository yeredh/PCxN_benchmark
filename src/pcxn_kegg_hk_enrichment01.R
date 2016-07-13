# 04/03/2016
#
# Test enrichment of housekeeping genes in KEGG human pathways 
# The idea is to select a KEGG pathway very significantly enriched 
# for housekeeping genes as a benchmark for PCxN

rm(list=ls())
# ==== Functions ====
SimpleEnrichmentTest <- function(background.gs,gs1,gs2){
  # Function to perform Fisher's exact test to see if the overlap
  # between the genes in class 1 (gs1) and category 1 (gs2) is 
  # significant give the total population of genes (background.gs)
  #
  #             Category 1    not Category 1  Total
  #     Class 1     n11             n12         n1.
  # not Class 1     n21             n22         n2.
  #       Total     n.1             n.2         n
  #
  # Args
  #   background.gs: vector with the gene names in the background as characters
  #   gs1: vector with the gene names in class 1 as characters
  #   gs2: vector with the gene names in category 1 as characters
  #
  # Returns
  # a list with the results of Fisher exact test for the alternatives
  #   - Two-Sided (Enrichment/Depletion)
  #   - Greater (Enrichment)
  #   - Less (Depletion)
  # adjust gene sets to background
  gs1=intersect(gs1,background.gs)
  gs2=intersect(gs2,background.gs)
  
  # get entries for 2x2 contingency table
  # [fixed by experiment]
  n=length(unique(background.gs))
  n1.=length(unique(gs1))
  n.1=length(unique(gs2))
  n11=length(intersect(gs1,gs2))
  # additional entries
  n12=n1.-n11
  n21=n.1-n11
  n22=n-n21-n11-n12
  # matrix for contingency table
  M <- matrix(c(n11,n12,n21,n22),ncol=2,nrow=2,byrow=T)
  # Fisher exact tests
  res <- list(two.sided=fisher.test(M,alternative="two.sided"),
              greater=fisher.test(M,alternative="greater"),
              less=fisher.test(M,alternative="less"))
  return(res)
}


# ===== KEGG Human Pathways =====
# downloaded using KEGGREST om 9/2014)
kegg_gs = readRDS("Hide Lab/PCxN/data/Gene Sets/KEGG_hsa.RDS")

# ==== Housekeeping genes ====
# Taken from Eisenberg et al. (2013)
hk_genes = readRDS("Irizarry/Project I/data/HK.genes.Eisenberg.2013.RDS")


# # ==== Microarray Annotation ====
# library(hgu133plus2.db)
# # get Entrez ids for all genes represented in the human array
# mcr_genes = select(hgu133plus2.db,
#                    keys=keys(hgu133plus2.db,keytype="PROBEID"),
#                    keytype="PROBEID",columns="ENTREZID")
# 
# entrez_ids = unique(mcr_genes$ENTREZID)
# 

# genes in microarray background
entrez_ids = readRDS("Hide Lab/PCxN/data/GPL570_genes.RDS")

# filter by genes present in background
# Housekeeping genes
hk_genes = hk_genes[hk_genes$ENTREZID %in% entrez_ids,]
# KEGG pathways
kegg_gs$gs = lapply(kegg_gs$gs,function(x){x[x %in% entrez_ids]})


# ==== Enrichment test ====
res = vector("list",length(kegg_gs$gs))
names(res) = kegg_gs$pathway_names

pb = txtProgressBar(min=0,max=length(res),initial=0,style=3)
for(k in 1:length(res)){
  res[[k]] = SimpleEnrichmentTest(background.gs=entrez_ids,
                                  gs1=kegg_gs$gs[[k]],
                                  gs2=hk_genes$ENTREZID)
  setTxtProgressBar(pb,k)
}
close(pb)

# adjust p-values for multiple comparison
kegg_pvals = data.frame(Pathway = names(res),
                        two.sided=p.adjust(sapply(res,function(x){x$two.sided$p.value}),"fdr"),
                        greater=p.adjust(sapply(res,function(x){x$greater$p.value}),"fdr"),
                        less=p.adjust(sapply(res,function(x){x$less$p.value}),"fdr"))


kegg_pvals = kegg_pvals[order(kegg_pvals$greater),]

# ==== Table ====
library(xtable)

# Top 15 KEGG pathways most enriched for housekeeping genes
tab =  xtable(head(kegg_pvals[,c(1,3)],15) ,digits=-4)
print(tab,include.rownames = F)


