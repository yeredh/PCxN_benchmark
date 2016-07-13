# 03/23/2016
#
# Update gene set definitions from MSigDB to v5.1
# - Downloaded GMT files for three gene set collections
#   (Canonical pathways, GO: Biological Process, and Hallmark gene sets)
# - Read GMT files using GSEABase
# - Saved annotation as a list of gene ids for each gene set in an RDS object
#

library(GSEABase)

# ==== Canonical Pathways =====
# read gene set definitions from GMT file
cp_gsc = getGmt("Hide Lab/PCxN/data/Gene Sets/c2.cp.v5.1.entrez.gmt",
                geneIdType=EntrezIdentifier(),
                collectionType=BroadCollection(category = "c2",subCategory = "CP"))

# get entrez gene ids
cp_gs_lst = lapply(cp_gsc,geneIds)
names(cp_gs_lst) = names(cp_gsc)
# save as genes sets as list with gene ids
saveRDS(cp_gs_lst,"Hide Lab/PCxN/data/Gene Sets/cp_gs_v5.1.RDS")

# ==== GO: Biological Process =====
# read gene set definitions from GMT file
gobp_gsc = getGmt("Hide Lab/PCxN/data/Gene Sets/c5.bp.v5.1.entrez.gmt",
                geneIdType=EntrezIdentifier(),
                collectionType=BroadCollection(category = "c5",subCategory = "BP"))

# get entrez gene ids
gobp_gs_lst = lapply(gobp_gsc,geneIds)
names(gobp_gs_lst) = names(gobp_gsc)
# save as genes sets as list with gene ids
saveRDS(gobp_gs_lst,"Hide Lab/PCxN/data/Gene Sets/gobp_gs_v5.1.RDS")

# ===== Hallmark Gene Sets =====
# read gene set definitions from GMT file
h_gsc = getGmt("Hide Lab/PCxN/data/Gene Sets/h.all.v5.1.entrez.gmt",
               geneIdType=EntrezIdentifier(),
               collectionType=BroadCollection(category = "h"))

# get entrez gene ids
h_gs_lst = lapply(h_gsc,geneIds)
names(h_gs_lst) = names(h_gsc)
# save as genes sets as list with gene ids
saveRDS(h_gs_lst,"Hide Lab/PCxN/data/Gene Sets/h_gs_v5.1.RDS")

