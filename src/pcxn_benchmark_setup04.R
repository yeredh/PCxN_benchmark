# 03/08/2016
#
# Sample genes from background to create a set of the same
# size as the Ribosome pathway. Then, split this random set
# into two non overlapping gene sets
# These gene sets will be the random gene sets in the no overlap
# cases for the PCxN benchmark.
#
# ODYSSEY INTERACTIVE SESSION
# 
# source new-modules.sh
# module load R/3.2.2-fasrc01
# export R_LIBS_USER=$HOME/apps/R/3.2.2-fasrc01:$R_LIBS_USE
# srun -p irizarry,serial_requeue --mem-per-cpu=4000 -t 0-0:30 --pty R


rm(list=ls())
options(stringsAsFactors = F)
pcxn_dir = "/net/hsphfs1/srv/export/hsphfs1/share_root/hide_lab/PCxN/"
# ===== KEGG Human Pathways =====
# downloaded using KEGGREST om 9/2014)
kegg_gs = readRDS(paste0(pcxn_dir,"data/KEGG_hsa.RDS"))
# kegg_gs = readRDS("Hide Lab/PCxN/data/Gene Sets/KEGG_hsa.RDS")

# ==== Microarray Annotation ====
# genes in microarray background
entrez_ids = readRDS(paste0(pcxn_dir,"data/GPL570_genes.RDS"))
# entrez_ids = readRDS("Hide Lab/PCxN/data/GPL570_genes.RDS")
# KEGG pathways
kegg_gs$gs = lapply(kegg_gs$gs,function(x){x[x %in% entrez_ids]})

# Get Ribosome pathway
ribosome_gs = unlist(kegg_gs$gs[grep("Ribosome - Homo sapiens (human)",kegg_gs$pathway_names,fixed=T)])

# number of non overlapping random sets
iter=1000 

# sample genes from the background (|ribosome pathway|), and
# split them into two non overlapping gene sets 
nonover_rndgs = vector("list",iter)
set.seed(13)
for(k in 1:iter){
  rnd_gs = sample(entrez_ids,length(ribosome_gs))
  nonover_rndgs[[k]]$gs01 = sample(rnd_gs,length(rnd_gs)/2)
  names(nonover_rndgs[[k]]$gs01) = NULL
  nonover_rndgs[[k]]$gs02 = setdiff(rnd_gs,nonover_rndgs[[k]]$gs01)
}

names(nonover_rndgs) = paste0("Random",1:iter)
# save non overlapping cases
saveRDS(nonover_rndgs,paste0(pcxn_dir,"output/Benchmark/random_non_over_gs.RDS"))

rm(list=ls())
