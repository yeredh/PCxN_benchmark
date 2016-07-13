# 03/08/2016
#
# Split the Ribosome pathway into two non overlapping gene sets
# by resampling. These gene sets will be the no overlap
# cases for the PCxN benchmark.
#
# ODYSSEY INTERACTIVE SESSION
# 
# source new-modules.sh
# module load R/3.2.2-fasrc01
# export R_LIBS_USER=$HOME/apps/R/3.2.2-fasrc01:$R_LIBS_USE
# srun -p irizarry,serial_requeue --mem-per-cpu=8000 -t 0-0:30 --pty R


rm(list=ls())
options(stringsAsFactors = F)
pcxn_dir = "/net/hsphfs1/srv/export/hsphfs1/share_root/hide_lab/PCxN/"
# ===== KEGG Human Pathways =====
# downloaded using KEGGREST om 9/2014)
kegg_gs = readRDS(paste0(pcxn_dir,"data/KEGG_hsa.RDS"))


# ==== Microarray Annotation ====
# genes in microarray background
entrez_ids = readRDS(paste0(pcxn_dir,"data/GPL570_genes.RDS"))

# KEGG pathways
kegg_gs$gs = lapply(kegg_gs$gs,function(x){x[x %in% entrez_ids]})

# Get Ribosome pathway
ribosome_gs = unlist(kegg_gs$gs[grep("Ribosome - Homo sapiens (human)",kegg_gs$pathway_names,fixed=T)])

# number of non overlapping sets
iter=1000 

# split the ribosome pathway into two non overlapping gene sets 
noover_gs = vector("list",iter)
set.seed(13)
for(k in 1:iter){
  noover_gs[[k]]$gs01 = sample(ribosome_gs,length(ribosome_gs)/2)
  names(noover_gs[[k]]$gs01) = NULL
  noover_gs[[k]]$gs02 = setdiff(ribosome_gs,noover_gs[[k]]$gs01)
}

names(noover_gs) = paste0("Ribosome_",1:iter)
# save non overlapping cases
saveRDS(noover_gs,paste0(pcxn_dir,"output/Benchmark/ribosome_non_over_gs.RDS"))

rm(list=ls())
