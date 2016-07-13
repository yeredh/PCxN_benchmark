# 04/11/2016
#
# Split the Ribosome pathway into two gene sets with a 
# predetermined number of shared genes by resampling. 
# These gene sets will be the overlap
# cases for the PCxN benchmark.

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
# downloaded using KEGGREST on 9/2014
# kegg_gs = readRDS("Hide Lab/PCxN/data/Gene Sets/KEGG_hsa.RDS")
kegg_gs = readRDS(paste0(pcxn_dir,"data/KEGG_hsa.RDS"))

# genes in microarray background
# entrez_ids = readRDS("Hide Lab/PCxN/data/GPL570_genes.RDS")
entrez_ids = readRDS(paste0(pcxn_dir,"data/GPL570_genes.RDS"))

# filter by genes present in array
# KEGG pathways
kegg_gs$gs = lapply(kegg_gs$gs,function(x){x[x %in% entrez_ids]})

# Get Ribosome pathway
ribosome_gs = kegg_gs$gs[[grep("Ribosome - Homo sapiens (human)",kegg_gs$pathway_names,fixed=T)]]

# ===== Overlap Cases =====
# overlap cases
# gs_lst = readRDS("Hide Lab/PCxN/data/pcxn_overlap_cases.RDS")
gs_lst = readRDS(paste0(pcxn_dir,"data/pcxn_overlap_cases.RDS"))

# number of overlapping sets
iter=1000 
set.seed(13)
# loop thru overlap cases
for(ic in 1:length(gs_lst)){
  # list to store entrez gene ids for Ribosome pathway
  over_gs =  vector("list",iter)
  for(k in 1:iter){
    # permute genes in the Ribosome pathway
    gs_rnd = sample(ribosome_gs)
    # split into two gene sets with predefined overlap
    over_gs[[k]]$gs01 = gs_rnd[gs_lst[[ic]]$s1]
    over_gs[[k]]$gs02 = gs_rnd[gs_lst[[ic]]$s2]
  }
  names(over_gs) =  paste0("Ribosome_over",ic,"_",1:iter)
  
  # save overlapping cases
  saveRDS(over_gs,paste0(pcxn_dir,"output/Benchmark/ribosome_over",ic,"_gs.RDS"))
}
