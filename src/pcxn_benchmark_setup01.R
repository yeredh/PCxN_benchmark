# 04/08/2016
#
# Subset and a copy the matrix with the gene expression background to include
# only the genes in the Ribosome pathway in order to save memory while performing
# the benchmark for the Ribosome pathway.
#
# ODYSSEY (INTERACTIVE SESSION)
# source new-modules.sh
# module load R/3.2.2-fasrc01
# export R_LIBS_USER=$HOME/apps/R/3.2.2-fasrc01:$R_LIBS_USE
# srun -p irizarry,serial_requeue --mem-per-cpu=18000 -t 0-0:30 --pty R

rm(list=ls())
options(stringsAsFactors = F)
pcxn_dir = "/net/hsphfs1/srv/export/hsphfs1/share_root/hide_lab/PCxN/"

# ==== GSE annotation ====
gse_annot = readRDS(paste0(pcxn_dir,"data/GSE_annotation.RDS"))
# get sample size per GSE series
gse_count = table(gse_annot$GSE)
gse_count = sort(gse_count,decreasing=T)
# keep series with at least 5 samples
gse_ids = names(gse_count[gse_count >= 15])

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


# ==== Expression Background ====
load(paste0(pcxn_dir,"data/GPL570.R.mat.RData"))


# check if all genes from the annotation are included
if(mean(ribosome_gs %in% rownames(GPL570.R.mat)) == 1){
  cat("All genes included in the expression ranks matrix  =) \n")
}

# subset and save expression rank matrix to save memory 
# (i.e. only keep ranks for genes in the Ribosome pathway)
saveRDS(GPL570.R.mat[rownames(GPL570.R.mat) %in% ribosome_gs,],
        paste0(pcxn_dir,"output/Benchmark/GPL570_Rb_mat.RDS"))

rm(list=ls())


