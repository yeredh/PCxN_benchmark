# 05/14/2016
#
# Function to estimate the jackknife statistics of the weighted average 
# of the experiment-level correlation coefficients
# - Bias
# - Variance
#
# ODYSSEY (INTERACTIVE)
#
# source new-modules.sh
# module load R/3.2.2-fasrc01
# export R_LIBS_USER=$HOME/apps/R/3.2.2-fasrc01:$R_LIBS_USE
# srun -p irizarry,serial_requeue --mem-per-cpu=12000 -t 0-3:30 --pty R


rm(list=ls())
options(stringsAsFactors = F)
pcxn_dir = "/net/hsphfs1/srv/export/hsphfs1/share_root/hide_lab/PCxN/"

# adjust p-values and correlation estimates,
# otherwise the functions to combine the p-values
# cannot handle values near 0 and values near 1
AdjustRmat = function(r_mat,eps=1E-16){
  res = t(apply(r_mat,1,function(r){
    r[r == 0] = eps
    return(r)
  }))
  return(res)
}


library(genefilter)
# ==== GSE annotation ====
gse_annot = readRDS(paste0(pcxn_dir,"data/GSE_annotation.RDS"))
# get sample size per GSE series
gse_count = table(gse_annot$GSE)
gse_count = sort(gse_count,decreasing=T)
# keep series with at least 5 samples
gse_ids = names(gse_count[gse_count >= 15])

# ==== Benchmark Settings ====
el1 = c("NoOverlap","Overlap")
el2 = c("Mean","Median")
el3 = c("Pearson","Spearman")

# nested loops to get all benchmark settings
for(io in 1:2){
  for(jo in 1:2){
    for(ko in 1:2){
      # benchmark setting
      cat("Benchmark:",paste(el1[io],el2[jo],el3[ko],sep="/"),"\n")
      # Ribosome Gene Sets
      source(paste0(pcxn_dir,"src/pcxn_benchmark_jack01a.R"))
      # Random Gene Sets
      source(paste0(pcxn_dir,"src/pcxn_benchmark_jack01b.R"))
    }
  }
}
cat("DONE \n\n")
