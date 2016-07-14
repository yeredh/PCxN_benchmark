# 05/14/2016
#
# Function to estimate the variance of the weighted average 
# of the experiment-level correlation coefficients
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


# function to estimate the variance for the
# HS correlation estimator
GetHSVar = function(r_mat,r_bar,n_vec,verbose=T){
  # Args
  #
  # r_mat: matrix with the experiment level correlations, to columns
  #        correspond to experiments and the rows to pathway pairs
  # r_bar: vector with the weighted average for all pathway pairs
  # n_vec: vector with the sample sizes for all experiments
  #
  # Returns
  #
  # list with (s2_r = s2_p + s2_e)
  # - s2_r: variance of the sample correlations
  # - s2_e: variance due to sampling error
  # - s2_p: variance in the population of correlations
  
  # variance of sample correlations
  if(verbose){cat("Estimating variance of sample correlations .\n")}
  s2_r = (sweep(r_mat,1,r_bar))^2 %*% n_vec/sum(n_vec)
  if(verbose){cat("Done 1/3 .\n")}
  # sampling error variance
  if(verbose){cat("Estimating variance of sampling error .\n")}
  s2_e = 1/(mean(n_vec)-1)*(1-r_bar)^2
  if(verbose){cat("Done 2/3 .\n")}
  # variance in population correlations
  if(verbose){cat("Estimating variance of correlation population .\n")}
  s2_p = apply(s2_r - s2_e,1,function(x){max(x,0)})
  if(verbose){cat("Done 3/3 .\n")}
  
  res=list(s2_r=c(s2_r),s2_e=c(s2_e),s2_p=c(s2_p))
  return(res)
}


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
      source(paste0(pcxn_dir,"src/pcxn_benchmark_var01a.R"))
      # Random Gene Sets
      source(paste0(pcxn_dir,"src/pcxn_benchmark_var01b.R"))
    }
  }
}
cat("DONE \n\n")
