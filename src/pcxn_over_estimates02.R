# 03/13/2016
#
#
# Get experiment-level estimates for the Overlap cases
# of the random gene sets.
# The experiments are identified by its GSE id,
# and we get the estimates for every split with a predetermined overlap
# of the random gene sets.
#
# The job array index is used to select the GSE series with the script,
# and another command line argument to select the overlap case 
# 
# The SLURM script, has two command line arguments passed to the R script.
# The first argument is used to select the overlap case [1-10] and the second argument
# is the job array index which is used to select the random gene set pair [1-1000]
# The SLURM script is submitted using a shell script to select all overlap cases
#
# RUN ON ODYSSEY
#
# (pcxn_over_estimates02a.sh) [1-1000]
# (pcxn_over_estimates02a_shell.sh)

rm(list=ls())
library(parallel)
library(matrixStats)
options(stringsAsFactors = F)
pcxn_dir = "/net/hsphfs1/srv/export/hsphfs1/share_root/hide_lab/PCxN/"

# ==== PCxN Functions ====
OverlapCoefficient <- function(x,y){
  # function to calculate the overlap coefficient between x and y
  # which is defined as the size of the intersection divided by the
  # size of the smaller set
  #
  # Args
  #   x: a vector
  #   y: a vector
  #
  # Returns
  #   the overlap coefficient, a number between 0 and 1
  
  length(intersect(x,y))/min(length(unique(x)),length(unique(y)))
}

GetSummary = function(dat,gs,sum_fun){
  # function to calculate the summary statistic for the pathway
  #
  # Args.
  #   dat: genes by samples matrix
  #   gs: vector with the names of the genes in the gene set
  #   sum_fun: function to calculate the summary
  #
  # Returns
  #   a 1 by samples vector with the summary statistic for the pathway
  
  if(length(gs) > 1){
    # calculate summary for pathways with more than 1 element
    return(sum_fun(dat[rownames(dat) %in% gs,]))
  }else{
    # return actual value for pathways with a single element
    return(dat[rownames(dat) %in% gs,])
  }
}

GetCorEstimates = function(x,y,METHOD = "spearman"){
  # function to estimate the correlation coefficient between x and y,
  # the corresponding t-statistic and p-value
  #
  # Args
  #   x: a vector with n observations
  #   y: a vector with n observations
  #   method: character to pick either the Pearson or Spearman correlation coefficient
  #
  # Returns
  #   a named vector with the correlation estimate, the sample size n, the t-statistic
  #   and its corresponding p-value
  
  # function to estimate t-statistic
  GetStatistic = function(r,n){r*sqrt((n-2)/(1-r^2))}
  
  # get sample size
  if(length(x) == length(y)){
    n <- length(x)
  }else{
    cat("x and y have different lengths! >=( \n")
    return(NA)
  }
  
  # get correlation estimate
  estimate <- cor(x,y,method = METHOD)
  # get t-statistic
  statistic <- GetStatistic(estimate,n)
  # get p-value for the two-sided test
  p.value <- 2*pt(-abs(statistic),n-2)
  
  res <- c(estimate,n,statistic,p.value)
  names(res) <- c("estimate","n","statistic","p.value")
  return(res)
}


# ==== GSE annotation ====
gse_annot = readRDS(paste0(pcxn_dir,"data/GSE_annotation.RDS"))
# get sample size per GSE series
gse_count = table(gse_annot$GSE)
gse_count = sort(gse_count,decreasing=T)
# keep series with at least 5 samples
gse_ids = names(gse_count[gse_count >= 15])

# ==== Expression Background ====
exprs_rnk = readRDS(paste0(pcxn_dir,"data/GPL570_R_mat.RDS"))

# read command line arguments
args = as.numeric(commandArgs(trailingOnly = T))
# the first argument is the overlap case [1-10]
# the second argument is the random gene set [1-1000]

# ==== Random Overlap Cases ====
cat("Overlap Case:",args[1],"\n")
gs_lst = readRDS(paste0(pcxn_dir,"output/Benchmark/random_over",args[1],"_gs.RDS"))


# argument to pick random gene set 
gs = gs_lst[[args[2]]]

# subset expression ranks
exprs_rnk = exprs_rnk[rownames(exprs_rnk) %in% unique(unlist(gs)),]

# R garbage collection
gc()

# pathway gene sets
gsA = gs$gs01
gsB = gs$gs02

# split into three disjoint sets:
# 1. genes shared by pathwatys A and B
gsAB = intersect(gsA,gsB)
# 2. genes unique to pathway A
gsAu = setdiff(gsA,gsB)
# 3. genes unique to pathway B
gsBu = setdiff(gsB,gsA)

# overlap coefficient
over_coef = OverlapCoefficient(gsA,gsB)

# number of cores
(nc = detectCores())

# loop through all experiments
number_of_experiments = length(gse_ids)
input = 1:number_of_experiments

# Function to process elements
ProcessElement = function(ic,METHOD,sum_fun){
  # select GSE series
  gse = gse_ids[ic]
  gsm_targets = gse_annot$GSM[gse_annot$GSE == gse]
  gsm_ind = colnames(exprs_rnk) %in% gsm_targets
  
  # get pathway summaries for unique genes in each pathway
  summaryAu = GetSummary(dat=exprs_rnk[,gsm_ind],gs=gsAu,sum_fun=sum_fun)
  summaryBu = GetSummary(dat=exprs_rnk[,gsm_ind],gs=gsBu,sum_fun=sum_fun)
  
  # get correlation between the summaries for the unique genes
  tmp = data.frame(Pathway.A=names(gs_lst)[args],Pathway.B=names(gs_lst)[args])
  tmp = c(tmp,GetCorEstimates(x=summaryAu,y=summaryBu,METHOD = METHOD))
  # get overlap coefficient
  tmp$Overlap.Coeff= over_coef
  
  setTxtProgressBar(pb,ic)
  return(tmp)
}


# ==== PCxN Estimates: Mean/Pearson ====
ProcessElementMnPn = function(ic){ProcessElement(ic,METHOD="pearson",sum_fun=colMeans)}

pb = txtProgressBar(min=0,max=number_of_experiments,style=3,initial=0)
cat("\n")
resMnPn = mclapply(input,ProcessElementMnPn,mc.cores=nc)
close(pb)
names(resMnPn) = gse_ids
# save results
saveRDS(resMnPn,paste0(pcxn_dir,"output/Benchmark/Overlap/Mean/Pearson/Case",args[1],"/GSE/RandomGS",args[2],"_over",args[1],".RDS"))
rm(resMnPn)


# ==== PCxN Estimates: Mean/Spearman ====
ProcessElementMnSp = function(ic){ProcessElement(ic,METHOD="spearman",sum_fun=colMeans)}

pb = txtProgressBar(min=0,max=number_of_experiments,style=3,initial=0)
cat("\n")
resMnSp = mclapply(input,ProcessElementMnSp,mc.cores=nc)
close(pb)
names(resMnSp) = gse_ids
# save results
saveRDS(resMnSp,paste0(pcxn_dir,"output/Benchmark/Overlap/Mean/Spearman/Case",args[1],"/GSE/RandomGS",args[2],"_over",args[1],".RDS"))
rm(resMnSp)

# ==== PCxN Estimates: Median/Pearson ====
ProcessElementMdPn = function(ic){ProcessElement(ic,METHOD="pearson",sum_fun=colMedians)}


pb = txtProgressBar(min=0,max=number_of_experiments,style=3,initial=0)
cat("\n")
resMdPn = mclapply(input,ProcessElementMdPn,mc.cores=nc)
close(pb)
names(resMdPn) = gse_ids

# save results
saveRDS(resMdPn,paste0(pcxn_dir,"output/Benchmark/Overlap/Median/Pearson/Case",args[1],"/GSE/RandomGS",args[2],"_over",args[1],".RDS"))
rm(resMdPn)

# ==== PCxN Estimates: Median/Spearman ====
ProcessElementMdSp = function(ic){ProcessElement(ic,METHOD="spearman",sum_fun=colMedians)}

pb = txtProgressBar(min=0,max=number_of_experiments,style=3,initial=0)
cat("\n")
resMdSp = mclapply(input,ProcessElementMdSp,mc.cores=nc)
close(pb)
names(resMdSp) = gse_ids
# save results
saveRDS(resMdSp,paste0(pcxn_dir,"output/Benchmark/Overlap/Median/Spearman/Case",args[1],"/GSE/RandomGS",args[2],"_over",args[1],".RDS"))

rm(list=ls())
