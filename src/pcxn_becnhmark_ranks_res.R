# 05/15/2016
#
# Script to rank the following settings of PCxN
# - Summary
#   > Mean
#   > Median
# - Correlation
#   > Pearson
#   > Spearman
# - p-value combination method
#   > Liptak
#   > Logit
#
# For the both the no overlap and the overlap cases by combining 
# the ranks based on the following
# - t-statistic: comparing ribosome estimates with random gene sets
#                we expect higher correlation estimates for ribosome gene sets
# the measures below are based on distinguishing true postives and true negatives
# based on the p-values from the ribosome gene sets and the random gene sets
# We assume that all p-values from a ribosome gene sets should be significant
# while the p-values from the random gene sets should not be significant
# - AUC
# - Accuracy
# - F1-Score
# additionally, we estimate the variance of weighted average of the experiment-level 
# correlation estimates using results from meta-analysis and the bias using 
# jackknife statistics
# - Variance
# - Bias

rm(list=ls())

# ==== Rankings ====
# ranks based on correlation estimates and p-values
pcxn_ranks_nover = readRDS("Hide Lab/PCxN/output/pcxn_benchmark_ranks_nover.RDS")
pcxn_ranks_over = readRDS("Hide Lab/PCxN/output/pcxn_benchmark_ranks_over.RDS")
# ranks based on the variance of the weighted average of experiment-level correlations
pcxn_ranks_nover_var = readRDS("Hide Lab/PCxN/output/pcxn_benchmark_ranks_nover_var.RDS")
pcxn_ranks_over_var = readRDS("Hide Lab/PCxN/output/pcxn_benchmark_ranks_over_var.RDS")
# ranks based on the bias of the weighted average of experiment-level correlations
pcxn_ranks_nover_bias = readRDS("Hide Lab/PCxN/output/pcxn_benchmark_ranks_nover_bias.RDS")
pcxn_ranks_over_bias = readRDS("Hide Lab/PCxN/output/pcxn_benchmark_ranks_over_bias.RDS")


res = rbind(
  pcxn_ranks_nover,
  pcxn_ranks_nover_bias,
  pcxn_ranks_nover_var
)

for(k in 1:10){
  res = res + rbind(
    pcxn_ranks_over[[k]],
    pcxn_ranks_over_bias[[k]],
    pcxn_ranks_over_var[[k]])
}

which.max(colSums(res))

