# 05/15/2016
#
# Script to rank the following settings of PCxN
# - Summary
#   > Mean
#   > Median
# - Correlation
#   > Pearson
#   > Spearman
# For the No Overlap case based on the median variance
# of the weighted average of correlation coefficients




rm(list=ls())
# ==== Functions =====
# helper function to translate R's line positions
# stolen from: 
# http://stackoverflow.com/questions/29125019/get-margin-line-locations-mgp-in-user-coordinates
line2user <- function(line, side) {
  lh <- par('cin')[2] * par('cex') * par('lheight')
  x_off <- diff(grconvertX(0:1, 'inches', 'user'))
  y_off <- diff(grconvertY(0:1, 'inches', 'user'))
  switch(side,
         `1` = par('usr')[3] - line * y_off * lh,
         `2` = par('usr')[1] - line * x_off * lh,
         `3` = par('usr')[4] + line * y_off * lh,
         `4` = par('usr')[2] + line * x_off * lh,
         stop("side must be 1, 2, 3, or 4", call.=FALSE))
}

# ==== Ribosome (No overlap) =====
# Mean/Pearson
ribMnPn = readRDS("Hide Lab/PCxN/output/Benchmark/RibosomeNoOverMnPnShrkVar.RDS")
# Mean/Spearman
ribMnSp = readRDS("Hide Lab/PCxN/output/Benchmark/RibosomeNoOverMnSpShrkVar.RDS")
# Median/Pearson
ribMdPn = readRDS("Hide Lab/PCxN/output/Benchmark/RibosomeNoOverMdPnShrkVar.RDS")
# Median/Spearman
ribMdSp = readRDS("Hide Lab/PCxN/output/Benchmark/RibosomeNoOverMdSpShrkVar.RDS")

# ==== Random (No overlap) ====
# Mean/Pearson
rndMnPn = readRDS("Hide Lab/PCxN/output/Benchmark/RandomGSNoOverMnPnShrkVar.RDS")
# Mean/Spearman
rndMnSp = readRDS("Hide Lab/PCxN/output/Benchmark/RandomGSNoOverMnSpShrkVar.RDS")
# Median/Pearson
rndMdPn = readRDS("Hide Lab/PCxN/output/Benchmark/RandomGSNoOverMdPnShrkVar.RDS")
# Median/Spearman
rndMdSp = readRDS("Hide Lab/PCxN/output/Benchmark/RandomGSNoOverMdSpShrkVar.RDS")


library(rafalib)
# ==== Boxplots Correlation (No overlap) =====
mypar(1,1)
png("Hide Lab/PCxN/doc/Benchmark/figures/rib_v_rnd_nonover_shrk_var_boxplots01.png",
    width=800,height=450,pointsize = 15)
boxplot(
  list(
    ribMnPn$s2_r,rndMnPn$s2_r, # Mean/Pearson
    ribMnSp$s2_r,rndMnSp$s2_r, # Mean/Spearman
    ribMdPn$s2_r,rndMdPn$s2_r, # Median/Pearson
    ribMdSp$s2_r,rndMdSp$s2_r  # Median/Spearman
  ),
  #main="Correlation Estimates (No Overlap)",
  names=rep(c("Ribosome","Random"),times=4)
)
title("Variance Estimates (No Overlap)",line=2)
text(x=1.5,y=line2user(line=0.5,side=3),"Mean/Pearson",xpd=T)
text(x=3.5,y=line2user(line=0.5,side=3),"Mean/Spearman",xpd=T)
text(x=5.5,y=line2user(line=0.5,side=3),"Median/Pearson",xpd=T)
text(x=7.5,y=line2user(line=0.5,side=3),"Median/Spearman",xpd=T)
dev.off()

mypar(1,1)
png("Hide Lab/PCxN/doc/Benchmark/figures/rib_v_rnd_nonover_shrk_var_boxplots02.png",
    width=800,height=450,pointsize = 15)
boxplot(list(
  c(ribMnPn$s2_r,rndMnPn$s2_r),
  c(ribMnSp$s2_r,rndMnSp$s2_r),
  c(ribMdPn$s2_r,rndMdPn$s2_r),
  c(ribMdSp$s2_r,rndMdSp$s2_r)
),
names=c("Mean/Pearson","Mean/Spearman","Median/Pearson","Median/Spearman")
)
title("Variance Estimates (No Overlap)",line=2)
dev.off()



# === Median Variance ====
s2r_vec = rep(NA,4)
names(s2r_vec) = c("Mean/Pearson","Mean/Spearman","Median/Pearson","Median/Spearman")


s2r_vec[1] = median(c(ribMnPn$s2_r,rndMnPn$s2_r))
s2r_vec[2] = median(c(ribMnSp$s2_r,rndMnSp$s2_r))
s2r_vec[3] = median(c(ribMdPn$s2_r,rndMdPn$s2_r))
s2r_vec[4] = median(c(ribMdSp$s2_r,rndMdSp$s2_r))

# rankings by median variance (small is better)
rank(s2r_vec*-1)


# data frame to store the method rankings

el1 = c("Mean","Median")
el2 = c("Pearson","Spearman")
el3 = c("Liptak","Logit")

scoreBoard = matrix(NA,ncol=8,nrow=1)
colnames(scoreBoard) = paste(rep(rep(el1,each=2),2),rep(rep(el2,2),2),rep(el3,each=4),sep="/")
rownames(scoreBoard) = c("Var")
scoreBoard = as.data.frame(scoreBoard)

# ranks based on median variance
scoreBoard[1,] = rank(rep(s2r_vec,2)*-1)

saveRDS(scoreBoard,"Hide Lab/PCxN/output/pcxn_benchmark_ranks_nover_var.RDS")
