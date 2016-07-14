# 05/15/2016
#
# Script to rank the following settings of PCxN
# - Summary
#   > Mean
#   > Median
# - Correlation
#   > Pearson
#   > Spearman
#
# For the no overlap case based on the median |bias|
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
ribMnPn = readRDS("Hide Lab/PCxN/output/Benchmark/RibosomeNoOverMnPnShrkJack.RDS")
# Mean/Spearman
ribMnSp = readRDS("Hide Lab/PCxN/output/Benchmark/RibosomeNoOverMnSpShrkJack.RDS")
# Median/Pearson
ribMdPn = readRDS("Hide Lab/PCxN/output/Benchmark/RibosomeNoOverMdPnShrkJack.RDS")
# Median/Spearman
ribMdSp = readRDS("Hide Lab/PCxN/output/Benchmark/RibosomeNoOverMdSpShrkJack.RDS")

# ==== Random (No overlap) ====
# Mean/Pearson
rndMnPn = readRDS("Hide Lab/PCxN/output/Benchmark/RandomGSNoOverMnPnShrkJack.RDS")
# Mean/Spearman
rndMnSp = readRDS("Hide Lab/PCxN/output/Benchmark/RandomGSNoOverMnSpShrkJack.RDS")
# Median/Pearson
rndMdPn = readRDS("Hide Lab/PCxN/output/Benchmark/RandomGSNoOverMdPnShrkJack.RDS")
# Median/Spearman
rndMdSp = readRDS("Hide Lab/PCxN/output/Benchmark/RandomGSNoOverMdSpShrkJack.RDS")

library(rafalib)
# ==== Boxplots Bias (No overlap) =====
mypar(1,1)
png("Hide Lab/PCxN/doc/Benchmark/figures/rib_v_rnd_nonover_shrk_bias_boxplots01.png",
    width=800,height=450,pointsize = 15)
boxplot(
  list(
    abs(ribMnPn$bias_jack),abs(rndMnPn$bias_jack), # Mean/Pearson
    abs(ribMnSp$bias_jack),abs(rndMnSp$bias_jack), # Mean/Spearman
    abs(ribMdPn$bias_jack),abs(rndMdPn$bias_jack), # Median/Pearson
    abs(ribMdSp$bias_jack),abs(rndMdSp$bias_jack)  # Median/Spearman
  ),
  #main="Correlation Estimates (No Overlap)",
  names=rep(c("Ribosome","Random"),times=4)
)
title("|Bias| Estimates (No Overlap)",line=2)
abline(h=0,col="red")
text(x=1.5,y=line2user(line=0.5,side=3),"Mean/Pearson",xpd=T)
text(x=3.5,y=line2user(line=0.5,side=3),"Mean/Spearman",xpd=T)
text(x=5.5,y=line2user(line=0.5,side=3),"Median/Pearson",xpd=T)
text(x=7.5,y=line2user(line=0.5,side=3),"Median/Spearman",xpd=T)
dev.off()


png("Hide Lab/PCxN/doc/Benchmark/figures/rib_v_rnd_nonover_shrk_bias_boxplots02.png",
    width=800,height=450,pointsize = 15)
boxplot(list(abs(c(ribMnPn$bias_jack,rndMnPn$bias_jack)),
             abs(c(ribMnSp$bias_jack,rndMnSp$bias_jack)),
             abs(c(ribMdPn$bias_jack,rndMdPn$bias_jack)),
             abs(c(ribMdSp$bias_jack,rndMdSp$bias_jack))),
        names=c("Mean/Pearson","Mean/Spearman","Median/Pearson","Median/Spearman")
)
title("|Bias| Estimates (No Overlap)",line=2)
abline(h=0,col="red")
dev.off()

# ==== Median |Bias| ====
bias_vec = rep(NA,4)
names(bias_vec) = c("Mean/Pearson","Mean/Spearman","Median/Pearson","Median/Spearman")


bias_vec[1] = median(abs(c(ribMnPn$bias_jack,rndMnPn$bias_jack)))
bias_vec[2] = median(abs(c(ribMnSp$bias_jack,rndMnSp$bias_jack)))
bias_vec[3] = median(abs(c(ribMdPn$bias_jack,rndMdPn$bias_jack)))
bias_vec[4] = median(abs(c(ribMdSp$bias_jack,rndMdSp$bias_jack)))

rank(bias_vec*-1)

# data frame to store the method rankings
el1 = c("Mean","Median")
el2 = c("Pearson","Spearman")
el3 = c("Liptak","Logit")

scoreBoard = matrix(NA,ncol=8,nrow=1)
colnames(scoreBoard) = paste(rep(rep(el1,each=2),2),rep(rep(el2,2),2),rep(el3,each=4),sep="/")
rownames(scoreBoard) = c("Bias")
scoreBoard = as.data.frame(scoreBoard)

# ranks based on median variance
scoreBoard[1,] = rank(rep(bias_vec,2)*-1)

saveRDS(scoreBoard,"Hide Lab/PCxN/output/pcxn_benchmark_ranks_nover_bias.RDS")



