# 05/14/2016
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
# For the overlap cases
#
# The metrics are the ranks based on the following
# - t-statistic: comparing ribosome estimates with random gene sets
#                we expect higher correlation estimates for ribosome gene sets
# the measures below are based on distinguishing true postives and true negatives
# based on the p-values from the ribosome gene sets and the random gene sets
# We assume that all p-values from a ribosome gene sets should be significant
# while the p-values from the random gene sets should not be significant
# - AUC
# - Accuracy
# - F1-Score


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

pvalPerformance = function(pv01,pv00,cut_off){
  # Computes performance measures based on a pvalue cut-off to distinguish significant and null cases
  #
  # Args:
  #   pv01: pvalues for the significant cases
  #   pv00: pvalues for the null case
  #   cut_off: cut-off to determine significance
  #
  # Returns:
  #   a named vector with the following 
  #     - Accuracy = (TP + TN)/(TP + FP + TN + FN)
  #     - Precision = TP/(TP + FP) 
  #     - Recall = TP/(TP + FN)
  #     - F1-score = 2(precision*recall)/(precision+recall)
  #     - False Positive Rate = FP/(FP+TN)
  
  pv <- c(pv01,pv00)
  y <- c(rep(1,length(pv01)),rep(0,length(pv00)))
  # compute confusion table
  conf_tab = table(Test=1*(pv<cut_off),Condition=y)
  if(dim(conf_tab)[1] == 1){
    if(rownames(conf_tab) == "1"){
      conf_tab = rbind(c(0,0),conf_tab)
    }else{
      conf_tab = rbind(conf_tab,c(0,0))
    }
    rownames(conf_tab) <- c("0","1")
  }
  # estimate accuracy=(TP + TN)/(TP + FP + TN + FN)
  acc = (conf_tab[1,1]+conf_tab[2,2])/sum(conf_tab)
  # estimate precision=TP/(TP + FP) 
  precision = conf_tab[2,2]/sum(conf_tab[2,])
  # estimate recall=TP/(TP + FN)
  recall = conf_tab[2,2]/sum(conf_tab[,2])
  # calculate F1 score
  F1 = (2*(precision*recall))/(precision+recall)
  # estimate false positive rate = FP/(FP+TN)
  FPR = conf_tab[2,1]/sum(conf_tab[,1])
  
  res <- c(acc,precision,recall,FPR,F1)
  names(res) <- c("Accuracy","Precision","Recall (TPR)","FPR","F1-score")
  return(res)
}

# list to store rankings for each overlap case
res = vector("list",10)
names(res) = paste0("over",1:10)

library(rafalib)
library(ROCR)
for(over in 1:10){
  # ==== Ribosome (Overlap) =====
  # Mean/Pearson
  ribMnPn = readRDS(paste0("Hide Lab/PCxN/output/Benchmark/RibosomeOver",over,"MnPnShrk.RDS"))
  # Mean/Spearman
  ribMnSp = readRDS(paste0("Hide Lab/PCxN/output/Benchmark/RibosomeOver",over,"MnSpShrk.RDS"))
  # Median/Pearson
  ribMdPn = readRDS(paste0("Hide Lab/PCxN/output/Benchmark/RibosomeOver",over,"MdPnShrk.RDS"))
  # Median/Spearman
  ribMdSp = readRDS(paste0("Hide Lab/PCxN/output/Benchmark/RibosomeOver",over,"MdSpShrk.RDS"))
  
  # ==== Random (Overlap) ====
  # Mean/Pearson
  rndMnPn = readRDS(paste0("Hide Lab/PCxN/output/Benchmark/RandomGSOver",over,"MnPnShrk.RDS"))
  # Mean/Spearman
  rndMnSp = readRDS(paste0("Hide Lab/PCxN/output/Benchmark/RandomGSOver",over,"MnSpShrk.RDS"))
  # Median/Pearson
  rndMdPn = readRDS(paste0("Hide Lab/PCxN/output/Benchmark/RandomGSOver",over,"MdPnShrk.RDS"))
  # Median/Spearman
  rndMdSp = readRDS(paste0("Hide Lab/PCxN/output/Benchmark/RandomGSOver",over,"MdSpShrk.RDS"))
  
  
  # ==== Boxplots Correlation (Overlap) =====
  mypar(1,1)
  png(paste0("Hide Lab/PCxN/doc/Benchmark/figures/rib_v_rnd_over",over,"_shrk_cor_boxplots.png"),
      width=800,height=450,pointsize = 15)
  boxplot(
    list(
      ribMnPn$PathCor,rndMnPn$PathCor, # Mean/Pearson
      ribMnSp$PathCor,rndMnSp$PathCor, # Mean/Spearman
      ribMdPn$PathCor,rndMdPn$PathCor, # Median/Pearson
      ribMdSp$PathCor,rndMdSp$PathCor  # Median/Spearman
    ),
    #main="Correlation Estimates (No Overlap)",
    names=rep(c("Ribosome","Random"),times=4)
  )
  title(paste0("Correlation Estimates (Overlap ",over,")"),line=2)
  abline(h=0,col="red")
  text(x=1.5,y=line2user(line=0.5,side=3),"Mean/Pearson",xpd=T)
  text(x=3.5,y=line2user(line=0.5,side=3),"Mean/Spearman",xpd=T)
  text(x=5.5,y=line2user(line=0.5,side=3),"Median/Pearson",xpd=T)
  text(x=7.5,y=line2user(line=0.5,side=3),"Median/Spearman",xpd=T)
  dev.off()
  
  # === t-test ====
  t_vec = rep(NA,4)
  names(t_vec) = c("Mean/Pearson","Mean/Spearman","Median/Pearson","Median/Spearman")
  
  
  t_vec[1] = t.test(ribMnPn$PathCor,rndMnPn$PathCor,"greater")$statistic
  t_vec[2] = t.test(ribMnSp$PathCor,rndMnSp$PathCor,"greater")$statistic
  t_vec[3] = t.test(ribMdPn$PathCor,rndMdPn$PathCor,"greater")$statistic
  t_vec[4] = t.test(ribMdSp$PathCor,rndMdSp$PathCor,"greater")$statistic
  
  # rankings by t-test statistic
  rank(t_vec)
  
  
  # ==== ROC Curves:  Liptak's Method (Sample Size) ====
  # labels for Ribosome and random gene set cases
  sig_labels = rep(c(1,0),each=1000)
  
  # create prediction object based on pvalues
  predLipMnPn = prediction(c(ribMnPn$LiptakSS < 0.05,rndMnPn$LiptakSS<0.05)*1,sig_labels)
  predLipMnSp = prediction(c(ribMnSp$LiptakSS < 0.05,rndMnSp$LiptakSS<0.05)*1,sig_labels)
  predLipMdPn = prediction(c(ribMdPn$LiptakSS < 0.05,rndMdPn$LiptakSS<0.05)*1,sig_labels)
  predLipMdSp = prediction(c(ribMdSp$LiptakSS < 0.05,rndMdSp$LiptakSS<0.05)*1,sig_labels)
  
  
  # calculate true positive rate and false positive rate for the ROC curve
  rocLipMnPn = performance(predLipMnPn,"tpr","fpr")
  rocLipMnSp = performance(predLipMnSp,"tpr","fpr")
  rocLipMdPn = performance(predLipMdPn,"tpr","fpr")
  rocLipMdSp = performance(predLipMdSp,"tpr","fpr")
  
  # calculate AUC
  aucLipMnPn = as.numeric(performance(predLipMnPn,"auc")@y.values)
  aucLipMnSp = as.numeric(performance(predLipMnSp,"auc")@y.values)
  aucLipMdPn = as.numeric(performance(predLipMdPn,"auc")@y.values)
  aucLipMdSp = as.numeric(performance(predLipMdSp,"auc")@y.values)
  
  LipAuc = c(
    aucLipMnPn,
    aucLipMnSp,
    aucLipMdPn,
    aucLipMdSp
  )
  names(LipAuc) = c("Mean/Pearson","Mean/Spearman","Median/Pearson","Median/Spearman")
  # rankings by AUC
  rank(LipAuc)
  
  
  # ROC curves
  png(paste0("Hide Lab/PCxN/doc/Benchmark/figures/roc_liptak_over",over,"_sep.png"),
      width=800,height=450,pointsize = 15)
  mypar(2,2)
  plot(rocLipMnPn,lty=1,col=1,main="Mean/Pearson")
  abline(c(0,1),col="gray",lty=2)
  plot(rocLipMnSp,lty=1,col=2,main="Mean/Spearman")
  abline(c(0,1),col="gray",lty=2)
  plot(rocLipMdPn,lty=1,col=3,main="Median/Pearson")
  abline(c(0,1),col="gray",lty=2)
  plot(rocLipMdSp,lty=1,col=4,main="Median/Spearman")
  abline(c(0,1),col="gray",lty=2)
  dev.off()
  
  
  # ROC curves
  png(paste0("Hide Lab/PCxN/doc/Benchmark/figures/roc_liptak_over",over,"_comb.png"),
      width=800,height=450,pointsize = 15)
  mypar(1,1)
  plot(rocLipMnPn,lty=1,col=1)
  plot(rocLipMnSp,lty=1,col=2,add=T)
  plot(rocLipMdPn,lty=1,col=3,add=T)
  plot(rocLipMdSp,lty=1,col=4,add=T)
  legend("bottomright",c("Mn/Pn","Mn/Sp","Md/Pn","Md/Sp"),pch=19,col=1:4,bty="n")
  dev.off()
  
  # ==== P-value Metrics:  Liptak's Method (Sample Size) ====
  ppLipMnPn = pvalPerformance(pv01=ribMnPn$LiptakSS,pv00=rndMnPn$LiptakSS,cut_off=0.05)
  ppLipMnSp = pvalPerformance(pv01=ribMnSp$LiptakSS,pv00=rndMnSp$LiptakSS,cut_off=0.05)
  ppLipMdPn = pvalPerformance(pv01=ribMdPn$LiptakSS,pv00=rndMdPn$LiptakSS,cut_off=0.05)
  ppLipMdSp = pvalPerformance(pv01=ribMdSp$LiptakSS,pv00=rndMdSp$LiptakSS,cut_off=0.05)
  
  # Accuracy
  LipAcc = c(
    ppLipMnPn[1],
    ppLipMnSp[1],
    ppLipMdPn[1],
    ppLipMdSp[1]
  )
  names(LipAcc) = c("Mean/Pearson","Mean/Spearman","Median/Pearson","Median/Spearman")
  # rankings by accuracy
  rank(LipAcc)
  
  # F1-score
  LipF1 = c(
    ppLipMnPn[5],
    ppLipMnSp[5],
    ppLipMdPn[5],
    ppLipMdSp[5]
  )
  names(LipF1) = c("Mean/Pearson","Mean/Spearman","Median/Pearson","Median/Spearman")
  # rankings by accuracy
  rank(LipF1)
  
  
  # ==== ROC Curves:  Logit Method ====
  # labels for Ribosome and random gene set cases
  sig_labels = rep(c(1,0),each=1000)
  
  # create prediction object based on pvalues
  predLogMnPn = prediction(c(ribMnPn$LogitP < 0.05,rndMnPn$LogitP<0.05)*1,sig_labels)
  predLogMnSp = prediction(c(ribMnSp$LogitP < 0.05,rndMnSp$LogitP<0.05)*1,sig_labels)
  predLogMdPn = prediction(c(ribMdPn$LogitP < 0.05,rndMdPn$LogitP<0.05)*1,sig_labels)
  predLogMdSp = prediction(c(ribMdSp$LogitP < 0.05,rndMdSp$LogitP<0.05)*1,sig_labels)
  
  
  # calculate true positive rate and false positive rate for the ROC curve
  rocLogMnPn = performance(predLogMnPn,"tpr","fpr")
  rocLogMnSp = performance(predLogMnSp,"tpr","fpr")
  rocLogMdPn = performance(predLogMdPn,"tpr","fpr")
  rocLogMdSp = performance(predLogMdSp,"tpr","fpr")
  
  # calculate AUC
  aucLogMnPn = as.numeric(performance(predLogMnPn,"auc")@y.values)
  aucLogMnSp = as.numeric(performance(predLogMnSp,"auc")@y.values)
  aucLogMdPn = as.numeric(performance(predLogMdPn,"auc")@y.values)
  aucLogMdSp = as.numeric(performance(predLogMdSp,"auc")@y.values)
  
  LogAuc = c(
    aucLogMnPn,
    aucLogMnSp,
    aucLogMdPn,
    aucLogMdSp
  )
  names(LogAuc) = c("Mean/Pearson","Mean/Spearman","Median/Pearson","Median/Spearman")
  # rankings by AUC
  rank(LogAuc)
  
  
  # ROC curves
  png(paste0("Hide Lab/PCxN/doc/Benchmark/figures/roc_logit_over",over,"_sep.png"),
      width=800,height=450,pointsize = 15)
  mypar(2,2)
  plot(rocLogMnPn,lty=1,col=1,main="Mean/Pearson")
  abline(c(0,1),col="gray",lty=2)
  plot(rocLogMnSp,lty=1,col=2,main="Mean/Spearman")
  abline(c(0,1),col="gray",lty=2)
  plot(rocLogMdPn,lty=1,col=3,main="Median/Pearson")
  abline(c(0,1),col="gray",lty=2)
  plot(rocLogMdSp,lty=1,col=4,main="Median/Spearman")
  abline(c(0,1),col="gray",lty=2)
  dev.off()
  
  # ROC curves
  png(paste0("Hide Lab/PCxN/doc/Benchmark/figures/roc_logit_over",over,"_comb.png"),
      width=800,height=450,pointsize = 15)
  mypar(1,1)
  plot(rocLogMnPn,lty=1,col=1)
  plot(rocLogMnSp,lty=1,col=2,add=T)
  plot(rocLogMdPn,lty=1,col=3,add=T)
  plot(rocLogMdSp,lty=1,col=4,add=T)
  legend("bottomright",c("Mn/Pn","Mn/Sp","Md/Pn","Md/Sp"),pch=19,col=1:4,bty="n")
  dev.off()
  
  
  
  # ==== P-value Metrics:  Logit Method ====
  ppLogMnPn = pvalPerformance(pv01=ribMnPn$LogitP,pv00=rndMnPn$LogitP,cut_off=0.05)
  ppLogMnSp = pvalPerformance(pv01=ribMnSp$LogitP,pv00=rndMnSp$LogitP,cut_off=0.05)
  ppLogMdPn = pvalPerformance(pv01=ribMdPn$LogitP,pv00=rndMdPn$LogitP,cut_off=0.05)
  ppLogMdSp = pvalPerformance(pv01=ribMdSp$LogitP,pv00=rndMdSp$LogitP,cut_off=0.05)
  
  # Accuracy
  LogAcc = c(
    ppLogMnPn[1],
    ppLogMnSp[1],
    ppLogMdPn[1],
    ppLogMdSp[1]
  )
  names(LogAcc) = c("Mean/Pearson","Mean/Spearman","Median/Pearson","Median/Spearman")
  # rankings by accuracy
  rank(LogAcc)
  
  # F1-score
  LogF1 = c(
    ppLogMnPn[5],
    ppLogMnSp[5],
    ppLogMdPn[5],
    ppLogMdSp[5]
  )
  names(LogF1) = c("Mean/Pearson","Mean/Spearman","Median/Pearson","Median/Spearman")
  # rankings by accuracy
  rank(LogF1)
  
  
  # ==== Overall Rankings ====
  
  # data frame to store the method rankings
  el1 = c("Mean","Median")
  el2 = c("Pearson","Spearman")
  el3 = c("Liptak","Logit")
  
  
  
  scoreBoard = matrix(NA,ncol=8,nrow=4)
  colnames(scoreBoard) = paste(rep(rep(el1,each=2),2),rep(rep(el2,2),2),rep(el3,each=4),sep="/")
  rownames(scoreBoard) = c("t.test","AUC","Accuracy","F1-score")
  scoreBoard = as.data.frame(scoreBoard)
  
  # ranks based on t-statistic
  scoreBoard[1,] = rank(rep(t_vec,2))
  # ranks based on AUC
  scoreBoard[2,] = rank(c(LipAuc,LogAuc))
  # ranks based on Accuracy
  scoreBoard[3,] = rank(c(LipAcc,LogAcc))
  # ranks based on F1-scores
  scoreBoard[4,] = rank(c(LipF1,LogF1))
  
  res[[over]] = scoreBoard
}

saveRDS(res,"Hide Lab/PCxN/output/pcxn_benchmark_ranks_over.RDS")
