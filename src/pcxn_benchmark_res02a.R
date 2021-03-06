# ==== Ribosome Gene Sets =====
if(el1[io] == "NoOverlap"){
  # output directory
  setwd(paste0(pcxn_dir,"output/Benchmark/",paste(el1[io],el2[jo],el3[ko],sep="/"),"/GSE"))
  # matrix for experiment-level correlation estimates
  r_mat = matrix(NA,ncol=length(gse_ids),nrow=1000)
  colnames(r_mat) = gse_ids
  # matrix for experiment-level p-values
  p_mat = matrix(NA,ncol=length(gse_ids),nrow=1000)
  colnames(p_mat) = gse_ids
  # vector to store experiment sample size
  n_vec = rep(NA,length(gse_ids))
  names(n_vec) = gse_ids

  pb = txtProgressBar(min=0,max=length(gse_ids),style=3,initial=0)
  # loop thru GSE series (experiments)
  for(k in seq_along(gse_ids)){
    # load experiment-level results
    tmp = readRDS(paste0(gse_ids[k],"_ribosome_shrk_nonover.RDS"))
    # correlation estimates
    r_mat[,k] = sapply(tmp,function(x){x$estimate})
    # p-values
    p_mat[,k] = sapply(tmp,function(x){x$p.value})
    # sample size
    n_vec[k] = unique(sapply(tmp,function(x){x$n}))
    setTxtProgressBar(pb,k)
  }
  close(pb)


  rownames(r_mat) = sapply(tmp,function(x){x$Pathway.A})
  rownames(p_mat) = sapply(tmp,function(x){x$Pathway.A})
  rm(tmp)

  # adjust p-values and correlation estimates
  p_mat = AdjustPmat(p_mat)
  r_mat = AdjustRmat(r_mat)

  # data frame to store results
  res = data.frame(PathCor=rep(NA,1000),LogitP=rep(NA,1000),LiptakSS=rep(NA,1000),LiptakES=rep(NA,1000))

  # get correlation weighted average
  res$PathCor = c((r_mat %*% n_vec)/sum(n_vec))

  # combined p-values
  # Logit Method
  res$LogitP = apply(p_mat,1,function(x){logitp(x)$p})
  # Liptak's Method (sample size)
  res$LiptakSS = apply(p_mat,1,function(x){sumz(p=x,weights=n_vec)$p})
  # Liptak's Method (effect size)
  for(i in 1:1000){
    res$LiptakES[i] = sumz(p=p_mat[i,],weights=abs(r_mat[i,]))$p
  }
  # save results
  fname = paste0(
    "Ribosome",
    ifelse(el1[io] == "Overlap","Over","NoOver"),
    ifelse(el2[jo] == "Mean","Mn","Md"),
    ifelse(el3[ko] == "Pearson","Pn","Sp"),
    "Shrk",
    ".RDS"
  )
  cat("Benchmark:",paste(el1[io],el2[jo],el3[ko],sep="/"),"\n")
  cat(fname,"\n")
  saveRDS(res,paste0(pcxn_dir,"output/Benchmark/",fname))
}

if(el1[io] == "Overlap"){
  # overlap cases
  for(over in 1:10){
    # output directory
    setwd(paste0(pcxn_dir,"output/Benchmark/",paste(el1[io],el2[jo],el3[ko],sep="/"),"/Case",over,"/GSE"))

    # matrix for experiment-level correlation estimates
    r_mat = matrix(NA,ncol=length(gse_ids),nrow=1000)
    colnames(r_mat) = gse_ids
    # matrix for experiment-level p-values
    p_mat = matrix(NA,ncol=length(gse_ids),nrow=1000)
    colnames(p_mat) = gse_ids
    # vector to store experiment sample size
    n_vec = rep(NA,length(gse_ids))
    names(n_vec) = gse_ids

    pb = txtProgressBar(min=0,max=length(gse_ids),style=3,initial=0)
    # loop thru GSE series (experiments)
    for(k in seq_along(gse_ids)){
      # load experiment-level results
      tmp = readRDS(paste0(gse_ids[k],"_ribosome_shrk_over",over,".RDS"))
      # correlation estimates
      r_mat[,k] = sapply(tmp,function(x){x$estimate})
      # p-values
      p_mat[,k] = sapply(tmp,function(x){x$p.value})
      # sample size
      n_vec[k] = unique(sapply(tmp,function(x){x$n}))
      setTxtProgressBar(pb,k)
    }
    close(pb)


    rownames(r_mat) = sapply(tmp,function(x){x$Pathway.A})
    rownames(p_mat) = sapply(tmp,function(x){x$Pathway.A})
    rm(tmp)

    # adjust p-values and correlation estimates
    p_mat = AdjustPmat(p_mat)
    r_mat = AdjustRmat(r_mat)


    # data frame to store results
    res = data.frame(PathCor=rep(NA,1000),LogitP=rep(NA,1000),LiptakSS=rep(NA,1000),LiptakES=rep(NA,1000))

    # get correlation weighted average
    res$PathCor = c((r_mat %*% n_vec)/sum(n_vec))

    # combined p-values
    # Logit Method
    res$LogitP = apply(p_mat,1,function(x){logitp(x)$p})
    # Liptak's Method (sample size)
    res$LiptakSS = apply(p_mat,1,function(x){sumz(p=x,weights=n_vec)$p})
    # Liptak's Method (effect size)
    for(i in 1:1000){
      res$LiptakES[i] = sumz(p=p_mat[i,],weights=abs(r_mat[i,]))$p
    }
    # save results
    fname = paste0(
      "Ribosome",
      ifelse(el1[io] == "Overlap","Over","NoOver"),
      over,
      ifelse(el2[jo] == "Mean","Mn","Md"),
      ifelse(el3[ko] == "Pearson","Pn","Sp"),
      "Shrk",
      ".RDS"
    )
    cat("Benchmark:",paste(el1[io],el2[jo],el3[ko],sep="/"),"\n")
    cat(fname,"\n")
    saveRDS(res,paste0(pcxn_dir,"output/Benchmark/",fname))
  }
}
