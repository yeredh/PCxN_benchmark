# ==== Ribosome Gene Sets =====
if(el1[io] == "NoOverlap"){
  # output directory
  setwd(paste0(pcxn_dir,"output/Benchmark/",paste(el1[io],el2[jo],el3[ko],sep="/"),"/GSE"))
  # matrix for experiment-level test statistics
  t_mat = matrix(NA,ncol=length(gse_ids),nrow=1000)
  colnames(t_mat) = gse_ids
  
  pb = txtProgressBar(min=0,max=length(gse_ids),style=3,initial=0)
  # loop thru GSE series (experiments)
  for(k in seq_along(gse_ids)){
    # load experiment-level results
    tmp = readRDS(paste0(gse_ids[k],"_ribosome_shrk_nonover.RDS"))
    # t-statistic
    t_mat[,k] = sapply(tmp,function(x){x$statistic})
    setTxtProgressBar(pb,k)
  }
  close(pb)
  
  
  rownames(t_mat) = sapply(tmp,function(x){x$Pathway.A})
  rm(tmp)
  
  # save matrix with t-statistics
  fname = paste0(
    "Ribosome",
    ifelse(el1[io] == "Overlap","Over","NoOver"),
    ifelse(el2[jo] == "Mean","Mn","Md"),
    ifelse(el3[ko] == "Pearson","Pn","Sp"),
    "Shrk",
    "Tstat",
    ".RDS"
  )
  cat("Benchmark:",paste(el1[io],el2[jo],el3[ko],sep="/"),"\n")
  cat(fname,"\n")
  saveRDS(t_mat,paste0(pcxn_dir,"output/Benchmark/",fname))
  
  
}

if(el1[io] == "Overlap"){
  # overlap cases
  for(over in 1:10){
    # output directory
    setwd(paste0(pcxn_dir,"output/Benchmark/",paste(el1[io],el2[jo],el3[ko],sep="/"),"/Case",over,"/GSE"))
    
    # matrix for experiment-level test statistics
    t_mat = matrix(NA,ncol=length(gse_ids),nrow=1000)
    colnames(t_mat) = gse_ids
    
    pb = txtProgressBar(min=0,max=length(gse_ids),style=3,initial=0)
    # loop thru GSE series (experiments)
    for(k in seq_along(gse_ids)){
      # load experiment-level results
      tmp = readRDS(paste0(gse_ids[k],"_ribosome_shrk_over",over,".RDS"))
      # t-statistic
      t_mat[,k] = sapply(tmp,function(x){x$statistic})
      setTxtProgressBar(pb,k)
    }
    close(pb)
    
    
    rownames(t_mat) = sapply(tmp,function(x){x$Pathway.A})
    rm(tmp)
    
    
    # save matrix with t-statistics
    fname = paste0(
      "Ribosome",
      ifelse(el1[io] == "Overlap","Over","NoOver"),
      over,
      ifelse(el2[jo] == "Mean","Mn","Md"),
      ifelse(el3[ko] == "Pearson","Pn","Sp"),
      "Shrk",
      "Tstat",
      ".RDS"
    )
    cat("Benchmark:",paste(el1[io],el2[jo],el3[ko],sep="/"),"\n")
    cat(fname,"\n")
    saveRDS(t_mat,paste0(pcxn_dir,"output/Benchmark/",fname))
  }
}
