# ==== Random Gene Sets ====
if(el1[io] == "NoOverlap"){
  # output directory
  setwd(paste0(pcxn_dir,"output/Benchmark/",paste(el1[io],el2[jo],el3[ko],sep="/"),"/GSE"))
  # matrix for experiment-level test statistics
  t_mat = matrix(NA,ncol=length(gse_ids),nrow=1000)
  colnames(t_mat) = gse_ids
  
  pb = txtProgressBar(min=0,max=1000,style=3,initial=0)
  # loop thru random gene sets
  for(k in 1:1000){
    # load experiment-level results
    tmp = readRDS(paste0("RandomGS",k,"_shrk_nonover.RDS"))
    # correlation estimates
    t_mat[k,] = sapply(tmp,function(x){x$statistic})
    setTxtProgressBar(pb,k)
  }
  close(pb)
  
  rownames(t_mat) = paste0("RandomGS",1:1000)
  rm(tmp)
  
  # save matrix with t-statistics
  fname = paste0(
    "RandomGS",
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
    
    pb = txtProgressBar(min=0,max=1000,style=3,initial=0)
    # loop thru GSE series (experiments)
    for(k in 1:1000){
      # load experiment-level results
      tmp = readRDS(paste0("RandomGS",k,"_shrk_over",over,".RDS"))
      # correlation estimates
      t_mat[k,] = sapply(tmp,function(x){x$statistic})
      setTxtProgressBar(pb,k)
    }
    
    close(pb)
    
    rownames(t_mat) = paste0("RandomGS",1:1000)
    rm(tmp)
    
    # save matrix with t-statistics
    fname = paste0(
      "RandomGS",
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

cat("Benchmark:",paste(el1[io],el2[jo],el3[ko],sep="/"),"\n")
