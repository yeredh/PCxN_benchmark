# ==== Random Gene Sets ====
if(el1[io] == "NoOverlap"){
  # output directory
  setwd(paste0(pcxn_dir,"output/Benchmark/",paste(el1[io],el2[jo],el3[ko],sep="/"),"/GSE"))
  # matrix for experiment-level correlation estimates
  r_mat = matrix(NA,ncol=length(gse_ids),nrow=1000)
  colnames(r_mat) = gse_ids
  
  pb = txtProgressBar(min=0,max=1000,style=3,initial=0)
  # loop thru random gene sets
  for(k in 1:1000){
    # load experiment-level results
    tmp = readRDS(paste0("RandomGS",k,"_shrk_nonover.RDS"))
    # correlation estimates
    r_mat[k,] = sapply(tmp,function(x){x$estimate})
    setTxtProgressBar(pb,k)
  }
  close(pb)
  # sample size
  n_vec = sapply(tmp,function(x){x$n})
  rownames(r_mat) = paste0("RandomGS",1:1000)
  rm(tmp)
  
  # adjust correlation estimates
  r_mat = AdjustRmat(r_mat)
  # get correlation weighted average
  r_bar = c((r_mat %*% n_vec)/sum(n_vec))
  # get variance estimates
  hs_var = GetHSVar(r_mat,r_bar,n_vec)
  
  # save results
  fname = paste0(
    "RandomGS",
    ifelse(el1[io] == "Overlap","Over","NoOver"),
    ifelse(el2[jo] == "Mean","Mn","Md"),
    ifelse(el3[ko] == "Pearson","Pn","Sp"),
    "Shrk",
    "Var",
    ".RDS"
  )
  cat("Benchmark:",paste(el1[io],el2[jo],el3[ko],sep="/"),"\n")
  cat(fname,"\n")
  saveRDS(hs_var,paste0(pcxn_dir,"output/Benchmark/",fname))
}


if(el1[io] == "Overlap"){
  # overlap cases
  for(over in 1:10){
    # output directory
    setwd(paste0(pcxn_dir,"output/Benchmark/",paste(el1[io],el2[jo],el3[ko],sep="/"),"/Case",over,"/GSE"))
    
    # matrix for experiment-level correlation estimates
    r_mat = matrix(NA,ncol=length(gse_ids),nrow=1000)
    colnames(r_mat) = gse_ids
    
    pb = txtProgressBar(min=0,max=1000,style=3,initial=0)
    # loop thru GSE series (experiments)
    for(k in 1:1000){
      # load experiment-level results
      tmp = readRDS(paste0("RandomGS",k,"_shrk_over",over,".RDS"))
      # correlation estimates
      r_mat[k,] = sapply(tmp,function(x){x$estimate})
      setTxtProgressBar(pb,k)
    }
    
    close(pb)
    
    # sample size
    n_vec = sapply(tmp,function(x){x$n})
    rownames(r_mat) = paste0("RandomGS",1:1000)
    rm(tmp)
    
    # adjust correlation estimates
    r_mat = AdjustRmat(r_mat)
    # get correlation weighted average
    r_bar = c((r_mat %*% n_vec)/sum(n_vec))
    # get variance estimates
    hs_var = GetHSVar(r_mat,r_bar,n_vec)
    # save results
    fname = paste0(
      "RandomGS",
      ifelse(el1[io] == "Overlap","Over","NoOver"),
      over,
      ifelse(el2[jo] == "Mean","Mn","Md"),
      ifelse(el3[ko] == "Pearson","Pn","Sp"),
      "Shrk",
      "Var",
      ".RDS"
    )
    cat("Benchmark:",paste(el1[io],el2[jo],el3[ko],sep="/"),"\n")
    cat(fname,"\n")
    saveRDS(hs_var,paste0(pcxn_dir,"output/Benchmark/",fname))
  }
}

cat("Benchmark:",paste(el1[io],el2[jo],el3[ko],sep="/"),"\n")
