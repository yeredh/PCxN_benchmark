# 04/08/2016
#
# Separate the n elements (in this case corresponding to the number of
# genes in the Ribosome pathway) into two sets with common elements
# Then, select 10 representative splits to represent various degrees of
# overlap (measured with the overlap coefficient) ranging from low
# overlap to high overlap

rm(list=ls())
options(stringsAsFactors = F)
# ==== Functions ====
# function to calculate overlap coefficient
OverlapCoeff = function(x,y){length(intersect(x,y))/min(c(length(x),length(y)))}

GetOverSets = function(n){
  # split n elements into two sets s1 and s2
  # with varying degrees of overlap
  # 
  # In the first step, the two sets s1 and s2 share all but one element. 
  # In each consecutive step, we shift the indexes of one of the sets 
  # to decrease the number of shared elements between s1 and s2 until
  # the last step when the two sets s1 and s2 do not have any elements in common.
  #
  # Input
  #
  # n: number of elements to split into two sets
  #
  # Returns
  #  a list with the indices for the elements in s1 and s2
  
  # list to store the indices for sets s1 and s2
  s_ind = vector("list",n)
  
  # left shift for s1 indices
  s1_Lshift = 1
  # counter for s1 left left shift
  s1_ic=0
  # right shift for s2 indices
  s2_Rshift=1
  # counter for s2 right shift
  s2_ic=1
  
  for(k in 1:n){
    # initialize shift counters
    s1_ic=s1_ic+1
    s2_ic=s2_ic+1
    
    # sliding window to separate the the n elements into two overlapping sets
    s_ind[[k]]$s1 = 1:(n-s1_Lshift)
    s_ind[[k]]$s2 = s2_Rshift:n
    
    # swift sliding window one direction (right/left) at a time
    if(s1_ic == 2){
      s1_Lshift = s1_Lshift + 1
      s1_ic = 0
    }
    if(s2_ic == 2){
      s2_Rshift=s2_Rshift+1
      s2_ic=0
    }
  }
  
  return(s_ind)
}


# ===== KEGG Human Pathways =====
# downloaded using KEGGREST on 9/2014
kegg_gs = readRDS("Hide Lab/PCxN/data/Gene Sets/KEGG_hsa.RDS")

# genes in microarray background
entrez_ids = readRDS("Hide Lab/PCxN/data/GPL570_genes.RDS")

# filter by genes present in array
# KEGG pathways
kegg_gs$gs = lapply(kegg_gs$gs,function(x){x[x %in% entrez_ids]})


# ===== Overlap Cases for Ribosome Pathway ====
n = length(kegg_gs$gs[[grep("Ribosome - Homo sapiens",kegg_gs$pathway_names,value=F)]])
# split genes into overlapping sets
gs_ind = GetOverSets(n)
# estimate overlap coefficient
over_vec = sapply(gs_ind,function(x){OverlapCoeff(x$s1,x$s2)})
summary(over_vec)

# a grid over different values for the overlap coefficient
over_grid = seq(0.05,0.95,length.out = 10)

# get the absolute value of the difference between the overlap for each 
# of the sets generated from the sliding window with the overlap grid
# to select 10 representative overlap cases
dist_mat = matrix(NA,nrow=length(over_vec),ncol=length(over_grid))
colnames(dist_mat) = over_grid
for(k in 1:length(over_vec)){
  dist_mat[k,] = abs(over_grid - over_vec[k])
}

# find the sets with the closest overlap to the grid
sel_ind = apply(dist_mat,2,which.min)

# for(j in 1:length(sel_ind)){
#   cat(paste0("\\subfigure[Case ",j," ($o_{AB}=",round(over_vec[sel_ind[j]],4),
# "$) ]{\\includegraphics[width=4.75cm]{figures/overlap_case",j,"_venn.png}} \n"))
# }

# save selected overlap cases
saveRDS(gs_ind[sel_ind],"Hide Lab/PCxN/data/pcxn_overlap_cases.RDS")

# ==== Figures ====
# overlap coefficient histogram for selected cases
png("Hide Lab/PCxN/doc/Benchmark/figures/overlap_cases_hist.png",
    width=600,height=450,pointsize = 15)
hist(over_vec[sel_ind],main="Overlap Cases",
     ylab="",xlab="Overlap Coefficient")
dev.off()


## Venn Diagrams
library(VennDiagram)
venn_counts = t(sapply(gs_ind[sel_ind],
                       function(x){
                         c(length(x$s1),length(x$s2),length(intersect(x$s1,x$s2)))
                         }))

for(i in 1:dim(venn_counts)[1]){
  png(paste0("Hide Lab/PCxN/doc/Benchmark/figures/overlap_case",i,"_venn.png"),pointsize=17,width=650,height=550)
  venn_plot = draw.pairwise.venn(venn_counts[i,1],venn_counts[i,2],venn_counts[i,3],
                                 c("A", "B"),cex=2.5,cat.cex=2.5)
  grid.draw(venn_plot)
  dev.off()
  grid.newpage()
}

