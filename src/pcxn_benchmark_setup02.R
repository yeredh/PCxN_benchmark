# 04/08/2016
#
# Create directories to save the results from the
# benchmark
#
# ODYSSEY (INTERACTIVE SESSION)
#
# source new-modules.sh
# module load R/3.2.2-fasrc01
# export R_LIBS_USER=$HOME/apps/R/3.2.2-fasrc01:$R_LIBS_USE
# srun -p irizarry,serial_requeue --mem-per-cpu=2000 -t 0-0:30 --pty R

rm(list=ls())
options(stringsAsFactors = F)
pcxn_dir = "/net/hsphfs1/srv/export/hsphfs1/share_root/hide_lab/PCxN/"

# ==== Benchmark Directories ====
# create directories to store output from benchmark

# No Overlap Cases
for(el1 in c("Mean/","Median/")){
  for(el2 in c("Pearson/","Spearman/")){
    new_dir = paste0(pcxn_dir,"output/Benchmark/NoOverlap/",el1,el2,"GSE")
    dir.create(new_dir,recursive = T)
  }
}

# Overlap cases
for(el1 in c("Mean/","Median/")){
  for(el2 in c("Pearson/","Spearman/")){
    for(el3 in paste0("Case",1:10)){
      new_dir = paste0(pcxn_dir,"output/Benchmark/Overlap/",el1,el2,el3,"/GSE")
      dir.create(new_dir,recursive = T)
    }
  }
}

rm(list=ls())

