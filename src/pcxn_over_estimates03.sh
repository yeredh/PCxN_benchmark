#!/bin/sh
#
# the parameter for overlap case {OVER} is passed from bash
#
#SBATCH -n 16 # Number of cores
#SBATCH -N 1 # All cores on one machine
#SBATCH -p serial_requeue,irizarry # Partition
#SBATCH --mem 4000 # Memory request
#SBATCH -t 0-00:03 # (D-HH:MM)
#SBATCH --mail-type=ALL # Type of email notification- BEGIN,END,FAIL,ALL
#SBATCH --mail-user=yered.h@gmail.com # Email to which notifications will be sent

source new-modules.sh
module load R/3.2.2-fasrc01
export R_LIBS_USER=$HOME/apps/R/3.2.2-fasrc01:$R_LIBS_USE

R CMD BATCH "--args ${OVER} ${SLURM_ARRAY_TASK_ID}" /net/hsphfs1/srv/export/hsphfs1/share_root/hide_lab/PCxN/src/pcxn_over_estimates03.R /net/hsphfs1/srv/export/hsphfs1/share_root/hide_lab/PCxN/log/pcxn_shrk_over_"${OVER}"_estimates_"${SLURM_ARRAY_TASK_ID}".Rout
