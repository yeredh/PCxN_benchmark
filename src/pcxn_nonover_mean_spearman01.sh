#!/bin/sh
#SBATCH -J pcxn_nonover_mean_spearman # A single job name for the array
#SBATCH -n 16 # Number of cores
#SBATCH -N 1 # All cores on one machine
#SBATCH -p serial_requeue,irizarry # Partition
#SBATCH --mem 50000 # Memory request
#SBATCH -t 0-00:15 # (D-HH:MM)
#SBATCH -o /net/hsphfs1/srv/export/hsphfs1/share_root/hide_lab/PCxN/log/pcxn_nonover_mean_spearman_%a.out # Standard output
#SBATCH -e /net/hsphfs1/srv/export/hsphfs1/share_root/hide_lab/PCxN/log/pcxn_nonover_mean_spearman_%a.err # Standard error
#SBATCH --mail-type=ALL # Type of email notification- BEGIN,END,FAIL,ALL
#SBATCH --mail-user=yered.h@gmail.com # Email to which notifications will be sent

source new-modules.sh
module load R/3.2.2-fasrc01
export R_LIBS_USER=$HOME/apps/R/3.2.2-fasrc01:$R_LIBS_USE

R CMD BATCH "--args ${SLURM_ARRAY_TASK_ID}" /net/hsphfs1/srv/export/hsphfs1/share_root/hide_lab/PCxN/src/pcxn_nonover_mean_spearman01.R /net/hsphfs1/srv/export/hsphfs1/share_root/hide_lab/PCxN/log/pcxn_nonover_mean_spearman_"${SLURM_ARRAY_TASK_ID}".Rout
