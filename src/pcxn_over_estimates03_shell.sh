# 04/19/2016
#
# BASH script to pass an argument to select the overlap case to get the estimates
# for the gene sets based on the Ribosome pathway
# The argument determines
# - job name
# - standard output file name
# - standard error file name
# - R output file name

for OVER in $(seq 9 10); do
  # select overlap case
  echo "Overlap Case ${OVER}"
  export OVER
  # set SLURM parameters and submit job array
  sbatch --job-name=pcxn_shrk_over${OVER}_estimates \
  -o /net/hsphfs1/srv/export/hsphfs1/share_root/hide_lab/PCxN/log/pcxn_shrk_over${OVER}_estimates%a.out \
  -e /net/hsphfs1/srv/export/hsphfs1/share_root/hide_lab/PCxN/log/pcxn_shrk_over${OVER}_estimates%a.err \
  --array=1-863 \
  pcxn_over_estimates03.sh
  # pause between each job submission
  sleep 1
done

#
# squeue -u ypitajuarez --format="%i %P %j %u %t %M %D %R"
# for JID in $(seq 864 1000); do
# echo "Cancelling ${JID}"
# scancel 60399310_${JID}
# sleep 0.5
# done
# sleep 10
# squeue -u ypitajuarez --format="%i %P %j %u %t %M %D %R"
#
