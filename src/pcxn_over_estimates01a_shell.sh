# 04/11/2016
#
# BASH script to pass an argument to select the overlap case to get the estimates
# for the gene sets based on the Ribosome pathway
# The argument determines
# - job name
# - standard output file name
# - standard error file name
# - R output file name

for OVER in $(seq 3 10); do
  #
  echo "Overlap Case ${OVER}"
  export OVER
  sbatch --job-name=pcxn_over${OVER}_estimates \
  -o /net/hsphfs1/srv/export/hsphfs1/share_root/hide_lab/PCxN/log/pcxn_over${OVER}_estimates%a.out \
  -e /net/hsphfs1/srv/export/hsphfs1/share_root/hide_lab/PCxN/log/pcxn_over${OVER}_estimates%a.err \
  --array=11-863 \
  pcxn_over_estimates01a.sh
  # pause between each job submission
  sleep 1
done
