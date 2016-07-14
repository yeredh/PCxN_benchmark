# 04/19/2016
#
# BASH script to pass an argument to select the overlap case to get the estimates
# for the random gene sets
# The argument determines
# - job name
# - standard output file name
# - standard error file name
# - R output file name

for OVER in $(seq 1 10); do
  # select overlap case
  echo "Overlap Case ${OVER}"
  export OVER
  # set SLURM parameters and submit job array
  sbatch --job-name=pcxn_shrk_over${OVER}_estimates_rnd \
  -o /net/hsphfs1/srv/export/hsphfs1/share_root/hide_lab/PCxN/log/pcxn_shrk_over${OVER}_estimates_rnd_%a.out \
  -e /net/hsphfs1/srv/export/hsphfs1/share_root/hide_lab/PCxN/log/pcxn_shrk_over${OVER}_estimates_rnd_%a.err \
  --array=1-4 \
  pcxn_over_estimates04.sh
  # pause between each job submission
  sleep 1
done
