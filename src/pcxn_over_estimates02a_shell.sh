
# BASH script to pass an argument to select the overlap case to get the estimates
# for the random gene sets
# The argument determines
# - job name
# - standard output file name
# - standard error file name
# - R output file name

for OVER in $(seq 8 10); do
  echo "Overlap Case ${OVER}"
  export OVER
  sbatch --job-name=pcxn_over${OVER}_estimates_rnd \
  -o /net/hsphfs1/srv/export/hsphfs1/share_root/hide_lab/PCxN/log/pcxn_over${OVER}_estimates_rnd_%a.out \
  -e /net/hsphfs1/srv/export/hsphfs1/share_root/hide_lab/PCxN/log/pcxn_over${OVER}_estimates_rnd_%a.err \
  --array=4-1000 \
  pcxn_over_estimates02a.sh
  # pause between each job submission
  sleep 1
done
