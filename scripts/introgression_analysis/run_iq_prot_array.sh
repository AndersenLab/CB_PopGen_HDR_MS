source activate iqtree3

list=$1
file=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$list")

dir=$(dirname "$file")
base=$(basename "$file")

cd "$dir" || exit 1

iqtree -s "$base" -m LG -bb 1000 -T 8
