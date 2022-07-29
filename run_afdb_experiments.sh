work_dir="results/afdb_experiments"
max_cpus=1    # cpus per experiment.
run_on_hpc="" # set true to use HPC. "" if not

mkdir -p "${work_dir}"

# MICROMINER SEARCH EXPERIMENTS

# do MicroMiner searches
search_dir="${work_dir}/search"
mkdir -p "${search_dir}"

if [ "$run_on_hpc" = true ] ; then
  echo "About to run on HPC. Copying necessary files."
  local_repo_dir=$(dirname "$0")
  echo "$local_repo_dir"
  rm -rf /work/sieg/delme/microminer_evaluation ;
  git clone ${USER}@localhost:/local/sieg/projekte/microminer_evaluation /work/sieg/delme/microminer_evaluation
fi

run_search() {
  local dataset="$1"
  result_dir="${search_dir}/${dataset}"
  mkdir -p "${result_dir}"
  python search.py --dataset "${work_dir}/${dataset}.tsv" \
    --cpus ${max_cpus} \
    --outdir "${result_dir}" \
    --mode "standard" \
    ${run_on_hpc:+--hpc} || { return 1; }
}

python create_dataset.py --dataset afdb -o "${work_dir}"
run_search "afdb"


# you can do do something like this to collect results
#head -n1 /scratch/sieg/microminer_distributed/20220817_160600z1o51oit/results/100D/resultStatistic.csv > /local/sieg/pdb_all.tsv && find /scratch/sieg/microminer_distributed/20220817_160600z1o51oit/results/ -name resultStatistic.csv -exec tail -n+2 -q {} >> /local/sieg/pdb_all.tsv +
