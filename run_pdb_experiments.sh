max_cpus=800    # cpus per experiment.
run_on_hpc=true # set true to use HPC. "" if not
mm_repr="ppi"  # representation mode for MicroMiner
work_dir="results/pdb_experiments_${mm_repr}"

mkdir -p "${work_dir}"

# MICROMINER SEARCH EXPERIMENTS

# do MicroMiner searches against the PDB with wild-type protein structures
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
  local repr="$2"
  result_dir="${search_dir}/${dataset}"
  mkdir -p "${result_dir}"
  python search.py --dataset "${work_dir}/${dataset}.tsv" \
    --cpus ${max_cpus} \
    --outdir "${result_dir}" \
    --representation "${repr}" \
    ${run_on_hpc:+--hpc} || { return 1; }
}

python create_dataset.py --dataset pdb -o "${work_dir}"
run_search "pdb" ${mm_repr}


# you can do do something like this to collect results
#head -n1 /scratch/sieg/microminer_distributed/20220817_160600z1o51oit/results/100D/resultStatistic.csv > /local/sieg/pdb_all.tsv && find /scratch/sieg/microminer_distributed/20220817_160600z1o51oit/results/ -name resultStatistic.csv -exec tail -n+2 -q {} >> /local/sieg/pdb_all.tsv +
