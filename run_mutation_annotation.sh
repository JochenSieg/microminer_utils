work_dir="results/mutation_annotation"
max_cpus=1
run_on_hpc="" # true is HPC should be used. "" if not

mkdir -p "${work_dir}"

# MICROMINER SEARCH EXPERIMENTS

# do MicroMiner searches against the PDB with wild-type protein structures
search_dir="${work_dir}/search"
mkdir -p "${search_dir}"

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

# annotation: find structures for the mutant in the PDB
python create_dataset.py \
  --dataset protherm thermomutdb platinum fireprotdb prothermdb skempi2 \
  -o "${work_dir}"
run_search "protherm" "monomer"
run_search "thermomutdb" "monomer"
run_search "platinum" "full_complex"
# mutation data set without known mutant structures:
run_search "prothermdb" "monomer"
run_search "skempi2" "full_complex"
run_search "fireprotdb" "monomer"
