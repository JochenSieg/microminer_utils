work_dir="results/mutation_benchmark"
max_cpus=1
run_on_hpc="" # true is HPC should be used. "" if not

mkdir -p "${work_dir}"

# MICROMINER SEARCH EXPERIMENTS

# do MicroMiner searches against the PDB with wild-type protein structures
search_dir="${work_dir}/search"
mkdir -p "${search_dir}"

run_search() {
  local dataset="$1"
  local backward="$2"
  local do_eval="$3"
  local repr="$4"
  result_dir="${search_dir}/${dataset}"${backward:+"_backward"}
  mkdir -p "${result_dir}"
  python search.py --dataset "${work_dir}/${dataset}${backward:+"_backward"}.tsv" \
    --cpus ${max_cpus} \
    --outdir "${result_dir}" \
    --representation "${repr}" \
    ${run_on_hpc:+--hpc} || { return 1; }

  if [ "$do_eval" = true ]; then
    report_dir="${result_dir}/report"
    mkdir -p "${report_dir}"
    echo "python eval_known_mutations.py --csv ${result_dir} --dataset ${dataset} --outdir ${report_dir} ${backward:+--backward}"
    python eval_known_mutations.py --csv "${result_dir}" --dataset "${dataset}" \
      --outdir "${report_dir}" ${backward:+--backward}
  fi
}

# benchmark: try to re-find known structure pairs of wild-type/mutant with MicroMiner
python create_dataset.py \
  --dataset protherm thermomutdb platinum shanthirabalan \
  -o "${work_dir}"
run_search "protherm" "" true "monomer"
run_search "thermomutdb" "" true "monomer"
run_search "shanthirabalan" "" true "monomer"
run_search "platinum" "" true "full_complex"
