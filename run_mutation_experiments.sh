work_dir="mutation_experiments"
max_cpus=3   # cpus per experiment. Total number of parallel jobs will be max_cpus * experiments
run_on_hpc=""   # true is HPC should be used. "" if not

mkdir -p "${work_dir}"
#cd "${work_dir}" || exit 1

# MICROMINER PAIR EXPERIMENTS

pair_dir="${work_dir}/pair"
mkdir -p "${pair_dir}"

run_pair () {
  local dataset="$1"
  local backward="$2"
  result_dir="${pair_dir}/${dataset}"${backward:+"_backward"}
  mkdir -p "${result_dir}"
  python pair.py --dataset "${work_dir}/${dataset}_pair${backward:+"_backward"}.tsv" \
                 --cpus ${max_cpus} \
                 --outdir "${result_dir}" \
                 ${run_on_hpc:+--hpc} || { return 1; }

  report_dir="${result_dir}/report"
  mkdir -p "${report_dir}"
  python eval_known_mutations.py --csv "${result_dir}" --dataset "${dataset}" \
                                 --outdir "${report_dir}" ${backward:+--backward}
}
#python create_dataset_mutation_pairs.py --dataset protherm thermomutdb platinum shanthirabalan -o "${work_dir}"
#run_pair "protherm" &
#run_pair "thermomutdb" &
#run_pair "platinum" &
#run_pair "shanthirabalan" &
#
#python create_dataset_mutation_pairs.py --dataset protherm thermomutdb platinum -o "${work_dir}" \
#                                        --backward
#run_pair "protherm" true &
#run_pair "thermomutdb" true &
#run_pair "platinum" true &

# MICROMINER SEARCH EXPERIMENTS

# do MicroMiner searches against the PDB with wild-type protein structures
search_dir="${work_dir}/search"
mkdir -p "${search_dir}"

run_search () {
  local dataset="$1"
  local backward="$2"
  local do_eval="$3"
  result_dir="${search_dir}/${dataset}"${backward:+"_backward"}
  mkdir -p "${result_dir}"
  python search.py --dataset "${work_dir}/${dataset}${backward:+"_backward"}.tsv" \
                   --cpus ${max_cpus} \
                   --outdir "${result_dir}" \
                   ${run_on_hpc:+--hpc} || { return 1; }

  if [ "$do_eval" = true ]; then
    report_dir="${result_dir}/report"
    mkdir -p "${report_dir}"
    python eval_known_mutations.py --csv "${result_dir}" --dataset "${dataset}" \
                                   --outdir "${report_dir}" ${backward:+--backward}
  fi
}

python create_dataset.py \
  --dataset protherm prothermdb thermomutdb skempi2 platinum shanthirabalan -o "${work_dir}"
run_search "protherm" "" true
run_search "thermomutdb" "" true
run_search "shanthirabalan" "" true
run_search "platinum" "" true
# mutation data set without known mutant structures:
#run_search "prothermdb" &
#run_search "skempi2" &
#
#python create_dataset.py \
#  --dataset protherm thermomutdb platinum -o "${work_dir}" --backward
#run_search "protherm" true true &
#run_search "thermomutdb" true true &
#run_search "platinum" true true &

