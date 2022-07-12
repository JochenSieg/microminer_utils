set -e

work_dir="fold_experiments"
max_cpus=6    # cpus per experiment. Total number of parallel jobs will be max_cpus * experiments
run_on_hpc="" # TODO currently not supported

mkdir -p "${work_dir}"
#cd "${work_dir}" || exit 1

run_search() {
  local dataset="$1"
  result_dir="${work_dir}/${dataset}"
  mkdir -p "${result_dir}"
  python search.py --dataset "${work_dir}/${dataset}.tsv" \
    --cpus "${max_cpus}" \
    --outdir "${result_dir}" \
    ${run_on_hpc:+--hpc} || { return 1; }
}

# MICROMINER SEARCH SCOPe (https://scop.berkeley.edu)
export SITE_SEARCH_DB="/local/sieg/projekte/microminer/test/scope40_index.db"

dataset="scope"

python create_dataset.py \
  --dataset "${dataset}" -o "${work_dir}"
run_search "${dataset}"

python create_dataset_structure_pairs.py --csv "${work_dir}/${dataset}" \
  -o "${work_dir}/${dataset}_microminer_pdb_pairs.tsv" \
  --pdb_mirror1 "${dataset}" \
  --pdb_mirror2 "${dataset}"

python tmalign.py --dataset "${work_dir}/${dataset}_microminer_pdb_pairs.tsv" \
  --cpus "${max_cpus}" \
  -o "${work_dir}/${dataset}/"   #TODO ist das Dir mit dataset-namen nicht zu viel
tmalign_tsv="${work_dir}/${dataset}_microminer_pdb_pairs_tmalign.tsv"
if [ ! -f "${tmalign_tsv}" ]; then
  echo "Failed to generate TMalign similarities!"
  exit 1
fi

## make all MicroMiner result CSVs to a single CSV file
#big_csv="${work_dir}/${dataset}/big_resultStatistic.csv"
#find "${work_dir}/${dataset}" -name "resultStatistic.csv" -print0 -quit | xargs -0 head -n1 \
#  >"${big_csv}" # writes header to new CSV file
#find "${work_dir}/${dataset}" -name "resultStatistic.csv" -exec awk 'FNR > 1' {} + >>"${big_csv}"
