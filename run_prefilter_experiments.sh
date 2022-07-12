# Evaluates the capabilities of the prefilter.

set -e

work_dir="prefilter_experiments"
max_cpus=5    # cpus per experiment. Total number of parallel jobs will be max_cpus * experiments
run_on_hpc="" # TODO currently not supported

mkdir -p "${work_dir}"
#cd "${work_dir}" || exit 1

run_prefilter() {
  local dataset="$1"
  result_dir="${work_dir}/${dataset}"
  mkdir -p "${result_dir}"
  python prefilter.py --dataset "${work_dir}/${dataset}.tsv" \
    --cpus "${max_cpus}" \
    --outdir "${result_dir}" \
    ${run_on_hpc:+--hpc} || { return 1; }
}

# MICROMINER SEARCH SCOPe (https://scop.berkeley.edu)
#export SITE_SEARCH_DB="/local/sieg/projekte/microminer/test/scope40_index.db"
#export SITE_SEARCH_DB="/local/sieg/projekte/microminer/test/intkmer_scope40_index.db"
export SITE_SEARCH_DB="/local/sieg/projekte/microminer/test/intkmer_scope40_k6_index.db"

dataset="scope"

python create_dataset.py \
  --dataset "${dataset}" -o "${work_dir}"
run_prefilter "${dataset}"

#python create_dataset_structure_pairs.py --csv "${work_dir}/${dataset}" \
#  -o "${work_dir}/${dataset}_microminer_pdb_pairs.tsv" \
#  --pdb_mirror1 "${dataset}" \
#  --pdb_mirror2 "${dataset}"
#
#python tmalign.py --dataset "${work_dir}/${dataset}_microminer_pdb_pairs.tsv" \
#  --cpus "${max_cpus}" \
#  -o "${work_dir}/" \
#  --skip_same # Skip same pairs of same ID and same filepath
#tmalign_tsv="${work_dir}/${dataset}_microminer_pdb_pairs_tmalign.tsv"
#if [ ! -f "${tmalign_tsv}" ]; then
#  echo "Failed to generate TMalign similarities!"
#  exit 1
#fi

python eval_scope.py --microminer "${work_dir}/scope/" --tmalign "${tmalign_tsv}" \
                     -o "${work_dir}/${dataset}/"
