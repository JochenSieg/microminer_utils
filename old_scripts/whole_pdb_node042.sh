#! /bin/bash
#$ -N mut_pdb
#$ -cwd
#$ -l hostname=node042
# -l notshared
#$ -q 64c.q
#$ -t 1-55
#$ -tc 55
# -l mem_total=350G
# -l mem_free=300G
#$ -o /scratch/sieg/mutscreen_whole_pdb/cluster_out
#$ -e /scratch/sieg/mutscreen_whole_pdb/cluster_out

#set -e

N=1

 # source the conda env
pdb_data_dir="/data/pdb/current/data/structures/all/pdb/" # filename: pdb2lzm.ent.gz
scripts_src_dir="/work/sieg/delme/tool_eval"
exe="/work/sieg/packages/MutScreen_0.1/MutScreen"
local_sienadb="/ssd_local/sieg/databases/siena_local.db"

# result paths
base_result_dir="/scratch/sieg/mutscreen_whole_pdb"

# make my own working dir because of sym links
wd="/local/sieg/tmp_jobs/${JOB_ID}"
mkdir -p "${wd}"
cd "${wd}" || exit
# TODO das ist ein hotfix:
ln -s /ssd_local/sieg/data/pdb/sym sym 2>/dev/null
local_result_dir="${wd}/results"
mkdir -p "${local_result_dir}"

#for pdb in ${pdb_data_dir}/*; do

pdb_files=("${pdb_data_dir}"/*.ent.gz)
nof_pdbfiles=${#pdb_files[@]}

NOF_JOBS=${SGE_TASK_LAST}
JOB_ID=${SGE_TASK_ID}
chunksize=$((nof_pdbfiles / NOF_JOBS))
chunksize=$((chunksize + 1))
min_idx=$((JOB_ID - 1))
min_idx=$((min_idx * chunksize))
max_idx=$((JOB_ID * chunksize))
max_idx=$((max_idx < nof_pdbfiles ? max_idx : nof_pdbfiles))

for ((i=${min_idx}; i<${max_idx}; i+=1)) do
  pdb="${pdb_files[${i}]}"
  name="$(basename ${pdb})"
  name=${name:3:4} # get substring in basename of the PDB id. Name like pdb2lzm.ent.gz
  this_result_dir="${local_result_dir}/${name}"
  mkdir -p "${this_result_dir}"

  { time ${exe} \
     search \
     -q "${pdb}" \
     -s ${local_sienadb} \
     -o "${this_result_dir}" \
     ; } \
     2> "${this_result_dir}"/"${name}".err 1> "${this_result_dir}"/"${name}".out
  echo "$?" > "${this_result_dir}"/status
done

# clean up
cp -r "${wd}" ${base_result_dir}
#cd ..
#rm -r "${wd}"
