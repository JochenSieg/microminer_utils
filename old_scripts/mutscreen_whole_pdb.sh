
# fixed paths
exe="/local/sieg/projekte/pps/bin/MutScreen_release"
pdb_data_dir="/data/pdb/current/data/structures/all/pdb/" # filename: pdb2lzm.ent.gz
sienadb="/local/sieg/databases/siena_local.db"
complexdb="/local/sieg/databases/pdb.db"

# result paths
base_result_dir="/scratch/sieg/mutscreen_whole_pdb"
result_dir="${base_result_dir}/results"
log_dir="${base_result_dir}/log"
mkdir -p ${result_dir}
mkdir -p ${log_dir}


for pdb in ${pdb_data_dir}/*; do

  name="$(basename ${pdb})"
  name=${name:3:4} # get substring in basename of the PDB id. Name like pdb2lzm.ent.gz
  this_result_dir="${result_dir}/${name}"
  mkdir -p ${this_result_dir}

  { time ${exe} \
     --input ${pdb} \
     --sienadb ${sienadb} \
     --complexdb ${complexdb} \
     -q \
     -o ${this_result_dir} \
     ; } \
     2> ${log_dir}/${name}.err 1> ${log_dir}/${name}.out

done




