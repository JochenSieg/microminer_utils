#! /bin/bash
node_index_dir_suffix="{user}/microminer_distributed_cache"
node_index_dir=""
if [ -d "/ssd_local" ]; then
    node_index_dir="/ssd_local/${{node_index_dir_suffix}}"
elif [ -d "/local" ]; then
    node_index_dir="/local/${{node_index_dir_suffix}}"
else
    echo "Node has no expected drive! Abort"
    exit 1
fi
mkdir -p ${{node_index_dir}}

rsync -ra {user}@{host}:{kmer_index} ${{node_index_dir}}  