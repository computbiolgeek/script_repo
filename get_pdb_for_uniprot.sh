#!/bin/bash

uniprot_id=$1
id_mapping=$2
prefix=$3

for pdb_id in `grep -w ${uniprot_id} ${id_mapping} | grep PDB | awk '{print $3}'`; do
    curl https://files.rcsb.org/view/${pdb_id}.pdb -o ${prefix}_${pdb_id}.pdb
done
