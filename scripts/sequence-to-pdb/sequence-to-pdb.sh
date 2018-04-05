#!/bin/bash

# Use x3dna (http://x3dna.org/) to generate pdbs from a DNA sequence

# Untar x3dna if needed
x3dna_version="2.3"
x3dna_tar=./tools/x3dna-v${x3dna_version}-linux-64bit.tar.gz
[ ! -e "$x3dna_tar" ] && echo "$x3dna_tar not found. Aborting." && exit 1

out_dir="./build"
x3dna_dir="./${out_dir}/x3dna-v${x3dna_version}"
[ ! -e "$x3dna_dir" ] && mkdir -p $out_dir && echo "Bulding x3dna..." && tar pzxvf "$x3dna_tar" -C "./${out_dir}" &> /dev/null
[ ! -e "$x3dna_dir" ] && echo "$x3dna_dir not found. Aborting." && exit 1

# Environment variable required to set by x3dna
export X3DNA="$x3dna_dir"

# The fiber tool builds the pdbs ('-b' flag gives type B DNA).
fiber=${x3dna_dir}/bin/fiber

# Input data file. One sequence per line.
sequences="./sequences.txt"
[ ! -e "$sequences" ] && echo "$sequences not found. Aborting." && exit 1

pdb_dir="${out_dir}/pdbs"
mkdir -p $pdb_dir
echo  "Building pdb files with sequences from $sequences"
for seq in `cat ${sequences}`; do
    $fiber -seq=${seq} -b ${pdb_dir}/${seq}.pdb || { echo "Something is probably wrong with sequence '$seq' found in $sequences" && break; };
done
echo "pdb files written to $pdb_dir"
