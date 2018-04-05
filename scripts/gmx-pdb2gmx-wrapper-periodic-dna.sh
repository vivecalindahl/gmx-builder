#!/bin/bash

# ==========================================================================
# Read input
# ==========================================================================
function print_usage
{
  echo "Usage: $0 -f <.pdb> -water <enum> -ff <string> [-h]"
}
args="$@"
old_gmx_is_ok=false;
while [ $# -gt 0 ]; do 
  case "$1" in    
    -h)       print_usage; exit 0;;
    -f)       pdb_in="$2"; shift;;
    -ff)      forcefield="$2"; shift;;
    -water)   water="$2"; shift;;
    -try)   old_gmx_is_ok=true; shift;;
    *)    echo "invalid option: $1" && exit 1;;
  esac 
  shift 
done

# Sanity and existence checks
[ -z "$pdb_in" ] && print_usage && exit
[ ! "${pdb_in##*.}" == "pdb" ] && echo "$pdb_in is not a pdb-file" && print_usage && exit 1 
[ ! -e "$pdb_in" ] && echo "$pdb_in does not exist" && exit 1 

[ -z "$forcefield" ] && print_usage && exit 
[ -z "$water" ] && print_usage && exit 


# 'gmx' command is assumed to be available but can be changed to any binary here.
gmx="gmx"
[ -z "$(which $gmx)" ]  && { $gmx; exit 1; }

# This script is intended for gmx version >= 2018.1  but it's not always strictly necessary.
gmx_major=$($gmx --version | sed  -n 's/GROMACS version: *\([0-9]*\)\.*\([0-9]*\).*/\1/p')
gmx_minor=$($gmx --version | sed  -n 's/GROMACS version: *\([0-9]*\)\.*\([0-9]*\).*/\2/p')

[ -z ${gmx_major} ] && { "Failed to read gmx version"; exit 1; }
min_major=2018
min_minor=1

old_gmx=false
if  [ ${gmx_major} -eq ${min_major} ]; then
    if   [ -z ${gmx_minor} ] || [ ${gmx_minor} -lt ${min_minor} ]; then
	old_gmx=true
    fi
elif [ ${gmx_major} -lt ${min_major} ]; then
    old_gmx=true
fi

if $old_gmx && ! $old_gmx_is_ok; then
    echo "gmx version ${gmx_major}.${gmx_minor} is less than the minimum recommended gmx version ${min_major}.${min_minor}."
    echo "Use '-try' flag to try anyway."
    exit 1
elif $old_gmx; then
    echo "Warning: gmx version ${gmx_major}.${gmx_minor} is less than the minimum recommended gmx version ${min_major}.${min_minor}."\
         "Proceeding anyway."
fi

# ==========================================================================
# The actual work:
# create a periodically connected molecule topology for DNA
# Output: topol.top and conf.gro files.
# ==========================================================================
# Notes:
# 1) Likely this is easy to generalize to e.g. RNA.
# 2) Tested for charmm27 and parmbsc1 (amber99bsc1.ff)
# 3) It's assumed there is nothing but the DNA molecule in the pdb.

# The strategy: copy the first residue of each chain to the end of each n-residue long chain,
# making a temporary new n+1 long chain. This way we can "trick" gmx pd2gmx to
# generate all the interaction terms we need, and then we delete what is not needed.

# The output
gro_out="conf.gro"
top_out="topol.top"

# A pdb to work with
pdb_work="./work.pdb"

# and a file for temporary output
tmp_log="./tmp.log"

# Contents of columns in pdb file format
chain_col=5;
resid_col=6;

# Add residue to the work pdb, for each chain.
awk -v chain_col=$chain_col -v resid_col=6 \
'{chain=$chain_col;  resid=$resid_col;
if (chain != prev_chain){if (str != ""){print str; str=""}; first_res=resid;};
if ($1 ~ /^ATOM/ && resid==first_res){if(str != ""){str=str"\n"$0}else{str=$0};};
prev_chain=chain; print}' \
$pdb_in  > $pdb_work

# Generate topology and config gro-file for the n+1 chain.

# We need a local copy of the forcefield that we can modify.
$gmx pdb2gmx -f $pdb_in -ff $forcefield -water $water &> $tmp_log
forcefield_dir=$(awk '/Opening force field file/{gsub(/.ff.*/,".ff", $5);print $5; exit}' ${tmp_log})
[ -z "$forcefield_dir" ] && echo "Could not extract force field location from pdb2gmx. See ${tmp_log}." && exit 1

ff_work="${forcefield}_work"
ff_work_dir="./${ff_work}.ff"
cp -rL $forcefield_dir $ff_work_dir

# Remove special treatment of the end termini in the .r2b file if there is one.
r2b="${ff_work_dir}/dna.r2b"
if [ -e "$r2b" ]; then
    mv $r2b ${r2b}"_tmp"
    awk '/;/{print; next}; {if (NF==5){print $1, $2, $2, $2, $2} else{print}}' ${r2b}"_tmp" > $r2b
    rm ${r2b}"_tmp"
fi

# Need to add a hack to the terminal data base file for pdb2gmx to accept the ends
# without modification. To avoid having to parse the interactive session of pdb2gmx,
# remove all .tdb entries other than those we define, so there is only one choice
# for each end (here we assume there is only the periodic molecule in the pdb)

# Make a custom tdb-file for the 5' end to force gmx pdb2gmx to allow end residue
# that have (temporarily) missing bonds.
# This hack (loophole?) will avoid ax gmx pdb2gmx fatal error. 
rm -f ${ff_work_dir}/*.n.tdb  # 5' end
rm -f ${ff_work_dir}/*.c.tdb  # 3' end

# For a gmx version < 2018.1 the pdb2gmx '-missing' flag was more strict requiring
# this hack to avoid a fatal error.
tdb5="${ff_work_dir}/dna.n.tdb" 
if $old_gmx; then 
    if [ "${forcefield}" == 'amber99bsc1' ]; then
	echo -e "[ hack ]\n[ replace ]\nC5' C5' CI 12.01 -0.0069" > $tdb5
    elif [ "${forcefield}" == 'charmm27' ]; then
	echo -e "[ hack ]\n[ replace ]\nC5' C5' CN8B 12.011 -0.08" > $tdb5
    else
	echo "Forcefield $forcefield not supported."; exit 1
    fi
else
    echo -e "[ none ]" > $tdb5
fi

# Put a dummy entry "none" on the 3' side.
tdb3="${ff_work_dir}/dna.c.tdb" 
echo -e "[ none ]" > $tdb3

# Generate .top and .gro files from the unmodified pdb.
# The .gro file will be compatible with the periodic .top file.
# The .top file will be used a reference when modifying the n+1 
# topology to keep track of the start and end (atoms) of each
# molecule chain.
gro_nonperiodic="${gro_out}"
top_nonperiodic="nonperiodic.top"
$gmx pdb2gmx -v -missing -f $pdb_in -ff $ff_work -water $water -ter -o  $gro_nonperiodic -p $top_nonperiodic  &> $tmp_log || \
    { echo "gmx pdb2gmx exited with an error. See ${tmp_log}."; exit 1; }

# List with the non-periodic topology (.itp) files of each chain
itp_nonperiodic_list=(`grep  "chain.*itp" ${top_nonperiodic} | awk '{print $NF}' | awk -F  \" '{print $2}'`)

# Extract the last atom index of the unmodified topology, for each chain.
max_atom_index_list=()
for itp in ${itp_nonperiodic_list[@]}; do
    name='atoms'
    max_atom_index=`awk -v  name=$name 'BEGIN{found=0;};/^ *\[ /{ if ($2 == name){found=1; getline;} else { found=0 }};{if (found && $1 != ";"  && $1 !~ /^ *#/ && NF>0){max_atom_index=$1};}END{print max_atom_index}' $itp`
    max_atom_index_list+=($max_atom_index)
done

# Generate topology (.top, .itp) of n+1 config. The toplogy files will be modified below.
# The configuration file is not needed.
top_periodic="periodic.top"
gro_dummy="dummy.gro"
$gmx pdb2gmx -v -missing -f $pdb_work -ff $ff_work -water $water -ter -o $gro_dummy -p $top_periodic  &> $tmp_log || \
    { echo "gmx pdb2gmx exited with an error. See ${tmp_log}."; exit 1; }

# The .itp files included in the to-be periodic .top file,  one per chain.
itp_periodic_list=(`grep  "chain.*itp" "${top_periodic}" | awk '{print $NF}' | awk -F  \" '{print $2}'`)

# Transform the n+1 topology into the periodic topology by
# modifying the .itp file corresponding to each molecule chain.
for chain in ${!itp_periodic_list[@]}; do
    itp_in=${itp_periodic_list[chain]}
    max_atom_index=${max_atom_index_list[chain]}

    # Overwrite the input
    itp_out=$itp_in
    
    # Array with the number of atoms listed for each type of interaction
    declare -A natoms
    natoms['bonds']=2
    natoms['angles']=3
    natoms['dihedrals']=4
    natoms['pairs']=2
    natoms['atoms']=1

    # Modify the topology.
    # For each gmx itp-directive (interaction type), 
    # 1) remove entries with atom indices only containing the added residue n+1
    # 2) in "mixed" entries containing the connections between residues n and n+1,
    #    remap atom indices of residue n+1 to the first residue, i.e:
    #    i --> i + max_atom_index, if i > max_atom_index,  where max_atom_index = max index of the chain.
    itp_work="./work.itp"

    cp $itp_in $itp_work
    for name in ${!natoms[@]}; do
	natoms=${natoms[$name]}
	# 1) Don't print interactions where the all involved atoms are in residue  n+1.
	# 2) Modify (remap) interactions where some involved atoms are in residue  n+1.
	# 3) Otherwise, print as is.
	awk -v max_atom_index=$max_atom_index -v name=$name -v natoms=$natoms 'BEGIN{found=0};
/^ *\[ /{ if ($2 == name){found=1; print; getline;} else { found=0 }};
{ if (found && $1 != ";"  && $1 !~ /^ *#/ && NF>0){ count=0; for (j=1; j<=natoms; j++){if ($j >max_atom_index){$j=($j-max_atom_index); count++}}; if (count<natoms){print}}
else{print} }' $itp_work > $itp_out

	cp $itp_out $itp_work
    done; 
done

# Rename .top file to the output name.
mv $top_periodic $top_out

# Merge top and itps into a single topology fil by replacing the 
# include statements .itp files with their contents.
for itp in ${itp_periodic_list[@]}; do
    # sed r command appends the itp file after the matching include statement
    sed -i "/^ *\#include *\"$itp/r  $itp" ${top_out};

    # delete the include statement
    sed  -i "s/^ *\#include *\"$itp.*//g"  ${top_out};
done

# Final tweaks.

# Fix the force field entry so it doesn't map to the modified work version.
sed -i "s/${ff_work}/${forcefield}/" ${top_out}

# If we found the forcefield directory somewhere else,
# the include statement should not contain a relative path.
nonlocal_ff=true
[[  "$forcefield_dir" =~ ^\..* ]] && nonlocal_ff=false
if $nonlocal_ff; then
    sed -i "s/.\/${forcefield}/${forcefield}/" ${top_out}
fi

# Add comment about how the topology file was generated at the top.
# Below, '1s' refers to the first position in the file and '~' is used as the sed delimiter
# since the variables may contain '/'.
sed -i "1s~^~; >>>>>>>>>> NOTE: This file was automatically generated using \"$0 $args\"~" ${top_out}

# Clean up.
# Remove all the work files that were generated except the output files (top_out, gro_out).
rm -rf ${gro_dummy} \
    ${itp_periodic_list[@]} \
    ${itp_nonperiodic_list[@]} ${top_nonperiodic} \
    *posre_DNA*.itp* \
    ${pdb_work} ${itp_work} ${tmp_log} \
    ${ff_work_dir}
