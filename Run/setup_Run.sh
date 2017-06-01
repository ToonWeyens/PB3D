#!/bin/bash
echo -e "This script sets up the required files and scripts to run PB3D simulations."

# check if directory name was provided
if [[ $# -eq 1 ]]; then
    sim_dir=${1}
else
    echo -e "Where do you want to install?"
    read sim_dir
fi

# check for current directory
if [[ $sim_dir -ef $(pwd) ]]; then
    echo -e "Directory ${sim_dir%/} is the current directory"
    exit
fi
echo -e "Installing in directory ${sim_dir%/}"

# check for existence
if [[ ! -d "${sim_dir%/}" ]]; then
    echo -e "Directory ${sim_dir%/} does not exist. Creating it"
    mkdir -p ${sim_dir%/}
fi
echo -e ""

# save directory and go to simulation directory
PB3D_dir=$(pwd)
cd $sim_dir

# copy README
cp --interactive $PB3D_dir/README* .

# copy scripts
cp --interactive $PB3D_dir/*.sh .

# copy simulation inputs
cp --interactive $PB3D_dir/input_* .

# copy array inputs
cp --interactive $PB3D_dir/array_input .

# copy 3-D perturbation data files
cp --interactive $PB3D_dir/pert_*.dat .

# copy scorep files
cp --interactive $PB3D_dir/scorep.fil .

# copy VisIt expressions
cp --interactive $PB3D_dir/visit_expressions.xmf .

# change the PB3D directory in the run script
# (from https://unix.stackexchange.com/a/66989 and https://askubuntu.com/a/76842)
sed -i "s,^\( *prog_dir= *\)[^ ]*\(.*\)*$,\1"${PB3D_dir%/}/../"\2," run.sh

# done
echo -e ""
echo -e "done"
