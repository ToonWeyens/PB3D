#!/bin/bash
# set bash source
DIR="${BASH_SOURCE%/*}"
if [[ ! -d "$DIR" ]]; then DIR="$PWD"; fi
. "$DIR/run_aux.sh"

# set machine ID
set_machine_ID

# set program name
prog_name=POST

# Display usage function
display_usage() { 
    echo -e "Usage:\n$0 [OPTS] PB3D_DIR NR_PROCS \n" 
    print_opts
    echo -e ""
    echo -e "    PB3D_DIR:  PB3D directory"
    echo -e ""
    echo -e "    NR_PROCS:  nr. of MPI processes"
    echo -e ""
} 

# Setting some variables
init_vars
PB3D_out_name="PB3D_out.h5"
input_name="input_POST"

# Catch options
catch_options $@

# Shift arguments to skip options
shift $n_opt_args

# Check number of input arguments
# (need at least 1 for # procs and 1 for PB3D directory)
if [[ $# -lt 2 ]]; then
    display_usage
    exit 1
fi

# Check existence of PB3D folder
if [[ -d "$1" ]] ; then
    if [[ ! -f "${1%/}/$PB3D_out_name" ]] ; then
        echo -e 'ERROR: PB3D output file does not exist.'
        display_usage
        exit 1
    fi
else
    echo -e 'ERROR: PB3D Directory does not exist.'
    display_usage
    exit 1
fi

# loop over all inputs
# (from http://www.cyberciti.biz/faq/unix-linux-iterate-over-a-variable-range-of-numbers-in-bash/)
for (( input_i=1; input_i<=$n_inputs; input_i++ )); do
    # setup output directory
    setup_output_dir ${1%/}

    # user output
    echo -e "Working in directory $out_full/"
    echo -e ""
    echo -e "input file:       " $input_name
    echo -e "PB3D output file: " ${1%/}/$PB3D_out_name
    echo -e ""
    echo -e "Copying files"
    echo -e ""

    # Copy inputs and the program
    cp $input_name $out_full
    [[ $use_out_loc = true ]] && cp ${1%/}/$PB3D_out_name $out
    cp ../$prog_name $out_full
    chmod +x $out_full/$prog_name
    
    # go to run directory
    cd $out_full
    
    # modify input file
    modify_input_file

    # set mpi command
    mpi_command="./$prog_name $input_name $PB3D_out_name ${@:3}"
    nr_procs=$2

    # run mpi command, depending on machine ID
    run_mpi_command

    # finish
    finish
done
