#!/bin/bash
# set bash source
DIR="${BASH_SOURCE%/*}"
if [[ ! -d "$DIR" ]]; then DIR="$PWD"; fi
. "$DIR/run_aux.sh"

# set machine ID
set_machine_ID

# set program name
prog_name=PB3D

# Display usage function
display_usage() { 
    echo -e "Usage:\n$0 [OPTS] CASE NR_PROCS \n" 
    print_opts
    echo -e ""
    echo -e "    CASE:       1 cbm18a"
    echo -e "                2 cbm18a_HELENA"
    echo -e ""
    echo -e "               11 cbm18a_ripple_1"
    echo -e "               12 cbm18a_ripple_2"
    echo -e "               13 cbm18a_ripple_5"
    echo -e "               14 cbm18a_ripple_18"
    echo -e "               15 cbm18a_ripple_flat_R"
    echo -e "               16 cbm18a_ripple_flat_Z"
    echo -e "               17 cbm18a_ripple_flat_RZ"
    echo -e "               18 cbm18a_ripple_inverse"
    echo -e ""
    echo -e "               21 cdxu"
    echo -e ""
    echo -e "               31 qps"
    echo -e ""
    echo -e "               41 Hmode_ped4"
    echo -e "               42 Hmode_ped4_HELENA"
    echo -e "               43 Hmode_ped2"
    echo -e "               44 Hmode_ped2_HELENA"
    echo -e "               45 Hmode_ped1.5"
    echo -e "               46 Hmode_ped1.5_HELENA"
    echo -e "               47 Hmode_ped1"
    echo -e "               48 Hmode_ped1_HELENA"
    echo -e "               49 Hmode_ped0.8"
    echo -e "               50 Hmode_ped0.8_HELENA"
    echo -e "               51 Hmode_ped0.7"
    echo -e "               52 Hmode_ped0.7_HELENA"
    echo -e "               53 Hmode_ped0.6"
    echo -e "               54 Hmode_ped0.6_HELENA"
    echo -e "               55 Hmode_ped0.5"
    echo -e "               56 Hmode_ped0.5_HELENA"
    echo -e ""
    echo -e "    NR_PROCS:  nr. of MPI processes"
    echo -e ""
} 

# Setting some variables
init_vars
#slepc_opt="-st_pc_factor_shift_type NONZERO -st_pc_type lu -st_pc_factor_mat_solver_package mumps -eps_monitor -eps_view"
slepc_opt="-st_pc_type lu -st_pc_factor_mat_solver_package mumps -eps_monitor -eps_ncv 100 -eps_mpd 100"
#slepc_opt="-st_pc_type jacobi -st_pc_factor_mat_solver_package mumps -eps_monitor -eps_ncv 100 -eps_mpd 100"

# Catch options
catch_options $@

# Shift arguments to skip options
shift $n_opt_args

# Check number of input arguments
# (need at least 1 for # procs and 1 for case)
if [[ $# -lt 2 ]]; then
    display_usage
    exit 1
fi

# Set parameters: input_name, eq_name
case $1 in
    1|2)
        input_name=cbm18a
    ;;
    11|12|13|14|15|16|17|18)
        input_name=cbm18a_ripple
    ;;
    21)
        input_name=cdxu
    ;;
    31)
        input_name=qps
    ;;
    41|42|43|44|45|46|47|48|49|50|51|52|53|54|55|56)
        input_name=Hmode
    ;;
    *)
        echo -e "ERROR: No case found"
        display_usage
        exit 1
    ;;
esac
input_name=input_$input_name
case $1 in
    1)
        eq_name=wout_cbm18a.nc
    ;;
    2)
        eq_name=cbm18a
    ;;
    11)
        eq_name=wout_cbm18a_ripple.nc
    ;;
    12)
        eq_name=wout_cbm18a_ripple_2.nc
    ;;
    13)
        eq_name=wout_cbm18a_ripple_5.nc
    ;;
    14)
        eq_name=wout_cbm18a_ripple_18.nc
    ;;
    15)
        eq_name=wout_cbm18a_ripple_flat_R.nc
    ;;
    16)
        eq_name=wout_cbm18a_ripple_flat_Z.nc
    ;;
    17)
        eq_name=wout_cbm18a_ripple_flat_RZ.nc
    ;;
    18)
        eq_name=wout_cbm18a_ripple_inverse.nc
    ;;
    21)
        eq_name=wout_cdxu.nc
    ;;
    31)
        eq_name=wout_qps.nc
    ;;
    41)
        eq_name=wout_Hmode_ped4.nc
    ;;
    42)
        eq_name=Hmode_ped4
    ;;
    43)
        eq_name=wout_Hmode_ped2.nc
    ;;
    44)
        eq_name=Hmode_ped2
    ;;
    45)
        eq_name=wout_Hmode_ped1.5.nc
    ;;
    46)
        eq_name=Hmode_ped1.5
    ;;
    47)
        eq_name=wout_Hmode_ped1.nc
    ;;
    48)
        eq_name=Hmode_ped1
    ;;
    49)
        eq_name=wout_Hmode_ped0.8.nc
    ;;
    50)
        eq_name=Hmode_ped0.8
    ;;
    51)
        eq_name=wout_Hmode_ped0.7.nc
    ;;
    52)
        eq_name=Hmode_ped0.7
    ;;
    53)
        eq_name=wout_Hmode_ped0.6.nc
    ;;
    54)
        eq_name=Hmode_ped0.6
    ;;
    55)
        eq_name=wout_Hmode_ped0.5.nc
    ;;
    56)
        eq_name=Hmode_ped0.5
    ;;
    *)
        echo -e "ERROR: No case found"
        display_usage
        exit 1
    ;;
esac

# loop over all inputs
# (from http://www.cyberciti.biz/faq/unix-linux-iterate-over-a-variable-range-of-numbers-in-bash/)
for (( input_i=1; input_i<=$n_inputs; input_i++ )); do
    # setup output directory
    setup_output_dir $(date +"%Y-%m-%d-%H-%M-%S")

    # user output
    echo -e "Working in directory $out_full/"
    echo -e ""
    echo -e "input file:       " $input_name
    echo -e "equilibrium file: " $eq_name
    echo -e ""
    echo -e "Copying files"
    echo -e ""

    # Copy inputs and the program
    cp $input_name $out_full
    cp $eq_name $out_full
    cp ../$prog_name $out_full
    chmod +x $out_full/$prog_name
    cd $out_full
    
    # modify input file
    modify_input_file

    # set mpi command
    mpi_command="./$prog_name $input_name $eq_name $slepc_opt ${@:3}"
    nr_procs=$2

    # run mpi command, depending on machine ID
    run_mpi_command
    
    # finish
    finish
done
