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
    echo -e "               41 Hmode_ped4.0"
    echo -e "               42 Hmode_ped4.0_HELENA"
    echo -e "               43 Hmode_ped2.0"
    echo -e "               44 Hmode_ped2.0_HELENA"
    echo -e "               45 Hmode_ped1.5"
    echo -e "               46 Hmode_ped1.5_HELENA"
    echo -e "               47 Hmode_ped1.0"
    echo -e "               48 Hmode_ped1.0_HELENA"
    echo -e "               49 Hmode_ped0.8"
    echo -e "               50 Hmode_ped0.8_HELENA"
    echo -e "               51 Hmode_ped0.7"
    echo -e "               52 Hmode_ped0.7_HELENA"
    echo -e "               53 Hmode_ped0.6"
    echo -e "               54 Hmode_ped0.6_HELENA"
    echo -e "               55 Hmode_ped0.5"
    echo -e "               56 Hmode_ped0.5_HELENA"
    echo -e ""
    echo -e "               61 Hmode_ped1.0_ripple16_0.5"
    echo -e "               62 Hmode_ped1.0_ripple16_1.0"
    echo -e "               63 Hmode_ped1.0_ripple16_1.5"
    echo -e "               64 Hmode_ped1.0_ripple16_2.0"
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
    41|42|43|44|45|46|47|48|49|50|51|52|53|54|55|56|61|62|63|64)
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
        eq_name=wout_Hmode_ped4.0.nc
    ;;
    42)
        eq_name=Hmode_ped4.0
    ;;
    43)
        eq_name=wout_Hmode_ped2.0.nc
    ;;
    44)
        eq_name=Hmode_ped2.0
    ;;
    45)
        eq_name=wout_Hmode_ped1.5.nc
    ;;
    46)
        eq_name=Hmode_ped1.5
    ;;
    47)
        eq_name=wout_Hmode_ped1.0.nc
    ;;
    48)
        eq_name=Hmode_ped1.0
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
    61)
        eq_name=wout_Hmode_ped1.0_ripple16_0.5.nc
    ;;
    62)
        eq_name=wout_Hmode_ped1.0_ripple16_1.0.nc
    ;;
    63)
        eq_name=wout_Hmode_ped1.0_ripple16_1.5.nc
    ;;
    64)
        eq_name=wout_Hmode_ped1.0_ripple16_2.0.nc
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
    # setup output name and home directory
    setup_out $(date +"%Y-%m-%d-%H-%M-%S")
    home=$(pwd)
    
    # set variables
    case $machine_ID in
        1)                                                                  # 1. XPS-L501X
            nr_procs=$2                                                     # take nr_procs from input
            out_full=$home/$out                                             # use subdirectory
        ;;
        2)                                                                  # 2. ITER
            nr_procs=1                                                      # only 1 MPI process
            n_nodes=1                                                       # use one node
            n_cores=4                                                       # 8 cores per node
            memory_factor=2                                                 # safety factor for requesting memory (needs to be integer)
            max_tot_mem=16000                                               # every core has 4GB, using 4
            mem_unit='mb'                                                   # MB
            out_full=/tmp/$out                                              # use temporary directory in local node
            out_full_loc=$home/$out                                         # afterwards copy to this directory
        ;;
    esac
    
    # user output
    echo -e "Working in directory '$out_full/'"
    echo -e ""
    echo -e "input file:       " $input_name
    echo -e "equilibrium file: " $eq_name
    echo -e ""
    
    # possibly modify copy of input file
    cp $input_name $input_name"_loc"
    modify_input_file $input_name"_loc"
    
    # set up local simulation script
    echo -e "Creating local simulation script in '"$prog_name"_loc.sh'"
    echo -e ""
    cat <<- END > $prog_name"_loc.sh"
        
        # create directory structure
        mkdir -p $out_full $out_full/Plots $out_full/Data $out_full/Scripts || {
            # failure
            echo "ERROR: unable to create directory $out_full/ and subdirectories"
            echo "Do you have the correct permissions?"
            exit 1
        }
        
        # Move/copy inputs and the program
        mv $out_full_loc/$input_name $out_full/$input_name 2> /dev/null
        cp $home/$eq_name $out_full
        cp $home/../$prog_name $out_full
        chmod +x $out_full/$prog_name
        
        # go to run directory
        cd $out_full
        
        # remove possible output and error file
        rm -f $prog_name.out $prog_name.err
        
        # set mpi command
        mpi_command="./$prog_name $input_name $eq_name $slepc_opt ${@:3} --minim_output"
        
        # run mpi command
        echo "$(echo $debug_opt mpirun | xargs) -np $nr_procs \$mpi_command" > command_$prog_name
        . command_$prog_name
END
    # cut 8 leading whitespaces and make executable
    sed -i 's/^.\{8\}//' $prog_name"_loc.sh"
    chmod u+x $prog_name"_loc.sh"
    
    if [[ $on_cluster = true ]]; then
        # user output
        echo -e "The output will be copied to '$out_full_loc' afterwards"
        echo -e ""
        
        # user output
        echo -e "Creating pbs script with"
        echo -e "    $n_nodes node(s) and $n_cores core(s) per node"
        echo -e "    $max_tot_mem$mem_unit memory in total"
        echo -e ""
        
        # initizalize pbs script
        init_pbs_script
        
        # create directory for output and errors
        mkdir -p $out_full_loc || {
            # failure
            echo "ERROR: unable to create directory $out_full_loc/"
            echo "Do you have the correct permissions?"
            exit 1
        }
        
        # move input file
        mv $input_name"_loc" $out_full_loc/$input_name
        
        # add commands to create symbolic link to output in home
        echo "" >> $prog_name".pbs"
        echo "# set local output and error and create symbolic links" >> $prog_name".pbs"
        echo "loc_out=\$PBS_O_HOME/$(echo $out | tr '/' '_').o\$(echo \$PBS_JOBID | cut -d'.' -f 1)" >> $prog_name".pbs"
        echo "loc_err=\$PBS_O_HOME/$(echo $out | tr '/' '_').e\$(echo \$PBS_JOBID | cut -d'.' -f 1)" >> $prog_name".pbs"
        echo "ln -s \$loc_out $out_full_loc/$prog_name.out" >> $prog_name".pbs"
        echo "ln -s \$loc_err $out_full_loc/$prog_name.err" >> $prog_name".pbs"
        
        # add local execution script to pbs script
        cat $prog_name"_loc.sh" >> $prog_name".pbs"
        
        # add commands to move output after finishing
        echo "" >> $prog_name".pbs"
        echo "# Done, copy files back" >> $prog_name".pbs"
        echo "cd $out_full_loc" >> $prog_name".pbs"
        echo "mv $out_full/* $out_full_loc/" >> $prog_name".pbs"
        echo "rmdir $out_full" >> $prog_name".pbs"
        echo "mv \$loc_out $out_full_loc/$prog_name.out" >> $prog_name".pbs"
        echo "mv \$loc_err $out_full_loc/$prog_name.err" >> $prog_name".pbs"
        echo "[[ -s $out_full_loc/$prog_name.err ]] || rm $out_full_loc/$prog_name.err" >> $prog_name".pbs"
        
        # complete pbs script
        echo "" >> $prog_name".pbs"
        echo "exit" >> $prog_name".pbs"
        
        # mv pbs script to output
        mv $prog_name".pbs" $out_full_loc/
        
        # user output
        echo -e "Submitting script:"
        echo -e "    \c"
        
        # submit
        qsub $out_full_loc/$prog_name."pbs"
    else
        # user output
        echo -e "Running the script"
        echo -e ""
        
        # run the local script
        . $prog_name"_loc.sh"
    fi
    
    # remove files
    rm $prog_name"_loc.sh"
    
    # finish
    finish
done
