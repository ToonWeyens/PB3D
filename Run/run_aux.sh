#!/bin/bash
# require number of arguments
# run this function as "[[ $# -ne A ]] && err_arg A"
err_arg() { 
    [[ $# -ge 1 ]] && nr_args_req=$1 || nr_args_req=some
    echo "ERROR: need $nr_args_req argument(s)"
    exit 1
}

# set machine ID
set_machine_ID() {
    machine_ID=0
    name=$(uname -n)
    uname -n | grep toon &> /dev/null && machine_ID=1
    uname -n | grep uranus &> /dev/null && machine_ID=2
    if (( machine_ID == 0 )); then
        echo -e "ERROR: invalid machine $(uname -n)"
        exit 1
    else
        echo -e "for user \"$(uname -n)\", machine ID: $machine_ID"
        echo ""
    fi
}

# set DrMemory location
set_drmemory_location() {
    case $machine_ID in
        1)                                                                      # laptop
            drmemory_location="/opt/DrMemory-Linux-1.10.1-3/bin64"
        ;;
        2)                                                                      # uranus
            drmemory_location="/home/tweyens/Programs/DrMemory-Linux-1.10.1-3/bin64"
        ;;
    esac
}

# initialize variables
init_vars() {
    debug_opt=""
    n_opt_args=0
    n_inputs=1
    use_out_loc=false
    use_n_nodes_loc=false
    use_input_mods=false
}

# print OPTS
print_opts() {
    echo -e "    OPTS:      -d use Dr. Memory debugging"
    echo -e "               -h display this help message"
    echo -e "               -i [NAME] specify file name for input modification"
    echo -e "                       * 1 line per element of array"
    echo -e "                       * per line a list of input parameters name and value(s)"
    echo -e "                       * format: name_1 = val_1, name_2 = val_2, ..."
    echo -e "                       * val_i can be an array, delimited by spaces"
    echo -e "               -n [NR] specify nr. of parallel computing nodes"
    echo -e "                       * only valid for clusters"
    echo -e "               -o [NAME] specify output name"
    echo -e "                       * default: date"
}

# Catch options
# input: needs to be run with argument "#@"!
catch_options() {
    local OPTIND opt
    while getopts :o:n:di:h opt; do
        case $opt in
            d)                                                                  # debug
                set_drmemory_location
                debug_opt="$drmemory_location/drmemory --"
                n_opt_args=$((n_opt_args+1))                                    # 1 argument
                echo -e "Using Dr. Memory for debugging"
                echo -e ""
            ;;
            h)                                                                  # help
                echo -e "launcher for program PB3D and POST"
                echo -e "written by Toon Weyens"
                echo -e ""
                display_usage
                exit 1
            ;;
            i)                                                                  # input file modifications
                use_input_mods=true
                n_opt_args=$((n_opt_args+2))                                    # 2 arguments
                IFS=$'\n' read -d '' -r -a input_mods < "${OPTARG%/}"
                n_inputs=${#input_mods[@]}
                echo -e "Using $n_inputs input file modification(s)"
                echo -e ""
            ;;
            n)                                                                  # nr. of nodes on server
                n_nodes_loc=${OPTARG%/}
                use_n_nodes_loc=true
                n_opt_args=$((n_opt_args+2))                                    # 2 arguments
                if (( machine_ID == 2 )); then                                  # uranus
                    echo -e "Using $n_nodes_loc nodes"
                    echo -e ""
                else
                    echo -e "WARNING: ignoring option -n because not in Uranus"
                fi
            ;;
            o)                                                                  # output file
                out_loc=${OPTARG%/}
                use_out_loc=true
                n_opt_args=$((n_opt_args+2))                                    # 2 arguments
            ;;
            \?)                                                                 # other
                echo -e "ERROR: invalid option -$OPTARG"
                display_usage
                exit 1
            ;;
        esac
    done
}

# Setup output directory
# input: - default directory of output
# Note: if the -i option is used, subfolder are created for every line of the input modification file.
setup_output_dir() {
    home_dir=$(pwd)
    
    [[ $# -ne 1 ]] && err_arg 1
    if [[ $use_out_loc = true ]]; then
        out=$out_loc
    else
        out=$1
    fi
    [[ $use_input_mods = true ]] && out=$out/$input_i
    
    out_full=$home_dir/$out
    
    mkdir -p $out_full $out_full/Plots $out_full/Data $out_full/Scripts || {
        # failure
        echo "ERROR: unable to create directory $out_full/ and subdirectories"
        echo "Do you have the correct permissions?"
        exit 1
    }
}

# setup pbs parameters
setup_pbs_params() {
    # set variables
    n_cores_max=16
    memory_factor=3                                                             # safety factor for requesting memory (needs to be integer)
    
    # get memory
    max_tot_mem=$(grep $input_name -e 'max_tot_mem_per_proc' | tr -dc '0-9')
    if [ -z "$max_tot_mem" ]; then
        # default for PB3D and POST
        max_tot_mem=6000
    else
        # multiply by nr. of procs.
        max_tot_mem=$(($nr_procs*$max_tot_mem))
    fi 
    
    # multiply by safety factor
    max_tot_mem=$(($max_tot_mem*$memory_factor))
    max_ind_mem=$((max_tot_mem/$nr_procs))
    mem_unit='mb'
    
    # set nodes and procs
    if [[ $use_n_nodes_loc = true ]]; then
        n_nodes=$n_nodes_loc
    else
        n_nodes=$(($nr_procs/$n_cores_max))
        [[ $nr_procs%$n_cores_max -gt 0 ]] && n_nodes=$(($n_nodes+1))
    fi
    n_cores=$(($nr_procs/$n_nodes))
    [[ $nr_procs%$n_nodes -gt 0 ]] && n_cores=$(($n_cores+1))
}

# setup pbs script
setup_pbs_script() {
    # print into script
    cat <<- END > $prog_name.pbs
        #!/bin/sh
        #PBS -N $out
        #PBS -o $out_full/$prog_name.o
        #PBS -e $out_full/$prog_name.e
        #PBS -l nodes=$n_nodes:ppn=$n_cores
        #PBS -l pvmem=$max_ind_mem$mem_unit
        #PBS -l pmem=$max_ind_mem$mem_unit
        #PBS -l walltime=04:00:00
        #PBS -m abe
        #PBS -M toon.weyens@gmail.com
        cd $out_full
        rm -f $prog_name.o $prog_name.e
        echo "Job statistics:"
        echo ""
        echo "id:   \$PBS_JOBID"
        echo "name: \$PBS_JOBNAME"
        echo "node: \$PBS_QUEUE"
        echo "Walltime: \$PBS_WALLTIME"
        echo ""
        echo "Computing nodes:"
        echo ""
        pbsdsh uname -n
        echo ""
        export PATH="/share/apps/openmpi-1.10.1/bin:/share/apps/gcc-5.2/bin:$PATH"
        export LD_LIBRARY_PATH="$HOME/lib:/share/apps/openmpi-1.10.1/lib:/share/apps/gcc-5.2/lib64:/share/apps/gcc-5.2/lib:$LD_LIBRARY_PATH"
        . /opt/torque/etc/openmpi-setup.sh
        echo "$debug_opt mpirun $mpi_command" > command_$prog_name
        $debug_opt mpirun $mpi_command
        exit
END
    
    # trim leading whitespaces and make executable
    sed -i 's/^ *//' $prog_name.pbs
    chmod u+x $prog_name.pbs
}

# modify input file if requested
modify_input_file() {
    if [[ $use_input_mods = true ]]; then
        echo -e "modifications for input file $input_i/$n_inputs:"
        input_mod_i=${input_mods[$input_i-1]}
        
        # split input line in ',' elements (from http://stackoverflow.com/a/918931)
        IFS=',' read -ra input_elem <<< "$input_mod_i"
        
        # set all the input modifications
        for j in "${input_elem[@]}"; do
            # first element is name
            var_name=$(echo $j | cut -d = -f 1 )
            
            # rest is value (no error checking!)
            val=$(echo $j | cut --complement -d = -f 1 )
            
            # change input file variable
            change_input_var "$input_name" "$var_name" "$val"
        done
    fi
    echo ""
}

# run mpi command, depending on machine ID
run_mpi_command() {
    case $machine_ID in
        1)                                                                          # laptop
            echo "$debug_opt mpirun -np $nr_procs $mpi_command" > command_$prog_name
            $debug_opt mpirun -np $nr_procs $mpi_command
        ;;
        2)                                                                          # uranus
            # setup pbs parameters
            setup_pbs_params
            
            # user output
            echo "creating pbs script with"
            echo "    $n_nodes node(s) and $n_cores core(s) per node"
            echo "    $max_ind_mem$mem_unit memory per process"
            echo ""
            echo "job:"
            echo -e "    \c"
            
            # setup pbs script
            setup_pbs_script
            
            # submit
            qsub $prog_name.pbs
        ;;
    esac
}

# change an input variable of a PB3D or POST input file.
# input: - input file
#        - variable name
#        - replacement value
change_input_var() {
    [[ $# -ne 3 ]] && err_arg 3
    # find the line of var
    var_line=$(grep -nr -m 1 $2 $1)
    # find the line number of var
    var_line_nr=$(echo $var_line | cut -d : -f 1)
    # replace whole line by new line
    sed -i "${var_line_nr}s/.*/    $2 = $3/" "$1"
    echo -e "    replace $2 by $3"

    # Note: in an old version, only integer variables were treated and they were replaced inline, using
    # # find the old value of var
    # var_old=$(echo $var_line | cut -d : -f 2 | tr -dc '0-9')
    # echo var_old=$var_old
    # # replace old value of var by new value
    # sed -i "${var_line_nr}s/${var_old}/${3}/" "$1"
}

# finish
finish() {
    # return
    echo ""
    echo "Leaving directory $out_full/"
    echo ""
    cd $home_dir
}
