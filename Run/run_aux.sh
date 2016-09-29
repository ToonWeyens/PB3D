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
    uname -n | grep toon &> /dev/null && machine_ID=1 && on_cluster=false
    uname -n | grep iter &> /dev/null && machine_ID=2 && on_cluster=true
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
                IFS=$'\n' input_mods=($(grep "^[^#;]" "${OPTARG%/}"))
                n_inputs=${#input_mods[@]}
                echo -e "Using $n_inputs input file modification(s)"
                echo -e ""
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

# Setup output name
# input: - default directory of output
# Note: if the -i option is used, subfolder are created for every line of the input modification file.
setup_out() {
    [[ $# -ne 1 ]] && err_arg 1
    if [[ $use_out_loc = true ]]; then
        out=$out_loc
    else
        out=$1
    fi
    [[ $use_input_mods = true ]] && out=$out/$input_i
}

# setup pbs script
init_pbs_script() {
    # print into script
    cat <<- END > $prog_name.pbs
        #!/bin/sh
        #PBS -N $(echo $out | tr '/' '_')
        #PBS -k oe
        #PBS -o $out_full_loc/$prog_name.out
        #PBS -e $out_full_loc/$prog_name.err
        #PBS -l nodes=$n_nodes:ppn=$n_cores
        #PBS -q batch
        #PBS -l vmem=$max_tot_mem$mem_unit
        #PBS -l mem=$max_tot_mem$mem_unit
        #PBS -l walltime=04:00:00
        #PBS -m a
        #PBS -M toon.weyens@gmail.com
        
        # set up general variables
        . /home/ITER/weyenst/load_MPICH3.1.3.sh
        
        # user output
        echo "Job statistics:"
        echo "    ID                \$PBS_JOBID"
        echo "    name              \$PBS_JOBNAME"
        echo "    environment       \$PBS_ENVIRONMENT"
        echo "    nodefile          \$PBS_NODEFILE"
        echo "    array ID          \$PBS_ARRAYID"
        echo "    procs             \$PBS_NP"
        echo "    queue             \$PBS_QUEUE"
        echo "    walltime          \$PBS_WALLTIME"
        echo "    submit directory  \$PBS_O_WORKDIR"
        echo "    host machine      \$PBS_O_HOST"
        echo "    procs per node    \$PBS_NUM_PPN"
        echo "    login name        \$PBS_O_LOGNAME"
        echo "    home              \$PBS_O_HOME"
        echo "" 
END
    # cut 8 leading whitespaces and make executable
    sed -i 's/^.\{8\}//' $prog_name".pbs"
    chmod u+x $prog_name".pbs"
}

# modify input file if requested
# input: - input file
modify_input_file() {
    [[ $# -ne 1 ]] && err_arg 1
     
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
            change_input_var "$1" "$var_name" "$val"
        done
    fi
    echo ""
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
    cd $home
}
