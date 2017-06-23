#!/bin/bash
# Structure from http://stackoverflow.com/a/13588485/3229162
# Makes use of auxiliary scripts
#   - gen_node_list.sh: to generate the correct nodes on the ITER cluster
#   - get_disk_space.sh: to create a list of broken nodes on the ITER cluster (has to be called manually)

# Main program
main() {
    # set bash source
    DIR="${BASH_SOURCE%/*}"
    if [[ ! -d "$DIR" ]]; then DIR="$PWD"; fi

    # set program ID
    set_prog_ID $1
    
    # shift first argument
    shift 1

    # set machine ID
    set_machine_ID

    # initialize variables
    init_vars

    # catch options
    catch_options $@

    # shift option arguments
    shift $n_opt_args

    # check number of input arguments
    # (need at least 1 for # procs and 1 for case)
    if [[ $# -lt 2 ]]; then
        echo "ERROR: Insufficient number of arguments"
        echo ""
        display_usage
        exit 1
    fi

    # Set parameters: input_name, eq_name
    set_input $1
    
    # set base
    base=$(pwd)
    
    # set program directory
    if [[ $on_cluster = true ]]; then
        prog_dir=$base/..
    else
        prog_dir=/opt/PB3D
    fi

    # loop over all inputs
    # (from http://www.cyberciti.biz/faq/unix-linux-iterate-over-a-variable-range-of-numbers-in-bash/)
    for (( input_i=1; input_i<=$n_inputs; input_i++ )); do
        # setup output name
        case $prog_ID in
            1)  # PB3D
                set_output $(date +"%Y-%m-%d-%H-%M-%S")
            ;;
            2)  # POST
                set_output ${1%/}
            ;;
        esac
        
        # set variables
        # Note  about difference  between out_full  and out_full_loc:  The input
        # file  is potentially  modified before  use. This  happens in  the base
        # folder. Afterwards, the run directory  is created on the local system,
        # and this  modified input file is  copied there. On a  normal computer,
        # this is also  where the simulations are done. On  server, however, the
        # actual  computation happens  in  another folder  on the  computational
        # node, and afterwards everything is copied back to the local node
        case $machine_ID in
            1)  # XPS-L501X
                nr_procs=$2                                                     # take nr_procs from input
                out_full=$base/$out                                             # use subdirectory
            ;;
            2)  # ITER
                ##########
                # select your queue here: compute, testqueue, ib or ib_gen8
                # + additional node name: batch,   testqueue, ib or ib_gen8
                ##########
                queue="compute"
                queue_name="batch"
                n_nodes=1                                                       # use one node
                node_list=$(./gen_node_list.sh -m $n_nodes $queue $(cat broken_nodes.txt 2> /dev/null))
                n_cores=${node_list#*=}
                nr_procs=$(( n_cores < $2 ? n_cores : $2 ))                     # take nr_procs from input, limited by number of cores
                memory_factor=2                                                 # safety factor for requesting memory (needs to be integer)
                case $queue in
                    compute)
                        max_tot_mem=30000
                    ;;
                    testqueue)
                        max_tot_mem=20000
                    ;;
                    ib)
                        max_tot_mem=20000
                    ;;
                    ib_gen8)
                        max_tot_mem=60000
                    ;;
                esac
                mem_unit='mb'                                                   # MB
                out_full=/tmp/$out                                              # use temporary directory in local node
            ;;
        esac
        
        out_full_loc=$base/$out                                                 # afterwards copy to this directory
        
        # save other options
        other_opts=${@:3}
        
        # user output
        echo -e "work directory:    $out_full/"
        echo -e ""
        
        # create local directory
        mkdir -p $out_full_loc || {
            # failure
            echo "ERROR: unable to create directory $out_full_loc/"
            echo "Do you have the correct permissions?"
            echo ""
            exit 1
        }
        
        # possibly modify copy of input file
        cp $input_name $out_full_loc/$input_name
        modify_input_file $out_full_loc/$input_name
        
        # set up local simulation script
        create_loc_script
        
        if [[ $on_cluster = true ]]; then
            # user output
            echo -e "The output will be copied to '$out_full_loc' afterwards"
            echo -e ""
            
            # user output
            echo -e "Creating pbs script in '$prog_name.pbs' with"
            echo -e "    $n_nodes node(s) and $n_cores core(s) per node"
            echo -e "    $max_tot_mem$mem_unit memory in total"
            echo -e ""
            
            # initizalize pbs script
            echo -e "    setting up part 1"
            setup_pbs_script_1
            
            # add local execution script to pbs script
            echo -e "    including '${prog_name}_loc.sh'"
            cat $prog_name"_loc.sh" >> $prog_name".pbs"
            
            # set up last part of pbs script
            echo -e "    setting up part 2"
            setup_pbs_script_2
            
            # mv pbs script to output
            mv $prog_name".pbs" $out_full_loc/
            
            # user output
            echo -e ""
            echo -e "Submitting script:"
            echo -e ""
            
            # submit
            qsub $out_full_loc/$prog_name."pbs"
            
            # write detailed job information
            qstat -f $(qselect -u weyenst | tail -n 1) > $out_full_loc/qstat
        else
            # user output
            echo -e "Running the script"
            echo -e ""
            
            # run the local script
            . $prog_name"_loc.sh"
        fi
        
        # remove files
        rm $base/$prog_name"_loc.sh"
    done
    
    # possibly copy input modification array to base
    [[ $use_input_mods = true ]] && cp $base/$input_mod_file $base/$out_base/array_input
    
    # finish
    echo -e ""
    echo -e "work directory:    $out_full/"
    echo -e "Finished"
    exit 0
}

# require number of arguments
# run this function as "[[ $# -ne A ]] && err_arg A"
err_arg() { 
    [[ $# -ge 1 ]] && nr_args_req=$1 || nr_args_req=some
    echo "ERROR: need $nr_args_req argument(s)"
    echo ""
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
        echo ""
        exit 1
    else
        echo -e "for user \"$(uname -n)\", machine ID: $machine_ID"
        echo ""
    fi
}

# set program ID and name
# input: - program name: 1. PB3D, 2. POST
# sets: prog_ID
#       prog_name
set_prog_ID() {
    if [[ $# -ne 1 ]]; then
        display_usage
    fi
    prog_name=$1
    prog_ID=0
    [[ "$prog_name" == "PB3D" ]] && prog_ID=1
    [[ "$prog_name" == "POST" ]] && prog_ID=2
    if (( prog_ID == 0 )); then
        echo "ERROR: invalid program name $prog_name"
        echo ""
        display_usage
    fi
}

# print OPTS
display_usage() {
    case $prog_ID in
        1)  # PB3D
            echo -e "Usage:\n$0 $prog_name [OPTS] CASE NR_PROCS \n" 
            echo -e "               -i [NAME] specify file name for input modification"
            echo -e "                       * 1 line per element of array"
            echo -e "                       * per line a list of input parameters name and value(s)"
            echo -e "                       * format: name_1 = val_1, name_2 = val_2, ..."
            echo -e "                       * val_i can be an array, delimited by spaces"
        ;;
        2)  # POST
            echo -e "Usage:\n$0 $prog_name [OPTS] PB3D_DIR NR_PROCS \n" 
        ;;
        *)  # unknown
            echo -e "Usage:\n$0 [PROG NAME] [OPTIONS] \n"
            echo -e "    PROG_NAME: PB3D"
            echo -e "               POST"
            echo -e ""
            echo -e "    OPTIONS:   depend on PROG_NAME"
            echo -e ""
            exit 1
        ;;
    esac
    echo -e "    OPTS:      -d use Dr. Memory debugging"
    echo -e "               -h display this help message"
    echo -e "               -o [NAME] specify output name"
    echo -e "                       * default: date"
    echo -e ""
    case $prog_ID in
        1)  # PB3D
            echo -e "    CASE:       1 cbm18a"
            echo -e "                2 cbm18a_HELENA"
            echo -e "                3 cbm18a_small"
            echo -e "                4 cbm18a_small_HELENA"
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
            echo -e "               41 D_ped4.0"
            echo -e "               42 D_ped4.0_HELENA"
            echo -e "               43 D_ped3.0"
            echo -e "               44 D_ped3.0_HELENA"
            echo -e "               45 D_ped2.0"
            echo -e "               46 D_ped2.0_HELENA"
            echo -e "               47 D_ped1.0"
            echo -e "               48 D_ped1.0_HELENA"
            echo -e "               49 D_ped0.9"
            echo -e "               50 D_ped0.9_HELENA"
            echo -e "               51 D_ped0.8"
            echo -e "               52 D_ped0.8_HELENA"
            echo -e "               53 D_ped0.7"
            echo -e "               54 D_ped0.7_HELENA"
            echo -e "               55 D_ped0.6"
            echo -e "               56 D_ped0.6_HELENA"
            echo -e "               57 D_ped0.5"
            echo -e "               58 D_ped0.5_HELENA"
            echo -e "               59 D_ped0.4"
            echo -e "               60 D_ped0.4_HELENA"
            echo -e "               61 D_ped0.3"
            echo -e "               62 D_ped0.3_HELENA"
            echo -e ""
            echo -e "               71 D_BS_ped4.0"
            echo -e "               72 D_BS_ped4.0_HELENA"
            echo -e "               73 D_BS_ped3.0"
            echo -e "               74 D_BS_ped3.0_HELENA"
            echo -e "               75 D_BS_ped2.0"
            echo -e "               76 D_BS_ped2.0_HELENA"
            echo -e "               77 D_BS_ped1.0"
            echo -e "               78 D_BS_ped1.0_HELENA"
            echo -e "               79 D_BS_ped0.9"
            echo -e "               80 D_BS_ped0.9_HELENA"
            echo -e "               81 D_BS_ped0.8"
            echo -e "               82 D_BS_ped0.8_HELENA"
            echo -e "               83 D_BS_ped0.7"
            echo -e "               84 D_BS_ped0.7_HELENA"
            echo -e "               85 D_BS_ped0.6"
            echo -e "               86 D_BS_ped0.6_HELENA"
            echo -e "               87 D_BS_ped0.5"
            echo -e "               88 D_BS_ped0.5_HELENA"
            echo -e "               89 D_BS_ped0.4"
            echo -e "               90 D_BS_ped0.4_HELENA"
            echo -e "               91 D_BS_ped0.3"
            echo -e "               92 D_BS_ped0.3_HELENA"
            echo -e ""
            echo -e "              141 D_ripple1.0_ped4.0"
            echo -e "              142 D_ripple1.0_ped4.0_HELENA"
            echo -e "              143 D_ripple1.0_ped3.0"
            echo -e "              144 D_ripple1.0_ped3.0_HELENA"
            echo -e "              145 D_ripple1.0_ped2.0"
            echo -e "              146 D_ripple1.0_ped2.0_HELENA"
            echo -e "              147 D_ripple1.0_ped1.0"
            echo -e "              148 D_ripple1.0_ped1.0_HELENA"
            echo -e "              149 D_ripple1.0_ped0.9"
            echo -e "              150 D_ripple1.0_ped0.9_HELENA"
            echo -e "              151 D_ripple1.0_ped0.8"
            echo -e "              152 D_ripple1.0_ped0.8_HELENA"
            echo -e "              153 D_ripple1.0_ped0.7"
            echo -e "              154 D_ripple1.0_ped0.7_HELENA"
            echo -e "              155 D_ripple1.0_ped0.6"
            echo -e "              156 D_ripple1.0_ped0.6_HELENA"
            echo -e "              157 D_ripple1.0_ped0.5"
            echo -e "              158 D_ripple1.0_ped0.5_HELENA"
            echo -e "              159 D_ripple1.0_ped0.4"
            echo -e "              160 D_ripple1.0_ped0.4_HELENA"
            echo -e "              161 D_ripple1.0_ped0.3"
            echo -e "              162 D_ripple1.0_ped0.3_HELENA"
            echo -e ""
            echo -e "              171 D_BS_ripple1.0_ped4.0"
            echo -e "              172 D_BS_ripple1.0_ped4.0_HELENA"
            echo -e "              173 D_BS_ripple1.0_ped3.0"
            echo -e "              174 D_BS_ripple1.0_ped3.0_HELENA"
            echo -e "              175 D_BS_ripple1.0_ped2.0"
            echo -e "              176 D_BS_ripple1.0_ped2.0_HELENA"
            echo -e "              177 D_BS_ripple1.0_ped1.0"
            echo -e "              178 D_BS_ripple1.0_ped1.0_HELENA"
            echo -e "              179 D_BS_ripple1.0_ped0.9"
            echo -e "              180 D_BS_ripple1.0_ped0.9_HELENA"
            echo -e "              181 D_BS_ripple1.0_ped0.8"
            echo -e "              182 D_BS_ripple1.0_ped0.8_HELENA"
            echo -e "              183 D_BS_ripple1.0_ped0.7"
            echo -e "              184 D_BS_ripple1.0_ped0.7_HELENA"
            echo -e "              185 D_BS_ripple1.0_ped0.6"
            echo -e "              186 D_BS_ripple1.0_ped0.6_HELENA"
            echo -e "              187 D_BS_ripple1.0_ped0.5"
            echo -e "              188 D_BS_ripple1.0_ped0.5_HELENA"
            echo -e "              189 D_BS_ripple1.0_ped0.4"
            echo -e "              190 D_BS_ripple1.0_ped0.4_HELENA"
            echo -e "              191 D_BS_ripple1.0_ped0.3"
            echo -e "              192 D_BS_ripple1.0_ped0.3_HELENA"
            echo -e ""
            echo -e "              201 Hmode_PHD_ped0.3_HELENA"
            echo -e "              202 Hmode_PHD_ped0.4_HELENA"
            echo -e "              203 Hmode_PHD_ped0.5_HELENA"
            echo -e "              204 Hmode_PHD_ped0.6_HELENA"
            echo -e "              205 Hmode_PHD_ped0.7_HELENA"
            echo -e "              206 Hmode_PHD_ped0.8_HELENA"
            echo -e "              207 Hmode_PHD_ped0.9_HELENA"
            echo -e "              208 Hmode_PHD_ped1.0_HELENA"
            echo -e "              209 Hmode_PHD_ped1.5_HELENA"
            echo -e "              210 Hmode_PHD_ped2.0_HELENA"
            echo -e "              211 Hmode_PHD_ped3.0_HELENA"
            echo -e "              212 Hmode_PHD_ped4.0_HELENA"
        ;;
        2)  # POST
            echo -e "    PB3D_DIR:  PB3D directory"
        ;;
    esac
    echo -e ""
    echo -e "    NR_PROCS:  nr. of MPI processes"
    echo -e ""
}

# initialize variables
# sets: debug_opt
#       n_opt_args
#       n_inputs
#       use_out_loc
#       use_input_mods
#       PB3D:   slepc_opt
init_vars() {
    debug_opt=""
    n_opt_args=0
    n_inputs=1
    use_out_loc=false
    use_input_mods=false
    case $prog_ID in
        1)  # PB3D
            #slepc_opt="-st_ksp_type preonly -st_pc_type lu -st_pc_factor_mat_solver_package mumps -eps_monitor -eps_ncv 32 -eps_view -log_view"
            slepc_opt="-st_ksp_type preonly -st_pc_type lu -st_pc_factor_mat_solver_package mumps -eps_monitor"
        ;;
        2)  # POST
            # nothing else
        ;;
    esac
}

# Catch options
# input: all arguments "#@"
# sets: d:  valgrind_location
#           debug_opt
#       o:  out_loc
#           use_out_loc
#       i:  use_input_mods  (only PB3D)
#           n_inputs        (only PB3D)
catch_options() {
    local OPTIND opt
    while getopts :o:n:di:h opt; do
        case $opt in
            h)  # help
                echo -e "launcher for program PB3D and POST"
                echo -e "written by Toon Weyens"
                echo -e ""
                display_usage
                exit 0
            ;;
            d)  # debug
                case $machine_ID in
                    1)  # XPS-L501X
                        valgrind_location="/usr/bin/"
                    ;;
                    2)  # ITER
                        valgrind_location="/usr/bin/"
                    ;;
                esac
                debug_opt="$valgrind_location/valgrind --"
                n_opt_args=$((n_opt_args+1))                                    # 1 argument
                echo -e "Using Dr. Memory for debugging"
                echo -e ""
            ;;
            i)  # input file modifications
                case $prog_ID in
                    1)  # PB3D
                        use_input_mods=true
                        n_opt_args=$((n_opt_args+2))                            # 2 arguments
                        input_mod_file=${OPTARG%/}
                        IFS=$'\n' input_mods=($(grep "^[^#;]" "$input_mod_file"))
                        n_inputs=${#input_mods[@]}
                        echo -e "Using $n_inputs input file modification(s)"
                        echo -e ""
                    ;;
                    2)  # POST
                        echo -e "ERROR: Input file modification(s) not possible for POST"
                        echo -e ""
                        display_usage
                        exit 1
                    ;;
                esac
            ;;
            o)  # output file
                out_loc=${OPTARG%/}
                use_out_loc=true
                n_opt_args=$((n_opt_args+2))                                    # 2 arguments
            ;;
            *)  # other
                echo -e "ERROR: invalid option -$OPTARG"
                echo -e ""
                display_usage
                exit 1
            ;;
        esac
    done
}

# Set input parameters
# input: $1
# sets: PB3D:   input_name
#               eq_name
#       POST:   input_name
#               PB3D_out_name
#               PB3D_out_full
set_input() {
    [[ $# -ne 1 ]] && err_arg 1
    case $prog_ID in
        1)  # PB3D
            case $1 in
                [1-4])
                    input_name=cbm18a
                ;;
                1[1-8])
                    input_name=cbm18a_ripple
                ;;
                21)
                    input_name=cdxu
                ;;
                4[1-9]|5[0-9]|6[0-2]|7[1-9]|8[0-9]|9[0-2]|301|501)
                    input_name=Hmode
                ;;
                13[1-9]|14[0-9]|15[0-9]|16[0-9]|17[0-9]|180|20[1-5]|30[2-7]|401|502|60[1-5]|70[1-3])
                    input_name=Hmode_ripple
                ;;
                *)
                    echo -e "ERROR: Case $1 not found"
                    echo -e ""
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
                3)
                    eq_name=wout_cbm18a_small.nc
                ;;
                4)
                    eq_name=cbm18a_small
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
                41)
                    eq_name=wout_D_ped4.0.nc
                ;;
                42)
                    eq_name=D_ped4.0
                ;;
                43)
                    eq_name=wout_D_ped3.0.nc
                ;;
                44)
                    eq_name=D_ped3.0
                ;;
                45)
                    eq_name=wout_D_ped2.0.nc
                ;;
                46)
                    eq_name=D_ped2.0
                ;;
                47)
                    eq_name=wout_D_ped1.0.nc
                ;;
                48)
                    eq_name=D_ped1.0
                ;;
                49)
                    eq_name=wout_D_ped0.9.nc
                ;;
                50)
                    eq_name=D_ped0.9
                ;;
                51)
                    eq_name=wout_D_ped0.8.nc
                ;;
                52)
                    eq_name=D_ped0.8
                ;;
                53)
                    eq_name=wout_D_ped0.7.nc
                ;;
                54)
                    eq_name=D_ped0.7
                ;;
                55)
                    eq_name=wout_D_ped0.6.nc
                ;;
                56)
                    eq_name=D_ped0.6
                ;;
                57)
                    eq_name=wout_D_ped0.5.nc
                ;;
                58)
                    eq_name=D_ped0.5
                ;;
                59)
                    eq_name=wout_D_ped0.4.nc
                ;;
                60)
                    eq_name=D_ped0.4
                ;;
                61)
                    eq_name=wout_D_ped0.3.nc
                ;;
                62)
                    eq_name=D_ped0.3
                ;;
                71)
                    eq_name=wout_D_BS_ped4.0.nc
                ;;
                72)
                    eq_name=D_BS_ped4.0
                ;;
                73)
                    eq_name=wout_D_BS_ped3.0.nc
                ;;
                74)
                    eq_name=D_BS_ped3.0
                ;;
                75)
                    eq_name=wout_D_BS_ped2.0.nc
                ;;
                76)
                    eq_name=D_BS_ped2.0
                ;;
                77)
                    eq_name=wout_D_BS_ped1.0.nc
                ;;
                78)
                    eq_name=D_BS_ped1.0
                ;;
                79)
                    eq_name=wout_D_BS_ped0.9.nc
                ;;
                80)
                    eq_name=D_BS_ped0.9
                ;;
                81)
                    eq_name=wout_D_BS_ped0.8.nc
                ;;
                82)
                    eq_name=D_BS_ped0.8
                ;;
                83)
                    eq_name=wout_D_BS_ped0.7.nc
                ;;
                84)
                    eq_name=D_BS_ped0.7
                ;;
                85)
                    eq_name=wout_D_BS_ped0.6.nc
                ;;
                86)
                    eq_name=D_BS_ped0.6
                ;;
                87)
                    eq_name=wout_D_BS_ped0.5.nc
                ;;
                88)
                    eq_name=D_BS_ped0.5
                ;;
                89)
                    eq_name=wout_D_BS_ped0.4.nc
                ;;
                90)
                    eq_name=D_BS_ped0.4
                ;;
                91)
                    eq_name=wout_D_BS_ped0.3.nc
                ;;
                92)
                    eq_name=D_BS_ped0.3
                ;;
                141)
                    eq_name=wout_D_ripple1.0_ped4.0.nc
                ;;
                142)
                    eq_name=D_ripple1.0_ped4.0
                ;;
                143)
                    eq_name=wout_D_ripple1.0_ped3.0.nc
                ;;
                144)
                    eq_name=D_ripple1.0_ped3.0
                ;;
                145)
                    eq_name=wout_D_ripple1.0_ped2.0.nc
                ;;
                146)
                    eq_name=D_ripple1.0_ped2.0
                ;;
                147)
                    eq_name=wout_D_ripple1.0_ped1.0.nc
                ;;
                148)
                    eq_name=D_ripple1.0_ped1.0
                ;;
                149)
                    eq_name=wout_D_ripple1.0_ped0.9.nc
                ;;
                150)
                    eq_name=D_ripple1.0_ped0.9
                ;;
                151)
                    eq_name=wout_D_ripple1.0_ped0.8.nc
                ;;
                152)
                    eq_name=D_ripple1.0_ped0.8
                ;;
                153)
                    eq_name=wout_D_ripple1.0_ped0.7.nc
                ;;
                154)
                    eq_name=D_ripple1.0_ped0.7
                ;;
                155)
                    eq_name=wout_D_ripple1.0_ped0.6.nc
                ;;
                156)
                    eq_name=D_ripple1.0_ped0.6
                ;;
                157)
                    eq_name=wout_D_ripple1.0_ped0.5.nc
                ;;
                158)
                    eq_name=D_ripple1.0_ped0.5
                ;;
                159)
                    eq_name=wout_D_ripple1.0_ped0.4.nc
                ;;
                160)
                    eq_name=D_ripple1.0_ped0.4
                ;;
                161)
                    eq_name=wout_D_ripple1.0_ped0.3.nc
                ;;
                162)
                    eq_name=D_ripple1.0_ped0.3
                ;;
                171)
                    eq_name=wout_D_BS_ripple1.0_ped4.0.nc
                ;;
                172)
                    eq_name=D_BS_ripple1.0_ped4.0
                ;;
                173)
                    eq_name=wout_D_BS_ripple1.0_ped3.0.nc
                ;;
                174)
                    eq_name=D_BS_ripple1.0_ped3.0
                ;;
                175)
                    eq_name=wout_D_BS_ripple1.0_ped2.0.nc
                ;;
                176)
                    eq_name=D_BS_ripple1.0_ped2.0
                ;;
                177)
                    eq_name=wout_D_BS_ripple1.0_ped1.0.nc
                ;;
                178)
                    eq_name=D_BS_ripple1.0_ped1.0
                ;;
                179)
                    eq_name=wout_D_BS_ripple1.0_ped0.9.nc
                ;;
                180)
                    eq_name=D_BS_ripple1.0_ped0.9
                ;;
                181)
                    eq_name=wout_D_BS_ripple1.0_ped0.8.nc
                ;;
                182)
                    eq_name=D_BS_ripple1.0_ped0.8
                ;;
                183)
                    eq_name=wout_D_BS_ripple1.0_ped0.7.nc
                ;;
                184)
                    eq_name=D_BS_ripple1.0_ped0.7
                ;;
                185)
                    eq_name=wout_D_BS_ripple1.0_ped0.6.nc
                ;;
                186)
                    eq_name=D_BS_ripple1.0_ped0.6
                ;;
                187)
                    eq_name=wout_D_BS_ripple1.0_ped0.5.nc
                ;;
                188)
                    eq_name=D_BS_ripple1.0_ped0.5
                ;;
                189)
                    eq_name=wout_D_BS_ripple1.0_ped0.4.nc
                ;;
                190)
                    eq_name=D_BS_ripple1.0_ped0.4
                ;;
                191)
                    eq_name=wout_D_BS_ripple1.0_ped0.3.nc
                ;;
                192)
                    eq_name=D_BS_ripple1.0_ped0.3
                ;;
                201)
                    eq_name=Hmode_PHD_ped0.3
                ;;
                202)
                    eq_name=Hmode_PHD_ped0.4
                ;;
                203)
                    eq_name=Hmode_PHD_ped0.5
                ;;
                204)
                    eq_name=Hmode_PHD_ped0.6
                ;;
                205)
                    eq_name=Hmode_PHD_ped0.7
                ;;
                206)
                    eq_name=Hmode_PHD_ped0.8
                ;;
                207)
                    eq_name=Hmode_PHD_ped0.9
                ;;
                208)
                    eq_name=Hmode_PHD_ped1.0
                ;;
                209)
                    eq_name=Hmode_PHD_ped1.5
                ;;
                210)
                    eq_name=Hmode_PHD_ped2.0
                ;;
                211)
                    eq_name=Hmode_PHD_ped3.0
                ;;
                212)
                    eq_name=Hmode_PHD_ped4.0
                ;;
                *)
                    echo -e "ERROR: Case $1 not found"
                    echo -e ""
                    display_usage
                    exit 1
                ;;
            esac
            
            echo -e "input file:        $input_name"
            echo -e "equilibrium file:  $eq_name"
        ;;
        2)  # POST
            PB3D_out_name="PB3D_out.h5"
            PB3D_out_full="${1%/}/$PB3D_out_name"
            input_name="input_POST"
            if [[ -d "$1" ]] ; then
                if [[ ! -f "$PB3D_out_full" ]] ; then
                    echo -e "ERROR: no PB3D output file found in '${1%/}/'"
                    echo -e ""
                    display_usage
                    exit 1
                fi
            else
                echo -e "ERROR: PB3D Directory '${1%/}/' does not exist"
                echo -e ""
                display_usage
                exit 1
            fi
             
            echo -e "input file:        '$input_name'"
            echo -e "PB3D output file:  '$PB3D_out_full'"
        ;;
    esac
}

# Setup output name
# input: default directory of output
# sets: out_base
#       out
# Note: if the -i option is used, subfolder are created for every line of the input modification file.
set_output() {
    [[ $# -ne 1 ]] && err_arg 1
    if [[ $use_out_loc = true ]]; then
        out_base=$out_loc
    else
        out_base=$1
    fi
    if [[ $use_input_mods = true ]]; then
        out=$out_base/$input_i
    else
        out=$out_base
    fi
}

# modify input file if requested
# input: - input file to potentially be modified
modify_input_file() {
    [[ $# -ne 1 ]] && err_arg 1
    if [[ $use_input_mods = true ]]; then
        echo -e "modifications for input file $input_i/$n_inputs:"
        local input_mod_i=${input_mods[$input_i-1]}
        
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
        echo ""
    fi
}

# change an input variable of a PB3D or POST input file.
# input: - input file
#        - variable name
#        - replacement value
change_input_var() {
    [[ $# -ne 3 ]] && err_arg 3
    # find the line of var
    var_line=$(grep -E -nr -m 1 "(^| )$2" $1)
    # find the line number of var
    var_line_nr=$(echo $var_line | cut -d : -f 1)
    # replace whole line by new line
    sed -i "${var_line_nr}s/.*/    $2 = $3/" "$1"
    echo -e "    replace $2 by $3"
}

# Create local run script
# sets: PB3D:   PB3D_loc.sh
#       POST:   POST_loc.sh
create_loc_script() {
    echo -e "Creating local simulation script in '"$prog_name"_loc.sh'"
    echo -e ""
    cat <<- END > $prog_name"_loc.sh"
        
        # create directory structure
        out_full=$out_full
        base=$base
        mkdir -p \$out_full \$out_full/Plots \$out_full/Data \$out_full/Scripts || {
            # failure
            echo "ERROR: unable to create directory $out_full/ and subdirectories"
            echo "Do you have the correct permissions?"
            echo ""
            exit 1
        }
        
        # copy inputs and the program
        rsync --progress -zvhL $out_full_loc/$input_name \$out_full/ 2> /dev/null
        rsync --progress -zvhL $out_full_loc/${prog_name}_out.txt \$out_full/ 2> /dev/null
        $(aux_copy_inputs)
        rsync --progress -zvhL $prog_dir/$prog_name \$out_full/
        chmod +x \$out_full/$prog_name
        
        # go to run directory
        cd \$out_full
        
        # remove possible output and error file
        rm -f $prog_name.out $prog_name.err
        
        # set mpi command
        MPI_command="$(init_MPI_command)"
        
        # run mpi command
        my_mpirun=$(which mpirun)
        echo "$(echo $debug_opt \$my_mpirun | xargs) -np $nr_procs \$MPI_command" > command_$prog_name
        . command_$prog_name
        
        # go back to base directory
        cd \$base
END
    # cut leading whitespaces and make executable
    sed -i 's/^[ \t]*//' $prog_name"_loc.sh"
    chmod u+x $prog_name"_loc.sh"
}

# Auxiliary procedure to copy inputs
aux_copy_inputs() {
    case $prog_ID in
        1)  # PB3D
            echo "rsync --progress -vhL --compress-level=9 $out_full_loc/PB3D_out.h5 \$out_full/ 2> /dev/null"
            echo "rsync --progress -zvhL \$base/$eq_name \$out_full/"
        ;;
        2)  # POST
            echo "ln -sf $base/$PB3D_out_full $out_full/ 2> /dev/null"
        ;;
    esac
}

# Initialize MPI command
### Note: --minim_ouptut WAS hard-coded, but it can be removed
# sets: MPI_command
init_MPI_command() {
    case $prog_ID in
        1)  # PB3D
            echo "./$prog_name $input_name $eq_name $slepc_opt $other_opts"
        ;;
        2)  # POST
            echo "./$prog_name $input_name $PB3D_out_name $other_opts"
        ;;
    esac
}

# setup first part of pbs script
# sets: PB3D:   PB3D.pbs
#       POST:   POST.pbs
### previously: PBS -l nodes=$n_nodes:ppn=$n_cores
setup_pbs_script_1() {
    # print into script
    cat <<- END > $prog_name.pbs
        #!/bin/sh
        #PBS -N $(echo $out | tr '/' '_')
        #PBS -k oe
        #PBS -o $out_full_loc/$prog_name.out
        #PBS -e $out_full_loc/$prog_name.err
        #PBS -l nodes=$node_list
        #PBS -q $queue_name
        #PBS -l vmem=$max_tot_mem$mem_unit
        #PBS -l mem=$max_tot_mem$mem_unit
        #PBS -l walltime=12:00:00
        #PBS -m a
        #PBS -M toon.weyens@gmail.com
        
        # set up general variables
        . /home/ITER/weyenst/load_MPICH3.1.3_nomods.sh
        
        # user output
        echo "Job statistics:"
        echo "    ID                \$PBS_JOBID"
        echo "    name              \$PBS_JOBNAME"
        echo "    environment       \$PBS_ENVIRONMENT"
        echo "    nodefile          \$PBS_NODEFILE"
        echo "    nodes             \$(cat \$PBS_NODEFILE | paste -sd ',' -)"
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
        
        # set memory limit
        ulimit -l unlimited
        
        # set local output and error and create symbolic links
        loc_out=\$PBS_O_HOME/$(echo $out | tr '/' '_').o\$(echo \$PBS_JOBID | cut -d'.' -f 1)
        loc_err=\$PBS_O_HOME/$(echo $out | tr '/' '_').e\$(echo \$PBS_JOBID | cut -d'.' -f 1)
        ln -sf \$loc_out $out_full_loc/$prog_name.out
        ln -sf \$loc_err $out_full_loc/$prog_name.err
END
    # cut leading whitespaces and make executable
    sed -i 's/^[ \t]*//' $prog_name".pbs"
    chmod u+x $prog_name".pbs"
}

# setup second part of  pbs script
# sets: PB3D:   PB3D.pbs
#       POST:   POST.pbs
setup_pbs_script_2() {
    # add commands to move output after finishing, but only if succesfully transfered
    echo "" >> $prog_name".pbs"
    echo "# Done, copy files back with rsync" >> $prog_name".pbs"
    echo "echo '# Copy results with rsync'" >> $prog_name".pbs"
    echo "find $out_full/ -name '*temp.h5' | xargs rm -f" >> $prog_name".pbs"
    echo "cd $out_full_loc" >> $prog_name".pbs"
    echo "rsync --remove-source-files --progress --exclude='PB3D_out.h5' -zvhr $out_full/* ." >> $prog_name".pbs"
    if (( prog_ID == 1 )); then
        echo "rsync --remove-source-files --progress -zvhr $out_full/* ." >> $prog_name".pbs"
    fi
    echo "[[ \$? -eq 0 ]] && rm -r $out_full" >> $prog_name".pbs"
    echo "mv \$loc_out $prog_name.out" >> $prog_name".pbs"
    echo "mv \$loc_err $prog_name.err" >> $prog_name".pbs"
    echo "[[ -s $prog_name.err ]] || rm $prog_name.err" >> $prog_name".pbs"
    
    # terminate pbs script
    echo "" >> $prog_name".pbs"
    echo "exit 0" >> $prog_name".pbs"
}

# Run the main program
main "$@"
