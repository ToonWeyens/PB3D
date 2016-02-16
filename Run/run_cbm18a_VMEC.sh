#!/bin/bash
# Display usage function
display_usage() { 
    echo -e "\nUsage:\n$0 [OPTS] NR_PROCS \n" 
    echo -e "    OPTS: -o specify output name"
    echo -e "          -d use Valgrind debugging"
    echo -e "          -s trace sources of errors in Valgrind\n"
    echo -e "          -l check leaks in Valgrind\n"
    } 
#
# Setting some variables
#slepc_opt="-st_pc_factor_shift_type NONZERO -st_pc_type lu -st_pc_factor_mat_solver_package mumps -eps_monitor -eps_view"
slepc_opt="-st_pc_factor_mat_solver_package mumps -eps_monitor -eps_ncv 100 -eps_mpd 100"
debug_opt=""
extra_debug_opt=""
n_opt_args=0
use_out_loc=false
#
# Catching options
while getopts "o:dsl" opt; do
    case $opt in
        o) 
            out_loc=${OPTARG%/}
            use_out_loc=true
            n_opt_args=$((n_opt_args+2))                                        # 2 arguments
        ;;
        d)
            debug_opt="valgrind"
            n_opt_args=$((n_opt_args+1))                                        # 1 argument
        ;;
        s)
            extra_debug_opt=$extra_debug_opt" --track-origins=yes"
            n_opt_args=$((n_opt_args+1))                                        # 1 argument
        ;;
        l)
            extra_debug_opt=$extra_debug_opt" --leak-check=full"
            n_opt_args=$((n_opt_args+1))                                        # 1 argument
        ;;
        \?)
            display_usage
            exit 1
        ;;
    esac
done
#
# Checking for number of input arguments (need at least 1 for # procs)
if [ "$#" -lt $((n_opt_args+1)) ]; then
    display_usage
    exit 1
fi
if [ -n "$debug_opt" ]; then
    echo "Using valgrind for debugging"
    if [ -n "$extra_debug_opt" ]; then
        echo "Also, tracking sources of errors"
    fi
else
    if [ -n "$extra_debug_opt" ]; then
        echo "Ignoring extra debugging options because no debugging"
        extra_debug_opt=""
    fi
fi
#
# Shift arguments to skip options
shift $n_opt_args
#
# Make new folder
if [ "$use_out_loc" = true ]; then
    out=$out_loc
else
    get_date=$(date +"%Y-%m-%d-%H-%M-%S")
    out=$get_date
fi
mkdir -p $out $out/Plots $out/Data $out/Scripts &&
{
# success
echo "Working in directory $out/"
echo ""
# Copy inputs and the program
cp input_cbm18a $out
cp wout_cbm18a.nc $out
cp ../PB3D $out
chmod +x $out/PB3D
cd $out
rm -f .lock_file*
echo "mpirun -np $1 $debug_opt $extra_debug_opt ./PB3D input_cbm18a wout_cbm18a.nc $slepc_opt ${@:2}" > command
mpirun -np $1 $debug_opt $extra_debug_opt ./PB3D input_cbm18a wout_cbm18a.nc $slepc_opt ${@:2}
cd ../
echo ""
echo "Leaving directory $out/"
echo ""
#
# failure
} || {
echo "Error: unable to create directory $out/"
echo "Do you have the correct permissions?"
}
