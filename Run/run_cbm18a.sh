# Display usage function
display_usage() { 
    echo -e "\nUsage:\n$0 [OPTS] NR_PROCS \n" 
    echo -e "    OPTS: -d use Valgrind debugging"
    echo -e "          -s trace sources of errors in Valgrind\n"
    } 
#
# Setting some variables
slepc_opt="-st_pc_factor_shift_type NONZERO -st_pc_type lu -st_pc_factor_mat_solver_package mumps"
debug_opt=""
extra_debug_opt=""
n_opt_args=0
#
# Catching options
while getopts "ds" opt; do
    case $opt in
        d)
            debug_opt="valgrind"
            n_opt_args=$((n_opt_args+1))
        ;;
        s)
            extra_debug_opt="--track-origins=yes"
            n_opt_args=$((n_opt_args+1))
        ;;
        \?)
            display_usage
            exit 1
        ;;
    esac
done
#
# Checking for number of input arguments
if [ "$#" != $((n_opt_args+1)) ]; then
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
        echo "Ignoring tracking of sources of errors because no debugging"
        extra_debug_opt=""
    fi
fi
#
# Running command
cp ../PB3D . && mpirun -np $1 $debug_opt $extra_debug_opt ./PB3D input_cbm18a cbm18a $slepc_opt
