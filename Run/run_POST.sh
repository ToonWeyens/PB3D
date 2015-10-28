#!/bin/bash
# Display usage function
display_usage() { 
    echo -e "\nUsage:\n$0 [OPTS] PB3D_DIR NR_PROCS \n" 
    echo -e "    OPTS: -o specify output name"
    echo -e "          -d use Valgrind debugging"
    echo -e "          -s trace sources of errors in Valgrind\n"
    echo -e "          -l check leaks in Valgrind\n"
    } 
#
# Setting some variables
debug_opt=""
extra_debug_opt=""
n_opt_args=0
use_out_loc=false
PB3D_out_name="PB3D_out.h5"
POST_in_name="input_POST"
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
if [ "$#" -lt $((n_opt_args+2)) ]; then
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
# Shift arguments to skip options
shift $n_opt_args
#
# Check existence of PB3D folder
#
if [[ -d "$1" ]] ; then
    if [[ ! -f "${1%/}/$PB3D_out_name" ]] ; then
        echo 'Error: PB3D output file does not exist.'
        display_usage
        exit 1
    fi
else
    echo 'Error: PB3D Directory does not exist.'
    display_usage
    exit 1
fi
#
# Make new output folder if requested
if [ "$use_out_loc" = true ]; then
    out=$out_loc
else
    out=${1%/}
fi
mkdir -p $out $out/Plots $out/Data $out/Scripts &&
{
# success
echo "Working in directory $out/"
echo ""
# Copy inputs and the program
cp $POST_in_name $out
cp ../POST $out
chmod +x $out/POST
if [ "$use_out_loc" = true ]; then
    cp ${1%/}/$PB3D_out_name $out
fi
cd $out
echo "mpirun -np $2 $debug_opt $extra_debug_opt ./POST $POST_in_name PB3D_out.h5 ${@:3}" > command_POST
mpirun -np $2 $debug_opt $extra_debug_opt ./POST $POST_in_name $PB3D_out_name  ${@:3}
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
