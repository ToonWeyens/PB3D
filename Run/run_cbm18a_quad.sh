#!/bin/bash
# Display usage function
!!!!!!!!!!!!!!!! NEEDS TO BE UPDATED WITH DRMEMORY AND OTHER THINGS !!!!!!!!!!!
display_usage() { 
    echo -e "\nUsage:\n$0 [OPTS] NR_PROCS \n" 
    echo -e "    OPTS: -o specify output name"
    echo -e "          -d use Valgrind debugging"
    echo -e "          -s trace sources of errors in Valgrind\n"
    echo -e "          -l check leaks in Valgrind\n"
    } 
#
# Setting some variables
#slepc_opt="-st_pc_factor_shift_type NONZERO -st_pc_type lu -st_pc_factor_mat_solver_package mumps -eps_monitor"
slepc_opt="-eps_monitor"
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
            debug_opt="valgrind --db-attach=yes"
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
cp cbm18a $out
cp ../PB3D_PREP $out
chmod +x $out/PB3D_PREP
cd $out
rm -f .lock_file*
echo "creating pbs script"
# Create pbs script
rm -f PB3D_PREP.pbs
cat > PB3D_PREP.pbs << END
#!/bin/sh
#PBS -N $out
#PBS -l nodes=3:ppn=$1
#PBS -l mem=30GB
#PBS -q parallel_16
#PBS -l walltime=04:00:00
#PBS -o $(pwd)/PB3D_PREP.o
#PBS -e $(pwd)/PB3D_PREP.e
export PATH="/share/apps/gcc/4.6.4/bin:/share/apps/openmpi/1.6.5/gcc-4.6.4/bin:$PATH"
export LD_LIBRARY_PATH="/share/apps/gcc/4.6.4/lib:/share/apps/gcc/4.6.4/lib64:/share/apps/openmpi/1.6.5/gcc-4.6.4/lib:usr/lib64:$LD_LIBRARY_PATH"
cd $(pwd)
echo "job is run on"
pbsdsh uname -n
. /opt/torque/etc/openmpi-setup.sh
echo "mpirun -np $1 $debug_opt $extra_debug_opt ./PB3D_PREP input_cbm18a cbm18a $slepc_opt ${@:2} --no_plots" > command
mpirun -np $1 $debug_opt $extra_debug_opt ./PB3D_PREP input_cbm18a cbm18a $slepc_opt ${@:2} --no_plots
exit
END
chmod u+x PB3D_PREP.pbs
echo "Submitting job to"
qsub PB3D_PREP.pbs
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
