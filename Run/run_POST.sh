#!/bin/bash
# set machine ID
machine_ID=0
name=$(uname -n)
uname -n | grep toon &> /dev/null && machine_ID=1
uname -n | grep uranus &> /dev/null && machine_ID=2
if (($machine_ID==0)); then
    echo -e "ERROR: invalid machine $(uname -n)"
    exit 1
else
    echo -e "for user \"$(uname -n)\", machine ID: $machine_ID"
    echo -e ""
fi
# Display usage function
display_usage() { 
    echo -e "\nUsage:\n$0 [OPTS] PB3D_DIR NR_PROCS \n" 
    echo -e "    OPTS:      -o [NAME] specify output name"
    echo -e "               -d use Dr. Memory debugging"
    echo -e "               -h help"
    echo -e ""
    echo -e "    PB3D_DIR:  PB3D directory"
    echo -e ""
    echo -e "    NR_PROCS:  nr. of MPI processes"
    echo -e ""
    } 
#
# Setting some variables
debug_opt=""
extra_debug_opt=""
n_opt_args=0
use_out_loc=false
n_procs=1
PB3D_out_name="PB3D_out.h5"
input_name="input_POST"
case $machine_ID in
    1)                                                                          # laptop
        drmemory_location="/opt/DrMemory-Linux-1.10.1-3/bin64"
    ;;
    2)                                                                          # uranus
        drmemory_location="/home/tweyens/Programs/DrMemory-Linux-1.10.1-3/bin64"
    ;;
esac
#
# Catching options
while getopts :o:n:dh opt; do
    case $opt in
        o)                                                                      # output file
            out_loc=${OPTARG%/}
            use_out_loc=true
            n_opt_args=$((n_opt_args+2))                                        # 2 arguments
        ;;
        n)                                                                      # nr. of nodes on server
            n_procs=${OPTARG%/}
            n_opt_args=$((n_opt_args+2))                                        # 2 arguments
            if (($machine_ID==2)); then                                         # uranus
                echo -e "using $n_procs processe(s) per node"
            else
                echo -e "WARNING: ignoring option -n because not on Uranus"
            fi
        ;;
        d)
            debug_opt="$drmemory_location/drmemory --"
            n_opt_args=$((n_opt_args+1))                                        # 1 argument
        ;;
        h)
            echo -e "program PB3D"
            echo -e "written by Toon Weyens"
            display_usage
            exit 1
        ;;
        \?)
            echo -e "ERROR: invalid option -$OPTARG"
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
#
# Shift arguments to skip options
shift $n_opt_args
# Check existence of PB3D folder
#
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
#
if [ -n "$debug_opt" ]; then
    echo "Using Dr. Memory for debugging"
fi
#
# Make new output folder if requested
if [ "$use_out_loc" = true ]; then
    out=$out_loc
else
    out=${1%/}
fi
out_full=$(pwd)/$out
mkdir -p $out_full $out_full/Plots $out_full/Data $out_full/Scripts &&
{
# success
echo -e "input file:       " $input_name
echo -e "PB3D output file: " ${1%/}/$PB3D_out_name
echo -e ""
echo -e "Working in directory $out_full/"
echo -e ""
# Copy inputs and the program
cp $input_name $out_full
cp ../POST $out_full
chmod +x $out_full/POST
if [ "$use_out_loc" = true ]; then
    cp ${1%/}/$PB3D_out_name $out
fi
cd $out_full
# Do actions depending on machine ID
case $machine_ID in
    1)                                                                          # laptop
        echo "$debug_opt mpirun -np $2 ./POST $input_name $PB3D_out_name ${@:3}" > command_POST
        $debug_opt mpirun -np $2 ./POST $input_name $PB3D_out_name ${@:3}
    ;;
    2)                                                                          # uranus
        # get memory
        max_tot_mem=$(grep $input_name -e 'max_tot_mem_per_proc' | tr -dc '0-9')
        if [ -z "$max_tot_mem" ]; then
            # PB3D default
            max_tot_mem=6000
        else
            # multiply by nr. of procs.
            max_tot_mem=$(($2*$max_tot_mem))
        fi 
        mem_unit='kb'
        # Create pbs script
        echo "creating pbs script"
        rm -f POST.pbs
        cat > POST.pbs << END
#!/bin/sh
#PBS -N $out
#PBS -o $out_full/POST.o
#PBS -e $out_full/POST.e
#PBS -l nodes=$2:ppn=$n_procs
#PBS -l pvmem=$max_tot_mem$mem_unit
#PBS -l walltime=04:00:00
#PBS -m abe
#PBS -M toon.weyens@gmail.com
cd $out_full
rm -f POST.o POST.e
echo "Job id:   \$PBS_JOBID"
echo "Job name: \$PBS_JOBNAME"
echo "Job node: \$PBS_QUEUE"
echo "Walltime: \$PBS_WALLTIME"
echo ""
export PATH="/share/apps/openmpi-1.10.1/bin:/share/apps/gcc-5.2/bin:$PATH"
export LD_LIBRARY_PATH="$HOME/lib:/share/apps/openmpi-1.10.1/lib:/share/apps/gcc-5.2/lib64:/share/apps/gcc-5.2/lib:$LD_LIBRARY_PATH"
pbsdsh uname -n
. /opt/torque/etc/openmpi-setup.sh
echo "$debug_opt mpirun -np $2 ./POST $input_name $PB3D_out_name ${@:3}" > command_POST
$debug_opt mpirun -np $2 ./POST $input_name $PB3D_out_name ${@:3}
exit
END
        chmod u+x POST.pbs
        qsub POST.pbs
    ;;
esac
cd ../
echo ""
echo "Leaving directory $out_full/"
echo ""
#
# failure
} || {
echo "Error: unable to create directory $out_full/"
echo "Do you have the correct permissions?"
}
