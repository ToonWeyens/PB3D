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
    echo -e "\nUsage:\n$0 [OPTS] CASE NR_PROCS \n" 
    echo -e "    OPTS:      -o [NAME] specify output name"
    echo -e "               -d use Dr. Memory debugging"
    echo -e "               -h help"
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
    echo -e ""
    echo -e "    NR_PROCS:  nr. of MPI processes"
    echo -e ""
    } 
#
# Setting some variables
#slepc_opt="-st_pc_factor_shift_type NONZERO -st_pc_type lu -st_pc_factor_mat_solver_package mumps -eps_monitor -eps_view"
slepc_opt="-st_pc_type lu -st_pc_factor_mat_solver_package mumps -eps_monitor -eps_ncv 100 -eps_mpd 100"
#slepc_opt="-st_pc_type jacobi -st_pc_factor_mat_solver_package mumps -eps_monitor -eps_ncv 100 -eps_mpd 100"
debug_opt=""
n_opt_args=0
use_out_loc=false
n_procs=1
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
        d)                                                                      # debug
            debug_opt="$drmemory_location/drmemory --"
            n_opt_args=$((n_opt_args+1))                                        # 1 argument
        ;;
        h)                                                                      # help
            echo -e "program PB3D"
            echo -e "written by Toon Weyens"
            display_usage
            exit 1
        ;;
        \?)                                                                     # other
            echo -e "ERROR: invalid option -$OPTARG"
            display_usage
            exit 1
        ;;
    esac
done
#
# Checking for number of input arguments (need at least 1 for # procs and 1 for case)
if [ "$#" -lt $((n_opt_args+2)) ]; then
    display_usage
    exit 1
fi
#
# Shift arguments to skip options
shift $n_opt_args
#
# Set case parameters: input_name
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
    41|42|43|44)
        input_name=Hmode
    ;;
    *)
        echo -e "ERROR: No case found"
        display_usage
        exit 1
    ;;
esac
# Set case parameters: eq_name
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
    *)
        echo -e "ERROR: No case found"
        display_usage
        exit 1
    ;;
esac
input_name=input_$input_name
#
# Shift arguments to skip case
shift 1
#
if [ -n "$debug_opt" ]; then
    echo "Using Dr. Memory for debugging"
fi
# Make new folder
if [ "$use_out_loc" = true ]; then
    out=$out_loc
else
    get_date=$(date +"%Y-%m-%d-%H-%M-%S")
    out=$get_date
fi
out_full=$(pwd)/$out
mkdir -p $out_full $out_full/Plots $out_full/Data $out_full/Scripts &&
{
# success
echo -e "input file:       " $input_name
echo -e "equilibrium file: " $eq_name
echo -e ""
echo -e "Working in directory $out_full/"
echo -e ""
# Copy inputs and the program
cp $input_name $out_full
cp $eq_name $out_full
cp ../PB3D $out_full
chmod +x $out_full/PB3D
cd $out_full
rm -f .lock_file*
# Do actions depending on machine ID
case $machine_ID in
    1)                                                                          # laptop
        echo "$debug_opt mpirun -np $1 ./PB3D $input_name $eq_name $slepc_opt ${@:2}" > command
        $debug_opt mpirun -np $1 ./PB3D $input_name $eq_name $slepc_opt ${@:2}
    ;;
    2)                                                                          # uranus
        # get memory
        max_tot_mem=$(grep $input_name -e 'max_tot_mem_per_proc' | tr -dc '0-9')
        if [ -z "$max_tot_mem" ]; then
            # PB3D default
            max_tot_mem=6000
        else
            # multiply by nr. of procs.
            max_tot_mem=$(($1*$max_tot_mem))
        fi 
        mem_unit='kb'
        # Create pbs script
        echo "creating pbs script"
        rm -f PB3D.pbs
        cat > PB3D.pbs << END
#!/bin/sh
#PBS -N $out
#PBS -o $out_full/PB3D.o
#PBS -e $out_full/PB3D.e
#PBS -l nodes=$1:ppn=$n_procs
#PBS -l pvmem=$max_tot_mem$mem_unit
#PBS -l walltime=04:00:00
#PBS -m abe
#PBS -M toon.weyens@gmail.com
cd $out_full
rm -f PB3D.o PB3D.e
echo "Job id:   \$PBS_JOBID"
echo "Job name: \$PBS_JOBNAME"
echo "Job node: \$PBS_QUEUE"
echo "Walltime: \$PBS_WALLTIME"
echo ""
export PATH="/share/apps/openmpi-1.10.1/bin:/share/apps/gcc-5.2/bin:$PATH"
export LD_LIBRARY_PATH="$HOME/lib:/share/apps/openmpi-1.10.1/lib:/share/apps/gcc-5.2/lib64:/share/apps/gcc-5.2/lib:$LD_LIBRARY_PATH"
pbsdsh uname -n
. /opt/torque/etc/openmpi-setup.sh
echo "$debug_opt mpirun -np $1 ./PB3D $input_name $eq_name $slepc_opt ${@:2}" > command
$debug_opt mpirun -np $1 ./PB3D $input_name $eq_name $slepc_opt ${@:2}
exit
END
        chmod u+x PB3D.pbs
        qsub PB3D.pbs
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
