#!/bin/bash
declare -a names=("0.3" "0.4" "0.5" "0.6" "0.7" "0.8" "0.9" "1.0" "1.5" "2.0" "3.0" "4.0")
declare -a ripples=("0.0" "1.0" "2.0" "3.0" "4.0" "5.0" "6.0")
for name in "${names[@]}"
do
    for ripple in "${ripples[@]}"
    do
        ln -sf /home/ITER/weyenst/Programs_MPICH3.1.3/VMEC_simulations/Hmode_PHD/ped${name}/VMEC/${ripple}/wout_Hmode_PHD_ripple${ripple}_ped${name}.nc /home/ITER/weyenst/Programs_MPICH3.1.3/PB3D/Run/wout_Hmode_PHD_ripple${ripple}_ped${name}.nc
    done
done

# OLD PHD equilibria
declare -a names=("0.3" "0.4" "0.5" "0.6" "0.7")
declare -a ripples=("0.001" "0.002" "0.003" "0.004" "0.005" "0.006" "0.007" "0.008" "0.009" "0.010")
for name in "${names[@]}"
do
    for ripple in "${ripples[@]}"
    do
        cp /home/ITER/weyenst/Archive/Run/Hmod_ripple_150/Hmode_ped${name}_ripple16/${ripple}/wout_Hmode_ped${name}_ripple16_${ripple}.nc /home/ITER/weyenst/Programs_MPICH3.1.3/PB3D/Run/wout_Hmode_OLD_ripple${ripple}_ped${name}.nc
    done
done

declare -a names=("0.3" "0.4" "0.5" "0.6" "0.7")
declare -a ripples=("0.000")
for name in "${names[@]}"
do
    for ripple in "${ripples[@]}"
    do
        cp /home/ITER/weyenst/Archive/Run/Hmod_ripple_150/Hmode_ped${name}_ripple16/${ripple}/wout_Hmode_ped${name}.nc /home/ITER/weyenst/Programs_MPICH3.1.3/PB3D/Run/wout_Hmode_OLD_ripple${ripple}_ped${name}.nc
    done
done

