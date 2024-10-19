#!/bin/bash

# This script explores the irregular behaviour seen for some of the
# parameters based on Esslinger

export GFORTRAN_UNBUFFERED_ALL="Y"

E_0=40
Kappa=8.1 
N_A=100000 
Omega_r=0.05
nmax=10

#RunType='time'

RunType='sweep'
SweepVariable='omegatil'
SweepDir='down'

OmegaTil=2.1
P=0.6

OmegaTil_init=1.9
OmegaTil_end=2.1
OmegaTil_no=2

Pump_init=$P
Pump_end=$P
Pump_no=1

# Linear stability parameters:
lmax=10
imposed_nmax=1001

. common-functions.sh

main()
{ 
    if [ $RunType == "sweep" ]; then
	if [ $SweepVariable == "omegatil" ]; then
	    NEWDIR="Esslinger_${RunType}_${SweepVariable}_P=${P}"
	else
	    NEWDIR="Esslinger_${RunType}_${SweepVariable}_OmegaTil=${OmegaTil}"
	fi
    else
	NEWDIR="Esslinger_${RunType}_Omtil=${OmegaTil}_P=${P}"
    fi
    other_params
    setup_new_directory
    create_and_run_script



    cd ..
}

main
