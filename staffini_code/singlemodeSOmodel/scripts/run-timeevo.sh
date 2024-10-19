#!/bin/bash

# This script explores the irregular behaviour seen for some of the
# parameters based on Esslinger


E_0=40
Kappa=8.1 
N_A=100000 
Omega_r=0.05
nmax=10

RunType='single'

# SET PARAMS FROM INPUT
OmegaTil=$1
P=$2

# Unused vars that still go in (frown)
SweepVariable='None'
SweepDir='None'
OmegaTil_init=$OmegaTil
OmegaTil_end=$OmegaTil
OmegaTil_no=1
Pump_init=$P
Pump_end=$P
Pump_no=1
lmax=10
imposed_nmax=1001

. common-functions.sh

main() {
    NEWDIR="Esslinger_${RunType}_OmTil=${OmegaTil}_P=${P}"

    other_params
    setup_new_directory
    create_and_run_script

    cd ..
}

main
