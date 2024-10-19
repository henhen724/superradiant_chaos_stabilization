#!/bin/bash

# This script explores the irregular behaviour seen for some of the
# parameters based on Esslinger.

# In particular, this script uses a table of values for P, and calculating
# omtil with the lin stab code, to get a collection of sweeps around the
# section of the phase diagram we know the tricritical point should be in.

# This script produces multiple folders, one for each initial pump input,
# containing each sweep.

# From input:
P=$1
Pump_init=$P
Pump_end=$P
Pump_no=5

# Other parameters:

E_0=40
Kappa=8.1 
N_A=100000 
Omega_r=0.05
nmax=10

RunType='tricrit'
SweepVariable='Pump'
SweepDir='updown'

# Omega is overwritten anyway:
OmegaTil=1.000
OmegaTil_init=$OmegaTil
OmegaTil_end=$OmegaTil
OmegaTil_no=1


# Linear stability parameters:
lmax=10
imposed_nmax=1001

. common-functions.sh

main()
{ 
    NEWDIR="Esslinger_Tricrit_Sweep_P=${P}"

    other_params
    setup_new_directory
    create_and_run_script

    cd ..
}

main

# NOTE TO SELF I added this later to automate running on personal
# but means if I uncomment this cannot use kelpie code
# gridfile="Pgrid.dat"
# while read line || [[ -n "$line" ]]; # short-circuit boolean*
# do
#    P=$line
#    Pump_init=$P
#    Pump_end=$P
#    Pump_no=5
#
#    main
# done < "$gridfile"
#
# *Tests both conditions (read condition *and* input to line variable) so that loop reads last line of file too.
