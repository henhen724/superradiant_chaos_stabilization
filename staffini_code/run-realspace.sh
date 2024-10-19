#!/bin/bash
#$ -S /bin/bash                 
#$ -cwd

#################################
# PARAMETERS FOR SIMULATION

OMP_NUM_THREADS=2

# This doesn't go into the code, it's just so I remember to friggin change the time label
tname=20

# Separate control of short range and long range contributions
shortrange=1
longrange=1
# Control whether to use flat profile [yes for flat, otherwise gaussian]:
flat=0
# Control harmonic trap on/off:
extpot=1

# Control whether we assume interaction w mirror image:
mirror=1

# Control whether we `correct' the effective potential
balancing=0

## DEBUGGING STUFF - use one or the other
# turns off interactions if 0:
interaction=1
# full debug -- this parameter is not actually anywhere in xmds code
debugging=0
# If still debugging non-bessel:
ISBESSEL=0

## GEOMETRY PARAMS
# so many gridpoints on a domain (-range, range)
grid=$3
range=$4

### PHYSICAL PARAMS
# pump str: on kelpie take from qsub code
fP=$1
# detuning: on kelpie take from qsub code
Delta0=$2

# atom position
zA=0
# atom cloud width
sigmaA=$5
# pump width
sigmaP=2

# confocality parameter
epsilon=50
# stength of interaction
E0=-50
# Oscillator length to beamwaist ratio
eta=1
# modes kept
modes=50
# atomic oscillator freq
omegaT=1
# cavity losses
kappamin=1000
# Atom-atom interaction -- this is outdated
U=0.001
# No of atoms
N=100000
# eta is lho/(wA/sqrt(2))
#eta=1

# Custom name to append to the end of the filename
#CUSTOM=""
CUSTOM="_cloud-${sigmaA}_pump-${sigmaP}-FEW-MODES-NEG-E"

#################################

. common-functions.sh

set_host_and_prog

# Single brakets are a shortcut for 'test', much more
# portable but doesn't parse multiple inputs which strings containing
# special characters may be interpreted as. Here, "1Drealspace"* is
# throws "too many argument error". We need instead [[]], which are a
# shortcut for 'new test' (bash and a few other specific)

echo $PROG
BASE=$(basename $PROG)

# If non-weak coupling
if [[ "$BASE" == "realspace-1D"* ]]
then
    DIM="1D"
else
    DIM="2D"
fi

if [ "$BASE" = *"-atoms-only" ]
then
    EOMS="atoms"
else
    EOMS="both"
fi

# If weak coupling
if [[ "$BASE" == "weak-coupling"* ]]
then
    WEAKCOUPLING="True"
    DIM="1D"
    DIMSUFFX="1D_"
else
    WEAKCOUPLING="False"
fi

#################################
# NAMING

if [ $ISBESSEL = 1 ]; then
    BESSLSUFFX="_R"
else
    BESSLSUFFX="_X"
fi

if [[ "$DIM" == "1D" ]]
then
    DIMSUFFX="_1D${BESSLSUFFX}"
else
    DIMSUFFX=""
fi

if [ $(($shortrange+$longrange)) = 2 ]; then
    INTSUFFX=""
elif [ $(($shortrange+$longrange)) = 1 ];then
    if [ $shortrange = 1 ]; then
	INTSUFFX="_SHORTRANGE"
    else
	INTSUFFX="_LONGRANGE"
    fi
else
    INTSUFFX="_NOINT"
fi


if [ $flat = 0 ]; then
    gaussian=1
    DISTRSUFFX="_GAUSSIAN"
else
    gaussian=0
    DISTRSUFFX="_FLAT"
fi

if [ $extpot = 1 ]; then
    POTSUFFX=""
else
    POTSUFFX="_NO-POT"
fi

if [[ $mirror = 0 ]] && [[ "$DIM" == "1D" ]];
then
    MIRRORSUFFX="_NO-MIRROR"
else
    MIRRORSUFFX=""
fi

if [[ $flat = 1 ]] && [[ $WEAKCOUPLING = "True" ]];
then
    DISTRSUFFX="_FLAT"
else
    DISTRSUFFX=""
fi

if [ $balancing = 1 ];
then
    POTSUFFX="${POTSUFFX}_NO-VEFF"
fi


#################################

main() {

    pad_down
    
    set_name_by_run_type
     
    # OTHER CUSTOM NAMES
    VERBOSENAME="zA=${zA}_epsilon=${epsilon}"

    SIMVALS="${grid}-grid_${range}-range_tmax=${tname}_"
    #SIMSETUP=$BESSLSUFFX$DISTRSUFFX$POTSUFFX$INTSUFFX$MIRRORSUFFX
    SIMSETUP=$DIMSUFFX$POTSUFFX$INTSUFFX$MIRRORSUFFX$DISTRSUFFX

    NEWDIR=$NEWDIR$SIMVALS$VERBOSENAME$SIMSETUP$CUSTOM
    #NEWDIR=$NEWDIR$SIMVALS$SIMSETUP$CUSTOM

    
    echo "Parameters set: fP=${fP}, Delta0=${Delta0}, zA=${zA}, epsilon=${epsilon}"
    echo $NEWDIR
    setup_dir_and_run

    pad_up
}

main
