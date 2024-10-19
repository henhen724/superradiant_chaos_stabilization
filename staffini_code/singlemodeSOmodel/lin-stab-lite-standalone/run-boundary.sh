#! /bin/bash

#################################################
# LINEAR STABILITY FILENAME SETUP
#################################################

filename() {
if [ $imposed_nmax == 1 ]; then 
    filename="boundary_${setup}_lmax=${lmax}_Dicke_coeffs.dat"
else
    filename="boundary_${setup}_lmax=${lmax}.dat"
fi 
}

##################################################
# MODEL SETUP
##################################################

# P interval and number of points:
Pmin=0
Pmax=0.7
nP=605

# System parameters for Esslinger setup:
setup="Esslinger"
E_0=40  # Note: E_0=N*(g**2)/(4*Delta_a)
Kappa=8.1 
omega_r=0.05
# N_A=100000 

# Sum over l for states runs up to:
lmax=10

#Truncation of the matrix yielding coeffs:
imposed_nmax=1001

##################################################
# OUTPUT SETUP
##################################################

# Set file names for output:
filename
echo $filename

cat <<EOF > PARAMS
$Pmin, $Pmax, $nP
$E_0, $Kappa, $omega_r
$lmax, $imposed_nmax
$filename
EOF

# Running code with parameters as specified above:
./find_transition < PARAMS
