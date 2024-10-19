setup_new_directory() {

    mkdir -p $NEWDIR
    cd $NEWDIR/

    ln -sf ../../main ./
    
    rm -f *.dat

}

create_and_run_script() {
    cat << EOF > $CONTROL
$RunType, $SweepVariable, $SweepDir
$nmax, $E_0, $Kappa, $N_A, $Omega_r
$Pump_init, $Pump_end, $Pump_no
$OmegaTil_init, $OmegaTil_end, $OmegaTil_no
$tmin, $tmax, $tstep
$convergence_threshold
$lmax, $imposed_nmax
EOF

./main

}

cleanup_data() {
    tar -czf "$NEWDIR-Evol.tar.gz" UP* DN*
    rm UP* DN*
}

other_params() {

CONTROL=params

tmin=0.0
tmax=5.0E6
tstep=2

convergence_threshold=0.002

}
