#!/bin/bash

pad_down() {
    echo "##############################"
    echo
}

pad_up() {
    echo
    echo "##############################"
}

# Creates PARAMS file
dump_params() {
    cat << EOF > params
# Variables #
fP=$fP    
Delta0=$Delta0    
sigmaA=$sigmaA
sigmaP=$sigmaP
zA=$zA  
epsilon=$epsilon  

# Constant parameters #
E0=$E0    
eta=$eta  
omegaT=$omegaT    
kappamin=$kappamin  
No of atoms=$N   
modes=$modes

# Simulation values #
grid points=$grid    
domain=(-${range},${range})

# Logical simulation values #
interaction=$interaction
shortrange=$shortrange
longrange=$longrange 
gaussian=$gaussian
flat=$flat
extpot=$extpot
mirror=$mirror
balancing=$balancing
EOF
}

# Sets hostname variable and program name variable
set_host_and_prog() {
    HOST=$(hostname)
    if [[ "$HOST" = "compute"* ]] || [[ "$HOST" = "kelpie.st-and.ac.uk" ]]
    then
	export LD_LIBRARY_PATH=/share/apps/gcc/4.6.2//lib/:/share/apps/gcc/4.6.2//lib64:{LD_LIBRARY_PATH}
    fi
    SOURCES=`find . -maxdepth 1 -name *.cc`
    if [ "$SOURCES" = "" ]
    then
	pad_down
	echo "ERROR: Please compile something first."
	pad_up
	exit
    fi
    PROG=${SOURCES/.cc/}
    
}

# Picks out the correct run type and sets the NEWDIR parameter, plus prints out appropriate messages
set_name_by_run_type() {

    NEWDIR="fP=${fP}_Delta0=${Delta0}_"

    # 1D, atoms and light field:
    if [ $WEAKCOUPLING = "False" ]
    then
	if [ $DIM = "1D" ] && [ $EOMS = "both" ]
	then
	    echo "RUNNING 1D VERSION OF CODE;"
	    NEWDIR="${NEWDIR}"
	    # 1D, atoms only:
	elif [ $DIM = "1D" ] && [ $EOMS = "atoms" ]
	then
	    if [ $epsilon -ne 0 ]
	    then
		echo "ERROR - You cannot run atoms-only EoM out of confocality - substitution breaks."
		exit
	    fi
	    echo "RUNNING 1D VERSION OF CODE FOR ATOMIC EOM ONLY;"
	    NEWDIR="ATOMS_${NEWDIR}"
	    # 2D, atoms only:
	elif [ $DIM = "2D" ] && [ $EOMS = "atoms" ]
	then
	    if [ $epsilon -ne 0 ]
	    then
		echo "ERROR - You should not run atoms-only EoM out of confocality."
		exit
	    fi
	    echo "RUNNING CODE FOR ATOMIC EOM ONLY;"
	    NEWDIR="ATOMS_${NEWDIR}"
	else
	    # 2D, atoms and light field: 
	    echo "RUNNING FULL CODE;"
	fi
    else
	echo "RUNNING WEAK COUPLING CODE IN 1D:"
	NEWDIR="WEAK-CPL_${NEWDIR}"
    fi
    # Debugging warnings:

    if [ $debugging -eq 1 ]
    then
	echo "WARNING - You're debugging something and the code won't be running in full."
	echo
	str="DEBUGGING-"
	NEWDIR=$str$NEWDIR
    fi
    
    if [ $interaction -eq 0 ]
    then
	echo "WARNING - You're running your simulation with interactions turned off.";

	str="TESTING-SHO-ONLY-"
	NEWDIR=$str$NEWDIR
    fi

    echo
    
}

# Copies over the appropriate Mathematica notebook template and renames it 
make_notebook() {

    if [[ "$HOST" = "mls-workstation" ]]
    then
	if [[ "$DIM" = "1D" ]]
	then
	    if [[ $WEAKCOUPLING = "True" ]]
	    then	
		cp ../mathematica/weak-coupling.nb .
		mv weak-coupling.nb ${outname}.nb
	    else
	    	cp ../mathematica/1Dlistplottingtemplate.nb .
		mv 1Dlistplottingtemplate.nb ${outname}.nb
	    fi
	else
	    cp ../mathematica/listplottingtemplate.nb .
	    mv listplottingtemplate.nb ${outname}.nb
	fi
	
	echo "(I made you a mathematica notebook too. I know you like that sort of thing.)"
    fi
    
}

# Tidies up outputs in folder by renaming them appropriately
rename_output() {

    outname="output_${NEWDIR}"
    outxil=`find . -maxdepth 1 -name *.xsil`
    outh5=`find . -maxdepth 1 -name *.h5`
    mv $outxil ${outname}.xsil
    mv $outh5 ${outname}.h5
  
    make_notebook
    
}

# Sets up new folder for results, clears if existing, dumps new parameters in it, and runs things
setup_dir_and_run() {

    mkdir -p $NEWDIR
    cd $NEWDIR/
    
    dump_params

    ln -sf ../${BASE} ./
    
    rm -f *.xsil *.h5

    ./${BASE} --fP $fP --Delta0 $Delta0 --zA $zA --sigmaA $sigmaA --sigmaP $sigmaP --epsilon $epsilon --E0 $E0 --modes $modes --omegaT $omegaT --kappamin $kappamin --U $U --N $N --eta $eta --interaction $interaction --grid $grid --range $range --shortrange $shortrange --longrange $longrange --gaussian $gaussian --flat $flat --extpot $extpot --mirror $mirror --balancing $balancing

    echo ""
    echo "That was ${NEWDIR}."

    rename_output

    cd ..

    echo
    echo "Looks like we survived this one."
}


##############################################
# EXPLICIT LINKING IF NEEDED
# mls-workstation
link_workstation() {
    
     /usr/bin/g++  -fstack-protector-strong -fopenmp -ftree-vectorize -march=native -mtune=native -O3 -ffast-math -funroll-all-loops -fomit-frame-pointer -falign-loops -fstrict-aliasing -momit-leaf-frame-pointer   -I/usr/lib/python2.7/dist-packages/xpdeint/includes -I/usr/include/atlas -I/usr/include/hdf5/serial -DHAVE_HDF5_HL -DHAVE_H5LEXISTS -D_LARGEFILE64_SOURCE -D_LARGEFILE_SOURCE -D_FORTIFY_SOURCE=2  ${PROG}.cc -c -o ${PROG}.cc.1.o  -lsz
    
    /usr/bin/g++ -Wl,-Bsymbolic-functions -Wl,-z,relro -Wl,-rpath -Wl,/usr/lib/x86_64-linux-gnu/hdf5/serial -fopenmp  ${PROG}.cc.1.o -o /home/mls/Code/multimode-mls/${PROG}     -Wl,-Bstatic  -lhdf5_hl -lhdf5 -lfftw3_omp -lfftw3 -Wl,-Bdynamic -L/usr/lib/x86_64-linux-gnu/hdf5/serial -L/usr/lib/x86_64-linux-gnu/hdf5/serial -lpthread -lz -ldl -lm   -lsz
    
   
}

# Kelpie
link_cluster() {

XPDEINT=/export/apps/xmds-2.2.2/xpdeint/
FFTW=/export/apps/xmds-2.2.2/fftw-3.3.4/
HDF5=/export/apps/xmds-2.2.2/hdf5-1.8.8/


/usr/bin/g++  -fstack-protector -fopenmp -ftree-vectorize -march=native -mtune=native -O3 -ffast-math -funroll-all-loops -fomit-frame-pointer -falign-loops -fstrict-aliasing -momit-leaf-frame-pointer   -I${XPDEINT}/includes -I/usr/include/atlas -I${FFTW}/include -I${HDF5}/include -DHAVE_HDF5_HL -DHAVE_H5LEXISTS -D_FILE_OFFSET_BITS=64 -D_LARGEFILE_SOURCE -D_LARGEFILE64_SOURCE -D_BSD_SOURCE -D_FORTIFY_SOURCE=2  ${PROG}.cc -c -o ${PROG}.cc.1.o


/usr/bin/g++ -Wl,-Bsymbolic-functions -Wl,-z,relro -Wl,-rpath -Wl,/usr/lib/i386-linux-gnu -fopenmp  ${PROG}.cc.1.o -o ${PROG} -Wl,-rpath,/usr/lib/atlas-base -Wl,-rpath,${FFTW}/lib    -Wl,-Bstatic  -lhdf5_hl -lhdf5 -lfftw3_omp -lfftw3 -Wl,-Bdynamic -L/usr/lib/atlas-base -L${FFTW}/lib -L${HDF5}/lib -L/usr/lib/i386-linux-gnu -L/usr/lib/i386-linux-gnu -lpthread -lz -ldl -lm  


    }

