#!/bin/bash

echo "Re-plot any existing? [y/n]"
read replotq

function plotit {
    echo $x
    
    if [[ $x == "WEAK"* ]]; then
	FILENAME="pngmaker-weak-coupling-output-X.m"
    else
	FILENAME="pngmaker-realspace-1D.m"
    fi

    cp ~/Code/multimode/mathematica/${FILENAME} .

    MathematicaScript -script ./${FILENAME}

    echo "================="
}

for x in *fP=*;
do
    if [[ $x != *png ]];
    then
	cd $x

	if [[ $replotq = "y" ]]; then
	    plotit
	else
	    filename="output_${x}.png"
	    echo $filename
	    if [ -f $filename ]; then
		echo "Png already plotted."
		echo "================="
	    else
		plotit
	    fi
	fi		
	cd ..
    fi
done

echo "Copy all png files just generated here? [y/n]"
read moveq

if [[ $moveq = "y" ]]; then
    for x in *fP=*;
    do
	if [[ $x != *png ]];
	then
	    cp $x/*png .
	fi
    done
fi
