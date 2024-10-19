#!/bin/bash

for x in *Delta0=-*.png;
do
    neg=$x
    pos="${x/Delta0=-/Delta0=}"

    #merged="${x/Delta0=-*tmax=/tmax=}"
    
    #Verbose option leaving delta in:
    merged="${pos/128-grid*tmax=/tmax=}"
    merged="${merged/Delta0/modDelta}"

    merged="${merged/_1D_*GAUSSIAN_/_}"
    merged="${merged/output/compare}"
    convert $neg $pos +append $merged
    echo "Combined $merged..."
done

echo "Delete the single png files? [y/n]"
read query

if [[ $query = "y" ]]; then
    rm -f ./output*png
fi
