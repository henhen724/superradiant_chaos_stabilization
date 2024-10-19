#!/bin/bash

########################################################
# Thank you for agreeing to do this! I've set it up as #
# simple as possible: f90 code needs compiling once,   #
# then you just run this for various omtil/subgrid.    #
########################################################

# This is the only stuff you should ever need to change:
# this value here
OmegaTil=1.900

# and picking a subgrid
subgrid="even"
#subgrid="odd"

#########################################################
# From here it's all stuff you don't gotta worry about  #
#########################################################

gridfile="Pgrid-${subgrid}.dat"
nmax=1

#read in the P grid & queue code
while read line || [[ -n "$line" ]];
do
    qsub run-timeevo.sh $OmegaTil $line
    #./run-timeevo.sh $OmegaTil $line
    nmax=$nmax+1
done < "$gridfile"

echo "All done- there should now be" $(echo "${nmax}-1" | bc) "jobs in the queue."

#Results management stuff
resultsgohere="OmTil=${OmegaTil}-${subgrid}-results"
echo "Once the codes have run please place all the results folders in this master folder: ${resultsgohere}"
if [[ -d "$resultsgohere" ]]
then
    echo "The folder $resultsgohere exists, just pointing it out to make sure we didn't already run these parameters."
else
    mkdir $resultsgohere
fi

#Easter egg
whichdo=$(( $RANDOM % 3 + 1 ))
case $whichdo in
    "1") echo "Thank you Kristín <3" ;;
    "2") echo "You're the best Kristín!";;
    "3") echo "★ You're a star Kristín! ★ ";;
esac
