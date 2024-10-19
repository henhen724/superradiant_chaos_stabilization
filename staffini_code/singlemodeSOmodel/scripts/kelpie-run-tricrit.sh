#!/bin/bash

gridfile="Pgrid-tricrit.dat"
nmax=1

while read line || [[ -n "$line" ]];
do
    qsub run-tricrit-code.sh "$line"
    nmax=$nmax+1
done < "$gridfile"

echo "All done- there should now be" $(echo "${nmax}-1" | bc) "jobs in the queue."
