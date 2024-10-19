#!/bin/bash
#$ -S /bin/bash                 
#$ -V -cwd

HOST=$(hostname)
SOURCES=`find . -maxdepth 1 -name *.cc`
PROG=${SOURCES/.cc/}

. common-functions.sh

main() {

    if [[ "$HOST" = "mls-workstation" ]]
    then
	echo "I recognise ${HOST} so I'm gonna do my thing."
	echo
	link_workstation
    elif [[ "$HOST" = "compute"* ]] || [[ "$HOST" = "kelpie.st-and.ac.uk" ]]
    then
	echo "I recognise ${HOST} so I'm gonna do my thing."
	echo
	link_cluster
    else
	echo "I don't know who you are, $HOST."
	exit 1
    fi
}

main
