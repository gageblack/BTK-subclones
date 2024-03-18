#!/bin/bash

PATIENT=$1
SAMPLE=$2

if [ ! $# -eq 2 ]; then
	echo "You need to give PATIENT and SAMPLE. Try again."
	exit
fi

echo "Loading modules..."
module load gcc/10.2.0 htslib/1.9 scbayes/1.0.0
cd $PATIENT
YML=${SAMPLE}.assign.yml
GEN=${SAMPLE}.genotype

echo "Running scAssign.."
scAssign -q $MIN_QUAL $YML $GEN > ${SAMPLE}.assigned

echo "Parsing the output..."
cut -f 2 ${SAMPLE}.assigned | sort | uniq -c > ${SAMPLE}.assigned.stats

cd ..
