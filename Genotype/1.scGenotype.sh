#!/bin/bash

PATIENT=$1
SAMPLE=$2
BAM=$3

if [ ! $# -eq 2 ]; then
	echo "You need to give PATIENT and SAMPLE. Try again."
	exit
fi

echo "Loading modules..."
module load gcc/9.2.0 htslib/1.9 scbayes/1.0.0

VCF=${PATIENT}/${PATIENT}.somatic.vcf
BARCODES=${PATIENT}/${SAMPLE}.barcodes_fixed.tsv

## Get and fix the barcode file if needed.
if [ ! -f $BARCODES ]; then
    echo "Fixing Barcodes"
    cut -d "-" -f 1 seurat/${SAMPLE}/genes_seurat/barcodes.tsv > $BARCODES
fi

#echo "Running scGenoytpe.."
scGenotype $BAM <(cat $BARCODES) $VCF > ${PATIENT}/${SAMPLE}.genotype

echo "Parsing the output..."
python 1b.parse_variants.py ${PATIENT}/${SAMPLE}.genotype

