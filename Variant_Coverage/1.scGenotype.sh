PATIENT=$1
SAMPLE=$2
VCF=$3
BARCODES=$4
BAM=$5

#if [ $PATIENT == "" ] || [ $SAMPLE == "" ]
if [ ! $# -eq 2 ]; then
	echo "You need to give a PATIENT and SAMPLE. Try again."
	exit
fi

echo "Running scGenotype for $SAMPLE from $PATIENT"

echo "Loading modules..."
module load gcc/10.2.0 htslib/1.9 scbayes/1.0.0

echo "Running scGenoytpe.."
scGenotype $BAM <(cat $BARCODES) $VCF > genotypes/${SAMPLE}.genotype

echo "Parsing the output..."
python 1a.parse_variants.py genotypes/${SAMPLE}.genotype
