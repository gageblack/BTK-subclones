import pandas as pd
import pysam
import sys
import time
from multiprocessing import Pool, cpu_count

def read_bed_file(bed_file):
    """Reads gene regions from a BED file into a DataFrame."""
    columns = ['chr', 'start', 'end', 'gene_name']
    df = pd.read_csv(bed_file, sep='\t', names=columns, header=None, dtype={'chr': str, 'start': int, 'end': int, 'gene_name': str})
    return df

def calculate_coverage_for_gene(args):
    """Calculates coverage for a single gene."""
    bam_file, gene_row, exon_df, min_coverage_fraction = args  # Unpack arguments

    gene_name = gene_row['gene_name']
    gene_start = gene_row['start']
    gene_end = gene_row['end']
    chr = gene_row['chr']
    result = []
    
    # Get Exons for the gene
    exons = exon_df[exon_df['gene_name'] == gene_name]
    # sort exons by start position
    exons = exons.sort_values(by='start', ascending=True)

    total_exonic_length = (exons['end'] - exons['start']).sum()

    with pysam.AlignmentFile(bam_file, "rb") as bam:
        # Fetch reads that overlap the gene region
        reads = bam.fetch(chr, gene_start, gene_end)

        for read in reads:
            # Initialize a variable to keep track of the read's coverage
            coverage = 0

            read_start = read.reference_start
            read_end = read.reference_end

            # skip reads where less than X% of the read overlaps with the gene
            if read_start < gene_start or read_end > gene_end:
                overlap_with_gene = abs(max(read_start, gene_start) - min(read_end, gene_end))
                if overlap_with_gene < (min_coverage_fraction * (read_end - read_start)):
                    continue

            for _, exon in exons.iterrows(): # The exons in the bed are in reverse order, so we iterate in reverse
                exon_start = exon['start']
                if exon_start > read_end:
                    break
                exon_end = exon['end'] 
                if exon_end < read_start:
                    continue               
                overlap_start = max(read_start, exon_start)
                overlap_end = min(read_end, exon_end)
                overlap = max(0, overlap_end - overlap_start)
                coverage += overlap
            
            if coverage == 0:
                continue

            # Calculate and print coverage percentage for the read
            if total_exonic_length > 0:  # Avoid division by zero
                coverage_percentage = (coverage / total_exonic_length) * 100
                result.append(coverage_percentage)
    return result

############################################

if len(sys.argv) != 6:
    print("Usage: python read_coverage_parallel.py <bam_file> <gene_file> <exon_file> <outfile_path> <min_coverage_fraction>")
    sys.exit(1)
bam_file = sys.argv[1]
gene_file = sys.argv[2]
exon_file = sys.argv[3]
outfile_path = sys.argv[4]
min_coverage_fraction = float(sys.argv[5])

print("Output file: ", outfile_path)

print("Reading gene and exon regions...")
genes_df = read_bed_file(gene_file)
exons_df = read_bed_file(exon_file)

print("Calculating coverage...")
start=time.time()
# Prepare arguments for parallel processing
tasks = [(bam_file, row, exons_df, min_coverage_fraction) for index, row in genes_df.iterrows()]
num_processes = cpu_count()

# Initialize pool with the desired number of workers
with Pool(processes=num_processes) as pool:  # Adjust the number of processes as necessary
    results = pool.map(calculate_coverage_for_gene, tasks)

# Write results to file
with open(outfile_path, 'w') as outfile:
    for result in results:
        for coverage_percentage in result:
            outfile.write(str(coverage_percentage)+"\n")

print("Time taken: ", round((time.time()-start)/60, 2), " minutes")