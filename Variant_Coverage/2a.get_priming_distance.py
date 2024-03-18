import pybedtools
import pandas as pd
import sys

# For each variant, figure out what gene it belongs to and which exons are in that gene.
# for 5', take the variant postition and subtract the start of the exon, and then add the length of any exon before it.
# for 3', take the variant position and subtract it from the end of the exon, and then add the length of any exon after it.

sample = sys.argv[1]

# Load the bed file
bed_file = pybedtools.BedTool('Variant_Coverage/meta/protein_coding_exons.bed')
genotype_file = 'Variant_Coverage/genotypes/'+sample+'.genotype.Variant_Coverage.tsv'
exon_df = pd.read_csv('Variant_Coverage/meta/protein_coding_exons.bed', sep='\t', header=None, names=['chrom', 'start', 'end', 'gene'])

# Read in protein coding gene info and pull out gene names.
pcg_info = pd.read_csv('Variant_Coverage/meta/protein_coding_genes.gtf', sep='\t', names=['chrom', 'source', 'feature', 'start', 'end', 'score', 'strand', 'frame', 'attribute'])
def extract_gene_name(attributes):
    """Extract the gene_name from the attributes column."""
    parts = attributes.split(';')
    gene_name_part = [part for part in parts if part.strip().startswith('gene_name')]
    if gene_name_part:
        # Extract the gene name, removing quotes
        return gene_name_part[0].split('"')[1]
    return None
pcg_info['gene_name'] = pcg_info['attribute'].apply(extract_gene_name)
pcg_info = pcg_info[['chrom', 'start', 'end', 'strand', 'gene_name']]

# Read in the TSV file to a pandas DataFrame
df = pd.read_csv(genotype_file, sep='\t')

# Split the variant column into separate columns
df[['chrom', 'position', 'ref_allele', 'alt_allele']] = df['Variant'].str.split(':', expand=True)
df['position'] = df['position'].astype(int)

# Add chr to the chromosome name if it is not already there
df['chrom'] = df['chrom'].apply(lambda x: x if x.startswith('chr') else 'chr' + x)

# Convert the variant positions into a BED format for intersection
df_bed = df[['chrom', 'position']].copy()
# rename the position column to start
df_bed.rename(columns={'position': 'end'}, inplace=True)
df_bed['start'] = df_bed['end'] - 1
df_bed = df_bed[['chrom', 'start', 'end']]

# Now convert the correctly formatted DataFrame to a BedTool object
variants_bed = pybedtools.BedTool.from_dataframe(df_bed)

# Intersect variants with the exon BED file to find overlapping exons
intersected = variants_bed.intersect(bed_file, wa=True, wb=True)

# Convert the intersection results back to a DataFrame
intersected_df = intersected.to_dataframe(names=['chrom', 'start', 'end', 'gene', 'exon_start', 'exon_end', 'exon_id'])

result = pd.DataFrame(columns=['chrom', 'position', 'gene', '5p_dist', '3p_dist'])

line_num = 0
for i, row in intersected_df.iterrows():
    if line_num % 1000 == 0:
        print(f'Processing line {line_num}')
    chrom = row['chrom']
    position = int(row['end'])
    gene_name = row['exon_id']
    strand = pcg_info[pcg_info['gene_name'] == gene_name]['strand'].values[0]
    before = position - int(row['exon_start'])
    after = int(row['exon_end']) - position
    # extract all rows in the bed file that have the same gene name
    gene_exons = exon_df[exon_df['gene'] == gene_name]
    for j, exon in gene_exons.iterrows():
        if int(exon['end']) < position:
            before += int(exon['end']) - int(exon['start'])
        if int(exon['start']) > position:
            after += int(exon['end']) - int(exon['start'])
    if strand == '+':
        result = result.append({'chrom': chrom, 'position': position, 'gene': gene_name, '5p_dist': before, '3p_dist': after}, ignore_index=True)
    elif strand == '-':
        result = result.append({'chrom': chrom, 'position': position, 'gene': gene_name, '5p_dist': after, '3p_dist': before}, ignore_index=True)
    else:
        print(f"ERROR: Strand unknown for {gene_name}, aborting!")
        quit()
    line_num += 1

result.to_csv('Variant_Coverage/priming_distances/'+sample+'.priming_distances.tsv', sep='\t', index=False)