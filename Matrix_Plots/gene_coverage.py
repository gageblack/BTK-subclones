import pandas as pd
import sys

patient = sys.argv[1]
sample = sys.argv[2]

print("Using seurat sparse matrix to make gene coverage matrix for "+sample)

## Read in sparse matrix as pandas dataframe
sparse_matrix = pd.read_csv('gene_depth_matricies/'+sample+'.sparse_matrix.tsv', index_col=0, sep='\t')

sparse_matrix = sparse_matrix.groupby(sparse_matrix.index).sum()

# Read in the vcf file and make a list of genes and positions
gene_list = []
position_list = []
vcf_file = open(patient+'/'+patient+'.subclones.ordered.vcf', 'r')
for line in vcf_file:
    if line.startswith('#'):
        continue
    fields = line.strip().split('\t')
    position_list.append(fields[0]+':'+fields[1])
    info=fields[7].split(';')
    annotation = fields[7].split(';')[41]
    gene = annotation.split('|')[3]
    gene_list.append(gene)

# Create a DataFrame with the same columns as the sparse_matrix
new_df = pd.DataFrame(columns=sparse_matrix.columns)

for gene in gene_list:
    if gene in sparse_matrix.index:
        # If the gene is in the sparse_matrix, add the corresponding row
        new_df = new_df.append(sparse_matrix.loc[gene])
    else:
        # If the gene is not in the sparse_matrix, create a new row with all values set to 0
        new_df = new_df.append(pd.Series(0, index=new_df.columns, name=gene))

new_df.index = position_list

# Save the new_df as a tsv file
new_df.to_csv('gene_depth_matricies/'+sample+'.gene_coverage.tsv', sep='\t')