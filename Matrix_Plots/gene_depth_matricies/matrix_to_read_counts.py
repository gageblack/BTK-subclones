import pandas as pd
import sys

sample = sys.argv[1]
matrix_dir="gene_depth_matricies/"

print("Loading sparse matrix for sample "+sample)
sparse_matrix = pd.read_csv(matrix_dir+sample+".sparse_matrix.tsv", delimiter='\t', index_col=0)

print("Calculating total reads per cell")
total_reads_per_cell = sparse_matrix.sum(axis=0)

print("Saving total reads per cell to file")
# Create a new dataframe with cell_barcode and total_num_reads columns
total_reads_df = pd.DataFrame({
    'Barcode': total_reads_per_cell.index,
    'total_num_reads': total_reads_per_cell.values
})

# save to file
total_reads_df.to_csv(matrix_dir+sample+".total_reads_per_cell.tsv", sep='\t', index=False)
