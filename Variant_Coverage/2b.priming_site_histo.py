import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import sys

directory = 'Variant_Coverage/'

sample = sys.argv[1] # Using 16150X1, 20583X3, 208020X6.
print(f"Making plot for sample: {sample}")

if len(sys.argv) > 2 and (sys.argv[2] == "5" or sys.argv[2] == "3"):
    priming_site = sys.argv[2]
else:
    print("No valid priming site given (5 or 3). Setting to 3\'")
    priming_site = "3"

# Read in the Cell number TSV file to a pandas DataFrame
cell_df = pd.read_csv(f'{directory}meta/num_cells_per_sample.tsv', sep='\t', header=0, index_col=0)
total_cells = cell_df.loc[sample, 'num_cells']
print(f"Total cells: {total_cells}")

priming_distances = pd.read_csv(f'{directory}priming_distances/{sample}.priming_distances.tsv', sep='\t', header=0)
priming_distances['position'] = priming_distances['position'].astype(int)
priming_distances['chrom'] = priming_distances['chrom'].astype(str)
# make 6 different subsets of priming_distances based on the priming distance.
priming_distances['bin'] = pd.cut(priming_distances[f'{priming_site}p_dist'], bins=[0, 201, 401, 601, 801, 1001,float('inf')], labels=['<=200bp', '201bp-400bp', '401bp-600bp', '601bp-800bp', '801bp-1000bp', '>1000bp'])

# Read in the Variant Coverage TSV file to a pandas DataFrame
file_path = f'{directory}genotypes/{sample}.genotype.Variant_Coverage.tsv'
df = pd.read_csv(file_path, sep='\t')
df[['chrom', 'position', 'ref_allele', 'alt_allele']] = df['Variant'].str.split(':', expand=True)
df['position'] = df['position'].astype(int)
df['chrom'] = df['chrom'].astype(str)
df['chrom'] = df['chrom'].apply(lambda x: x if x.startswith('chr') else 'chr' + x)
df['cell_percentage'] = (df['Num_Cells'] / total_cells) * 100
df['alt_percentage'] = (df['Num_Alt_Cells'] / total_cells) * 100

bin_size = 1  # bin size of 1%
bin_edges = np.linspace(0, 100, 101) # Create a bin for 0-100% cells for the x axis

label_size = 25
tick_size = 24
# Prepare the figure for subplots
fig, axes = plt.subplots(1, 6, figsize=(32, 7), sharey=True)
# Add the legend to your plot

# Loop through each bin and generate a plot
for index, label in enumerate(priming_distances['bin'].cat.categories):
    # Merge the subset for the current bin with the coverage data
    bin_df = priming_distances[priming_distances['bin'] == label]
    merged_df = pd.merge(df, bin_df[['chrom', 'position', 'bin']], on=['chrom', 'position'])
    cell_matrix = merged_df['cell_percentage'].values[:, None] >= bin_edges # For each variant, determine if the cell_percentage is greater than or equal to the bin. Makes a 2D array of comparisons with True or False as values.
    binned_counts = np.sum(cell_matrix, axis=0) # Sums the counts per bin, representing the number of variants with at least X% of cells.
    binned_counts = binned_counts / len(merged_df) * (100/bin_size) # Get the percentage of variants for each bin.

    alt_cell_matrix = merged_df['alt_percentage'].values[:, None] >= bin_edges # For each variant, determine if the cell_percentage is greater than or equal to the bin. Makes a 2D array of comparisons with True or False as values.
    alt_binned_counts = np.sum(alt_cell_matrix, axis=0) # Sums the counts per bin, representing the number of variants with at least X% of cells.
    alt_binned_counts = alt_binned_counts / len(merged_df) * (100/bin_size) # Get the percentage of variants for each bin.


    # Plotting
    axes[index].bar(bin_edges[1:], binned_counts[1:], width=bin_size, color='darkgray', edgecolor='black', linewidth=0.5)
    axes[index].bar(bin_edges[1:], alt_binned_counts[1:], width=bin_size, color='red', alpha=0.5, edgecolor='black', linewidth=0.5)
    axes[index].set_title(f'{label}', fontsize=label_size)
    axes[index].set_xlabel('% of Cells', fontsize=label_size)
    if index == 0:  # Only add y-label to the first subplot to avoid repetition
        axes[index].set_ylabel('% Variants Covered', fontsize=label_size)
        axes[index].set_yticks([5,10,15,20,25,30,35])
        axes[index].set_yticklabels(['5%', '10%', '15%', '20%', '25%', '30%', '35%'], fontsize=tick_size)
    
    axes[index].set_xticks([1, 20, 40, 60, 80, 100], fontsize=tick_size)
    axes[index].set_xticklabels(['1%', '20%', '40%', '60%', '80%', '100%'], fontsize=tick_size)
    axes[index].set_ylim(0, 35)
    axes[index].set_xlim(0, 100)

plt.tight_layout()
plt.subplots_adjust(left=0.05) 
plt.show()