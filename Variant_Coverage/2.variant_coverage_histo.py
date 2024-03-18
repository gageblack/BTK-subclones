import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import sys

directory = 'Variant_Coverage/genotypes/'

# Read in the TSV file to a pandas DataFrame
cell_df = pd.read_csv('meta/num_cells_per_sample.tsv', sep='\t', header=0, index_col=0)
sample = sys.argv[1] # Using 16150X1, 20583X3, 208020X6.
total_cells = cell_df.loc[sample, 'num_cells']
bin_size = 1  # bin size of 1%

# if there is a second argument, use it as the total number of cells
if len(sys.argv) > 2 and sys.argv[2] == "log":
    plot_type = sys.argv[2]
else:
    plot_type = "percentage"

print(f"Making plot for sample: {sample}")
print(f"Total cells: {total_cells}")

file_path = f'{directory}{sample}.genotype.Variant_Coverage.tsv'
df = pd.read_csv(file_path, sep='\t')

# Efficiently count the number of 0s using numpy sum on boolean array
num_0s = np.sum(df['Num_Cells'] == 0)

# Count the number of variants for each percentage bin
max_cells = int(total_cells * (bin_size / 100))  # number of cells per bin
bins = np.arange(0, total_cells + max_cells, max_cells)  # create bins
bin_labels = (bins / total_cells) * 100  # convert bins to percentages

# Use numpy broadcasting to create a 2D array. For each variant, determine if the number of cells is greater than or equal to the bin
cell_matrix = df['Num_Cells'].values[:, None] >= bins

# Efficiently sum over the columns to get the counts per bin
binned_counts = np.sum(cell_matrix, axis=0)

# Use numpy broadcasting to create a 2D array of comparisons
alt_cell_matrix = df['Num_Alt_Cells'].values[:, None] >= bins

# Efficiently sum over the columns to get the counts per bin
alt_binned_counts = np.sum(alt_cell_matrix, axis=0)
x_tick_positions = [0, 20, 40, 60, 80, 100]

if plot_type == "percentage":
    # Get the percentage of variants 
    binned_counts = binned_counts / len(df) * (100/bin_size)
    alt_binned_counts = alt_binned_counts / len(df) * (100/bin_size)
    x_tick_positions = [1, 20, 40, 60, 80, 100]
    y_tick_positions = [2,4,6,8,10,12,14,16,18,20]
    y_tick_labels = [f"{x}%" for x in y_tick_positions]

x_tick_labels = [f"{x}%" for x in x_tick_positions]

label_size = 19.5
tick_size = 17
print("Plotting bar chart...")
plt.figure(figsize=(5, 6))
if plot_type == "log":
    plt.bar(bin_labels, binned_counts, align='center', color='darkgray', edgecolor='black', linewidth=0.75, width=bin_size)
    plt.bar(bin_labels, alt_binned_counts, align='center', color='red', edgecolor='black', linewidth=0.75, width=bin_size, alpha=0.5)
    plt.ylabel(f'# variants covered by at least X% cells', fontsize=label_size)
    plt.yscale('log')
else:
    plt.bar(bin_labels[1:], binned_counts[1:], align='center', color='darkgray', edgecolor='black', linewidth=0.75, width=bin_size)
    plt.bar(bin_labels[1:], alt_binned_counts[1:], align='center', color='red', edgecolor='black', linewidth=0.75, width=bin_size, alpha=0.5)
    plt.ylabel(f'% variants covered by at least X% cells', fontsize=label_size)
    plt.xticks(ticks=x_tick_positions, labels=x_tick_labels)  # Set x-ticks and format them as percentages
    plt.yticks(ticks=y_tick_positions, labels=y_tick_labels)

plt.xticks(ticks=x_tick_positions, labels=x_tick_labels)  # Set x-ticks and format them as percentages
plt.xlabel(f'% of cells ({total_cells} total cells)', fontsize=label_size)
plt.xticks(fontsize=tick_size)
plt.yticks(fontsize=tick_size)
plt.xlim(-1, 101)
plt.tight_layout()

if plot_type == "log":
    plt.savefig('Variant_Coverage_Plots/'+sample+'.variant_coverage.log.png')
else:
    plt.savefig('Variant_Coverage_Plots/'+sample+'.variant_coverage.percentage.png')
