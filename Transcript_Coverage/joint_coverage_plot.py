print("Loading python libraries...")
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
import seaborn as sns
import numpy as np
import os

def load_counts(counts_file):
    # Initialize counts as a zero array with the desired size
    counts = np.zeros(101, dtype=int)

    with open(counts_file, 'r') as f:
        # Read lines and convert to integer coverage values
        coverage_values = np.array([int(float(line.strip())) for line in f])

    # Determine the maximum coverage to avoid unnecessary operations
    max_coverage = min(coverage_values.max(), counts.size)

    # Update counts for each coverage value efficiently
    for coverage in coverage_values:
        if coverage > 0:  # Check to avoid indexing error for coverage=0
            counts[:min(coverage, max_coverage)] += 1

    # Calculate the final counts as percentages
    num_lines = len(coverage_values)
    final_counts_np = counts / num_lines * 100

    return final_counts_np.tolist()

def load_final_counts(counts_file):
    final_counts = []
    with open(counts_file, 'r') as f:
        for line in f:
            final_counts.append(float(line.strip()))
    return final_counts

files = [f for f in os.listdir('Read_Coverages/') if f.endswith('.final_counts.txt')]

print('Plotting...')
# Plot settings
label_size = 21
tick_size = 19
colors = sns.color_palette()

fig, ax = plt.subplots(figsize=(7, 6))

for file in files:
    counts_file = 'Read_Coverages/' + file
    run = ''
    line_color = ''

    if file[0:5] == "20583":
        run = "Sequel II (N=4)"
        line_color = colors[0]
    elif file[0:5] == "20802":
        run = "Revio (N=8)"
        line_color = colors[1]
    elif file[0:5] == "21183":
        run = "Revio (N=8)"
        line_color = colors[1]
    elif file[0:5] == "16150":
        run = "Short-Read Sequencing (N=8)"
        line_color = "black"
    else:
        print("Error: Unknown File", file)
        exit(1)

    final_counts = load_final_counts(counts_file)
    x_values = list(range(0, len(final_counts)))  # Create x-values for positions
    ax.plot(x_values, final_counts, alpha=1, color=line_color)

ax.set_xlabel('% of Transcript', fontsize=label_size)
ax.set_ylabel('% Reads', fontsize=label_size)
ax.tick_params(axis='y', labelsize=tick_size)
ax.tick_params(axis='x', labelsize=tick_size)

handles = [
    mlines.Line2D([], [], color=colors[0],  linewidth=2, label='Sequel II (N=4)'),
    mlines.Line2D([], [], color=colors[1], linewidth=2, label='Revio (N=8)'),
    mlines.Line2D([], [], color="black", linewidth=2, label='Short-Read Seq (N=8)')
]

# Use the handles in the legend
ax.legend(handles=handles, loc='upper right', fontsize=tick_size-4)#, title=r'Sequencing Run', title_fontsize=tick_size-4)

plt.title('Read Coverage of Transcript', fontsize=label_size)
plt.tight_layout()
plt.savefig('joint_coverage_plot.png')

