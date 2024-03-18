import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator, ScalarFormatter
import seaborn as sns
import pandas as pd

df = pd.read_csv('table1_stats.csv')

numerical_columns = ['Total HiFi Reads', 'Mean HiFi Read Length', 'Total Bases', 
                     'Segmented Reads', 'Average Length of S-reads', '# of cells', 
                     'Mean Reads per Cell', 'Median UMIs per cell']

for col in numerical_columns:
    if col in df.columns:
        df[col] = df[col].str.replace(',', '').astype(float)

# Setting the style for the plot
sns.set_context("paper")  # Larger font size for publication readiness
sns.set(style="whitegrid")

# Creating a figure and a set of subplots
fig, axes = plt.subplots(3, 1, figsize=(5.3, 4.5), sharex=True)

# Plotting each metric for individual samples
sns.barplot(x='Sample ID', y='Total HiFi Reads', hue='Sequencer', data=df, ax=axes[0],dodge=False)
sns.barplot(x='Sample ID', y='Segmented Reads', hue='Sequencer', data=df, ax=axes[1],dodge=False)
sns.barplot(x='Sample ID', y='Mean Reads per Cell', hue='Sequencer', data=df, ax=axes[2],dodge=False)

# Removing the legends for the first two plots and adjusting the legend for the last plot
axes[0].legend(title='Sequencer', bbox_to_anchor=(1.05, 1), loc='upper left')
axes[1].get_legend().remove()
axes[2].get_legend().remove()

# Removing the X-axis tick labels for all subplots
for ax in axes:
    ax.set_xticklabels([])

# Setting titles and labels
for i, ax in enumerate(axes):
    ax.set_title(['Total HiFi Reads', 'Segmented Reads', 'Mean Reads per Cell'][i])
    ax.set_ylabel('Count')
    ax.set_xlabel('')
    ax.yaxis.set_major_locator(MaxNLocator(nbins=4, integer=True))
    ax.yaxis.set_major_formatter(ScalarFormatter(useMathText=True))
    ax.ticklabel_format(axis='y', style='scientific', scilimits=(0,0))
axes[-1].set_xlabel('Sample')

plt.tight_layout()
plt.savefig('metrics_figure.png')
