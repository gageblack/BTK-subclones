import sys
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import math

def filter_and_sort_assignment(assignment, subclone_order, qual_threshold, sort_method = "alt_count", filter_unassigned=False):
    """Filter the assignment DataFrame as specified by the quality and filter values.
    Also, sort the DataFrame by the subclone order and then by quality score.
    """

    # Make "." 0 in quality column and convert to float. Then filter by QUAL > qual_threshold.
    assignment = assignment.copy()
    assignment['QUAL'] = assignment['QUAL'].replace(".", 0)
    assignment['QUAL'] = assignment['QUAL'].astype(float)
   
    # If the QUAL <= qual_threshold, then change the ASIG value to UNASSIGN.
    assignment.loc[assignment['QUAL'] <= qual_threshold, 'ASIG'] = 'UNASSIGN'

    # Filter out the UNASSIGN cells if filter_unassigned is True
    if filter_unassigned:
        assignment = assignment[(assignment['ASIG'] != 'UNASSIGN')]

    if sort_method == 'total_num_reads':
        read_counts = pd.read_csv(working_dir+"gene_depth_matricies/"+sample+".total_reads_per_cell.tsv", delimiter='\t', index_col=0)
        new_assignment = pd.merge(assignment, read_counts, left_on='Barcode', right_index=True)
        new_assignment['total_num_reads'] = new_assignment['total_num_reads'].fillna(0)

    if sort_method == 'alt_count':
        alt_counts = pd.read_csv(working_dir+"scaf_files/"+sample+".scaf.tsv", delimiter='\t')
        if MAKE_DRIVER_PLOT:
            cll_vars = get_driver_variants(patient)
            alt_counts = alt_counts[alt_counts["variant"].str.split(":").str[:2].str.join(":").isin(cll_vars)]
  
        # For each column, count how many rows are > 0
        alt_counts = alt_counts.astype(bool).sum(axis=0)
        alt_counts = alt_counts.to_frame()
        alt_counts.columns = ["alt_count"]
        alt_counts.index.name = "Barcode"
        new_assignment = pd.merge(assignment, alt_counts, left_on='Barcode', right_index=True)
    
    if sort_method == 'num_alts':
        alt_counts = pd.read_csv(working_dir+"scaf_files/"+sample+".scaf.tsv", delimiter='\t')
        alt_counts = alt_counts.set_index("variant")

        final_barcode_order = []

        # For each subclone in subclone_order, make a new dataframe that is a subset of assignment
        for clone in subclone_order:
            clone_assignment = assignment[assignment["ASIG"] == clone]
            # subset alt_counts to include the barcodes in clone_assignment
            clone_alt_counts = alt_counts[clone_assignment["Barcode"]]
            variant_coverage = pd.DataFrame({
                "counts":clone_alt_counts.astype(bool).sum(axis=1),
                "variant":clone_alt_counts.index})

            sorted_variant_coverage = variant_coverage.sort_values(by="counts", ascending=False)
            column_order = clone_alt_counts.loc[sorted_variant_coverage["variant"].tolist()].sort_values(by=sorted_variant_coverage["variant"].tolist(),ascending=False, axis=1).columns

            final_barcode_order = final_barcode_order + column_order.tolist()
        
        # make a new dataframe with a "Barcode" column and a "postition" column that is an index of 0 to len(column_order)
        column_order_df = pd.DataFrame({"Barcode":final_barcode_order, "num_alts":range(0,len(final_barcode_order))})
        column_order_df = column_order_df.set_index("Barcode") 

        new_assignment = pd.merge(assignment, column_order_df, left_on='Barcode', right_index=True)
        new_assignment = new_assignment.sort_values(by=sort_method, ascending=True)

        return new_assignment

    # Sort the DataFrame by the subclone order and then by quality value.
    new_assignment['ASIG'] = pd.Categorical(new_assignment['ASIG'], categories=subclone_order, ordered=True)
    new_assignment = new_assignment.sort_values(by=['ASIG', sort_method], ascending=[True, False])

    return new_assignment

def process_genotype_data(assignment, alt_data, ref_data):
    """Combine the Alt and Ref count DFs into one. 
    Alt values will be positive, Ref values will be negative.
    only the filtered cells.
    """
    # Subset the data matrix to include only the filtered cells
    selected_columns = ["variant"] + list(assignment["Barcode"])
    filtered_alt_data = alt_data[selected_columns]
    filtered_ref_data = ref_data[selected_columns]

    ## Make the first column the index and then convert all values to numeric
    filtered_alt_data = filtered_alt_data.set_index("variant")
    filtered_alt_data = filtered_alt_data.apply(pd.to_numeric, errors='coerce')
    filtered_ref_data = filtered_ref_data.set_index("variant")
    filtered_ref_data = filtered_ref_data.apply(pd.to_numeric, errors='coerce')

    # Create a copy of filtered_alt_data
    merged_data = filtered_alt_data.copy()

    # Apply the condition to set values to -filtered_ref_data
    condition = (filtered_alt_data == 0) & (filtered_ref_data > 0)
    merged_data[condition] = filtered_ref_data[condition].applymap(lambda x: -x)

    return merged_data

def get_x_coords(assignment, subclone_order):
    x_coords = []
    subclones_assigned = []
    for i in range(len(subclone_order)):
        ## Get the position in the index of the first cell that has the subclone label
        for j in range(len(assignment)):
            if assignment.iloc[j]["ASIG"] == subclone_order[i]:
                index = j
                if len(subclones_assigned) == 0:
                    subclones_assigned.append(subclone_order[i])
                elif index != x_coords[-1]:
                    subclones_assigned.append(subclone_order[i])
                
                x_coords.append(index+0.5) ##Add 0.5 to center the line between the subclone clusters
                break
    x_coords.append(len(assignment)+0.5)
    return list(set(x_coords)), subclones_assigned

def get_y_coords(patient, geno_df, subclone_order):
    parsed_vcf_file = open("parsed_vcfs/"+patient+".subclones.parsed.tsv", "r")
    ## Make a dictionary that holds the cluster assignments of each variant
    var_clus = dict()
    variant_cluster_id = []
    parsed_vcf_file.readline()
    for line in parsed_vcf_file:
        fields = line.strip().split()
        genotype_index_list = [x.split(":")[0] + ":" + x.split(":")[1] for x in geno_df.index]

        if fields[0]+":"+fields[1] in genotype_index_list: ## Only add the variant to the dictionary if it is in the genotype matrix. Some variants are not in the matrix because they were removed for the scBayes step.
            var_clus[fields[0]+":"+fields[1]] = "C"+fields[5]
            variant_cluster_id.append(fields[5])

    y_coords = []
    subclones_assigned = []
    for i in range(len(subclone_order)):
        ## Get the position in the index of the first cell that has the subclone label
        for j in range(len(genotype_index_list)):
            if var_clus[genotype_index_list[j]] == subclone_order[i]:
                index = j
                if len(subclones_assigned) == 0:
                    subclones_assigned.append(subclone_order[i])
                elif index != y_coords[-1]:
                    subclones_assigned.append(subclone_order[i])
                
                y_coords.append(index+0.5) ##Add 0.5 to center the line between the subclone clusters
                break
    y_coords.append(len(genotype_index_list)+0.5)
    return list(set(y_coords)), subclones_assigned

def normalize(ref_values, max_val):
    ref_values = ref_values.copy()
    ref_values[ref_values > max_val] = max_val
    array_max = np.max(ref_values)
    ref_values = ref_values.astype(float) / array_max
    return ref_values

def get_driver_variants(patient):
    ## Make a list of the cll-relevant variants found in the parsed CLL VCF file
    parsed_cll_vcf_file = open("parsed_vcfs/"+patient+".subclones.cll.parsed.tsv", "r")
    cll_vars = []
    parsed_cll_vcf_file.readline()
    for line in parsed_cll_vcf_file:
        fields = line.strip().split()
        cll_vars.append(fields[0]+":"+fields[1])
    return cll_vars

def make_plot(data, cells, assignment_data, patient, subclone_order, subclone_colors, COVERAGE_PLOT=False, mainlabel="", filename="", scale_size = False, markersize=50, plot_width=20, plot_height=6):
    # Get the variant index of the driver variants
    new_df = pd.DataFrame()
    if COVERAGE_PLOT:
        new_df["variant"] = data.index
    else:
        new_df["variant"] = data.index.str.split(":").str[:2].str.join(":")
    new_df["variant_index"] = list(range(1, len(new_df.index) + 1))
    new_df.set_index("variant", inplace=True)
    
    ## Set the plot size based on the DF size.
    if scale_size:
        plot_height = math.ceil(math.sqrt(data.shape[0])) # If manual: 6
        plot_width = math.ceil(math.sqrt(data.shape[1])/1.5) # If manual: 20
        markersize = 1000/(plot_width+plot_height) # If manual: 50
        if plot_height < 3:
            plot_height = 3
        if markersize > 100:
            markersize = 100
    plot_height = 2
    plot_width = 6
    markersize = 10

    label_size = 15
    tick_size = 12
    
    matrix = data.iloc[:, data.columns.isin(cells)]
    matrix = pd.concat([data.index.to_frame(), matrix], axis=1)
    
    # Change row and column names to just a number, 1 through the number of variants and cell barcodes.
    matrix.columns = ["variant"] + list(range(1, len(matrix.columns)))
    matrix["variant"] = list(range(1, len(matrix.index) + 1))
    
    # Transform the matrix from wide format to long format.
    df_heatmap = pd.melt(matrix, id_vars="variant")
    
    # Make all the values in the matrix numeric
    df_heatmap["variant"] = pd.to_numeric(df_heatmap["variant"])
    df_heatmap["variable"] = pd.to_numeric(df_heatmap["variable"])
    df_heatmap["value"] = pd.to_numeric(df_heatmap["value"])
    
    # Create a Figure and Axes
    fig, ax = plt.subplots(figsize=(plot_width, plot_height))

    # Define and create a custom marker path as a rectangle (adjust width and height as needed)
    '''width = 0.05  # Width of the rectangle
    height = 1.5  # Height of the rectangle (will be centered around 0)
    marker_path = [(0, -height), (0, height), (width, height), (width, -height), (0, -height)]
    custom_marker = mpl.markers.MarkerStyle(marker_path)'''
    if COVERAGE_PLOT:
        # Filter rows based on value
        gene_AF = df_heatmap[df_heatmap["value"] > 0]
        
        # Create custom colormaps
        gene_colormap = mpl.cm.Blues(np.linspace(0,1,20))
        gene_colormap = mpl.colors.ListedColormap(gene_colormap[5:17,:-1])
        
        # Normalize the counts so it is scaled from 0 to 1
        scaled_gene_values = normalize(gene_AF["value"], 6)

        # Create the scatter plots on the Axes
        ax.scatter(gene_AF["variable"], gene_AF["variant"], marker="s", s=markersize, c=scaled_gene_values, cmap=gene_colormap, label="Gene AF < 0")

    else:
        print("New DataFrame shape: ", data.shape)
        if MAKE_DRIVER_PLOT:
            print("Reverse order of Indicies of the variants:")
            ## Print the new_df with the last column and first column last
            print(new_df.iloc[::-1])

        # Filter rows based on value
        no_coverage = df_heatmap[df_heatmap["value"] == 0]
        alt_AF = df_heatmap[df_heatmap["value"] > 0]
        ref_AF = df_heatmap[df_heatmap["value"] < 0]

        # Create custom colormaps
        ref_colormap = mpl.cm.Greens(np.linspace(0,1,20))
        ref_colormap = mpl.colors.ListedColormap(ref_colormap[5:17,:-1])
        alt_colormap = mpl.cm.Reds(np.linspace(0,1,20))
        alt_colormap = mpl.colors.ListedColormap(alt_colormap[5:17,:-1])
        
        # Normalize the counts so it is scaled from 0 to 1
        scaled_ref_values = normalize(-ref_AF["value"], 4)
        scaled_alt_values = normalize(alt_AF["value"], 4)

        # Create the scatter plots on the Axes
        ax.scatter(ref_AF["variable"], ref_AF["variant"], s=markersize, c=scaled_ref_values, cmap=ref_colormap, label="Ref AF < 0", marker="s")
        ax.scatter(alt_AF["variable"], alt_AF["variant"], s=markersize, c=scaled_alt_values, cmap=alt_colormap, label="Alt AF > 0", marker="s")

    # Create the verticle lines between cell subclone clusters
    x_line_positions, subclones_with_cells = get_x_coords(assignment_data, subclone_order)
    if len(subclones_with_cells) == 0:
        subclones_with_cells.append("No cells remaining")
    #print("Subclones that have cells assigned: ", subclones_with_cells)
    y_coords = np.arange(0.5, len(data.index) + 1.5) 
    for x in x_line_positions:
        x_coords = np.full(len(y_coords), x)
        ax.plot(x_coords, y_coords, color="black", linewidth=1.5)

    # Create the horizontal lines between variant subclone clusters
    y_line_positions, subclones_with_variants = get_y_coords(patient, data, subclone_order)
    #print("Subclones that have variants assigned: ", subclones_with_variants)
    x_coords = np.arange(0.5, len(matrix.columns) + 0.5)
    for y in y_line_positions:
        y_coords = np.full(len(x_coords), y)
        ax.plot(x_coords, y_coords, color="black", linewidth=1.5)

    x_line_positions = sorted(x_line_positions)
    y_line_positions = sorted(y_line_positions)
    
    # Create the vertical color bars next to the Y-axis
    start = y_line_positions[0]
    clu_index = 0
    for i in range(1, len(y_line_positions)):
        end = y_line_positions[i]
        y_coords = np.arange(start+.5, end+1)
        x_coords = np.full(len(y_coords), -2)#-x_line_positions[-1]/50)
        clu = subclones_with_variants[clu_index]
        ax.plot(x_coords, y_coords, color=subclone_colors[clu], label=clu, linewidth=4)
        start = end
        clu_index += 1
    
    # Create the horizontal bar below the X-axis
    x_line_positions = sorted(x_line_positions)
    start = x_line_positions[0]
    clu_index = 0
    for i in range(1, len(x_line_positions)):
        end = x_line_positions[i]
        x_coords = np.arange(start+2, end+1)
        y_coords = np.full(len(x_coords), 0)
        clu = subclones_with_cells[clu_index]
        ax.plot(x_coords, y_coords, color=subclone_colors[clu], label=clu, linewidth=4)
        start = end
        clu_index += 1

    ## Make a new dictionary for legend that only includes the subclones that have cells or variants assigned
    subclone_legend_colors = {k: subclone_colors[k] for k in subclone_colors if k in subclones_with_cells or k in subclones_with_variants}
    # Create the legend
    legend_handles = [plt.Rectangle((0, 0), 1, 1, color=color, label=label)
                    for label, color in subclone_legend_colors.items()]
    
    ax.set_xlabel('Cells', fontsize=label_size)
    ax.set_ylabel('Variants', fontsize=label_size)
    plt.xticks(fontsize = tick_size)
    plt.yticks([])
    plt.xlim(-5, len(matrix.columns))
    plt.title("Cell genotypes at subclone-defining variants", fontsize=label_size)
    plt.tight_layout()

    # Show the plot or save to a file
    if filename:
        fig.savefig(filename, dpi=500)
    else:
        plt.show()
    return

############################################################################################################
## Parameters to change
## You will need to run the make_scaf.py if you havent already.

patient = sys.argv[1]
sample = sys.argv[2]
output_opt = int(sys.argv[3]) # 1 = filter unassigned cells, 0 = do not filter unassigned cells

print("Making plots for patient:",patient," sample:",sample," output option:",output_opt)
working_dir = "Matrix_Plots/"

## Settings ##
min_assignment_quality = 0 # Minimum scBayes QUAL score. Any cell below is labeled unassigned
SORT_METHOD = "alt_count" # "alt_count", "total_num_reads", or "QUAL"
FILTER_UNASSIGNED = False # Remove any cell labeled UNASSIGN
MAKE_DRIVER_PLOT=False
NEEDS_VARIANT = True # For a cell to be included, it needs at least 1 alt read
ONLY_PERCENTAGE = False # Only include the top X% of cells in each subclone cluster
TOP_PERCENTAGE = 1.0 # The top X% of cells to include in each subclone, default is all cells
HIGH_COV_VARS = False # Only include variants with high coverage
MAKE_COVERAGE_PLOT = False # Make a coverage plot

mainlabel = patient+": "+sample+" genotype. QUAL > "+str(min_assignment_quality)

## The order you want the subclones to appear in the heatmap (Should match VCF order)##
## Subclone colors are also being specified here so that the order of colors is the same across patients. 
## Keep the order of the colors the same, make the "CX" order match subclone_order.
if patient == "patient1":
    subclone_order = ["C0", "C3", "C2", "C1", "C4", "C5","UNASSIGN"] # 111330019
    subclone_colors = {
        "C0": "blue",
        "C3": "green",
        "C2": "red",
        "C1": "black",
        "C4": "orange",
        "C5": "cyan",
        "UNASSIGN": "lightgray",
    }
elif patient == "patient2":
    subclone_order = ["C2", "C0", "C1", "C3","UNASSIGN"] # 111330031
    subclone_colors = {
        "C2": "blue",
        "C0": "green",
        "C1": "red",
        "C3": "black",
        "UNASSIGN": "lightgray",
    }
elif patient == "patient3":
    subclone_order = ["C3", "C0", "C2", "C1", "C4","UNASSIGN"] # 251882
    subclone_colors = {
        "C3": "blue",
        "C0": "green",
        "C2": "red",
        "C1": "black",
        "C4": "orange",
        "UNASSIGN": "lightgray",
    }
elif patient == "patient4":
    subclone_order = ["C2", "C0", "C4", "C1", "C3","UNASSIGN"] # 111330026
    subclone_colors = {
        "C0": "red",
        "C1": "green",
        "C2": "blue",
        "C3": "black",
        "C4": "orange",
        "UNASSIGN": "lightgray",
    }
elif patient == "patient5":
    subclone_order = ["C2", "C4", "C1", "C3", "C0", "UNASSIGN"] 
    subclone_colors = {
        "C0": "red",
        "C1": "green",
        "C2": "blue",
        "C3": "black",
        "C4": "orange",
        "UNASSIGN": "lightgray",
    }
elif patient == "patient6":
    subclone_order = ["C1", "C3", "C0", "C4", "C2", "UNASSIGN"] 
    subclone_colors = {
        "C0": "red",
        "C1": "green",
        "C2": "orange",
        "C3": "blue",
        "C4": "black",
        "UNASSIGN": "lightgray",
    }
else:
    print("Figure out what is wrong with your patient subclone order!")
    quit()
 
## Set the output options ##
if output_opt == 1: ## All cells, sorted by total reads per cell
    SORT_METHOD = "total_num_reads"
    FILTER_UNASSIGNED = True
    plot_file_name = sample+".by_total_reads.png"
    min_assignment_quality = 2
    TOP_PERCENTAGE = 0.9
    ONLY_PERCENTAGE = True
    HIGH_COV_VARS = True
elif output_opt == 2: ## All cells, sorted by alt count
    SORT_METHOD = "alt_count"
    plot_file_name = sample+".by_alt_count.png"
    FILTER_UNASSIGNED = True
elif output_opt == 3: ## Take the top 50% of cells in each subclone cluster, sorted by total reads per cell
    SORT_METHOD = "total_num_reads"
    ONLY_PERCENTAGE = True
    NEEDS_VARIANT = True
    TOP_PERCENTAGE = 0.50
    plot_file_name = sample+".by_total_reads.top50p.png"
elif output_opt == 4: ## Take the top 50% of cells in each subclone cluster, sorted by alt count
    SORT_METHOD = "alt_count"
    ONLY_PERCENTAGE = True
    NEEDS_VARIANT = True
    TOP_PERCENTAGE = 0.50
    plot_file_name = sample+".by_alt_count.top50p.png"
elif output_opt == 5: ## Make a plot of the driver variants
    MAKE_DRIVER_PLOT=True
    NEEDS_VARIANT = True
    min_assignment_quality = 0
    plot_file_name = sample+".drivers.png"
elif output_opt == 6: ## All cells, sorted by the number of alt alleles in the cell
    SORT_METHOD = "num_alts"
    plot_file_name = sample+".by_num_alt_order.png"
else:
    print("Specify an output_opt!")
    quit()

print("Loading in data frames...")
alt_data = pd.read_csv(working_dir+"scaf_files/"+sample+".scaf.tsv", delimiter='\t')
ref_data = pd.read_csv(working_dir+"scref_files/"+sample+".scaf.tsv", delimiter='\t')
assignment = pd.read_csv("../Genotype"+patient+"/"+sample+".assigned", delimiter='\t')
coverage_data = pd.read_csv(working_dir+"gene_depth_matricies/"+sample+".gene_coverage.tsv", delimiter='\t', index_col=0)

# Filter and manipulate the assignment DataFrame. Sort by ASIG, then by QUAL
print("Filtering and sorting assignments...")
filtered_assignment = filter_and_sort_assignment(assignment, subclone_order, min_assignment_quality, SORT_METHOD, FILTER_UNASSIGNED)

# Using the filtered assignment DataFrame, filter and manipulate the genotype data 
# to include only the filtered cells, and where ref reads are set to -1.
print("Processing genotype data...")
merged_data = process_genotype_data(filtered_assignment, alt_data, ref_data)

if MAKE_DRIVER_PLOT:
    cll_vars = get_driver_variants(patient)
    merged_data = merged_data[merged_data.index.str.split(":").str[:2].str.join(":").isin(cll_vars)]
    coverage_data = coverage_data[coverage_data.index.isin(cll_vars)]

if HIGH_COV_VARS:
    # Keep only variants with high coverage
    merged_data = merged_data.loc[(merged_data == 0).sum(axis=1) < len(merged_data.columns)*0.95]
    # Remove rows from coverage_data that are not in merged_data
    coverage_data = coverage_data[coverage_data.index.isin(merged_data.index)]

# remove columns that have all 0s from merged_data
merged_data = merged_data.loc[:, (merged_data != 0).any(axis=0)]

if NEEDS_VARIANT:
    # remove columns that do not have any postive values from merged_data
    merged_data = merged_data.loc[:, (merged_data > 0).any(axis=0)]

if ONLY_PERCENTAGE:
    # Remove rows from filtered_assignment if the barcode value is not a column name in merged_data
    updated_assignment = filtered_assignment[filtered_assignment['Barcode'].isin(merged_data.columns)]

    # Filter new_assignment to only include the first X% of cells in each subclone cluster
    updated_assignment = updated_assignment.groupby('ASIG').apply(lambda x: x.head(int(len(x) * TOP_PERCENTAGE))).reset_index(drop=True)

    # remove columns from merged_data if the barcode value is not in updated_assignment
    merged_data = merged_data[updated_assignment['Barcode']]

# Remove rows from filtered_assignment if the barcode value is not a column name in merged_data
filtered_assignment = filtered_assignment[filtered_assignment['Barcode'].isin(merged_data.columns)]

print("Making variant plot...")
make_plot(merged_data, merged_data.columns, filtered_assignment, patient, subclone_order, 
        subclone_colors, mainlabel=mainlabel, 
        filename=working_dir+"plots/"+plot_file_name)

if MAKE_COVERAGE_PLOT:
    ## Make the coverage plot ##
    selected_columns = list(filtered_assignment["Barcode"])
    filtered_coverage_data = coverage_data[selected_columns]
    print("Making gene coverage plot...")
    make_plot(filtered_coverage_data, filtered_coverage_data.columns, filtered_assignment, patient, 
                subclone_order, subclone_colors, COVERAGE_PLOT=True, mainlabel=mainlabel, 
                filename=working_dir+"plots/"+patient+"/"+plot_file_name[:-4]+".coverage.png")