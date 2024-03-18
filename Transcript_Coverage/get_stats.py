import numpy as np
import os

files = [f for f in os.listdir('Read_Coverages/') if f.endswith('.0.9.txt')]

print("Sample\tMedian\t>80%")
for file in files:
    file_path = 'Read_Coverages/' + file
    sample_name = file.split('.')[0]

    values = np.loadtxt(file_path)
    median_value = np.median(values)
    values_gt_80 = values[values > 80]
    percentage_gt_80 = (len(values_gt_80) / len(values)) * 100
    print(f"{sample_name}\t{median_value}\t{percentage_gt_80}")