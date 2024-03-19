import sys
import pandas as pd

# To run: give a patient id and a sample id as command line input.

###########################################################
## Read the the genotype matrix, split each value by the colon, and divide the last value by the first value to make the new matrix.
###########################################################

assigned=False
assignments = pd.DataFrame()
patient=sys.argv[1]
sample=sys.argv[2]

genotype_file = open(patient+"/"+sample+".genotype", "r")
scaf_file = open("scaf_files/"+sample+".scaf.tsv", "w")
scref_file = open("scref_files/"+sample+".scaf.tsv", "w")
for line in genotype_file:
    fields = line.strip().split()
    if fields[0] == "variant":
        scaf_file.write(line)
        scref_file.write(line)
    else:
        new_af_row = []
        new_ref_row = []
        for i in range(1,len(fields)):
            if fields[i] == "0:0:0":
                new_af_row.append(0)
                new_ref_row.append(0)
            else:
                ## Get Allele Counts
                new_af_row.append(float(fields[i].split(":")[2]))
                new_ref_row.append(float(fields[i].split(":")[1]))
        scaf_file.write(fields[0]+"\t"+"\t".join([str(i) for i in new_af_row])+"\n")
        scref_file.write(fields[0]+"\t"+"\t".join([str(i) for i in new_ref_row])+"\n")


scaf_file.close()
scref_file.close()
genotype_file.close()
