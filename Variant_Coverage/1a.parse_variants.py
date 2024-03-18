import sys
import statistics as st

infile = open(sys.argv[1], "r")
coverage_file = open(sys.argv[1]+".Variant_Coverage.tsv","w")

header = infile.readline().strip().split()
coverage_file.write("Variant\tNum_Cells\tNum_Alt_Cells\tTotal_Coverage\tRef_Count\tAlt_Count\tMedian_Coverage\tMedian_Ref\tMedian_Alt\n")
for line in infile:
	cell_count = 0
	alt_cells = 0
	total_coverages = []
	ref_counts = []
	alt_counts = []
	afs = set()

	fields = line.strip().split()
	for i in range(1,len(fields)):
		if fields[i] != "0:0:0":
			info = fields[i].split(":")
			total = int(info[0])
			ref = int(info[1])
			alt = int(info[2])
			
			afs.add(alt/total)
			total_coverages.append(total)
			ref_counts.append(ref)
			alt_counts.append(alt)
			cell_count = cell_count+1
			if alt > 0:
				alt_cells = alt_cells+1

	if len(total_coverages) == 0:
		median_total_coverages = 0
	else:
		median_total_coverages = st.median(total_coverages)
	if len(ref_counts) == 0:
		median_ref_counts = 0
	else:
		median_ref_counts = st.median(ref_counts)
	if len(alt_counts) == 0:
		median_alt_counts = 0
	else:
		median_alt_counts = st.median(alt_counts)

	coverage_file.write(fields[0]+"\t"+str(cell_count)+"\t"+str(alt_cells)+"\t"+
		     str(sum(total_coverages))+"\t"+str(sum(ref_counts))+"\t"+str(sum(alt_counts))+"\t"+
			 str(median_total_coverages)+"\t"+str(median_ref_counts)+"\t"+str(median_alt_counts)+
			 "\n")

infile.close()
coverage_file.close()


