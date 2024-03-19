library(Matrix)

# Get command-line arguments
args <- commandArgs(trailingOnly = TRUE)

# Check if at least one argument (the input directory) is provided
if (length(args) < 1) {
  cat("Usage: Rscript make_matrix.R <input_directory>\n")
  quit(status = 1)
}

# Extract the input directory name
sample_id <- args[1]
print(paste("Making sparse matrix for:",sample_id))

print("Loading in data...")
input_directory <- paste("seurat/",sample_id,"/genes_seurat", sep = "")

# Read the sparse matrix from the Matrix Market format
sparse_matrix <- readMM(paste(input_directory,"matrix.mtx", sep="/"))

# Read the barcodes and features (gene names)
barcodes <- read.table(paste(input_directory,"barcodes.tsv", sep="/"), header = FALSE, stringsAsFactors = FALSE, sep = "\t")
features <- read.table(paste(input_directory,"genes.tsv", sep="/"), header = FALSE, stringsAsFactors = FALSE, sep = "\t")

# Split features$V1 by ":" and keep the second part to only include gene name.
features$V1 <- sapply(strsplit(features$V1, ":"), function(x) x[2])

# split barcodes$V1 by "-" and keep the first part to only include the sample id.
barcodes$V1 <- sapply(strsplit(barcodes$V1, "-"), function(x) x[1])

# Assign row and column names to the sparse matrix
rownames(sparse_matrix) <- features$V1
colnames(sparse_matrix) <- barcodes$V1

print("Writing the sparse matrix to file...")
write.table(as.matrix(sparse_matrix), file = paste("gene_depth_matricies/", sample_id,".sparse_matrix.tsv", sep = ""), sep = "\t", quote = FALSE)

print("Done!")