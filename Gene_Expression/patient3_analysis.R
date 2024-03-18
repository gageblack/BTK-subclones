library(Seurat)
library(tidyverse)

patient = "patient3"
sample1 = "20583X1"
sample2 = "20583X2"

## Read in pre-treatment sample and make seurat object ##
matrix <- ReadMtx(paste(sample1,"/genes_seurat/matrix.mtx", sep = ""),
                  features = paste(sample1,"/genes_seurat/genes.tsv", sep = ""),
                  cells = paste(sample1,"/genes_seurat/barcodes.tsv", sep = ""),
                  feature.column = 2)
seurat_df_pre <- CreateSeuratObject(counts = matrix, project="BTK_pre", 
                                    min.cells = 1, min.features = 1)

## Read in pre-treatment sample and make seurat object ##
matrix <- ReadMtx(paste(sample2,"/genes_seurat/matrix.mtx", sep = ""),
                  features = paste(sample2,"/genes_seurat/genes.tsv", sep = ""),
                  cells = paste(sample2,"/genes_seurat/barcodes.tsv", sep = ""),
                  feature.column = 2)
seurat_df_post <- CreateSeuratObject(counts = matrix, project="BTK_post", 
                                     min.cells = 1, min.features = 1)

# Run SCTransform on both samples.
seurat_df_pre_sct <- SCTransform(seurat_df_pre, method = "glmGamPoi", vst.flavor = "v2")
seurat_df_post_sct <- SCTransform(seurat_df_post, method = "glmGamPoi", vst.flavor = "v2")

#################################################################################
# Process pre-treatment sample
seurat_df_pre_sct <- RunPCA(seurat_df_pre_sct, verbose = FALSE)
seurat_df_pre_sct <- RunUMAP(seurat_df_pre_sct, dims = 1:20, verbose = FALSE)
seurat_df_pre_sct <- FindNeighbors(seurat_df_pre_sct, dims = 1:20, verbose = FALSE)
seurat_df_pre_sct <- FindClusters(seurat_df_pre_sct, verbose = FALSE)
umap_plot_pre = DimPlot(seurat_df_pre_sct, label = TRUE) + NoLegend()
ggsave(paste("umaps/",patient,"/",patient,".pre.sct.umap.pdf",sep = ""), umap_plot_pre, width=8, height=6, useDingbats=FALSE)

## Read in subclone assignments, and write cell coordinates.
pre_sc_assignments = read.csv(paste("Genotype/",patient,"/",sample1,".assigned", sep=""), sep = "\t")
pre_sc_assignments$Barcode = paste(pre_sc_assignments$Barcode, "-1", sep="")

write.csv(seurat_df_pre_sct@reductions$umap@cell.embeddings, paste("umap_coordinate_files/",patient,".pre.umapcoordinate.csv",sep=""))
umap_pre_coordinate = read.csv(paste("umap_coordinate_files/",patient,".pre.umapcoordinate.csv",sep=""))
rownames(umap_pre_coordinate) = umap_pre_coordinate$X

sc0 = filter(pre_sc_assignments, ASIG == "C0")$Barcode
sc1 = filter(pre_sc_assignments, ASIG == "C1")$Barcode
sc2 = filter(pre_sc_assignments, ASIG == "C2")$Barcode
sc3 = filter(pre_sc_assignments, ASIG == "C3")$Barcode
sc4 = filter(pre_sc_assignments, ASIG == "C4")$Barcode
unassign = filter(pre_sc_assignments, ASIG == "UNASSIGN")$Barcode
normal = filter(pre_sc_assignments, ASIG == "normal")$Barcode

# Make umap plot coloring cells by subclone.
pdf(paste("umaps/",patient,"/",patient,".pre.sct.by_subclone.pdf",sep = ""))
plot(umap_pre_coordinate$UMAP_1,umap_pre_coordinate$UMAP_2,pch=16, cex=0.2,col="gray",xlab="UMAP1",ylab="UMAP2")#,xlim=c(-10,7),ylim=c(-8,9))  ## X and Y lim only if yo uwant to focus on somewhere.
points(umap_pre_coordinate[rownames(umap_pre_coordinate) %in% sc0,][,2],umap_pre_coordinate[rownames(umap_pre_coordinate) %in% sc0,][,3], pch=16, cex=0.7,col="green")
points(umap_pre_coordinate[rownames(umap_pre_coordinate) %in% sc1,][,2],umap_pre_coordinate[rownames(umap_pre_coordinate) %in% sc1,][,3], pch=16, cex=0.7,col="black")
points(umap_pre_coordinate[rownames(umap_pre_coordinate) %in% sc2,][,2],umap_pre_coordinate[rownames(umap_pre_coordinate) %in% sc2,][,3], pch=16, cex=0.7,col="red")
points(umap_pre_coordinate[rownames(umap_pre_coordinate) %in% sc3,][,2],umap_pre_coordinate[rownames(umap_pre_coordinate) %in% sc3,][,3], pch=16, cex=0.7,col="blue")
points(umap_pre_coordinate[rownames(umap_pre_coordinate) %in% sc4,][,2],umap_pre_coordinate[rownames(umap_pre_coordinate) %in% sc4,][,3], pch=16, cex=0.7,col="orange")
dev.off()

pre_diff_genes = FindMarkers(seurat_df_pre_sct,ident.1=c(1,2,5),ident.2=c(3,4,0)) ##Ident 1 is the green/orange clones, Ident 2 is the original blue clone
pre_diff_genes$gene = rownames(pre_diff_genes)
cll_drivers = read_tsv("../../meta/CLL_Drivers_2023.tsv")$GENE_NAME
pre_driver_diff = subset(pre_diff_genes, gene %in% cll_drivers)

#################################################################################
# Process post-treatment sample
seurat_df_post_sct <- RunPCA(seurat_df_post_sct, verbose = FALSE)
seurat_df_post_sct <- RunUMAP(seurat_df_post_sct, dims = 1:20, verbose = FALSE)
seurat_df_post_sct <- FindNeighbors(seurat_df_post_sct, dims = 1:20, verbose = FALSE)
seurat_df_post_sct <- FindClusters(seurat_df_post_sct, verbose = FALSE)
umap_plot_post = DimPlot(seurat_df_post_sct, label = TRUE) + NoLegend()
ggsave(paste("umaps/",patient,"/",patient,".post.sct.umap.pdf",sep = ""), umap_plot_post, width=8, height=6, useDingbats=FALSE)

# Pull in subclone assignment info
post_sc_assignments = read.csv(paste("Genotype/",patient,"/",sample2,".assigned", sep=""), sep = "\t")
post_sc_assignments$Barcode = paste(post_sc_assignments$Barcode, "-1", sep="")

write.csv(seurat_df_post_sct@reductions$umap@cell.embeddings, paste("umap_coordinate_files/",patient,".post.umapcoordinate.csv",sep=""))
umap_post_coordinate = read.csv(paste("umap_coordinate_files/",patient,".post.umapcoordinate.csv",sep=""))
rownames(umap_post_coordinate) = umap_post_coordinate$X

# Assign the cells of each subclone to a variable.
sc0 = filter(post_sc_assignments, ASIG == "C0")$Barcode
sc1 = filter(post_sc_assignments, ASIG == "C1")$Barcode
sc2 = filter(post_sc_assignments, ASIG == "C2")$Barcode
sc3 = filter(post_sc_assignments, ASIG == "C3")$Barcode
sc4 = filter(post_sc_assignments, ASIG == "C4")$Barcode
unassign = filter(post_sc_assignments, ASIG == "UNASSIGN")$Barcode
normal = filter(post_sc_assignments, ASIG == "normal")$Barcode

# Make umap plot coloring cells by subclone.
pdf(paste("umaps/",patient,"/",patient,".post.sct.by_subclone.pdf",sep = ""))
plot(umap_post_coordinate$UMAP_1,umap_post_coordinate$UMAP_2,pch=16, cex=0.2,col="gray",xlab="UMAP1",ylab="UMAP2")#,xlim=c(-10,7),ylim=c(-8,9))  ## X and Y lim only if yo uwant to focus on somewhere.
points(umap_post_coordinate[rownames(umap_post_coordinate) %in% sc0,][,2],umap_post_coordinate[rownames(umap_post_coordinate) %in% sc0,][,3], pch=16, cex=0.6,col="green")
points(umap_post_coordinate[rownames(umap_post_coordinate) %in% sc1,][,2],umap_post_coordinate[rownames(umap_post_coordinate) %in% sc1,][,3], pch=16, cex=0.6,col="black")
points(umap_post_coordinate[rownames(umap_post_coordinate) %in% sc2,][,2],umap_post_coordinate[rownames(umap_post_coordinate) %in% sc2,][,3], pch=16, cex=0.6,col="red")
points(umap_post_coordinate[rownames(umap_post_coordinate) %in% sc3,][,2],umap_post_coordinate[rownames(umap_post_coordinate) %in% sc3,][,3], pch=16, cex=0.6,col="blue")
points(umap_post_coordinate[rownames(umap_post_coordinate) %in% sc4,][,2],umap_post_coordinate[rownames(umap_post_coordinate) %in% sc4,][,3], pch=16, cex=0.6,col="orange")
dev.off()

post_diff_genes = FindMarkers(seurat_df_post_sct,ident.1=c(7),ident.2=NULL)
post_diff_genes$gene = rownames(post_diff_genes)
cll_drivers = read_tsv("../../meta/CLL_Drivers_2023.tsv")$GENE_NAME
post_driver_diff = subset(post_diff_genes, gene %in% cll_drivers)
