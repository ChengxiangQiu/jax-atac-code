
### Demo-1: Performing dimension reduction on a subset of the kidney trajectory

print("Loading basic data analysis packages")
suppressMessages(library(Seurat))
suppressMessages(library(BPCells))

print("Loading packages for regular data analysis, e.g. dplyr")
suppressMessages(library(Matrix))
suppressMessages(library(dplyr))
suppressMessages(library(reshape2))

print("Loading packages for plotting, e.g. ggplot2")
suppressMessages(library(ggplot2))
suppressMessages(library(plotly))
suppressMessages(library(htmlwidgets))
suppressMessages(library(gridExtra))
suppressMessages(library(gplots))
suppressMessages(library(patchwork))
suppressMessages(library(viridis))
suppressMessages(library(RColorBrewer))

###########################################################
### Step-1: read the peak x cell matrix and cell meta table

### The demo data can be downloaded from:
### https://shendure-web.gs.washington.edu/content/members/cxqiu/public/backup/jax_atac/download/demo
### (1) demo_input.h5; peak x cell matrix
### (2) demo_cell_metadata.csv; cell meta table

setwd("/path/to/your/directory")
peak_mat = open_matrix_10x_hdf5("demo_input.h5")
use_meta = read.csv("demo_cell_metadata.csv", row.names = 1)

####################################################
### Step-2: performing dimension reduction using LSI

tf = peak_mat %>%
    multiply_cols(1 / Matrix::colSums(peak_mat))

scale_factor=100000
tf_idf = log1p(tf * scale_factor) %>% 
    multiply_rows(log1p(1 / Matrix::rowMeans(peak_mat)))

dims = 2:50
svd_atac = irlba::irlba(tf_idf, nv=max(dims))
pca_mat = multiply_cols(svd_atac$v, svd_atac$d)
rownames(pca_mat) = colnames(peak_mat)
colnames(pca_mat) = paste0("PC_", 1:ncol(pca_mat))
write.csv(pca_mat, "pca_coor.csv", quote=F)

pca_l2_mat = Seurat:::L2Norm(pca_mat)

set.seed(12341512)
umap_mat = uwot::umap(pca_l2_mat[,dims], n_components = 3, metric = "cosine", 
                      n_neighbors = 30, min_dist = 0.3)

rownames(umap_mat) = colnames(peak_mat)
colnames(umap_mat) = paste0("UMAP_", 1:ncol(umap_mat))
write.csv(umap_mat, "umap_coor.csv", quote=F)

######################################
### Step-3: performing cell clustering

seurat_obj = Seurat::CreateSeuratObject(peak_mat[1:100,])
seurat_obj[['umap']] = Seurat::CreateDimReducObject(embeddings=umap_mat, key='UMAP_', assay='RNA')
seurat_obj = seurat_obj %>% 
    Seurat::FindNeighbors(reduction="umap", nn.eps=0, dims=1:3, k.param=20) %>%
    Seurat::FindClusters(n.start=20, resolution=0.1) %>%
    Seurat::FindClusters(n.start=20, resolution=0.3) %>%
    Seurat::FindClusters(n.start=20, resolution=0.5) %>%
    Seurat::FindClusters(n.start=20, resolution=0.8)

write.csv(data.frame(seurat_obj[[]]), "clustering_result.csv", quote=F)

############################
### Step-4: ploting 3D umaps
### Of note, this is a subset of the full dataset (only 10%), so the resulting 3D UMAP may not closely match the version presented in the paper.

pd = read.csv("demo_cell_metadata.csv", row.names=1)
umap = read.csv("umap_coor.csv", row.names=1)
sum(rownames(umap) == rownames(pd))

pd$UMAP_1 = umap[,1]
pd$UMAP_2 = umap[,2]
pd$UMAP_3 = umap[,3]

plot_ly(pd, x=~UMAP_1, y=~UMAP_2, z=~UMAP_3, size = I(30), color = ~celltype_L3)










