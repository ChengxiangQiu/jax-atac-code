
########################################
### Clustering and cell type annotations

### Support data can be downloaded from:
### https://shendure-web.gs.washington.edu/content/members/cxqiu/public/backup/jax_atac/download/

### The final dataset (n = 3,937,903 cells) with cell lineage, cell class, and cell type annotations is available at:
### https://shendure-web.gs.washington.edu/content/members/cxqiu/public/backup/jax_atac/download/atac.pd.rds

### Please contact Chengxiang (CX) Qiu for any questions!
### cxqiu@uw.edu or chengxiang.qiu@dartmouth.edu

############################################################
### Step-1: removing doublets and creating the final dataset

work_path = ""
web_path = "https://shendure-web.gs.washington.edu/content/members/cxqiu/public/backup/jax_atac/download"
source("help_code/utils.R")

peak_mat = open_matrix_dir(paste0(work_path, "/atac_all"))
### 376,574 peaks x 4,344,905 cells

use_meta = readRDS(paste0(work_path, "/atac_all.pd.rds"))

use_meta_exclude = use_meta$cluster_PCA_1 %in% paste0("cluster_", c(6,8)) | 
    use_meta$amulet_pval < 0.05
print(table(use_meta_exclude))
### 3937903 retained vs. 407002 removed; 9.4% doublets

use_meta = use_meta[!use_meta_exclude,]
peak_mat = peak_mat[,!use_meta_exclude]
print(dim(peak_mat))

### exclude peaks overlapping ENCODE blacklist regions or falling on sex chromosomes
exclude = as.vector(read.table(paste0(web_path, "/exclude_peak_sample_matrix.list"))$V1)
exclude = gsub("-","_",exclude)

### exclude peaks which express less in 0.1% of cells
peak_sum = Matrix::rowSums(peak_mat)
peak_include = names(peak_sum)[peak_sum >= ncol(peak_mat) * 0.001]
peak_include = peak_include[!peak_include %in% exclude]

peak_mat = peak_mat[rownames(peak_mat) %in% peak_include,]
print(dim(peak_mat))
### 295143 peaks X  3937903 cells

### First, performing LSI based on AH's approach
tf = peak_mat %>%
    multiply_cols(1 / Matrix::colSums(peak_mat))

scale_factor=100000
tf_idf = log1p(tf * scale_factor) %>% 
    multiply_rows(log1p(1 / Matrix::rowMeans(peak_mat)))

dims = 2:50
svd_atac = irlba::irlba(tf_idf, nv=max(dims))
pca_mat = multiply_cols(svd_atac$v, svd_atac$d)

rownames(pca_mat) = colnames(peak_mat)
saveRDS(pca_mat, paste0(work_path, "/atac.pca.rds"))

pca_l2_mat = Seurat:::L2Norm(pca_mat)
saveRDS(pca_l2_mat, paste0(work_path, "/atac.pca_l2.rds"))

write.csv(round(pca_l2_mat[,dims], 7), paste0(work_path, "/atac.pca_l2.csv"))

set.seed(12341512)
umap_mat = uwot::umap(pca_l2_mat[,dims], n_components = 3, metric = "cosine", 
                      n_neighbors = 50, min_dist = 0.1)
umap_mat[1:4,]

rownames(umap_mat) = colnames(peak_mat)
saveRDS(umap_mat, paste0(work_path, "/atac.umap.rds"))

library(future)
options(future.globals.maxSize = 50000 * 1024^2)

seurat_obj = Seurat::CreateSeuratObject(peak_mat[1:100,])
seurat_obj[['umap']] = Seurat::CreateDimReducObject(embeddings=umap_mat, key='UMAP_', assay='RNA')
seurat_obj = seurat_obj %>% 
    Seurat::FindNeighbors(reduction="umap", nn.eps=0, dims=1:3, k.param=20) %>%
    Seurat::FindClusters(n.start=20, resolution=0.1)

saveRDS(data.frame(seurat_obj[[]]), paste0(work_path, "/atac.clustering_UMAP.rds"))

use_meta = readRDS(paste0(work_path, "/atac_all.pd.rds"))
umap_mat = readRDS(paste0(work_path, "/atac.umap.rds"))
cluster_res = readRDS(paste0(work_path, "/atac.clustering_UMAP.rds"))
print(sum(rownames(umap_mat) == rownames(cluster_res)))
rownames(umap_mat) = rownames(cluster_res) = unlist(lapply(rownames(umap_mat), function(x) strsplit(x,"[-]")[[1]][1])) 

rownames(use_meta) = as.vector(use_meta$cell_id)
use_meta = use_meta[rownames(umap_mat),]
print(sum(rownames(umap_mat) == rownames(use_meta)))
rownames(use_meta) = NULL

use_meta = use_meta[,c("cell_id", "sample", "day", "sex",  "reads", "reads_peak", "reads_tss", 
                       "pct_reads_in_peaks", "pct_reads_in_tss", "blacklist_ratio",
                       "amulet_pval", "doublet_likelihood")]
use_meta$UMAP_1 = as.vector(umap_mat[,1])
use_meta$UMAP_2 = as.vector(umap_mat[,2])
use_meta$UMAP_3 = as.vector(umap_mat[,3])

use_meta$global_cluster = paste0("cluster_", as.numeric(cluster_res$RNA_snn_res.0.1))
use_meta$global_cluster = factor(use_meta$global_cluster, levels = paste0("cluster_", names(table(as.numeric(cluster_res$RNA_snn_res.0.1)))))

atac_scglue = readRDS(paste0(work_path, "/atac.rds"))
use_meta$major_trajectory_RNA = use_meta$celltype_update_RNA = NULL
use_meta_x = use_meta %>% left_join(atac_scglue[,c("cell_id", "major_trajectory", "celltype_update")])
use_meta$major_trajectory_RNA = as.vector(use_meta_x$major_trajectory)
use_meta$celltype_update_RNA = as.vector(use_meta_x$celltype_update)

pca_cluster = read.csv(paste0(work_path, "/atac.clustering_PCA.csv"), header=T, as.is=T)

use_meta$cluster_PCA_res_1 = paste0("cluster_", as.vector(pca_cluster$leiden_res_1))
use_meta$cluster_PCA_res_2 = paste0("cluster_", as.vector(pca_cluster$leiden_res_2))

saveRDS(use_meta, paste0(work_path, "/atac.pd.rds"))
### n = 3,937,903 cells left in the final dataset




