
#####################################
### Detecting doublets using scrublet

### Support data can be downloaded from:
### https://shendure-web.gs.washington.edu/content/members/cxqiu/public/backup/jax_atac/download/

### Please contact Chengxiang (CX) Qiu for any questions!
### cxqiu@uw.edu or chengxiang.qiu@dartmouth.edu


#############################################
### Step-1: detecting doublets using Scrublet

work_path = ""
web_path = "https://shendure-web.gs.washington.edu/content/members/cxqiu/public/backup/jax_atac/download"
source("help_code/utils.R")

pd_all = readRDS(paste0(work_path, "/pd_all.rds"))
rownames(pd_all) = as.vector(pd_all$cell_id)
pd_all = pd_all[pd_all$keep,]
print(dim(pd_all))
### n = 4,344,905 cells, before filtering doublets

### first, randomly split the whole dataset into 40 batches (roughly 100K cells per batch)
pd_all$split_group = sample(1:40, nrow(pd_all), replace = T)
pd_x = pd_all[,c("cell_id", "split_group", "sample_id", "sample")]
rownames(pd_x) = NULL
saveRDS(pd_x, paste0(work_path, "/pd_split_batch.rds"))


##############################################################
### Step-2: performing scrublet on individual batch (n = 1:40)

args = commandArgs(trailingOnly=TRUE)
kk = as.numeric(args[1])

work_path = ""
web_path = "https://shendure-web.gs.washington.edu/content/members/cxqiu/public/backup/jax_atac/download"
source("help_code/utils.R")

source("help_code/scrublet.R")
source("help_code/lsi.R")

### exclude peaks overlapping ENCODE blacklist regions or falling on sex chromosomes
exclude = as.vector(read.table(paste0(web_path, "/exclude_peak_sample_matrix.list"))$V1)
exclude = gsub("-","_",exclude)

pd_all = readRDS(paste0(work_path, "/pd_split_batch.rds"))
rownames(pd_all) = as.vector(pd_all$cell_id)
pd_all = pd_all[pd_all$split_group == kk,]
print(kk)

pd_x = unique(pd_all[,c("sample_id","sample")])
sample_list = as.vector(pd_x$sample)

peak_matrix = NULL
for(sample_i in sample_list){
    x = pd_x$sample_id[pd_x$sample == sample_i]
    
    counts = Read10X(paste0(work_path, sample_i), gene.column = 1) 
    colnames(counts) = paste0(x, "_", colnames(counts))
    
    counts_sub = counts[!rownames(counts) %in% exclude,
                        colnames(counts) %in% as.vector(pd_all$cell_id)]
    
    peak_matrix = cbind(peak_matrix, counts_sub)
}

### binary matrix and remove low express cells and peaks
### filter out peaks which are detected in less than 0.1% cells
### filter out cells which are deteced less than 100 peaks
peak_matrix@x[peak_matrix@x > 0] = 1
peak_matrix = filter_features(peak_matrix, cells=ncol(peak_matrix) * 0.001)

use_meta = pd_all[colnames(peak_matrix),]

print(dim(use_meta))
print(dim(peak_matrix))

saveRDS(use_meta, paste0(work_path, "/scrublet_", kk, "_pd.rds"))

scrub_res = atac_scrublet(peak_matrix)
saveRDS(scrub_res, paste0(work_path, "/scrublet_", kk, ".rds"))
### done ###


######################################
### Step-3: merging all the batch 1:40

work_path = ""

res_all = NULL
for(kk in 1:40){
    print(kk)
    res_i = readRDS(paste0(work_path, "/scrublet_", kk, ".rds"))
    res_all = rbind(res_all, res_i)
}
saveRDS(res_all, paste0(work_path, "/scrublet_result.rds"))

pd_all = readRDS(paste0(work_path, "/pd_all.rds"))
res = res_all[!res_all$simulated_doublet,]
res$cell_id = as.vector(res$cell)
pd_x = pd_all %>% left_join(res[,c("cell_id","doublet_score","doublet_likelihood")])
pd_all$doublet_score = as.vector(pd_x$doublet_score)
pd_all$doublet_likelihood = as.vector(pd_x$doublet_likelihood)
saveRDS(pd_all, paste0(work_path, "/pd_all.rds"))


###########################################################
### Step-4: two clusters enriched by doublets were detected

work_path = ""
web_path = "https://shendure-web.gs.washington.edu/content/members/cxqiu/public/backup/jax_atac/download"
source("help_code/utils.R")

### Of note, the atac_all.h5ad was created by Scanpy; please check the script "help_code/create_scanpy_object.py"

if (!file.exists(paste0(work_path, "/atac_all"))) {
    peak_mat = open_matrix_anndata_hdf5(paste0(work_path, "/atac_all.h5ad")) %>% 
        convert_matrix_type(type = "uint32_t") %>%
        write_matrix_dir(paste0(work_path, "/atac_all"))
} else {
    peak_mat = open_matrix_dir(paste0(work_path, "/atac_all"))
}
peak_mat
### 376,574 peaks x 4,344,905 cells

### exclude peaks overlapping ENCODE blacklist regions or falling on sex chromosomes
exclude = as.vector(read.table(paste0(web_path, "/exclude_peak_sample_matrix.list"))$V1)
exclude = gsub("-","_",exclude)

### exclude peaks which express less in 0.1% of cells
peak_sum = Matrix::rowSums(peak_mat)
peak_include = names(peak_sum)[peak_sum >= ncol(peak_mat) * 0.001]
peak_include = peak_include[!peak_include %in% exclude]

peak_mat = peak_mat[rownames(peak_mat) %in% peak_include,]

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
saveRDS(pca_mat, paste0(work_path, "/atac_all.pca.rds"))

pca_l2_mat = Seurat:::L2Norm(pca_mat)
saveRDS(pca_l2_mat, paste0(work_path, "/atac_all.pca_l2.rds"))

write.csv(round(pca_l2_mat[,dims], 3), paste0(work_path, "/atac_all.pca_l2.csv"))

set.seed(12341512)
umap_mat = uwot::umap(pca_l2_mat[,dims], n_components = 3, metric = "cosine", 
                      n_neighbors = 50, min_dist = 0.1)
umap_mat[1:4,]

rownames(umap_mat) = colnames(peak_mat)
saveRDS(umap_mat, paste0(work_path, "/atac_all.umap.rds"))

library(future)
options(future.globals.maxSize = 50000 * 1024^2)

seurat_obj = Seurat::CreateSeuratObject(peak_mat[1:100,])
seurat_obj[['umap']] = Seurat::CreateDimReducObject(embeddings=umap_mat, key='UMAP_', assay='RNA')
seurat_obj = seurat_obj %>% 
    Seurat::FindNeighbors(reduction="umap", nn.eps=0, dims=1:3, k.param=20) %>%
    Seurat::FindClusters(n.start=20, resolution=0.1) %>%
    Seurat::FindClusters(n.start=20, resolution=0.3) %>%
    Seurat::FindClusters(n.start=20, resolution=0.5)

saveRDS(data.frame(seurat_obj[[]]), paste0(work_path, "/atac_all.clustering_UMAP.rds"))

use_meta = read.csv(paste0(work_path, "/atac_all.obs.csv"), row.names=1, header=T)
print(sum(rownames(use_meta) == rownames(umap_mat)))
print(sum(rownames(use_meta) == colnames(seurat_obj)))

use_meta$UMAP_1 = as.vector(umap_mat[,1])
use_meta$UMAP_2 = as.vector(umap_mat[,2])
use_meta$UMAP_3 = as.vector(umap_mat[,3])

use_meta$cluster_UMAP_0.1 = paste0("cluster_", as.vector(seurat_obj$RNA_snn_res.0.1))
use_meta$cluster_UMAP_0.3 = paste0("cluster_", as.vector(seurat_obj$RNA_snn_res.0.3))
use_meta$cluster_UMAP_0.5 = paste0("cluster_", as.vector(seurat_obj$RNA_snn_res.0.5))

rownames(use_meta) = NULL
saveRDS(use_meta, paste0(work_path, "/atac_all.pd.rds"))

