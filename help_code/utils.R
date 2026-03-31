
print("Loading Monocle3 and Seurat")
suppressMessages(library(monocle3))
suppressMessages(library(Seurat))

print("Loading packages for regular data analysis, e.g. dplyr")
suppressMessages(library(Matrix))
suppressMessages(library(dplyr))
suppressMessages(library(reshape2))
suppressMessages(library(scales))

print("Loading packages for plotting, e.g. ggplot2")
suppressMessages(library(ggplot2))
suppressMessages(library(plotly))
suppressMessages(library(htmlwidgets))
suppressMessages(library(gridExtra))
suppressMessages(library(gplots))
suppressMessages(library(patchwork))
suppressMessages(library(viridis))
suppressMessages(library(RColorBrewer))

print("Loading packages for atac-seq analysis")
suppressMessages(library(BPCells))
suppressMessages(library(Signac))
suppressMessages(library(GenomeInfoDb))
suppressMessages(library(EnsDb.Mmusculus.v79))
suppressMessages(library(genomation))

Colors=rev(brewer.pal(11,"Spectral"))
Colors=colorRampPalette(Colors)(120)

web_path = "https://shendure-web.gs.washington.edu/content/members/cxqiu/public/backup/jax_atac/download"

### read mouse gene ID
mouse_gene_v12 <- read.table(paste0(web_path, "/mouse.v12.geneID.txt"), header = T, as.is = T)
rownames(mouse_gene_v12) = as.vector(mouse_gene_v12$gene_ID)

### global setting
# options(scipen = 999)

### adding some help functions
### my_plot_cells()
### estimateSex()
### my_sparse_prcomp_irlba()
### my_cluster_cells()
### split_view()
source("./help_code.R")
source("./color_code_jax_atac.R")

### basic setting for ggplot2
plot_tmp = theme_classic(base_size = 12) +
    theme(legend.position="none") +
    theme(plot.title = element_text(hjust = 0.5)) +
    theme(axis.text.x = element_text(color="black"), axis.text.y = element_text(color="black")) 

#####################################################
### Function: doing regular analysis using Seurat ###
#####################################################

doClusterSeurat <- function(obj, 
                            nfeatures = 2500, 
                            resolution = 1, 
                            k.filter = 200, 
                            doClustering = TRUE, 
                            n_dim = 30, 
                            min.dist = 0.3, 
                            n.components = 2){
  
  if(length(table(obj$group))!=1){
    
    obj.list <- SplitObject(object = obj, split.by = "group")
    
    for (i in 1:length(x = obj.list)) {
      obj.list[[i]] <- NormalizeData(object = obj.list[[i]], verbose = FALSE)
      obj.list[[i]] <- FindVariableFeatures(object = obj.list[[i]], 
                                            selection.method = "vst", nfeatures = nfeatures, verbose = FALSE)
    }
    
    reference.list <- obj.list[names(table(obj$group))]
    obj.anchors <- FindIntegrationAnchors(object.list = reference.list, dims = 1:n_dim, k.filter = k.filter)
    
    obj.integrated <- IntegrateData(anchorset = obj.anchors, dims = 1:n_dim)
    
    DefaultAssay(object = obj.integrated) <- "integrated"
    
    obj <- obj.integrated 
    
  } else {
    
    obj <- NormalizeData(obj, normalization.method = "LogNormalize", scale.factor = 10000)
    obj <- FindVariableFeatures(obj, selection.method = "vst", nfeatures = nfeatures)

  }
  
  obj <- ScaleData(object = obj, verbose = FALSE)
  obj <- RunPCA(object = obj, npcs = n_dim, verbose = FALSE)
    
  if (doClustering == TRUE){
    obj <- FindNeighbors(object = obj, dims = 1:n_dim, reduction = "pca")
    obj <- FindClusters(object = obj, resolution = resolution)
  }
  
  obj <- RunUMAP(object = obj, reduction = "pca", dims = 1:n_dim, min.dist = min.dist, n.components = n.components)

  return(obj)
}

###################################################################
### Function: transition data object between monocle and seurat ###
###################################################################

doObjectTransform <- function(x, transform_to = NULL){
  if (transform_to == "seurat"){
    
    count = exprs(x)
    meta.data = data.frame(pData(x))
    
    obj = CreateSeuratObject(count, meta.data = meta.data)
    obj = NormalizeData(obj, normalization.method = "LogNormalize", scale.factor = 10000)
    
    return(obj)
    
  } else {
    
    x = NormalizeData(x, normalization.method = "LogNormalize", scale.factor = 10000)
    x = FindVariableFeatures(x, selection.method = "vst", nfeatures = 2500)
    genes_include = VariableFeatures(x)
    
    count = GetAssayData(object = x, slot = "counts")
    meta.data = data.frame(x[[]])
    mouse_gene_sub = mouse_gene[rownames(count),]
    rownames(mouse_gene_sub) = rownames(count)
    
    cds = new_cell_data_set(count,
                            cell_metadata = meta.data,
                            gene_metadata = mouse_gene_sub)
    cds = preprocess_cds(cds, use_genes = genes_include)
    
    return(cds)
  }
}

####################################################
### Function: identify partitions using monocle3 ###
####################################################

doIdentifyPartitions = function(umap_emb, cell_meta){
    
    nn_control = monocle3:::set_nn_control(mode=3,
                                           nn_control=list(),
                                           nn_control_default=monocle3:::get_global_variable('nn_control_annoy_euclidean'),
                                           nn_index=NULL,
                                           k=20,
                                           verbose=FALSE)
    
    cluster_result = monocle3:::leiden_clustering(data=as.matrix(umap_emb),
                                                  pd=data.frame(cell_meta),
                                                  weight=FALSE,
                                                  nn_index=NULL,
                                                  k=20,
                                                  nn_control=nn_control,
                                                  num_iter=2,
                                                  resolution_parameter=NULL,
                                                  random_seed=42,
                                                  verbose=FALSE)
    
    cluster_graph_res = monocle3:::compute_partitions(cluster_result$g,
                                                      cluster_result$optim_res,
                                                      qval_thresh = 0.05, verbose=FALSE)
    
    partitions = igraph::components(cluster_graph_res$cluster_g)$membership[cluster_result$optim_res$membership]
    clusters = igraph::membership(cluster_result$optim_res)
    return(list(as.vector(partitions), as.vector(clusters)))
}




