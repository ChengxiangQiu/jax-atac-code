
#############################################################################
### Performing NNLS to exam correlation scATAC-seq with scRNA-seq annotations

### Support data can be downloaded from:
### https://shendure-web.gs.washington.edu/content/members/cxqiu/public/backup/jax_atac/download/

### Please contact Chengxiang (CX) Qiu for any questions!
### cxqiu@uw.edu or chengxiang.qiu@dartmouth.edu

work_path = ""
web_path = "https://shendure-web.gs.washington.edu/content/members/cxqiu/public/backup/jax_atac/download"
source("help_code/utils.R")

pd_atac = readRDS(paste0(work_path, "/atac.pd.rds"))
pd_rna = readRDS(paste0(work_path, "/rna.pd.rds"))
pd_rna = pd_rna[!pd_rna$day %in% c("E8.5","E8.75","E9.0","E9.25","E9.5","E9.75"),]

mouse_gene_sub = mouse_gene_v12[(mouse_gene_v12$gene_type %in% c('protein_coding')) & 
                                 mouse_gene_v12$chr %in% paste0("chr", c(1:19, "X", "Y")),]


##############################
### Step-1: Preparing datasets

### RNA ################

rna_dat = list()

embryo_list = names(table(as.vector(pd_rna$embryo_id)))

for(i in embryo_list){
    print(paste0(i, "/", length(embryo_list)))
    gene_count = readRDS(paste0(work_path_2, "/embryo/", i, "_gene_count.rds"))
    pd = pd_rna[pd_rna$embryo_id == i,]; rownames(pd) = as.vector(pd$cell_id)
    pd = pd[colnames(gene_count),]
    
    gene_count = t(t(gene_count) / Matrix::colSums(gene_count)) * 100000
    gene_count = gene_count[rownames(gene_count) %in% mouse_gene_sub$gene_ID,]
    gene_count@x = log(gene_count@x+1)
    
    celltype_list = names(table(pd$celltype_update))
    for(celltype_i in celltype_list){
        print(celltype_i)
        x = Matrix::rowSums(gene_count[, pd$celltype_update == celltype_i, drop=FALSE])
        if(celltype_i %in% names(rna_dat)){
            rna_dat[[celltype_i]] = rna_dat[[celltype_i]] + as.vector(x)
        } else {
            rna_dat[[celltype_i]] = as.vector(x)
        }
    }
}

saveRDS(rna_dat, paste0(work_path, "/NNLS/rna_dat.rds"))
saveRDS(rownames(gene_count), paste0(work_path, "/NNLS/rna_gene_ID.rds"))


### RNA (sub-celltype of lateral plate mesoderm) ################

pd_sub = pd_rna[!is.na(pd_rna$lateral_plate_mesoderm_sub_clustering),]
embryo_list = names(table(as.vector(pd_sub$embryo_id)))

gene_count = NULL
for(i in embryo_list){
    print(paste0(i, "/", length(embryo_list)))
    gene_count_tmp = readRDS(paste0(work_path, "/", i, "_gene_count.rds"))
    gene_count_tmp = gene_count_tmp[, colnames(gene_count_tmp) %in% pd_sub$cell_id, drop=FALSE]
    gene_count_tmp = t(t(gene_count_tmp) / Matrix::colSums(gene_count_tmp)) * 100000
    gene_count_tmp = gene_count_tmp[rownames(gene_count_tmp) %in% mouse_gene_sub$gene_ID,]
    gene_count_tmp@x = log(gene_count_tmp@x+1)
    gene_count = cbind(gene_count, gene_count_tmp)
}
gene_count = gene_count[,as.vector(pd_sub$cell_id)]

rna_dat = list()
celltype_list = names(table(pd_sub$lateral_plate_mesoderm_sub_clustering))
for(celltype_i in celltype_list){
    print(celltype_i)
    rna_dat[[celltype_i]] = 
        Matrix::rowSums(gene_count[, pd_sub$lateral_plate_mesoderm_sub_clustering == celltype_i, drop=FALSE])
}

saveRDS(rna_dat, paste0(work_path, "/NNLS/rna_dat_2.rds"))


### ATAC ################

atac_dat = list()

for(i in c(1:4)){
    print(i)
    gene_count = readRDS(paste0(work_path, "/gene_count_", i, ".rds"))
    pd = pd_atac; rownames(pd) = as.vector(pd$cell_id)
    pd = pd[colnames(gene_count),]
    
    gene_count = t(t(gene_count) / Matrix::colSums(gene_count)) * 100000
    gene_count = gene_count[rownames(gene_count) %in% mouse_gene_sub$gene_ID,]
    gene_count@x = log(gene_count@x+1)
    
    celltype_list = names(table(pd$celltype_L3))
    for(celltype_i in celltype_list){
        print(celltype_i)
        x = Matrix::rowSums(gene_count[, pd$celltype_L3 == celltype_i, drop=FALSE])
        if(celltype_i %in% names(atac_dat)){
            atac_dat[[celltype_i]] = atac_dat[[celltype_i]] + as.vector(x)
        } else {
            atac_dat[[celltype_i]] = as.vector(x)
        }
    }
}

saveRDS(atac_dat, paste0(work_path, "/NNLS/atac_dat.rds"))
saveRDS(rownames(gene_count), paste0(work_path, "/NNLS/atac_gene_ID.rds"))


### calculate how many cells are present in for each gene in either datasets ############

rna_cell_count = NULL

embryo_list = names(table(as.vector(pd_rna$embryo_id)))

for(i in embryo_list){
    print(paste0(i, "/", length(embryo_list)))
    gene_count = readRDS(paste0(work_path_2, "/", i, "_gene_count.rds"))
    gene_count = gene_count[rownames(gene_count) %in% mouse_gene_sub$gene_ID,]
    gene_count@x[gene_count@x > 0] = 1
    
    if(is.null(rna_cell_count)){
        rna_cell_count = Matrix::rowSums(gene_count)
    } else {
        rna_cell_count = rna_cell_count + Matrix::rowSums(gene_count)
    }
}

saveRDS(rna_cell_count, paste0(work_path, "/NNLS/rna_cell_count.rds"))


atac_cell_count = NULL

for(i in c(1:4)){
    print(i)
    gene_count = readRDS(paste0(work_path, "/gene_count_", i, ".rds"))
    gene_count = gene_count[rownames(gene_count) %in% mouse_gene_sub$gene_ID,]
    gene_count@x[gene_count@x > 0] = 1
    
    if(is.null(atac_cell_count)){
        atac_cell_count = Matrix::rowSums(gene_count)
    } else {
        atac_cell_count = atac_cell_count + Matrix::rowSums(gene_count)
    }
}

saveRDS(atac_cell_count, paste0(work_path, "/NNLS/atac_cell_count.rds"))




################################################
### Step-2: doing NNLS on the level-1 cell types


rna_dat = readRDS(paste0(work_path, "/NNLS/rna_dat.rds"))
atac_dat = readRDS(paste0(work_path, "/NNLS/atac_dat.rds"))

### matched cell type annotations between RNA-seq and ATAC-seq
celltype_match = read.table(paste0(web_path, "/celltype_match.txt"), header=T, as.is=T, sep="\t")

atac_celltype_list = table(pd_atac$celltype_L1)
atac_celltype_list = names(sort(atac_celltype_list, decreasing = T))

celltype_group = c("Adipocytes" = "Adipocytes",
                   "Cardiomyocytes" = "Muscle_cells",
                   "Definitive_erythroid" = "Erythroid_cells",
                   "Ependymal_cells" = "Epithelial_cells",
                   "Eye_and_other" = "Neuroectoderm",
                   "Intermediate_neuronal_progenitors" = "Neuroectoderm",
                   "Lung_and_airway" = "Epithelial_cells",
                   "Megakaryocytes" = "White_blood_cells",
                   "Neural_crest_PNS_glia" = "Neural_crest_PNS_glia",
                   "Neuroectoderm_and_glia" = "Neuroectoderm",
                   "Oligodendrocytes" = "Oligodendrocytes",
                   "T_cells" = "White_blood_cells",
                   "White_blood_cells" = "White_blood_cells",
                   "B_cells" = "White_blood_cells",
                   "CNS_neurons" = "Neuroectoderm",
                   "Endothelium" = "Endothelium",
                   "Epithelial_cells" = "Epithelial_cells",
                   "Hepatocytes" = "Hepatocytes",
                   "Intestine" = "Epithelial_cells",
                   "Mast_cells" = "White_blood_cells",
                   "Mesoderm" = "Mesoderm",
                   "Muscle_cells" = "Muscle_cells",
                   "Neural_crest_PNS_neurons" = "Neuroectoderm",
                   "Olfactory_sensory_neurons" = "Olfactory_sensory_neurons",
                   "Primitive_erythroid" = "Erythroid_cells")
celltype_group = data.frame(RNA_seq_L1 = names(celltype_group),
                            ATAC_seq_L1 = as.vector(celltype_group))

rna_matrix = NULL
atac_matrix = NULL
for(celltype_i in atac_celltype_list){
    print(celltype_i)
    atac_celltype_include = names(table(celltype_match$ATAC_seq_L3[celltype_match$ATAC_seq_L1 == celltype_i]))
    atac_celltype_include = atac_celltype_include[atac_celltype_include != "Missing"]
    
    tmp = as.vector(celltype_group$RNA_seq_L1[celltype_group$ATAC_seq_L1 == celltype_i])
    rna_celltype_include = names(table(celltype_match$RNA_seq_L2[celltype_match$RNA_seq_L1 %in% tmp]))
    rna_celltype_include = rna_celltype_include[rna_celltype_include != "Missing"]
    rna_celltype_include = rna_celltype_include[rna_celltype_include %in% names(rna_dat)]
    
    if(celltype_i == "Mesoderm"){
        rna_celltype_include = c(rna_celltype_include, "Lateral plate and intermediate mesoderm")
    }
    
    if(celltype_i == "NMPs_and_spinal_cord_progenitors"){
        rna_celltype_include = c("NMPs and spinal cord progenitors")
    }
    
    if(celltype_i == "Neuroectoderm"){
        rna_celltype_include = rna_celltype_include[rna_celltype_include != "NMPs and spinal cord progenitors"]
    }
    
    atac_vector = rep(0, length(atac_dat[[1]]))
    for(ii in atac_celltype_include){
        atac_vector = atac_vector + atac_dat[[ii]]
    }
    atac_vector = atac_vector/sum(pd_atac$celltype_L3 %in% atac_celltype_include)
    atac_matrix = cbind(atac_matrix, atac_vector)
    
    rna_vector = rep(0, length(rna_dat[[1]]))
    for(ii in rna_celltype_include){
        rna_vector = rna_vector + rna_dat[[ii]]
    }
    rna_vector = rna_vector/sum(pd_rna$celltype_update %in% rna_celltype_include)
    rna_matrix = cbind(rna_matrix, rna_vector)
}

rna_gene_ID = readRDS(paste0(work_path, "/NNLS/rna_gene_ID.rds"))
atac_gene_ID = readRDS(paste0(work_path, "/NNLS/atac_gene_ID.rds"))

rownames(rna_matrix) = rna_gene_ID; colnames(rna_matrix) = paste0("RNA_", atac_celltype_list)
rownames(atac_matrix) = atac_gene_ID; colnames(atac_matrix) = paste0("ATAC_", atac_celltype_list)
atac_matrix = atac_matrix[rna_gene_ID,]

rna_cell_count = readRDS(paste0(work_path, "/NNLS/rna_cell_count.rds"))
atac_cell_count = readRDS(paste0(work_path, "/NNLS/atac_cell_count.rds"))

rna_gene_include = names(rna_cell_count)[rna_cell_count > (nrow(pd_rna) * 0.01)]
atac_gene_include = names(atac_cell_count)[atac_cell_count > (nrow(pd_atac) * 0.01)]
gene_include = intersect(rna_gene_include, atac_gene_include)

rna_matrix = rna_matrix[gene_include,]
atac_matrix = atac_matrix[gene_include,]


library("gplots")
library("reshape2")
library(viridis)
library(dplyr)
library(tidyr)
source("help_code/NNLS_correlation.R")

conn = correlation_analysis_bidirection(as.matrix(rna_matrix), as.matrix(atac_matrix), fold.change = 1.5, top_gene_num = 1500, spec_gene_num = 1500)
conn$beta = 2*(conn$beta_1+0.01)*(conn$beta_2+0.01)
saveRDS(conn, paste0(work_path, "/NNLS/conn_L1.rds"))

dat <- conn[,c("source","target","beta")]
dat <- dcast(dat, source~target)
rownames(dat) <- dat[,1]; dat <- dat[,-1]

dat = dat[paste0("RNA_", atac_celltype_list),
          paste0("ATAC_", atac_celltype_list)]

pdf("conn_L1.pdf",12,12)
heatmap.2(as.matrix(dat), 
          col=viridis, 
          scale="row", 
          Rowv = FALSE, 
          Colv = FALSE, 
          key=T, 
          density.info="none", 
          trace="none", 
          cexRow = 0.5, 
          cexCol = 0.5,
          margins = c(15,15))
dev.off()




################################################
### Step-3: doing NNLS on the level-2 cell types


rna_dat = readRDS(paste0(work_path, "/NNLS/rna_dat.rds"))
rna_dat_2 = readRDS(paste0(work_path, "/NNLS/rna_dat_2.rds"))
rna_dat = c(rna_dat, rna_dat_2)
atac_dat = readRDS(paste0(work_path, "/NNLS/atac_dat.rds"))

pd_rna_1 = pd_rna[is.na(pd_rna$lateral_plate_mesoderm_sub_clustering),]
pd_rna_2 = pd_rna[!is.na(pd_rna$lateral_plate_mesoderm_sub_clustering),]
pd_rna_2$celltype_update = as.vector(pd_rna_2$lateral_plate_mesoderm_sub_clustering)
pd_rna_x = rbind(pd_rna_1, pd_rna_2); rownames(pd_rna_x) = as.vector(pd_rna_x$cell_id)
pd_rna_x = pd_rna_x[as.vector(pd_rna$cell_id),]

celltype_match = read.table(paste0(web_path, "/celltype_match.txt"), header=T, as.is=T, sep="\t")

atac_celltype_list = names(table(celltype_match$ATAC_seq_L2))
atac_celltype_list = atac_celltype_list[atac_celltype_list != "Missing"]

rna_matrix = NULL
atac_matrix = NULL
for(celltype_i in atac_celltype_list){
    print(celltype_i)
    atac_celltype_include = names(table(celltype_match$ATAC_seq_L3[celltype_match$ATAC_seq_L2 == celltype_i]))
    atac_celltype_include = atac_celltype_include[atac_celltype_include != "Missing"]
    rna_celltype_include = names(table(celltype_match$RNA_seq_L2[celltype_match$ATAC_seq_L2 == celltype_i]))
    rna_celltype_include = rna_celltype_include[rna_celltype_include != "Missing"]
    
    atac_vector = rep(0, length(atac_dat[[1]]))
    for(ii in atac_celltype_include){
        atac_vector = atac_vector + atac_dat[[ii]]
    }
    atac_vector = atac_vector/sum(pd_atac$celltype_L3 %in% atac_celltype_include)
    atac_matrix = cbind(atac_matrix, atac_vector)
    
    rna_vector = rep(0, length(rna_dat[[1]]))
    for(ii in rna_celltype_include){
        rna_vector = rna_vector + rna_dat[[ii]]
    }
    rna_vector = rna_vector/sum(pd_rna_x$celltype_update %in% rna_celltype_include)
    rna_matrix = cbind(rna_matrix, rna_vector)
}

rna_gene_ID = readRDS(paste0(work_path, "/NNLS/rna_gene_ID.rds"))
atac_gene_ID = readRDS(paste0(work_path, "/NNLS/atac_gene_ID.rds"))

rownames(rna_matrix) = rna_gene_ID; colnames(rna_matrix) = paste0("RNA_", atac_celltype_list)
rownames(atac_matrix) = atac_gene_ID; colnames(atac_matrix) = paste0("ATAC_", atac_celltype_list)
atac_matrix = atac_matrix[rna_gene_ID,]

rna_cell_count = readRDS(paste0(work_path, "/NNLS/rna_cell_count.rds"))
atac_cell_count = readRDS(paste0(work_path, "/NNLS/atac_cell_count.rds"))

rna_gene_include = names(rna_cell_count)[rna_cell_count > (nrow(pd_rna) * 0.01)]
atac_gene_include = names(atac_cell_count)[atac_cell_count > (nrow(pd_atac) * 0.01)]
gene_include = intersect(rna_gene_include, atac_gene_include)

rna_matrix = rna_matrix[gene_include,]
atac_matrix = atac_matrix[gene_include,]

rna_matrix = rna_matrix[,!colnames(rna_matrix) %in% c("RNA_Erythroid_progenitors_L2")]

library("gplots")
library("reshape2")
library(viridis)
library(dplyr)
library(tidyr)
source("help_code/NNLS_correlation.R")

conn = correlation_analysis_bidirection(as.matrix(rna_matrix), as.matrix(atac_matrix), fold.change = 1.5, top_gene_num = 1500, spec_gene_num = 1500)
conn$beta = 2*(conn$beta_1+0.01)*(conn$beta_2+0.01)
saveRDS(conn, paste0(work_path, "/NNLS/conn_L2.rds"))

dat <- conn[,c("source","target","beta")]
dat <- dcast(dat, source~target)
rownames(dat) <- dat[,1]; dat <- dat[,-1]

pdf("conn_L2.pdf",12,12)
heatmap.2(as.matrix(dat), 
          col=viridis, 
          scale="row", 
          Rowv = FALSE, 
          Colv = FALSE, 
          key=T, 
          density.info="none", 
          trace="none", 
          cexRow = 0.5, 
          cexCol = 0.5,
          margins = c(15,15))
dev.off()




###################################################################################################
### Step-4: doing NNLS on the level-3 cell types, after manually grouping cells into several groups


rna_dat = readRDS(paste0(work_path, "/NNLS/rna_dat.rds"))
rna_dat_2 = readRDS(paste0(work_path, "/NNLS/rna_dat_2.rds"))
rna_dat = c(rna_dat, rna_dat_2)
atac_dat = readRDS(paste0(work_path, "/NNLS/atac_dat.rds"))

pd_rna_1 = pd_rna[is.na(pd_rna$lateral_plate_mesoderm_sub_clustering),]
pd_rna_2 = pd_rna[!is.na(pd_rna$lateral_plate_mesoderm_sub_clustering),]
pd_rna_2$celltype_update = as.vector(pd_rna_2$lateral_plate_mesoderm_sub_clustering)
pd_rna_x = rbind(pd_rna_1, pd_rna_2); rownames(pd_rna_x) = as.vector(pd_rna_x$cell_id)
pd_rna_x = pd_rna_x[as.vector(pd_rna$cell_id),]

celltype_match = read.table(paste0(web_path, "/celltype_match.txt"), header=T, as.is=T, sep="\t")

celltype_group = c("Adipocytes" = "Mesoderm",
                   "Epithelial_cells" = "Epithelium_and_endoderm",
                   "Hepatocytes" = "Epithelium_and_endoderm",
                   "Neural_crest_PNS_glia" = "Neuroectoderm",
                   "NMPs_and_spinal_cord_progenitors" = "Neuroectoderm",
                   "Oligodendrocytes" = "Neuroectoderm",
                   "Endothelium" = "Endothelium",
                   "Erythroid_cells" = "Erythroid",
                   "Mesoderm" = "Mesoderm",
                   "Muscle_cells" = "Mesoderm",
                   "Neuroectoderm" = "Neuroectoderm",
                   "Olfactory_sensory_neurons" = "Neuroectoderm",
                   "White_blood_cells" = "White_blood_cells")
celltype_group = data.frame(ATAC_seq_L1 = names(celltype_group),
                            group = as.vector(celltype_group))
group_list = names(table(celltype_group$group))

for(group_i in group_list){
    
    print(group_i)
    atac_celltype_list = names(table(celltype_match$ATAC_seq_L3[celltype_match$ATAC_seq_L1 %in% celltype_group$ATAC_seq_L1[celltype_group$group == group_i]]))
    atac_celltype_list = atac_celltype_list[atac_celltype_list != "Missing"]
    
    rna_matrix = NULL
    atac_matrix = NULL
    for(celltype_i in atac_celltype_list){
        print(celltype_i)
        
        atac_vector = atac_dat[[celltype_i]]
        atac_vector = atac_vector/sum(pd_atac$celltype_L3 == celltype_i)
        atac_matrix = cbind(atac_matrix, atac_vector)
        
        rna_celltype_include = names(table(celltype_match$RNA_seq_L2[celltype_match$ATAC_seq_L3 == celltype_i]))
        rna_celltype_include = rna_celltype_include[rna_celltype_include != "Missing"]
    
        rna_vector = rep(0, length(rna_dat[[1]]))
        for(ii in rna_celltype_include){
            rna_vector = rna_vector + rna_dat[[ii]]
        }
        rna_vector = rna_vector/sum(pd_rna_x$celltype_update %in% rna_celltype_include)
        rna_matrix = cbind(rna_matrix, rna_vector)
    }
    
    rna_gene_ID = readRDS(paste0(work_path, "/NNLS/rna_gene_ID.rds"))
    atac_gene_ID = readRDS(paste0(work_path, "/NNLS/atac_gene_ID.rds"))
    
    rownames(rna_matrix) = rna_gene_ID; colnames(rna_matrix) = paste0("RNA_", atac_celltype_list)
    rownames(atac_matrix) = atac_gene_ID; colnames(atac_matrix) = paste0("ATAC_", atac_celltype_list)
    atac_matrix = atac_matrix[rna_gene_ID,]
    
    rna_cell_count = readRDS(paste0(work_path, "/NNLS/rna_cell_count.rds"))
    atac_cell_count = readRDS(paste0(work_path, "/NNLS/atac_cell_count.rds"))
    
    rna_gene_include = names(rna_cell_count)[rna_cell_count > (nrow(pd_rna) * 0.01)]
    atac_gene_include = names(atac_cell_count)[atac_cell_count > (nrow(pd_atac) * 0.01)]
    gene_include = intersect(rna_gene_include, atac_gene_include)
    
    rna_matrix = rna_matrix[gene_include,]
    atac_matrix = atac_matrix[gene_include,]
    
    rna_matrix = rna_matrix[,!colnames(rna_matrix) %in% c("RNA_Erythroid_progenitors", "RNA_Ureteric_bud_stalk_cells")]
    
    library("gplots")
    library("reshape2")
    library(viridis)
    library(dplyr)
    library(tidyr)
    source("help_code/NNLS_correlation.R")
    
    conn = correlation_analysis_bidirection(as.matrix(rna_matrix), as.matrix(atac_matrix), fold.change = 1.5, top_gene_num = 1500, spec_gene_num = 1500)
    conn$beta = 2*(conn$beta_1+0.01)*(conn$beta_2+0.01)
    saveRDS(conn, paste0(work_path, "/NNLS/conn_L3_", group_i, ".rds"))
    
    dat <- conn[,c("source","target","beta")]
    dat <- dcast(dat, source~target)
    rownames(dat) <- dat[,1]; dat <- dat[,-1]
    
    pdf(paste0(work_path, "/conn_L3_", group_i, ".pdf"),12,12)
    heatmap.2(as.matrix(dat), 
              col=viridis, 
              scale="row", 
              Rowv = FALSE, 
              Colv = FALSE, 
              key=T, 
              density.info="none", 
              trace="none", 
              cexRow = 0.5, 
              cexCol = 0.5,
              margins = c(15,15))
    dev.off()
}





