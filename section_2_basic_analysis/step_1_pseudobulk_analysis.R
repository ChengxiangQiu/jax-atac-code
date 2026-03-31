
##################################
### Performing pseudobulk analysis

### Support data can be downloaded from:
### https://shendure-web.gs.washington.edu/content/members/cxqiu/public/backup/jax_atac/download/

### Please contact Chengxiang (CX) Qiu for any questions!
### cxqiu@uw.edu or chengxiang.qiu@dartmouth.edu


#######################################
### Step-1: Pseudo-bulk analysis in PCA

work_path = ""
web_path = "https://shendure-web.gs.washington.edu/content/members/cxqiu/public/backup/jax_atac/download"
source("help_code/utils.R")

### read the whole dataset using BPCells
peak_mat = open_matrix_dir(paste0(work_path, "/atac_all"))
colnames(peak_mat) = unlist(lapply(colnames(peak_mat), function(x) strsplit(x,"[-]")[[1]][1]))
### 376,574 peaks x 4,344,905 cells

### read the cell meta data, after removing doublets
use_meta = readRDS(paste0(web_path, "/atac.pd.rds"))

use_meta_sub = use_meta %>% group_by(day) %>% slice_sample(n = 100000)

day_list = names(table(use_meta$day))

peak_mat_aggr = NULL
for(day_i in day_list){
    print(day_i)
    peak_mat_aggr = cbind(peak_mat_aggr, 
                          Matrix::rowSums(peak_mat[,colnames(peak_mat) %in% as.vector(use_meta_sub$cell_id[use_meta_sub$day == day_i]), drop=FALSE]))
}

### exclude peaks overlapping ENCODE blacklist regions or falling on sex chromosomes
exclude = as.vector(read.table(paste0(web_path, "/exclude_peak_sample_matrix.list"))$V1)
exclude = gsub("-","_",exclude)

peak_mat_aggr = peak_mat_aggr[!rownames(peak_mat_aggr) %in% exclude,]
colnames(peak_mat_aggr) = day_list
peak_mat_aggr = Matrix(peak_mat_aggr, sparse = TRUE)

obj = CreateSeuratObject(peak_mat_aggr)
obj = NormalizeData(obj, normalization.method = "LogNormalize", scale.factor = 10000)
obj = FindVariableFeatures(obj, selection.method = "vst", nfeatures = 10000)
gene_use = VariableFeatures(obj)

dat = t(t(peak_mat_aggr) / Matrix::colSums(peak_mat_aggr)) * 100000
dat@x = log(dat@x+1)

dat = GetAssayData(obj, slot = "data")

num_dim = 10
scaling = TRUE
set.seed(2016)
irlba_res = my_sparse_prcomp_irlba(Matrix::t(dat[gene_use,]), 
                                   n = min(num_dim, min(dim(dat)) - 1), 
                                   center = scaling, 
                                   scale. = scaling)
preproc_res = irlba_res$x
row.names(preproc_res) = colnames(dat)

prop_var_expl = irlba_res$sdev^2/sum(irlba_res$sdev^2)
print(prop_var_expl) ### 76.8%, 12.3%, 6.1%

df = data.frame(PC_1 = preproc_res[,1],
                PC_2 = preproc_res[,2],
                PC_3 = preproc_res[,3],
                day = colnames(dat))
df$day = factor(df$day, levels = day_list)

pdf("ATAC_BPCells_pseudobulk.pdf",5,6)
df %>%
    ggplot() +
    geom_point(aes(x = PC_2, y = PC_1), color = "grey80", size=5) +
    geom_point(aes(x = PC_2, y = PC_1, color = day), size=4.5) +
    theme_classic(base_size = 12) + 
    scale_color_manual(values = day_color_plate) + 
    labs(x = "PC_2 (12.3%)", y = "PC_1 (76.8%)") + 
    theme(legend.position="none") +
    theme(axis.text.x = element_text(color="black"), axis.text.y = element_text(color="black")) 
dev.off()

saveRDS(list(peak_mat_aggr = peak_mat_aggr,
             dat = dat,
             gene_use = gene_use,
             df = df), paste0(work_path, "/atac.pseudobulk.rds"))




##################################################
### Step-2: Correlation between timepoints in ATAC

dat_x = readRDS(paste0(work_path, "/atac.pseudobulk.rds"))
peak_mat_aggr = dat_x[["peak_mat_aggr"]]
dat = dat_x[["dat"]]
gene_use = dat_x[["gene_use"]]
df = dat_x[["df"]]

library("gplots")
library(RColorBrewer)
Colors=rev(brewer.pal(11,"Spectral"))
Colors=colorRampPalette(Colors)(120)

x = cor(as.matrix(peak_mat_aggr), method = "spearman")

pdf("ATAC_peak_timepoint_corr.pdf",12,12)
heatmap.2(as.matrix(x), 
          col=Colors, 
          scale="none", 
          Rowv = FALSE, 
          Colv = FALSE, 
          key=T, 
          density.info="none", 
          trace="none", 
          cexRow = 0.5, 
          cexCol = 0.5,
          margins = c(15,15))
dev.off()




##########################################################################################
### Step-3: Plotting the cell composition of each major trajectories as a function of time

pd_atac = readRDS(paste0(web_path, "/atac.pd.rds"))
cell_num = read.csv(paste0(web_path, "/cell_num_prediction.csv"))
x = as.vector(cell_num$day)
x[cell_num$day == "E14.333"] = "E14.375"
cell_num$day = as.vector(x)

df = pd_atac %>%
    group_by(day, celltype_L1) %>%
    tally() %>%
    dplyr::rename(cell_num = n)
df_sub = pd_atac %>%
    group_by(day) %>%
    tally() %>%
    dplyr::rename(cell_num_total = n)
df = df %>%
    left_join(df_sub, by = "day") %>%
    mutate(percentage = cell_num/cell_num_total) %>%
    left_join(cell_num[,c("day","cell_num_pred")], by = "day") %>%
    mutate(cell_num_pred_log2 = log2(cell_num_pred) * percentage)
df$day = factor(df$day, levels = names(day_color_plate_2))
df$celltype_L1 = factor(df$celltype_L1, levels = names(celltype_L1_color_plate))

cell_num_x = c("E10.5" = 10.10,
               "E11.5" = 22.98,
               "E12.5" = 45.03,
               "E13.5" = 60.59,
               "E14.375" = 131.00,
               "E15.5" = 216.79,
               "E16.5" = 353.17,
               "E17.5" = 515.85,
               "E18.5" = 584.78,
               "P0" = 671.50)

df_y = data.frame(day = as.vector(names(cell_num_x)),
                  log2_y = log2(cell_num_x * 1000000))

df_y$day = factor(df_y$day, levels = names(day_color_plate))

cell_num$day = factor(cell_num$day, levels = names(day_color_plate))

p = ggplot() +
    geom_bar(data = df, aes(x = day, y = cell_num_pred_log2, group = celltype_L1, fill = celltype_L1), stat="identity", width = 1) +
    geom_point(data=df_y, aes(x=day, y=log2_y), shape = 21, colour = "black", fill = "white", size = 2, stroke = 1.5, alpha = 0.8) +
    labs(x = "", y = "") +
    theme_classic(base_size = 12) +
    scale_fill_manual(values = celltype_L1_color_plate) +
    theme(legend.position="none") +
    theme(axis.text.x = element_text(color="black", angle = 90, hjust = 1), axis.text.y = element_text(color="black"))

pdf(paste0("Cell_composition_over_time.pdf"), 7, 5)
print(p)
dev.off()





#######################################################################################################################################
### Step-4: Compositions across 12 Level-1 cell lineages were compared between scRNA-seq and scATAC-seq using the Chi-squared statistic

pd_atac = readRDS(paste0(web_path, "/atac.pd.rds"))
pd_rna = readRDS(paste0(web_path, "/rna.pd.rds"))

### 1) keeping 36 matched timepoints
### 2) excluding Testis_and_adrenal from RNA-seq
### 3) merging Neuroectoderm and NMP in ATAC-seq

day = as.vector(pd_rna$day)
day[pd_rna$day == "E14.333"] = "E14.375"
pd_rna$day = as.vector(day)
pd_rna = pd_rna[pd_rna$day %in% names(table(pd_atac$day)),]

pd_rna = pd_rna[pd_rna$major_trajectory != "Testis_and_adrenal",]

tmp = as.vector(pd_atac$celltype_L1)
tmp[pd_atac$celltype_L1 == "NMPs_and_spinal_cord_progenitors"] = "Neuroectoderm"
pd_atac$celltype_L1 = as.vector(tmp)

celltype_match = data.frame(major_trajectory = c("Adipocytes", "Cardiomyocytes", "Definitive_erythroid", "Ependymal_cells", "Eye_and_other", "Intermediate_neuronal_progenitors", "Lung_and_airway", "Megakaryocytes", "Muscle_cells", "Neural_crest_PNS_neurons", "Olfactory_sensory_neurons", "Primitive_erythroid", "B_cells", "CNS_neurons", "Endothelium", "Epithelium", "Hepatocytes", "Intestine", "Mast_cells", "Mesoderm", "Neural_crest_PNS_glia", "Neuroectoderm_and_glia", "Oligodendrocytes", "T_cells", "White_blood_cells"),
                            celltype_L1 = c("Adipocytes", "Muscle_cells", "Erythroid_cells", "Epithelial_cells", "Neuroectoderm", "Neuroectoderm", "Epithelial_cells", "White_blood_cells", "Muscle_cells", "Neuroectoderm", "Olfactory_sensory_neurons", "Erythroid_cells", "White_blood_cells", "Neuroectoderm", "Endothelium", "Epithelial_cells", "Hepatocytes", "Epithelial_cells", "White_blood_cells", "Mesoderm", "Neural_crest_PNS_glia", "Neuroectoderm", "Oligodendrocytes", "White_blood_cells", "White_blood_cells"))

pd_rna = pd_rna %>% left_join(celltype_match, by = "major_trajectory")

### 4) group days to five groups
day = as.vector(pd_rna$day)
day[pd_rna$day %in% c("E10.0","E10.25","E10.5","E10.75")] = "day_group_1"
day[pd_rna$day %in% c("E11.0","E11.25","E11.5","E11.75")] = "day_group_2"
day[pd_rna$day %in% c("E12.0","E12.25","E12.5","E12.75")] = "day_group_3"
day[pd_rna$day %in% c("E13.0","E13.25","E13.5","E13.75")] = "day_group_4"
day[pd_rna$day %in% c("E14.0","E14.25","E14.375","E14.75")] = "day_group_5"
day[pd_rna$day %in% c("E15.0","E15.25","E15.5","E15.75")] = "day_group_6"
day[pd_rna$day %in% c("E16.0","E16.25","E16.5","E16.75")] = "day_group_7"
day[pd_rna$day %in% c("E17.0","E17.25","E17.5")] = "day_group_8"
day[pd_rna$day %in% c("E18.0","E18.25","E18.5","E18.75")] = "day_group_9"
day[pd_rna$day %in% c("P0")] = "day_group_10"
pd_rna$day_group = as.vector(day)

pd_rna$day_group = factor(pd_rna$day_group, levels = paste0("day_group_", 1:10))

day = as.vector(pd_atac$day)
day[pd_atac$day %in% c("E10.0","E10.25","E10.5","E10.75")] = "day_group_1"
day[pd_atac$day %in% c("E11.0","E11.25","E11.5","E11.75")] = "day_group_2"
day[pd_atac$day %in% c("E12.0","E12.25","E12.5","E12.75")] = "day_group_3"
day[pd_atac$day %in% c("E13.0","E13.25","E13.5","E13.75")] = "day_group_4"
day[pd_atac$day %in% c("E14.0","E14.25","E14.375","E14.75")] = "day_group_5"
day[pd_atac$day %in% c("E15.0","E15.25","E15.5","E15.75")] = "day_group_6"
day[pd_atac$day %in% c("E16.0","E16.25","E16.5","E16.75")] = "day_group_7"
day[pd_atac$day %in% c("E17.0","E17.25","E17.5")] = "day_group_8"
day[pd_atac$day %in% c("E18.0","E18.25","E18.5","E18.75")] = "day_group_9"
day[pd_atac$day %in% c("P0")] = "day_group_10"
pd_atac$day_group = as.vector(day)

pd_atac$day_group = factor(pd_atac$day_group, levels = paste0("day_group_", 1:10))


### calculating the proportions
df_rna = pd_rna %>% group_by(celltype_L1, day_group) %>% tally() %>%
    dcast(celltype_L1 ~ day_group, fill = 0)
rownames(df_rna) = df_rna[,1]
df_rna = as.matrix(df_rna[,-1])
df_rna = df_rna + 1
df_rna_norm = t(t(df_rna)/apply(df_rna, 2, sum))

df_atac = pd_atac %>% group_by(celltype_L1, day_group) %>% tally() %>%
    dcast(celltype_L1 ~ day_group, fill = 0)
rownames(df_atac) = df_atac[,1]
df_atac = as.matrix(df_atac[,-1])
df_atac = df_atac + 1
df_atac_norm = t(t(df_atac)/apply(df_atac, 2, sum))

res = matrix(NA, ncol(df_rna_norm), ncol(df_atac_norm))
for(i in 1:ncol(df_rna_norm)){
    for(j in 1:ncol(df_atac_norm)){
        res[i,j] = chisq.test(rbind(df_rna_norm[,i], df_atac_norm[,j]))$statistic
    }
}


library("gplots")
library(viridis)
library(RColorBrewer)
Colors=rev(brewer.pal(11,"Spectral"))
Colors=colorRampPalette(Colors)(120)

pdf(paste0("RNA_ATAC_celltype_composition_Chi_square_statistic.pdf"),12,12)
heatmap.2(as.matrix(res), 
         # col=viridis(100), 
          col = Colors,
          scale="none", 
          Rowv = FALSE, 
          Colv = FALSE, 
          key=T, 
          density.info="none", 
          trace="none", 
          cexRow = 0.5, 
          cexCol = 0.5,
          margins = c(15,15))
dev.off()







