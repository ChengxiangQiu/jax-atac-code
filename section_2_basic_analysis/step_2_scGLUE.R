
###################################################################
### Performing integration of scRNA-seq and scATAC-seq using scGLUE

### Support data can be downloaded from:
### https://shendure-web.gs.washington.edu/content/members/cxqiu/public/backup/jax_atac/download/

### Please contact Chengxiang (CX) Qiu for any questions!
### cxqiu@uw.edu or chengxiang.qiu@dartmouth.edu

############################################################################################
### Step-1: creating dataset from rna-seq, after downsampling 15K cells for individual stage

work_path = ""
web_path = "https://shendure-web.gs.washington.edu/content/members/cxqiu/public/backup/jax_atac/download"
source("help_code/utils.R")

pd_atac = readRDS(paste0(web_path, "/atac.pd.rds"))
day_list = names(table(pd_atac$day))
day_list = c(day_list, "E14.333")

pd_rna = readRDS(paste0(web_path, "/rna.pd.rds"))

pd_rna_downsample = pd_rna %>% filter(day %in% day_list) %>% group_by(day) %>% sample_n(15000)

pd_sub = pd_rna[pd_rna$cell_id %in% pd_rna_downsample$cell_id,]

rownames(pd_sub) = as.vector(pd_sub$cell_id)

fd = mouse_gene_v12[(mouse_gene_v12$gene_type %in% c('protein_coding', 'pseudogene', 'lincRNA')) & 
    mouse_gene_v12$chr %in% paste0("chr", c(1:19, "M")),]

embryo_list = names(table(pd_sub$embryo_id))
print(embryo_list)
gene_count = NULL
for(j in embryo_list){
    print(paste0(j,"/",length(embryo_list)))
    count_j = readRDS(paste0(work_path, "/", j, "_gene_count.rds"))
    keep = colnames(count_j) %in% rownames(pd_sub)
    gene_count = cbind(gene_count, count_j[rownames(count_j) %in% rownames(fd),keep])
}

df_cell = pd_sub[colnames(gene_count),]
df_gene = mouse_gene_v12[rownames(gene_count),]

Matrix::writeMM(t(gene_count), file = paste0(work_path, "/rna_count.mtx"))
write.csv(df_cell, paste0(work_path, "/rna_obs.csv"))
write.csv(df_gene, paste0(work_path, "/rna_var.csv"))


##########################################
### Step-2: creating dataset from atac-seq

work_path = ""
web_path = "https://shendure-web.gs.washington.edu/content/members/cxqiu/public/backup/jax_atac/download"
source("help_code/utils.R")

peak_mat = open_matrix_dir(paste0(work_path, "/atac_all"))
colnames(peak_mat) = unlist(lapply(colnames(peak_mat), function(x) strsplit(x,"[-]")[[1]][1])) 
### 376,574 peaks x 4,344,905 cells

use_meta = readRDS(paste0(web_path, "/atac.pd.rds"))

use_meta_x = use_meta %>% group_by(day) %>% sample_n(15000)
use_meta_sub = use_meta[use_meta$cell_id %in% use_meta_x$cell_id,]
rownames(use_meta_sub) = as.vector(use_meta_sub$cell_id)
use_meta_sub = use_meta_sub[,c("cell_id", "sample", "day", "sex", "reads", "reads_peak",
                               "global_cluster", "major_trajectory")]

peak_mat = peak_mat[,colnames(peak_mat) %in% rownames(use_meta_sub)]
use_meta_sub = use_meta_sub[colnames(peak_mat),]

### exclude peaks overlapping ENCODE blacklist regions or falling on sex chromosomes
exclude = as.vector(read.table(paste0(web_path, "/exclude_peak_sample_matrix.list"))$V1)
exclude = gsub("-","_",exclude)

### exclude peaks which express less in 0.1% of cells
peak_sum = Matrix::rowSums(peak_mat)
peak_include = names(peak_sum)[peak_sum >= ncol(peak_mat) * 0.001]
peak_include = peak_include[!peak_include %in% exclude]

peak_mat = peak_mat[rownames(peak_mat) %in% peak_include,]
print(dim(peak_mat))
### 310828 peaks X  540000 cells

fd = data.frame(peak_ID = rownames(peak_mat),
                chr = unlist(lapply(rownames(peak_mat), function(x) strsplit(x,"[_]")[[1]][1])),
                start = unlist(lapply(rownames(peak_mat), function(x) strsplit(x,"[_]")[[1]][2])),
                end = unlist(lapply(rownames(peak_mat), function(x) strsplit(x,"[_]")[[1]][3])), stringsAsFactors = FALSE)
rownames(fd) = as.vector(fd$peak_ID)

Matrix::writeMM(t(as(peak_mat, "dgCMatrix")), file = paste0(work_path, "/atac_count.mtx"))
write.csv(use_meta_sub, paste0(work_path, "/atac_obs.csv"))
write.csv(fd, paste0(work_path, "/atac_var.csv"))


######################
### Step-3: run scGLUE

python help_code/run_scglue.py


#############################################################################
### Step-4: after running scGLUE, performing UMAP on the shared feature space

rna = read.csv(paste0(work_path, "/rna_obs.csv"), row.names=1)
rna_coor = read.csv(paste0(work_path, "/rna_X_glue.csv"), header=F)

atac = read.csv(paste0(work_path, "/atac_obs.csv"), row.names=1)
atac_coor = read.csv(paste0(work_path, "/atac_X_glue.csv"), header=F)

combine_coor = rbind(rna_coor, atac_coor)

set.seed(2016)
umap_model = uwot::umap(as.matrix(combine_coor), 
                            n_components = 3,
                            n_neighbors = 30,
                            min_dist = 0.3,
                            metric = "cosine",
                            fast_sgd = FALSE,
                            nn_method = "annoy",
                            ret_model = TRUE,
                            n_threads = 1,
                            verbose = TRUE)

set.seed(2016)
umap_coor = uwot::umap_transform(as.matrix(combine_coor),
                                     umap_model)
umap_coor = data.frame(umap_coor)
names(umap_coor) = paste0("UMAP_", 1:3)

rna_select = rna[,c("cell_id", "day", "major_trajectory", 
                    "celltype_update", "embryo_id")]
rna_select$group = "RNA"

atac$embryo_id = atac$sample
atac$major_trajectory = as.vector(atac$major_trajectory)
atac$celltype_update = as.vector(atac$global_cluster)
atac$group = "ATAC"
atac_select = atac[,c("cell_id", "day", "major_trajectory", 
                      "celltype_update", "embryo_id", "group")]

df = rbind(rna_select, atac_select)
df = cbind(df, umap_coor)

saveRDS(df, paste0(work_path, "/combined_df_cell.rds"))

save_path = ""

fig = plot_ly(df[sample(1:nrow(df), 200000),], x=~UMAP_1, y=~UMAP_2, z=~UMAP_3, size = I(30), color = ~group) %>% 
    layout(scene = list(xaxis=list(title = list(text ='UMAP_1', font = t1), tickfont = t2),
                        yaxis=list(title = list(text ='UMAP_2', font = t1), tickfont = t2),
                        zaxis=list(title = list(text ='UMAP_3', font = t1), tickfont = t2)))
saveWidget(fig, paste0(save_path, "/RNA_ATAC_scglue_group.html"), selfcontained = FALSE, libdir = "tmp")

df$major_trajectory = paste0(df$group, "_", df$major_trajectory)
fig = plot_ly(df[sample(1:nrow(df), 200000),], x=~UMAP_1, y=~UMAP_2, z=~UMAP_3, size = I(30), color = ~major_trajectory)
saveWidget(fig, paste0(save_path, "/RNA_ATAC_scglue_major_trajectory.html"), selfcontained = FALSE, libdir = "tmp")

fig = plot_ly(df[sample(1:nrow(df), 200000),], x=~UMAP_1, y=~UMAP_2, z=~UMAP_3, size = I(30), color = ~celltype_update)
saveWidget(fig, paste0(save_path, "/RNA_ATAC_scglue_celltype_update.html"), selfcontained = FALSE, libdir = "tmp")

day_color_plate_2 = day_color_plate
day_color_plate = c(day_color_plate_2, "E14.333" = "#FEFAB7")
day_color_plate_1 = day_color_plate
names(day_color_plate_1) = paste0("RNA_", names(day_color_plate_1))
day_color_plate_2 = day_color_plate
names(day_color_plate_2) = paste0("ATAC_", names(day_color_plate_2))
day_color_plate = c(day_color_plate_1, day_color_plate_2)

df$group_day = paste0(df$group, "_", df$day)
fig = plot_ly(df[sample(1:nrow(df), 200000),], x=~UMAP_1, y=~UMAP_2, z=~UMAP_3, size = I(30), color = ~group_day, colors = day_color_plate)
saveWidget(fig, paste0(save_path, "/RNA_ATAC_scglue_day.html"), selfcontained = FALSE, libdir = "tmp")


##########################################################################################################
### Step-5: tranfering cell type labels for ATACseq data based on kNN in the coembedding learned by scGLUE

work_path = ""
web_path = "https://shendure-web.gs.washington.edu/content/members/cxqiu/public/backup/jax_atac/download"
source("help_code/utils.R")


rna = read.csv(paste0(work_path, "/rna_obs.csv"), row.names=1)
rna_coor = read.csv(paste0(work_path, "/rna_X_glue.csv"), header=F)
rna_coor = as.matrix(rna_coor)
rownames(rna_coor) = as.vector(rna$cell_id)

atac = read.csv(paste0(work_path, "/atac_obs.csv"), row.names=1)
atac_coor = read.csv(paste0(work_path, "/atac_X_glue.csv"), header=F)
atac_coor = as.matrix(atac_coor)
rownames(atac_coor) = as.vector(atac$cell_id)

k.param = 10; nn.method = "rann"; nn.eps = 0; annoy.metric = "euclidean"
nn.ranked = Seurat:::NNHelper(
    data = rna_coor,
    query = atac_coor,
    k = k.param,
    method = nn.method,
    searchtype = "standard",
    eps = nn.eps,
    metric = annoy.metric)
nn.ranked = Indices(object = nn.ranked)
nn_matrix = nn.ranked
saveRDS(nn_matrix, paste0(work_path, "/nn_matrix.rds"))

resultA = NULL
for(i in 1:k.param){
    print(i)
    resultA = cbind(resultA, as.vector(rna$major_trajectory)[as.vector(nn_matrix[,i])])
}
resultB = lapply(1:nrow(resultA), function(x){
    return(names(sort(table(resultA[x,]), decreasing = T))[1])
})
atac$major_trajectory = as.vector(unlist(resultB))

resultA = NULL
for(i in 1:k.param){
    print(i)
    resultA = cbind(resultA, as.vector(rna$celltype_update)[as.vector(nn_matrix[,i])])
}
resultB = lapply(1:nrow(resultA), function(x){
    return(names(sort(table(resultA[x,]), decreasing = T))[1])
})
atac$celltype_update = as.vector(unlist(resultB))

saveRDS(atac, paste0(work_path, "/atac.rds"))


#################################################################
### Step-6: checking the timepoints matching between RNA and ATAC

result = NULL
for(i in 1:k.param){
    print(i)
    result = cbind(result, as.vector(rna$day)[as.vector(nn_matrix[,i])])
}

df = data.frame(atac_day = rep(as.vector(atac$day), k.param),
                rna_day = c(result))

day_list = c("E10.0", "E10.25", 
             "E10.5", "E10.75", "E11.0", "E11.25", "E11.5", "E11.75", "E12.0", 
             "E12.25", "E12.5", "E12.75", "E13.0", "E13.25", "E13.5", "E13.75", 
             "E14.0", "E14.25", "E14.333", "E14.375", "E14.75", "E15.0", "E15.25", "E15.5", 
             "E15.75", "E16.0", "E16.25", "E16.5", "E16.75", "E17.0", "E17.25", 
             "E17.5", "E17.75", "E18.0", "E18.25", "E18.5", "E18.75", "P0")

df$atac_day = factor(df$atac_day, levels = day_list[day_list %in% df$atac_day])
df$rna_day = factor(df$rna_day, levels = day_list[day_list %in% df$rna_day])

saveRDS(df, paste0(work_path, "/df_overlap_day.rds"))

library(ggridges)
rna_day = as.vector(df$rna_day)
rna_day[df$rna_day == "P0"] = "E19"
df$rna_day = as.numeric(as.vector(gsub("E","",rna_day)))

breaks <- seq(10, 19, by = 1)
pdf("day_match.pdf", 8, 5)
ggplot(df, aes(x = rna_day, y = atac_day, fill = atac_day)) +
    geom_density_ridges(bandwidth = 0.6) +
    theme_classic(base_size = 10) +
    scale_fill_viridis(discrete = TRUE) +
    scale_x_continuous(
        breaks = breaks,
        labels = paste0("E", breaks)
    ) +
    theme(legend.position = "none") + 
    coord_flip() +
    labs(
        x = "Staging bins of scRNA-seq nearest neighbors",
        y = "scATAC-seq staging bin"
    ) +
    theme(
        axis.text.x = element_text(color="black", angle = 90, vjust = 1, hjust=1),
        axis.text.y = element_text(color="black")
    )
dev.off()



