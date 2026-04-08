
###################################################
### Performing UMAP on left windows after filtering

### Support data can be downloaded from:
### https://shendure-web.gs.washington.edu/content/members/cxqiu/public/backup/jax_atac/download/

### Please contact Chengxiang (CX) Qiu for any questions!
### cxqiu@uw.edu or chengxiang.qiu@dartmouth.edu


######################################
### Of note:
### 1) n = 2,790,521 windows (100bp), which passed Z > 3 for both predicted score and specificity, and merging across cell types
### 2) n = 1,465,915 windows (100bp), which passed liftover in at least 120 mammals
### 3) n = 547,317 windows (100bp), which passed median correlation > 0.6 across mammals


###################################################
### Step-1: save dataset with the subset of windows

work_path = ""
web_path = "https://shendure-web.gs.washington.edu/content/members/cxqiu/public/backup/jax_atac/download"
source("help_code/utils.R")

window_list = readRDS(paste0(work_path, "/corr_Mus_musculus/window_list.rds"))
window_list_uniq = readRDS(paste0(work_path, "/corr_Mus_musculus/window_list_uniq.rds"))
window_list_uniq_x = as.vector(window_list_uniq$window_ID)
### n = 547317, which are passed median corr > 0.6 acorss mammals

dat = read.table(paste0(work_path, "/tmp_Mus_musculus/dat_sub.txt"))
pd = read.table(paste0(work_path, "/tmp_Mus_musculus/peak_list_sub.txt"))
### n = 1,465,915 windows

colnames(pd) = c("chr", "start", "end")
pd$window_id = paste0(pd$chr, "_", pd$start, "_", pd$end)
rownames(dat) = as.vector(pd$window_id)
dat = dat[as.vector(window_list_uniq$window_id),]
rownames(dat) = window_list_uniq_x

write.table(dat, paste0(work_path, "/umap_Mus_musculus/dat.txt"), row.names = F, col.names = F, sep = "\t", quote = F)


###########################
### Step-2: performing UMAP

from annoy import AnnoyIndex
import numpy as np
import pandas as pd
import umap
import os, sys

work_path = "XXX"
mat = np.loadtxt(f"{work_path}/umap_Mus_musculus/dat.txt")

#### performing UMAP
reducer = umap.UMAP(
    n_neighbors=30,
    min_dist=0.1,
    metric="euclidean",
    n_components=2,
    random_state=42)
embedding = reducer.fit_transform(mat)
pd.DataFrame(embedding).to_csv(f"{work_path}/umap_Mus_musculus/umap_embedding_N30.txt", sep="\t", header=False, index=False)

reducer = umap.UMAP(
    n_neighbors=30,
    min_dist=0.1,
    metric="euclidean",
    n_components=3,
    random_state=42)
embedding = reducer.fit_transform(mat)
pd.DataFrame(embedding).to_csv(f"{work_path}/umap_Mus_musculus/umap_embedding_N30_3D.txt", sep="\t", header=False, index=False)

#### creating Knn graph using annoy (fast way)
dist_metric = 'euclidean'
k = 30
ncell = mat.shape[0]
npc = mat.shape[1]
annoy_index = AnnoyIndex(npc, metric=dist_metric)

for i in range(ncell):
    annoy_index.add_item(i, list(mat[i,:]))

annoy_index.build(15) ### bigger number will make the result more accurate

knn = []
for iCell in range(ncell):
    knn.append(annoy_index.get_nns_by_item(iCell, k + 1)[1:])

knn = np.array(knn, dtype=int)
np.savetxt(f"{work_path}/umap_Mus_musculus/kNN_30.txt", knn, delimiter=",", fmt='%s')



############################
### Step-3: plotting 2D umap

work_path = ""
web_path = "https://shendure-web.gs.washington.edu/content/members/cxqiu/public/backup/jax_atac/download"
source("help_code/utils.R")

window_list = readRDS(paste0(work_path, "/corr_Mus_musculus/window_list.rds"))
window_list_uniq = readRDS(paste0(work_path, "/corr_Mus_musculus/window_list_uniq.rds"))
window_list_uniq_x = as.vector(window_list_uniq$window_ID)
### n = 547317, which are passed median corr > 0.6 acorss mammals

emb = read.table(paste0(work_path, "/umap_Mus_musculus/umap_embedding_N30.txt"))
pd = data.frame(window_ID = window_list_uniq_x, UMAP_1 = emb[,1], UMAP_2 = emb[,2])
rownames(emb) = as.vector(pd$window_ID)

pd = pd %>% left_join(window_list_uniq[,c("window_ID","promoter")], by = "window_ID")

nn_matrix = read.table(paste0(work_path, "/umap_Mus_musculus/kNN_30.txt"), sep=",")
nn_matrix = as.matrix(nn_matrix)
nn_matrix = nn_matrix + 1 ### python index to R index

df_x = NULL
for(i in 1:ncol(nn_matrix)){
    df_x = cbind(df_x, pd$promoter[nn_matrix[,i]])
}
pd$promoter_enrich = apply(df_x == "yes", 1, sum)/ncol(nn_matrix)

p1 = ggplot() +
    geom_point(data = pd, aes(UMAP_1, UMAP_2, color = promoter_enrich), size = 0.03) +
    labs(x="UMAP 1", y="UMAP 2", title="Promoter enrichment (547,317 windows)") +
    theme_classic(base_size = 10) +
    theme(plot.title = element_text(hjust = 0.5)) +
    theme(legend.position="none") +
    theme(axis.text.x = element_text(color="black"), axis.text.y = element_text(color="black")) +
    scale_color_viridis_c(option = "viridis")

pd_big = pd %>% left_join(window_list[,c("window_ID","celltype")], by = "window_ID")
pd_big_coor = pd_big %>% group_by(celltype) %>% summarize(mean_UMAP_1 = mean(UMAP_1), mean_UMAP_2 = mean(UMAP_2))
p2 = ggplot() +
    geom_point(data = pd_big[sample(1:nrow(pd_big)),], aes(UMAP_1, UMAP_2, color = celltype) , size = 0.03) +
    geom_text(data = pd_big_coor, aes(mean_UMAP_1, mean_UMAP_2, label = celltype), size = 2) + 
    labs(x="UMAP 1", y="UMAP 2", title="Cell types (547,317 windows)") +
    theme_classic(base_size = 10) +
    theme(plot.title = element_text(hjust = 0.5)) +
    theme(legend.position="none") +
    theme(axis.text.x = element_text(color="black"), axis.text.y = element_text(color="black")) +
    scale_color_manual(values=celltype_L2_color_plate)

ggsave(paste0("Mouse_umap_promoter_celltype.png"), p1 + p2, height = 5, width = 10, dpi = 300)



#############################
### Step-4: leiden clustering

library(igraph)

edge_list = data.frame(A = rep(1:nrow(nn_matrix), ncol(nn_matrix)), B = c(nn_matrix))
g = graph_from_edgelist(as.matrix(edge_list), directed = FALSE)
g = simplify(g)
cl_louvain = cluster_louvain(g, resolution = 1)

pd$louvain_cluster = paste0("cluster_", as.vector(membership(cl_louvain)))

fig = plot_ly(pd, x=~UMAP_3D_1, y=~UMAP_3D_2, z=~UMAP_3D_3, size = I(30), color = ~louvain_cluster)
saveWidget(fig, paste0(save_path, "/14_crested/Mus_musculus_3D_UMAP_louvain_res_1.html"), selfcontained = FALSE, libdir = "tmp")

pd_big_num = pd_big %>% left_join(pd[,c("window_ID","louvain_cluster")], by = "window_ID") %>%
    group_by(louvain_cluster, celltype) %>% tally() %>%
    left_join(pd_big %>% left_join(pd[,c("window_ID","louvain_cluster")], by = "window_ID") %>% group_by(louvain_cluster) %>% tally() %>% rename(total_n = n), by = "louvain_cluster") %>%
    mutate(frac = n/total_n) %>% filter(frac > 0.1)

anno = rep(NA, nrow(pd))

anno[pd$louvain_cluster %in% paste0("cluster_", c(3,4,7,10,14,16,17,19,20,21))] = "Neuroectoderm_promoters"

anno[pd$louvain_cluster %in% paste0("cluster_", c(8))] = "Eye"
anno[pd$louvain_cluster %in% paste0("cluster_", c(9,27))] = "Intermediate_neuronal_progenitors"

anno[pd$louvain_cluster %in% paste0("cluster_", c(2))] = "Adipocytes"

anno[pd$louvain_cluster %in% paste0("cluster_", c(12))] = "Mesoderm"
anno[pd$louvain_cluster %in% paste0("cluster_", c(23))] = "Skeletal_muscle_cells"
anno[pd$louvain_cluster %in% paste0("cluster_", c(26))] = "Cardiomyocytes"

anno[pd$louvain_cluster %in% paste0("cluster_", c(5,18))] = "Epithelial_cells"

anno[pd$louvain_cluster %in% paste0("cluster_", c(15,28))] = "Endothelium"

anno[pd$louvain_cluster %in% paste0("cluster_", c(22))] = "White_blood_cells"

anno[pd$louvain_cluster %in% paste0("cluster_", c(25))] = "Erythroid_cells"

anno[pd$louvain_cluster %in% paste0("cluster_", c(1,13,24))] = "Glia"

anno[pd$louvain_cluster %in% paste0("cluster_", c(6,11))] = "Olfactory_neurons"

pd$anno_L1 = as.vector(anno)
fig = plot_ly(pd, x=~UMAP_3D_1, y=~UMAP_3D_2, z=~UMAP_3D_3, size = I(30), color = ~anno_L1)
saveWidget(fig, paste0(save_path, "/Mus_musculus_3D_UMAP_anno_L1.html"), selfcontained = FALSE, libdir = "tmp")

saveRDS(pd, paste0(work_path, "/umap_Mus_musculus/pd.rds"))

### create a heatmap for comparing cell type vs. major cluster 

df = pd %>% left_join(window_list[,c("window_ID","celltype")], by = "window_ID") %>%
    group_by(celltype, anno_L1) %>% tally() %>% dcast(celltype ~ anno_L1)
df[is.na(df)] = 0
rownames(df) = df[,1]
df = df[,-1]

df_x = t(t(df)/apply(df, 2, sum))
df_x = df_x/apply(df_x, 1, sum)

library("gplots")
library(RColorBrewer)
library(viridis)
Colors=rev(brewer.pal(11,"Spectral"))
Colors=colorRampPalette(Colors)(120)
pdf("heatmap.pdf", 8, 5)
heatmap.2(as.matrix(df_x), 
          col=viridis(100), 
          scale="none", 
          Rowv = F, 
          Colv = F, 
          key=T, 
          density.info="none", 
          trace="none", 
          cexRow = 1, 
          cexCol = 1,
          margins = c(5,5))
dev.off()




####################################################
### Step-5: checking distance to TSS for each window

window_list_uniq = readRDS(paste0(work_path, "/corr_Mus_musculus/window_list_uniq.rds"))
write.table(window_list_uniq[,c("chr", "start", "end", "window_ID")],
            paste0(work_path, "/umap_Mus_musculus/window_list.bed"), row.names=F, col.names=F, sep="\t", quote=F)


### bedtools intersect -a window_list.bed -b $DATAPATH/TSS_250kb.bed -wa -wb > window_list_overlap_TSS_250kb.bed

>>> R
dat = read.table(paste0(work_path,
                        "/umap_Mus_musculus/window_list_overlap_TSS_250kb.bed"))
colnames(dat) = c("chr","start","end","window_ID","gene_chr","gene_start","gene_end","gene_strand","transcript_ID")
dat$TSS = round((dat$gene_start + dat$gene_end)/2)
dat$distance_TSS = round((dat$start + dat$end)/2 - (dat$gene_start + dat$gene_end)/2)
dat$abs_distance_TSS = abs(dat$distance_TSS)
dat_min = dat %>% group_by(window_ID) %>% slice_min(order_by = abs_distance_TSS, n = 1, with_ties = F)

pd = readRDS(paste0(work_path, "/umap_Mus_musculus/pd.rds"))
pd = pd %>% left_join(dat_min[,c("window_ID","chr","start","end","distance_TSS","TSS","transcript_ID")], by = "window_ID")

mouse_gene = read.table("mouse.v12.geneID.transcriptID.txt", header = T, as.is = T)
pd = pd %>% left_join(mouse_gene[,c("transcript_ID", "gene_short_name")], by = "transcript_ID")

pd_out = pd[,c("anno_L2","window_ID","chr","start","end","TSS","distance_TSS","transcript_ID","gene_short_name")]
pd_out = pd_out[pd_out$anno_L2 != 'Promoters' & !is.na(pd_out$TSS),]
write.csv(pd_out, "window_list_nearest_TSS_gene.csv", quote=F)

pd$distance_TSS[is.na(pd$distance_TSS)] = 250000
pd$anno_L2 = factor(pd$anno_L2, levels = rev(col_names))
pd$abs_distance_TSS = log2(abs(pd$distance_TSS)+1)

p1 = pd %>% group_by(anno_L2) %>% tally() %>% mutate(log2_num = log2(n)) %>%
    ggplot(aes(x = anno_L2, y = log2_num, fill = anno_L2)) +
    geom_bar(stat="identity") +
    labs(x="", y="Log2 #", title="") +
    theme_classic(base_size = 10) +
    theme(plot.title = element_text(hjust = 0.5)) +
    theme(legend.position="none") +
    theme(axis.text.x = element_text(color="black"), axis.text.y = element_text(color="black")) +
    scale_fill_manual(values=window_cluster_color_plate) +
    coord_flip()

p2 = pd %>% ggplot() +
    geom_boxplot(aes(x = anno_L2, y = abs_distance_TSS, fill = anno_L2)) +
    labs(x="", y="Log2 (abs distance to nearest TSS)", title="") +
    theme_classic(base_size = 10) +
    theme(plot.title = element_text(hjust = 0.5)) +
    theme(legend.position="none") +
    theme(axis.text.x = element_text(color="black"), axis.text.y = element_text(color="black")) +
    scale_fill_manual(values=window_cluster_color_plate) +
    coord_flip()

ggsave("window_cluster_dist_num.pdf", p1 + p2, height = 10, width = 10)


#################################################
### Step-6: merging windows within each cell type

pd = readRDS(paste0(work_path, "/umap_Mus_musculus/pd.rds"))
window_list_uniq = readRDS(paste0(work_path, "/corr_Mus_musculus/window_list_uniq.rds"))
anno_L2 = names(table(pd$anno_L2))
for(anno in anno_L2){
    print(anno)
    x = window_list_uniq[pd$anno_L2 == anno,]
    write.table(x[,c("chr","start","end","window_ID")], 
                paste0(work_path, "/umap_Mus_musculus/window_merge/", anno, ".bed"), row.names=F, col.names=F, sep="\t", quote=F)
}

>>>
FILELIST=(Adipocyte_cells Adipocyte_cells_Cyp2e1 Cardiomyocytes Skeletal_muscle_cells Mesoderm Lateral_plate_and_intermediate_mesoderm White_blood_cells B_cells T_cells Brain_capillary_endothelial_cells Endocardial_cells Glomerular_endothelial_cells Liver_sinusoidal_endothelial_cells Lymphatic_vessel_endothelial_cells Endothelium Erythroid_cells Promoters Neuroectoderm_and_glia CNS_neurons Neural_crest_PNS_neurons Olfactory_neurons Corticofugal_neurons Glia Olfactory_ensheathing_cells Melanocyte_cells Oligodendrocytes Intermediate_neuronal_progenitors Eye Lung_and_airway Gut_epithelial_cells Hepatocytes Kidney Epithelial_cells)
for FILE in "${FILELIST[@]}"; do
    echo ${FILE}
    bedtools sort -i "$DATAPATH/${FILE}.bed" | bedtools merge > "$DATAPATH/${FILE}.merge.bed"
    bedtools intersect -a "$DATAPATH/${FILE}.bed" -b "$DATAPATH/${FILE}.merge.bed" -wa -wb > "$DATAPATH/${FILE}.merge.overlap.bed"
done

>>> how mnay # of merged windows in each cell type
merge_num = NULL
windows_merge = NULL
for(anno in anno_L2){
    print(anno)
    x = read.table(paste0(work_path, "/umap_Mus_musculus/window_merge/", anno, ".merge.bed"))
    merge_num = rbind(merge_num, data.frame(anno_L2 = anno, num = nrow(x)))
    x$anno_L2 = anno
    windows_merge = rbind(windows_merge, x)
}
merge_num$anno_L2 = factor(merge_num$anno_L2, levels = rev(col_names))
p3 = merge_num %>% mutate(log2_num = log2(num)) %>%
    ggplot(aes(x = anno_L2, y = log2_num, fill = anno_L2)) +
    geom_bar(stat="identity") +
    labs(x="", y="Log2 # of merged windows", title="") +
    theme_classic(base_size = 10) +
    theme(plot.title = element_text(hjust = 0.5)) +
    theme(legend.position="none") +
    theme(axis.text.x = element_text(color="black"), axis.text.y = element_text(color="black")) +
    scale_fill_manual(values=window_cluster_color_plate) +
    coord_flip()

ggsave("window_cluster_dist_num.pdf", p1 + p2 + p3, height = 10, width = 15)

x = pd %>% group_by(anno_L2) %>% tally() %>% rename(num_window = n) %>%
    left_join(merge_num, by = "anno_L2") %>% rename(num_merged_window = num)
x$anno_L2 = factor(merge_num$anno_L2, levels = col_names)
write.csv(x, "window_num.csv", quote=F)

windows_merge = windows_merge[windows_merge$anno_L2 != "Promoters",]
write.table(windows_merge, "window_merged.txt", row.names=F, col.names=F, sep="\t", quote=F)


