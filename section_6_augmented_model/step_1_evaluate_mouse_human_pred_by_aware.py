

##########################################################################################
### evaluating prediction on orthologs between human and mouse using evolution-aware model

### Support data can be downloaded from:
### https://shendure-web.gs.washington.edu/content/members/cxqiu/public/backup/jax_atac/download/

### Please contact Chengxiang (CX) Qiu for any questions!
### cxqiu@uw.edu or chengxiang.qiu@dartmouth.edu



###############################
### Step-1: prediction on human

import sys, os
import anndata as ad
import crested
import numpy as np
import pandas as pd
import gzip
import keras
import pysam
from collections import defaultdict

work_path = ""
model_id = ""
mamm = "Homo_sapiens"

model_path = f"evolution-aware-model.keras"
model = keras.models.load_model(model_path, compile=False)

window_include = set()
with open(f"{web_path}/orthologs_seq_include.all.txt") as f:
    for line in f:
        l = line.rstrip().split('\t')
        if l[1] == mamm:
            window_include.add(l[0])

window_list = []
window_coor = []
with gzip.open(f"{work_path}/genome/{mamm}/candidate_windows_sub.txt.gz", "rt") as file:
    for line in file:
        chr, start, end, window_id = line.rstrip().split("\t")
        if window_id in window_include:
            window_list.append(window_id)
            window_coor.append((chr, int(start), int(end)))

predictions = []
for iter in range(1,11):
    fasta = pysam.FastaFile(f"{work_path}/genome/{mamm}/fasta_trf/{mamm}_trf_{iter}.fa")
    candidate_sequences = []
    for region in window_coor:
        candidate_sequences.append(fasta.fetch(region[0], region[1], region[2]))
    predictions.append(crested.tl.predict(input=candidate_sequences, model=model))

predictions_arr = np.stack(predictions, axis = 0)
predictions_mean = np.mean(predictions_arr, axis = 0)

with gzip.open(f"{work_path}/{model_id}/human_mouse_aware_augmented/{mamm}_predictions.multispecies.txt.gz", 'wt') as f:
    np.savetxt(f, predictions_mean, fmt="%.3f")

with open(f"{work_path}/{model_id}/human_mouse_aware_augmented/{mamm}_window_list.txt", "w") as f:
    for item in window_list:
        f.write(f"{item}\n")


###############################
### Step-2: prediction on mouse


import sys, os
import anndata as ad
import crested
import numpy as np
import pandas as pd
import gzip
import keras
import pysam
from collections import defaultdict


mamm = "Mus_musculus"

model_path = f"evolution-aware-model.keras"
model = keras.models.load_model(model_path, compile=False)

window_include = set()
with open(f"{web_path}/orthologs_seq_include.all.txt") as f:
    for line in f:
        l = line.rstrip().split('\t')
        if l[1] == "Homo_sapiens":
            window_include.add(l[0])

window_include_2 = set()
with gzip.open(f"{work_path}/genome/Homo_sapiens/candidate_windows_sub.txt.gz", "rt") as file:
    for line in file:
        chr, start, end, window_id = line.rstrip().split("\t")
        if window_id in window_include:
            window_include_2.add(window_id)

window_list = []
window_coor = []
with open(f"{work_path}/window_include.txt") as file:
    for line in file:
        l = line.rstrip().split("\t")
        if l[0] in window_include_2:
            window_list.append(l[0])
            window_coor.append((l[3], int(l[4]) - 1007, int(l[5]) + 1007))

predictions = []
for iter in range(1,11):
    fasta = pysam.FastaFile(f"{work_path}/genome/{mamm}/fasta_trf/{mamm}_trf_{iter}.fa")
    candidate_sequences = []
    for region in window_coor:
        candidate_sequences.append(fasta.fetch(region[0], region[1], region[2]))
    predictions.append(crested.tl.predict(input=candidate_sequences, model=model))

predictions_arr = np.stack(predictions, axis = 0)
predictions_mean = np.mean(predictions_arr, axis = 0)

with gzip.open(f"{work_path}/{model_id}/human_mouse_aware_augmented/{mamm}_predictions.multispecies.txt.gz", 'wt') as f:
    np.savetxt(f, predictions_mean, fmt="%.3f")

with open(f"{work_path}/{model_id}/human_mouse_aware_augmented/{mamm}_window_list.x.txt", "w") as f:
    for item in window_list:
        f.write(f"{item}\n")



###################################
### Step-3: making the scatter plot

work_path = ""
web_path = "https://shendure-web.gs.washington.edu/content/members/cxqiu/public/backup/jax_atac/download"
source("help_code/utils.R")

celltype_list = c("Adipocyte_cells","Adipocyte_cells_Cyp2e1","B_cells","Brain_capillary_endothelial_cells","CNS_neurons","Cardiomyocytes","Corticofugal_neurons","Endocardial_cells","Endothelium","Epithelial_cells","Erythroid_cells","Eye","Glia","Glomerular_endothelial_cells","Gut_epithelial_cells","Hepatocytes","Intermediate_neuronal_progenitors","Kidney","Lateral_plate_and_intermediate_mesoderm","Liver_sinusoidal_endothelial_cells","Lung_and_airway","Lymphatic_vessel_endothelial_cells","Melanocyte_cells","Mesoderm","Neural_crest_PNS_neurons","Neuroectoderm_and_glia","Olfactory_ensheathing_cells","Olfactory_neurons","Oligodendrocytes","Skeletal_muscle_cells","T_cells","White_blood_cells")

window_list = read.table(paste0(web_path, "/window_include.txt"))
window_list = window_list[,c(1,3)]
colnames(window_list) = c("window_ID", "cell_class")

dat_mouse = read.table(paste0(work_path, "/Mus_musculus_predictions.aware.txt.gz"))
dat_id_mouse = read.table(paste0(work_path, "/Mus_musculus_window_list.txt"))

dat_human = read.table(paste0(work_path, "/Homo_sapiens_predictions.aware.txt.gz"))
dat_id_human = read.table(paste0(work_path, "/Homo_sapiens_window_list.txt"))

dat_mouse = as.matrix(dat_mouse)
rownames(dat_mouse) = as.vector(dat_id_mouse$V1)
colnames(dat_mouse) = celltype_list

df_mouse = melt(dat_mouse)
colnames(df_mouse) = c("window_ID", "cell_class", "mouse_pred")

dat_human = as.matrix(dat_human)
rownames(dat_human) = as.vector(dat_id_human$V1)
colnames(dat_human) = celltype_list

df_human = melt(dat_human)
colnames(df_human) = c("window_ID", "cell_class", "human_pred")

df_1 = window_list %>% left_join(df_mouse, by = c("window_ID", "cell_class")) %>% filter(!is.na(mouse_pred))
df_2 = window_list %>% left_join(df_human, by = c("window_ID", "cell_class")) %>% filter(!is.na(human_pred))

df = df_1 %>% left_join(df_2, by = c("window_ID", "cell_class"))

df$log_mouse_pred = log1p(df$mouse_pred)
df$log_human_pred = log1p(df$human_pred)

df$model = "aware"

p = ggplot(data = df, aes(x=log_mouse_pred, y=log_human_pred)) +
    geom_hex(bins = 70) +
    geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dashed", linewidth = 0.8) +
    scale_x_continuous(limits = c(-0.3, 5), expand = expansion(mult = 0, add = 0)) +
    scale_y_continuous(limits = c(-0.3, 5), expand = expansion(mult = 0, add = 0)) +
    labs(x="Log(mouse prediction + 1)", y="Log(human prediction + 1)", title="") +
    scale_fill_continuous(type = "viridis", trans = "log10") +
    theme_classic(base_size = 10) +
    theme(legend.position="none") +
    theme(plot.title = element_text(hjust = 0.5)) +
    theme(axis.text.x = element_text(color="black"), axis.text.y = element_text(color="black"))
ggsave(paste0("scatter_mouse_human_pred_by_aware.pdf"), p, height = 5, width = 5)

sum(df$log_human_pred > df$log_mouse_pred)/nrow(df)
### 30%

cor.test(df$log_human_pred, df$log_mouse_pred)
### 0.50


### n = 288982 pairs
df$log2fc = log2(df$human_pred/df$mouse_pred)

summary(df$log2fc)

   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
-7.1814 -1.4421 -0.5861 -0.6434  0.1716  7.0895

p = ggplot(data = df, aes(x = log2fc)) +
    geom_density(alpha = 0.5, fill = "#cb6751") +
    labs(x="Log2 fold-change of human pred vs. mouse pred", y="Density", title="") +
    theme_classic(base_size = 10) +
    theme(legend.position="none") +
    theme(axis.text.x = element_text(color="black"), axis.text.y = element_text(color="black")) +
    coord_cartesian(xlim = c(-8, 8)) 
ggsave(paste0("foldchange_mouse_human_pred_by_aware.pdf"), p, height = 3, width = 5)








work_path = ""
web_path = "https://shendure-web.gs.washington.edu/content/members/cxqiu/public/backup/jax_atac/download"
source("help_code/utils.R")

celltype_list = c("Adipocyte_cells","Adipocyte_cells_Cyp2e1","B_cells","Brain_capillary_endothelial_cells","CNS_neurons","Cardiomyocytes","Corticofugal_neurons","Endocardial_cells","Endothelium","Epithelial_cells","Erythroid_cells","Eye","Glia","Glomerular_endothelial_cells","Gut_epithelial_cells","Hepatocytes","Intermediate_neuronal_progenitors","Kidney","Lateral_plate_and_intermediate_mesoderm","Liver_sinusoidal_endothelial_cells","Lung_and_airway","Lymphatic_vessel_endothelial_cells","Melanocyte_cells","Mesoderm","Neural_crest_PNS_neurons","Neuroectoderm_and_glia","Olfactory_ensheathing_cells","Olfactory_neurons","Oligodendrocytes","Skeletal_muscle_cells","T_cells","White_blood_cells")

window_list = read.table(paste0(web_path, "/window_include.txt"))
window_list = window_list[,c(1,3)]
colnames(window_list) = c("window_ID", "cell_class")

dat_mouse = read.table(paste0(work_path, "/Mus_musculus_predictions.multispecies.txt.gz"))
dat_id_mouse = read.table(paste0(work_path, "/Mus_musculus_window_list.txt"))

dat_human = read.table(paste0(work_path, "/Homo_sapiens_predictions.multispecies.txt.gz"))
dat_id_human = read.table(paste0(work_path, "/Homo_sapiens_window_list.txt"))

dat_mouse = as.matrix(dat_mouse)
rownames(dat_mouse) = as.vector(dat_id_mouse$V1)
colnames(dat_mouse) = celltype_list

df_mouse = melt(dat_mouse)
colnames(df_mouse) = c("window_ID", "cell_class", "mouse_pred")

dat_human = as.matrix(dat_human)
rownames(dat_human) = as.vector(dat_id_human$V1)
colnames(dat_human) = celltype_list

df_human = melt(dat_human)
colnames(df_human) = c("window_ID", "cell_class", "human_pred")

df_1 = window_list %>% left_join(df_mouse, by = c("window_ID", "cell_class")) %>% filter(!is.na(mouse_pred))
df_2 = window_list %>% left_join(df_human, by = c("window_ID", "cell_class")) %>% filter(!is.na(human_pred))

df = df_1 %>% left_join(df_2, by = c("window_ID", "cell_class"))

df$log_mouse_pred = log1p(df$mouse_pred)
df$log_human_pred = log1p(df$human_pred)

df$model = "multispecies"

p = ggplot(data = df, aes(x=log_mouse_pred, y=log_human_pred)) +
    geom_hex(bins = 70) +
    geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dashed", linewidth = 0.8) +
    scale_x_continuous(limits = c(-0.3, 6), expand = expansion(mult = 0, add = 0)) +
    scale_y_continuous(limits = c(-0.3, 6), expand = expansion(mult = 0, add = 0)) +
    labs(x="Log(mouse prediction + 1)", y="Log(human prediction + 1)", title="") +
    scale_fill_continuous(type = "viridis", trans = "log10") +
    theme_classic(base_size = 10) +
    #theme(legend.position="none") +
    theme(plot.title = element_text(hjust = 0.5)) +
    theme(axis.text.x = element_text(color="black"), axis.text.y = element_text(color="black"))
ggsave(paste0("scatter_mouse_human_pred_by_multispecies.pdf"), p, height = 5, width = 5)

sum(df$log_human_pred > df$log_mouse_pred)/nrow(df)
### 38%

cor.test(df$log_human_pred, df$log_mouse_pred)
### 0.67


### n = 288982 pairs
df$log2fc = log2(df$human_pred/df$mouse_pred)

summary(df$log2fc[!is.na(df$log2fc) & !is.infinite(df$log2fc)])

    Min.  1st Qu.   Median     Mean  3rd Qu.     Max.
-13.5049  -0.8072  -0.2267  -0.2794   0.2788  11.621

p = ggplot(data = df[!is.na(df$log2fc) & !is.infinite(df$log2fc),], aes(x = log2fc)) +
    geom_density(alpha = 0.5, fill = "blue") +
    labs(x="Log2 fold-change of human pred vs. mouse pred", y="Density", title="") +
    theme_classic(base_size = 10) +
    theme(legend.position="none") +
    theme(axis.text.x = element_text(color="black"), axis.text.y = element_text(color="black")) +
    coord_cartesian(xlim = c(-8, 8)) 
ggsave(paste0("foldchange_mouse_human_pred_by_multispecies.pdf"), p, height = 3, width = 5)






















