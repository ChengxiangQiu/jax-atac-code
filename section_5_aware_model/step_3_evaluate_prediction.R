
#################################################################################
### Evaluating predition performance based on held-out set or genome-wide windows

### Support data can be downloaded from:
### https://shendure-web.gs.washington.edu/content/members/cxqiu/public/backup/jax_atac/download/

### Please contact Chengxiang (CX) Qiu for any questions!
### cxqiu@uw.edu or chengxiang.qiu@dartmouth.edu


############################################################
### Step-1: save the prediction outputs to R (not necessary)

work_path = ""
web_path = "https://shendure-web.gs.washington.edu/content/members/cxqiu/public/backup/jax_atac/download"
source("help_code/utils.R")

library(GenomicRanges)

model_id = ""

mamm = "Mus_musculus"

celltype_list = c("Adipocyte_cells","Adipocyte_cells_Cyp2e1","B_cells","Brain_capillary_endothelial_cells","CNS_neurons","Cardiomyocytes","Corticofugal_neurons","Endocardial_cells","Endothelium","Epithelial_cells","Erythroid_cells","Eye","Glia","Glomerular_endothelial_cells","Gut_epithelial_cells","Hepatocytes","Intermediate_neuronal_progenitors","Kidney","Lateral_plate_and_intermediate_mesoderm","Liver_sinusoidal_endothelial_cells","Lung_and_airway","Lymphatic_vessel_endothelial_cells","Melanocyte_cells","Mesoderm","Neural_crest_PNS_neurons","Neuroectoderm_and_glia","Olfactory_ensheathing_cells","Olfactory_neurons","Oligodendrocytes","Skeletal_muscle_cells","T_cells","White_blood_cells")

dat = read.table(paste0(work_path, "/14_crested/", model_id, "/prediction_mammals/prediction_", mamm, "_trf/dat.txt.gz"))
dat_loc = read.table(paste0(work_path, "/14_crested/", model_id, "/prediction_mammals/prediction_", mamm, "_trf/dat_loc.txt.gz"))

colnames(dat) = celltype_list
colnames(dat_loc) = c("chr", "start")
dat_loc$end = dat_loc$start + 100

saveRDS(dat, paste0(work_path, "/14_crested/", model_id, "/prediction_mammals/prediction_", mamm, "_trf/dat.txt.rds"))
saveRDS(dat_loc, paste0(work_path, "/14_crested/", model_id, "/prediction_mammals/prediction_", mamm, "_trf/dat_loc.txt.rds"))



#######################################################
### Step-2: based on the held-out set (TRF replacement)

import sys, os
import anndata as ad
import crested
import numpy as np
import matplotlib
import gzip
import keras
import pysam

model_id = ""
mamm = "Mus_musculus"

work_path = ""
adata = ad.read_h5ad(os.path.join(work_path, model_id, "data_window_cluster_top3K.h5ad"))
print(adata.var["split"].value_counts())

# load a trained model
model_path = f"{web_path}/CREsted_model/evolution_aware_model.keras"
model = keras.models.load_model(model_path, compile=False)

test_df = adata.var[adata.var["split"] == "test"]
all_regions = list(zip(test_df["chr"], test_df["start"], test_df["end"]))

dat_batch = []
for i in range(1, 11):
    print(i)
    fasta = pysam.FastaFile(f"{work_path}/genome/{mamm}/fasta_trf/{mamm}_trf_{i}.fa")
    candidate_sequences = [fasta.fetch(chrom, start, end) for chrom, start, end in all_regions]
    fasta.close()
    predictions = crested.tl.predict(input=candidate_sequences, model=model)
    dat_batch.append(predictions)

dat_batch_arr  = np.stack(dat_batch, axis=0)   # shape: (10, n_regions, n_outputs)
dat_batch_mean = np.mean(dat_batch_arr, axis=0) # shape: (n_regions, n_outputs)

np.savetxt(f"{work_path}/window_cluster_3K.pre.trf_replacement.txt", dat_batch_mean, fmt="%.3f")



###############################
### Step-3: Plotting the result

work_path = ""
web_path = "https://shendure-web.gs.washington.edu/content/members/cxqiu/public/backup/jax_atac/download"
source("help_code/utils.R")

model_id = "xxx"

celltype_list = c("Adipocyte_cells","Adipocyte_cells_Cyp2e1","B_cells","Brain_capillary_endothelial_cells","CNS_neurons","Cardiomyocytes","Corticofugal_neurons","Endocardial_cells","Endothelium","Epithelial_cells","Erythroid_cells","Eye","Glia","Glomerular_endothelial_cells","Gut_epithelial_cells","Hepatocytes","Intermediate_neuronal_progenitors","Kidney","Lateral_plate_and_intermediate_mesoderm","Liver_sinusoidal_endothelial_cells","Lung_and_airway","Lymphatic_vessel_endothelial_cells","Melanocyte_cells","Mesoderm","Neural_crest_PNS_neurons","Neuroectoderm_and_glia","Olfactory_ensheathing_cells","Olfactory_neurons","Oligodendrocytes","Skeletal_muscle_cells","T_cells","White_blood_cells")

obs = read.table(paste0(work_path, "/14_crested/", model_id, "/window_cluster/window_cluster_3K.obs.txt"), header=T, row.names=1, as.is=T)
obs = as.matrix(obs)

pre = read.table(paste0(work_path, "/14_crested/", model_id, "/window_cluster/window_cluster_3K.pre.trf_replacement.txt"), as.is=T)
pre = as.matrix(pre)

dat = read.table(paste0(work_path, "/14_crested/", model_id, "/window_cluster/window_cluster_3K.var.txt"), header=T, row.names=1, as.is=T)

keep = dat$split == "test"
obs = obs[keep,]
dat = dat[keep,]

obs = log1p(obs)
pre = log1p(pre)

df = data.frame(peak_id = rep(rownames(dat), ncol(obs)),
                obs = c(obs), pre = c(pre), 
                celltype_name = rep(colnames(obs), each = nrow(obs)))
df$celltype_name = factor(df$celltype_name, levels = celltype_list)

p = ggplot(data = df, aes(x=obs, y=pre)) +
    geom_hex(bins = 70) +
    geom_abline(intercept = 2, slope = 1, color = "red") +
    geom_abline(intercept = 3, slope = 1, color = "red") +
    geom_abline(intercept = -3, slope = 1, color = "red") +
    scale_x_continuous(limits = c(-0.3, 5.25), expand = expansion(mult = 0, add = 0)) +
    scale_y_continuous(limits = c(-0.3, 5.25), expand = expansion(mult = 0, add = 0)) +
#    labs(x="Log (observation + 1)", y="Log (prediction + 1)", title="") +
    labs(x="Log (observation + 1)", y="Log (prediction + 1)", title=paste0("Corr = ", round(cor(df$obs, df$pre), 2))) +
    scale_fill_continuous(type = "viridis", trans = "log10") +
    theme_classic(base_size = 10) +
    theme(legend.position="none") +
    theme(plot.title = element_text(hjust = 0.5)) +
    theme(axis.text.x = element_text(color="black"), axis.text.y = element_text(color="black")) 
ggsave(paste0("held_out_set.pdf"), p, height = 5, width = 5)

print(cor(df$obs, df$pre))




###############################################
### Step-4, based on the 100-bp windows on chr9


work_path = ""
web_path = "https://shendure-web.gs.washington.edu/content/members/cxqiu/public/backup/jax_atac/download"
source("help_code/utils.R")


mamm = "Mus_musculus"

chr = "chr9"

dat = read.table(paste0(work_path, "/prediction_", mamm, "/dat_", chr, "_obs.txt.gz"))
colnames(dat) = paste0("celltype_", c(1, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 2, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 3, 30, 31, 32, 33, 34, 35, 36, 4, 5, 6, 7, 8, 9))
peak_list = read.table(paste0(work_path, "/prediction_", mamm, "/dat_", chr, ".bed"))

celltype_convert = read.table(paste0(web_path, "/celltype_list_convert.txt"))
celltype_list = c("Adipocyte_cells","Adipocyte_cells_Cyp2e1","B_cells","Brain_capillary_endothelial_cells","CNS_neurons","Cardiomyocytes","Corticofugal_neurons","Endocardial_cells","Endothelium","Epithelial_cells","Erythroid_cells","Eye","Glia","Glomerular_endothelial_cells","Gut_epithelial_cells","Hepatocytes","Intermediate_neuronal_progenitors","Kidney","Lateral_plate_and_intermediate_mesoderm","Liver_sinusoidal_endothelial_cells","Lung_and_airway","Lymphatic_vessel_endothelial_cells","Melanocyte_cells","Mesoderm","Neural_crest_PNS_neurons","Neuroectoderm_and_glia","Olfactory_ensheathing_cells","Olfactory_neurons","Oligodendrocytes","Skeletal_muscle_cells","T_cells","White_blood_cells")

dat_obs = NULL
for(i in celltype_list){
    print(i)
    celltype_include = celltype_convert %>% filter(V1 == i) %>% pull(V3)
    if(length(celltype_include) > 1){
        dat_obs = cbind(dat_obs, apply(dat[,colnames(dat) %in% celltype_include], 1, max))
    } else {
        dat_obs = cbind(dat_obs, dat[,colnames(dat) == celltype_include])
    }
}

dat_pre = readRDS(paste0(work_path, "/prediction_", mamm, "_trf/dat.txt.rds"))
dat_loc = readRDS(paste0(work_path, "/prediction_", mamm, "_trf/dat_loc.txt.rds"))
keep = dat_loc$chr == chr
dat_loc = dat_loc[keep,]
dat_pre = dat_pre[keep,]


df = data.frame(pre = c(as.matrix(dat_pre)),
                obs = c(as.matrix(dat_obs)),
                chr = rep(dat_loc$chr, ncol(dat_pre)),
                start = rep(dat_loc$start, ncol(dat_pre)),
                end = rep(dat_loc$end, ncol(dat_pre)))
df$window_id = paste0(df$chr, "_", df$start, "_", df$end)
df$log_pre = log(df$pre + 1)
df$log_obs = log(df$obs + 1)
print(cor(df$log_pre, df$log_obs)) 

### chr9: 0.5505248

p = ggplot(data = df, aes(x=log_obs, y=log_pre)) +
    geom_hex(bins = 70) +
    geom_abline(intercept = 2, slope = 1, color = "red") +
    geom_abline(intercept = 3, slope = 1, color = "red") +
    geom_abline(intercept = -3, slope = 1, color = "red") +
    scale_x_continuous(limits = c(-0.3, 5.25), expand = expansion(mult = 0, add = 0)) +
    scale_y_continuous(limits = c(-0.3, 5.25), expand = expansion(mult = 0, add = 0)) +
#    labs(x="Log (observation + 1)", y="Log (prediction + 1)", title="") +
    labs(x="Log (observation + 1)", y="Log (prediction + 1)", title=paste0("Corr = ", round(cor(df$log_obs, df$log_pre), 2))) +
    scale_fill_continuous(type = "viridis", trans = "log10") +
    theme_classic(base_size = 10) +
#    theme(legend.position="none") +
    theme(plot.title = element_text(hjust = 0.5)) +
    theme(axis.text.x = element_text(color="black"), axis.text.y = element_text(color="black")) 
ggsave(paste0("chr9_windows_retrained.pdf"), p, height = 5, width = 5)



##############################################
### Step-5, Plotting TRF - Promoter enrichment


window_overlap_promoter = read.table(paste0(work_path, "/prediction_Mus_musculus/windows_overlap_with_TRF_promoters/chr9_windows_overlap_promoter.bed"))
colnames(window_overlap_promoter) = c("chr", "start", "end")
window_overlap_promoter$window_id = paste0(window_overlap_promoter$chr, "_", window_overlap_promoter$start, "_", window_overlap_promoter$end)

window_overlap_trf = read.table(paste0(work_path, "/prediction_Mus_musculus/windows_overlap_with_TRF_promoters/chr9_windows_overlap_trf.bed"))
colnames(window_overlap_trf) = c("chr", "start", "end")
window_overlap_trf$window_id = paste0(window_overlap_trf$chr, "_", window_overlap_trf$start, "_", window_overlap_trf$end)

df$overlap_promoter = df$window_id %in% window_overlap_promoter$window_id
df$overlap_trf = df$window_id %in% window_overlap_trf$window_id

df_long <- bind_rows(
  df %>% filter(overlap_promoter) %>% mutate(category = "promoter"),
  df %>% filter(overlap_trf) %>% mutate(category = "trf"),
  df %>% filter(!overlap_promoter & !overlap_trf) %>% mutate(category = "neither")
)

hex_df <- df_long %>%
  mutate(
    hx = cut(log_obs, 70),
    hy = cut(log_pre, 70)
  ) %>%
  group_by(hx, hy) %>%
  summarise(
    total = n(),
    promoter_n = sum(category == "promoter"),
    trf_n = sum(category == "trf"),
    neither_n = sum(category == "neither"),
    .groups = "drop"
  ) %>%
  mutate(
    hx_mid = (as.numeric(sub("\\((.+),(.+)\\]", "\\1", hx)) + as.numeric(sub("\\((.+),(.+)\\]", "\\2", hx))) / 2,
    hy_mid = (as.numeric(sub("\\((.+),(.+)\\]", "\\1", hy)) + as.numeric(sub("\\((.+),(.+)\\]", "\\2", hy))) / 2,
    score = (trf_n - promoter_n) / total
  )

p <- ggplot(hex_df, aes(x = hx_mid, y = hy_mid, color = score)) +
  geom_point(size = 1) +
  scale_color_gradientn(
    colors = Colors,
    limits = c(-1, 1),
    name = "TRF vs. Promoters"
  ) +
  geom_abline(intercept = 2, slope = 1, color = "red") +
  geom_abline(intercept = 3, slope = 1, color = "red") +
  geom_abline(intercept = -3, slope = 1, color = "red") +
  scale_x_continuous(limits = c(-0.3, 5.25), expand = expansion(mult = 0, add = 0)) +
  scale_y_continuous(limits = c(-0.3, 5.25), expand = expansion(mult = 0, add = 0)) +
  labs(
    x = "Log (observation + 1)",
    y = "Log (prediction + 1)",
    title = ""
  ) +
  theme_classic(base_size = 10) +
  theme(legend.position="none") +
  theme(
    plot.title = element_text(hjust = 0.5),
    axis.text.x = element_text(color="black"),
    axis.text.y = element_text(color="black")
  )

ggsave(paste0("chr9_windows_naive_overlap_promoter_trf.pdf"), p, height = 5, width = 5)










