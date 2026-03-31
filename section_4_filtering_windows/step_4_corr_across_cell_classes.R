
###############################################################
### Correlation coefficient between windows across three groups

### Support data can be downloaded from:
### https://shendure-web.gs.washington.edu/content/members/cxqiu/public/backup/jax_atac/download/

### Please contact Chengxiang (CX) Qiu for any questions!
### cxqiu@uw.edu or chengxiang.qiu@dartmouth.edu


###################################################
### Step-1: save dataset with the subset of windows

work_path = ""
web_path = "https://shendure-web.gs.washington.edu/content/members/cxqiu/public/backup/jax_atac/download"
source("help_code/utils.R")

data = read.table(paste0(work_path, "/14_crested/celltype_L2_cut_norm/prediction_mammals/umap_Mus_musculus/dat.txt"), sep = "\t")
window_list_uniq = readRDS(paste0(work_path, "/14_crested/celltype_L2_cut_norm/prediction_mammals/corr_Mus_musculus/window_list_uniq.rds"))
rownames(data) = window_list_uniq$window_id = paste0(window_list_uniq$chr, "_", window_list_uniq$start, "_", window_list_uniq$end)

candidate_windows = read.table(paste0(work_path, "/14_crested/mouse_fake_track_14/candidate_region_exclude_blacklist.bed"))
colnames(candidate_windows) = c("chr", "start", "end", "celltype")
candidate_windows$window_id = paste0(candidate_windows$chr, "_", candidate_windows$start, "_", candidate_windows$end)

set.seed(1234)
candidate_windows_sub = candidate_windows %>% group_by(celltype) %>% slice_sample(n = 500)
dat_sub = data[as.vector(candidate_windows_sub$window_id),]

mat_cor = cor(t(as.matrix(dat_sub)))
mat_cor[lower.tri(mat_cor, diag = TRUE)] = NA
df_cor = melt(mat_cor, na.rm = TRUE)

celltype_lookup <- candidate_windows_sub$celltype
names(celltype_lookup) <- candidate_windows_sub$window_id

df_cor <- df_cor %>%
  mutate(
    celltype1 = celltype_lookup[Var1],
    celltype2 = celltype_lookup[Var2],
    category = if_else(celltype1 == celltype2, "Same_cluster", "Diff_cluster")
  )

df_cor %>%
  group_by(category) %>%
  summarize(
    mean_value = mean(value),
    median_value = median(value),
    .groups = "drop"
  )

### adding correlation with other species
result_corr = readRDS(paste0(work_path, "/14_crested/celltype_L2_cut_norm/prediction_mammals/corr_Mus_musculus/prediction/result_corr.rds"))
window_ID_include = window_list_uniq %>% filter(window_id %in% as.vector(candidate_windows_sub$window_id)) %>% pull(window_ID)
result_corr = result_corr[result_corr$window_ID %in% window_ID_include,]

df = rbind(df_cor %>% select(corr = value, category), result_corr %>% select(corr) %>% mutate(category = "Diff_species"))

p = ggplot(df, aes(x = corr, fill = category)) +
    geom_density(alpha = 0.5) +
    scale_fill_manual(values = c("Diff_cluster" = "#cb6751", "Same_cluster" = "#7aa457", "Diff_species" = "#9e6ebd")) +
    labs(x="Corr", y="Density", title="") +
    theme_classic(base_size = 10) +
    theme(plot.title = element_text(hjust = 0.5)) +
    theme(legend.position="none") +
    theme(axis.text.x = element_text(color="black"), axis.text.y = element_text(color="black"))

ggsave("corr_across_celltypes.pdf", p, height = 5, width = 5)

# Diff_cluster: red
# Same_cluster: green
# Diff_species: purple






