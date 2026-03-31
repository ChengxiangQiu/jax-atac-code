


#############################################################
### Step-2, overlapping celltype-specific genes and enhancers

work_path=/net/shendure/vol2/projects/cxqiu/work/jax/atac_seq/novaseq/14_crested
model_id=mouse_fake_track_14
mkdir -p "$work_path"/"$model_id"/prediction_mammals/prediction_Mus_musculus_trf/compare_gene_exp
bedtools intersect \
-a "$work_path"/mouse_fake_track_12/compare_gene_exp/gene_loc.protein_coding.bed \
-b "$work_path"/"$model_id"/prediction_mammals/prediction_Mus_musculus_trf/window_list_10Mb.core_region.phred_score.bed \
-wa -wb \
> "$work_path"/"$model_id"/prediction_mammals/prediction_Mus_musculus_trf/compare_gene_exp/celltype_peaks_overlap_gene.bed



work_path = "/net/shendure/vol2/projects/cxqiu/work/jax/atac_seq/novaseq"
source("~/work/scripts/utils.R")

model_id = "mouse_fake_track_14"

dat = read.table(paste0(work_path, "/14_crested/", model_id, "/prediction_mammals/prediction_Mus_musculus_trf/compare_gene_exp/celltype_peaks_overlap_gene.bed"))
dat = dat[,c(1, 4, 5, 7, 8, 10, 11, 12, 13)]
colnames(dat) = c("chr", "strand", "gene_ID", "gene_short_name", "gene_TSS", "enhancer_start", "enhancer_end", "enhancer_score", "enhancer_celltype")

dat$enhancer_ID = paste0(dat$chr, "_", dat$enhancer_start, "_", dat$enhancer_end)
dat_1 = dat[dat$strand == '+',]
dat_1$distance = round((dat_1$enhancer_start + dat_1$enhancer_end)/2) - dat_1$gene_TSS
dat_2 = dat[dat$strand == '-',]
dat_2$distance = dat_2$gene_TSS - round((dat_2$enhancer_start + dat_2$enhancer_end)/2)
dat = rbind(dat_1, dat_2)
dat$abs_distance = abs(dat$distance)

dat_fc = readRDS(paste0(work_path, "/14_crested/mouse_fake_track_12/compare_gene_exp/gene_exp_log2fc_median.rds"))
dat_fc = data.frame(gene_ID = rep(rownames(dat_fc), ncol(dat_fc)),
                    log2fc = c(as.matrix(dat_fc)),
                    gene_celltype = rep(colnames(dat_fc), each = nrow(dat_fc)))

dat_median = dat %>% left_join(dat_fc, by = "gene_ID", relationship = "many-to-many")
saveRDS(dat_median, paste0(work_path, "/14_crested/", model_id, "/prediction_mammals/prediction_Mus_musculus_trf/compare_gene_exp/dat_median.rds"))



#############################################################
### Step-3, Plotting the mean/median log2FC across cell types

work_path = "/net/shendure/vol2/projects/cxqiu/work/jax/atac_seq/novaseq"
source("~/work/scripts/utils.R")

model_id = "mouse_fake_track_14"

celltype_list=c("Adipocyte_cells","Adipocyte_cells_Cyp2e1","B_cells","Brain_capillary_endothelial_cells","CNS_neurons","Cardiomyocytes","Corticofugal_neurons","Endocardial_cells","Endothelium","Epithelial_cells","Erythroid_cells","Eye","Glia","Glomerular_endothelial_cells","Gut_epithelial_cells","Hepatocytes","Intermediate_neuronal_progenitors","Kidney","Lateral_plate_and_intermediate_mesoderm","Liver_sinusoidal_endothelial_cells","Lung_and_airway","Lymphatic_vessel_endothelial_cells","Melanocyte_cells","Mesoderm","Neural_crest_PNS_neurons","Neuroectoderm_and_glia","Olfactory_ensheathing_cells","Olfactory_neurons","Oligodendrocytes","Skeletal_muscle_cells","T_cells","White_blood_cells")

dat_all = readRDS(paste0(work_path, "/14_crested/", model_id, "/prediction_mammals/prediction_Mus_musculus_trf/compare_gene_exp/dat_median.rds"))

dat_all <- dat_all %>%
  mutate(
    bin_label = cut(
      enhancer_score,
      breaks = c(0, 30, 32, 34, 36, 38, 40, 50, Inf),
      right = FALSE,
      labels = c(
        "<30", "[30, 32)", "[32, 34)", "[34, 36)",
        "[36, 38)", "[38, 40)", "[40, 50)", ">=50"
      )
    ),
    bin_label_big = cut(
      enhancer_score,
      breaks = c(0, 30, 40, Inf),
      right = FALSE,
      labels = c(
        "<30", "[30, 40)", ">=40"
      )
    )
  )

dat = dat_all[dat_all$enhancer_celltype == dat_all$gene_celltype,
                c("bin_label", "enhancer_celltype", "abs_distance", "log2fc")]

window_list = seq(0, 96, 2)
dat_x = NULL
for(i in 1:length(window_list)){
    print(i)
    df = dat %>% filter(abs_distance >= window_list[i]*1000, abs_distance < (window_list[i]*1000+5000)) %>%
        group_by(bin_label) %>%
        summarize(mean_log2fc = mean(log2fc), .groups = "drop") %>%
        mutate(distance = paste0(window_list[i], "-", window_list[i]+5))
    dat_x = rbind(dat_x, df)
}

window_list = data.frame(a = c(0, 25, 50), b = c(25, 50, 101))
dat_y = NULL
for(i in 1:nrow(window_list)){
    print(i)
    df = dat_all %>% filter(abs_distance >= window_list$a[i]*1000, abs_distance < window_list$b[i]*1000) %>%
        group_by(bin_label_big, enhancer_celltype, gene_celltype) %>%
        summarize(mean_log2fc = mean(log2fc), .groups = "drop") %>%
        mutate(distance = paste0(window_list$a[i], "-", window_list$b[i]))
    dat_y = rbind(dat_y, df)
}

saveRDS(list(dat_x, dat_y), paste0(work_path, "/14_crested/", model_id, "/prediction_mammals/prediction_Mus_musculus_trf/compare_gene_exp/result_median.rds"))



#########################################################
### Step-4, Permutation on peak - cell type relationships

work_path = "/net/shendure/vol2/projects/cxqiu/work/jax/atac_seq/novaseq"
source("~/work/scripts/utils.R")

model_id = "mouse_fake_track_14"

args = commandArgs(trailingOnly=TRUE)
perm = as.numeric(args[1])
print(perm)

dat_all = readRDS(paste0(work_path, "/14_crested/", model_id, "/prediction_mammals/prediction_Mus_musculus_trf/compare_gene_exp/dat_median.rds"))

dat_all <- dat_all %>%
  mutate(
    bin_label = cut(
      enhancer_score,
      breaks = c(0, 30, 32, 34, 36, 38, 40, 50, Inf),
      right = FALSE,
      labels = c(
        "<30", "[30, 32)", "[32, 34)", "[34, 36)",
        "[36, 38)", "[38, 40)", "[40, 50)", ">=50"
      )
    ),
    bin_label_big = cut(
      enhancer_score,
      breaks = c(0, 30, 40, Inf),
      right = FALSE,
      labels = c(
        "<30", "[30, 40)", ">=40"
      )
    )
  )

enhancer = unique(dat_all[,c("enhancer_ID", "enhancer_celltype")])
enhancer_x = enhancer
enhancer_x$enhancer_celltype = enhancer$enhancer_celltype[sample(1:nrow(enhancer))]

dat_all$enhancer_celltype = NULL
dat_all = dat_all %>% left_join(enhancer_x, by = "enhancer_ID", relationship = "many-to-many")

dat = dat_all[dat_all$enhancer_celltype == dat_all$gene_celltype,
                c("bin_label", "enhancer_celltype", "abs_distance", "log2fc")]

window_list = seq(0, 96, 2)
dat_x = NULL
for(i in 1:length(window_list)){
    print(i)
    df = dat %>% filter(abs_distance >= window_list[i]*1000, abs_distance < (window_list[i]*1000+5000)) %>%
        group_by(bin_label) %>%
        summarize(mean_log2fc = mean(log2fc), .groups = "drop") %>%
        mutate(distance = paste0(window_list[i], "-", window_list[i]+5))
    dat_x = rbind(dat_x, df)
}
dat_x$perm = perm

window_list = data.frame(a = c(0, 25, 50), b = c(25, 50, 101))
dat_y = NULL
for(i in 1:nrow(window_list)){
    print(i)
    df = dat_all %>% filter(abs_distance >= window_list$a[i]*1000, abs_distance < window_list$b[i]*1000) %>%
        group_by(bin_label_big, enhancer_celltype, gene_celltype) %>%
        summarize(mean_log2fc = mean(log2fc), .groups = "drop") %>%
        mutate(distance = paste0(window_list$a[i], "-", window_list$b[i]))
    dat_y = rbind(dat_y, df)
}
dat_y$perm = perm

saveRDS(list(dat_x, dat_y), paste0(work_path, "/14_crested/", model_id, "/prediction_mammals/prediction_Mus_musculus_trf/compare_gene_exp/permutations/result_perm_", perm, ".rds"))




#############################################
### Step-5: Plotting the results (adding permutation)

work_path = "/net/shendure/vol2/projects/cxqiu/work/jax/atac_seq/novaseq"
source("~/work/scripts/utils.R")

model_id = "mouse_fake_track_14"

celltype_list=c("Adipocyte_cells","Adipocyte_cells_Cyp2e1","B_cells","Brain_capillary_endothelial_cells","CNS_neurons","Cardiomyocytes","Corticofugal_neurons","Endocardial_cells","Endothelium","Epithelial_cells","Erythroid_cells","Eye","Glia","Glomerular_endothelial_cells","Gut_epithelial_cells","Hepatocytes","Intermediate_neuronal_progenitors","Kidney","Lateral_plate_and_intermediate_mesoderm","Liver_sinusoidal_endothelial_cells","Lung_and_airway","Lymphatic_vessel_endothelial_cells","Melanocyte_cells","Mesoderm","Neural_crest_PNS_neurons","Neuroectoderm_and_glia","Olfactory_ensheathing_cells","Olfactory_neurons","Oligodendrocytes","Skeletal_muscle_cells","T_cells","White_blood_cells")


### merging all cell types

dat_perm = NULL
for(perm in 1:100){
    print(perm)
    tmp = readRDS(paste0(work_path, "/14_crested/", model_id, "/prediction_mammals/prediction_Mus_musculus_trf/compare_gene_exp/permutations/result_perm_", perm, ".rds"))
    dat_perm = rbind(dat_perm, tmp[[1]])
}

dat_y = dat_perm %>% group_by(bin_label, distance) %>%
    summarize(mean_perm = mean(mean_log2fc))

dat_x = readRDS(paste0(work_path, "/14_crested/", model_id, "/prediction_mammals/prediction_Mus_musculus_trf/compare_gene_exp/result_median.rds"))

dat_x = dat_x[[1]] %>% left_join(dat_y, by = c("bin_label", "distance"))

dat_x$delta_exp = dat_x$mean_log2fc - dat_x$mean_perm

window_list = seq(0, 96, 2)
window_list_x = NULL
for(i in 1:length(window_list)){
    window_list_x = c(window_list_x, paste0(window_list[i], "-", window_list[i]+5))
}

dat_x$distance = factor(dat_x$distance, levels = window_list_x)

p <- ggplot(dat_x, aes(x = distance, y = bin_label, color = delta_exp)) +
  geom_point(size = 4) +
  scale_color_gradientn(
    colors = Colors,
    limits = c(min(dat_x$delta_exp), max(dat_x$delta_exp)),
    name = "Delta Gene Exp"
  ) +
  labs(
    x = "Distance from gene TSS to peak (kb)",
    y = "Predicted accessibility",
    title = ""
  ) +
  theme_classic(base_size = 10) +
#  theme(legend.position="none") +
  theme(
    plot.title = element_text(hjust = 0.5),
    axis.text.x = element_text(color="black", angle = 45, hjust = 1),
    axis.text.y = element_text(color="black")
  )

ggsave("~/share/Access_Exp_scatter.protein_coding.pdf", p, height = 5, width = 9)


distance_value = data.frame(distance = names(table(dat_x$distance)), distance_index = 1:length(table(dat_x$distance)))
dat_x = dat_x %>% left_join(distance_value, by = "distance")

bin_label_color = c("#440154", "#46327e", "#365c8d", "#277f8e", "#1fa187", "#4ac16d", "#a0da39", "#e65d2f")
names(bin_label_color) = names(table(dat_x$bin_label))

p = ggplot() +
    geom_point(data = dat_x, aes(x = distance_index, y = delta_exp, color = bin_label)) +
    geom_smooth(data = dat_x, aes(x = distance_index, y = delta_exp, group = bin_label, color = bin_label), method = 'loess', formula = 'y ~ x', se = FALSE) +
    labs(
    x = "Distance from gene TSS to peak (kb)",
    y = "Gene expression",
    title = ""
    ) +
  theme_classic(base_size = 10) +
  scale_color_manual(values=bin_label_color) +
  theme(legend.position="none") +
  theme(
    plot.title = element_text(hjust = 0.5),
    axis.text.x = element_text(color="black", angle = 45, hjust = 1),
    axis.text.y = element_text(color="black")
  ) + scale_x_continuous(breaks = 1:49)

ggsave("~/share/Access_Exp_scatter.protein_coding_x.pdf", p, height = 5, width = 6)


########## split by cell types

dat_perm = NULL
for(perm in 1:100){
    print(perm)
    tmp = readRDS(paste0(work_path, "/14_crested/", model_id, "/prediction_mammals/prediction_Mus_musculus_trf/compare_gene_exp/permutations/result_perm_", perm, ".rds"))
    dat_perm = rbind(dat_perm, tmp[[2]])
}

dat_y = dat_perm %>% group_by(bin_label_big, distance, enhancer_celltype, gene_celltype) %>%
    summarize(mean_perm = mean(mean_log2fc))

dat_x = readRDS(paste0(work_path, "/14_crested/", model_id, "/prediction_mammals/prediction_Mus_musculus_trf/compare_gene_exp/result_median.rds"))

dat_x = dat_x[[2]] %>% left_join(dat_y, by = c("bin_label_big", "distance", "enhancer_celltype", "gene_celltype"))

dat_x$delta_exp = dat_x$mean_log2fc - dat_x$mean_perm

dat_x_scaled <- dat_x %>% filter(enhancer_celltype == gene_celltype) %>% group_by(enhancer_celltype) %>%
  mutate(delta_exp_scaled = scales::rescale(delta_exp))

dat_x_scaled$enhancer_celltype = factor(dat_x_scaled$enhancer_celltype, levels = celltype_list)
dat_x_scaled$bin_label_big = factor(dat_x_scaled$bin_label_big, levels = c("<30", "[30, 40)", ">=40"))
dat_x_scaled$distance = factor(dat_x_scaled$distance, levels = c("0-25", "25-50", "50-101"))

#dat_x_scaled = dat_x %>% filter(enhancer_celltype == gene_celltype)
#dat_x_scaled$delta_exp_scaled = dat_x_scaled$delta_exp

p <- ggplot(dat_x_scaled, aes(x = distance, y = bin_label_big, fill = delta_exp_scaled)) +
  geom_tile() +
  scale_fill_gradientn(
    colors = Colors,
    limits = c(min(dat_x_scaled$delta_exp_scaled), max(dat_x_scaled$delta_exp_scaled)),
    name = "Delta Gene Exp"
  ) +
  labs(
    x = "Distance from gene TSS to peak (kb)",
    y = "Predicted accessibility",
    title = ""
  ) +
  theme_classic(base_size = 10) +
  theme(legend.position="none") +
  theme(
    plot.title = element_text(hjust = 0.5),
    axis.text.x = element_text(color="black", angle = 45, hjust = 1),
    axis.text.y = element_text(color="black")
  ) +
  facet_wrap(~ enhancer_celltype, ncol = 8, nrow = 4)

ggsave("~/share/Access_Exp_scatter.protein_coding.split_celltype_x.pdf", p, height = 6, width = 10)





### spliting cell types (32 x 32)

celltype_list_reorder = c("Erythroid_cells","T_cells","B_cells","White_blood_cells","Corticofugal_neurons","Olfactory_neurons","Adipocyte_cells_Cyp2e1","Adipocyte_cells","Hepatocytes","Gut_epithelial_cells","Cardiomyocytes","Lung_and_airway","Endocardial_cells","Endothelium","Lymphatic_vessel_endothelial_cells","Liver_sinusoidal_endothelial_cells","Glomerular_endothelial_cells","Brain_capillary_endothelial_cells","Oligodendrocytes","Glia","Melanocyte_cells","Olfactory_ensheathing_cells","Intermediate_neuronal_progenitors","Eye","Neural_crest_PNS_neurons","CNS_neurons","Neuroectoderm_and_glia","Skeletal_muscle_cells","Kidney","Epithelial_cells","Lateral_plate_and_intermediate_mesoderm","Mesoderm")

dat_perm = NULL
for(perm in 1:100){
    print(perm)
    tmp = readRDS(paste0(work_path, "/14_crested/", model_id, "/prediction_mammals/prediction_Mus_musculus_trf/compare_gene_exp/permutations/result_perm_", perm, ".rds"))
    dat_perm = rbind(dat_perm, tmp[[2]])
}

dat_y = dat_perm %>% group_by(bin_label_big, distance, enhancer_celltype, gene_celltype) %>%
    summarize(mean_perm = mean(mean_log2fc))

dat_x = readRDS(paste0(work_path, "/14_crested/", model_id, "/prediction_mammals/prediction_Mus_musculus_trf/compare_gene_exp/result_median.rds"))

dat_x = dat_x[[2]] %>% left_join(dat_y, by = c("bin_label_big", "distance", "enhancer_celltype", "gene_celltype"))

dat_x$delta_exp = dat_x$mean_log2fc - dat_x$mean_perm

dat_x$x_axis = paste0(dat_x$gene_celltype, "_", dat_x$distance)
dat_x$y_axis = paste0(dat_x$enhancer_celltype, "_", dat_x$bin_label_big)
dat_x$x_axis = factor(dat_x$x_axis, levels = paste0(rep(celltype_list_reorder, each = 3), "_", rep(c("0-25", "25-50", "50-101"), 32)))
dat_x$y_axis = factor(dat_x$y_axis, levels = paste0(rep(rev(celltype_list_reorder), each = 3), "_", rep(c("<30", "[30, 40)", ">=40"), 32)))


p <- ggplot(dat_x, aes(x = x_axis, y = y_axis, color = delta_exp)) +
  geom_point(size = 2) +
  scale_color_gradientn(
    colors = c("grey95", "red"),
    limits = c(min(dat_x$delta_exp), max(dat_x$delta_exp)),
    name = "Delta Gene Exp"
  ) +
  labs(
    x = "Distance from gene TSS to peak (kb)",
    y = "Predicted accessibility",
    title = ""
  ) +
  theme_classic(base_size = 10) +
#  theme(legend.position="none") +
  theme(
    plot.title = element_text(hjust = 0.5),
    axis.text.x = element_text(color="black", angle = 45, hjust = 1),
    axis.text.y = element_text(color="black")
  )

ggsave("~/share/Access_Exp_scatter.protein_coding.split_celltype.pdf", p, height = 12, width = 12)



p <- ggplot(dat_x %>% filter(distance == "0-25"), aes(x = x_axis, y = y_axis, color = delta_exp)) +
  geom_point(size = 1) +
  scale_color_gradientn(
    colors = c("grey95", "red"),
    limits = c(min(dat_x$delta_exp), max(dat_x$delta_exp)),
    name = "Delta Gene Exp"
  ) +
  labs(
    x = "Distance from gene TSS to peak (kb)",
    y = "Predicted accessibility",
    title = ""
  ) +
  theme_classic(base_size = 10) +
#  theme(legend.position="none") +
  theme(
    plot.title = element_text(hjust = 0.5),
    axis.text.x = element_text(color="black", angle = 45, hjust = 1),
    axis.text.y = element_text(color="black")
  )

ggsave("~/share/Access_Exp_scatter.protein_coding.split_celltype_subset.pdf", p, height = 10, width = 7)










