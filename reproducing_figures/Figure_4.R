
#################################################################################################################
### Figure 4. Evolutionary coherence separates distal enhancer grammar from promoter & repeat-associated signals.

####################################################################################################################################################################
### Fig. 4b: For 2.79M x 100-bp mouse genomic windows predicted to be highly and specifically accessible, we performed liftover to the genomes of 240 other mammals.

web_path = "https://shendure-web.gs.washington.edu/content/members/cxqiu/public/backup/jax_atac/download/reproducing_figures"

library(ggplot2)
library(dplyr)

df = read.csv(paste0(web_path, "/Fig4b_num_species_liftover.csv"))

p = ggplot() +
    geom_histogram(data = df, aes(count), bins = 50) +
    theme_classic(base_size = 10) +
    geom_vline(xintercept = 120) +
    labs(x = "# of mammals were successfully lifted over", y = "# of candidate windows (100 bp)")

ggsave("Fig4b_num_species_liftover.pdf", p, width = 6, height = 4.5)


#################################################################################################################################################################################################################################################################
### Fig. 4c: For retained windows, we applied the mouse-based, evolution-naive CREsted model to syntenic regions of other mammals, and calculated the pairwise Pearson's r between the resulting 36-value predicted accessibility vectors for mouse vs. ortholog.

web_path = "https://shendure-web.gs.washington.edu/content/members/cxqiu/public/backup/jax_atac/download/reproducing_figures"

library(ggplot2)
library(dplyr)

df = read.csv(paste0(web_path, "/Fig4c_median_corr_across_species.csv"))

p = ggplot() +
    geom_histogram(data = df, aes(x = median_corr), bins = 100) +
    labs(x="Median correlation coefficient across mammals", y="# of candidate windows (100 bp)") +
    theme_classic(base_size = 10) +
    theme(plot.title = element_text(hjust = 0.5)) +
    theme(legend.position="none") +
    theme(axis.text.x = element_text(color="black"), axis.text.y = element_text(color="black")) 

ggsave("Fig4c_median_corr_across_species.pdf", p, width = 10, height = 5)


####################################################################################################################################
### Fig. 4d: The predicted score matrix for the remaining 547,317 windows x 36 Level-2 cell classes was subjected to UMAP embedding.

web_path = "https://shendure-web.gs.washington.edu/content/members/cxqiu/public/backup/jax_atac/download/reproducing_figures"

library(ggplot2)
library(dplyr)

window_cluster_color_plate = c("Adipocyte_cells_Cyp2e1" = "#7f3e39",
                               "Adipocyte_cells" = "#5e7fbf",
                               "Brain_capillary_endothelial_cells" = "#b25525",
                               "Endocardial_cells" = "#c6a851",
                               "Endothelium" = "#00a34e",
                               "Glomerular_endothelial_cells" = "#7a66d6",
                               "Liver_sinusoidal_endothelial_cells" = "#748132",
                               "Lymphatic_vessel_endothelial_cells" = "#ca47a7",
                               "Epithelial_cells" = "#af9fb6",
                               "Gut_epithelial_cells" = "#ff007a",
                               "Kidney" = "#65c17d",
                               "Lung_and_airway" = "#02b0d1",
                               "Erythroid_cells" = "#dc453e",
                               "Hepatocytes" = "#185700",
                               "Lateral_plate_and_intermediate_mesoderm" = "#d78860",
                               "Mesoderm" = "#bb46c5",
                               "Cardiomyocytes" = "#643e8c",
                               "Skeletal_muscle_cells" = "#ffa1f5",
                               "Promoters" = "grey80",
                               "Melanocyte_cells" = "#88b439",
                               "Glia" = "#fff167",
                               "Olfactory_ensheathing_cells" = "#b46aa2",
                               "CNS_neurons" = "#e5c000",
                               "Eye" = "#00d450",
                               "Intermediate_neuronal_progenitors" = "#2e0ab7",
                               "Neural_crest_PNS_neurons" = "#b5ce92",
                               "Neuroectoderm_and_glia" = "#f96100",
                               "Corticofugal_neurons" = "#8d6a32",
                               "Olfactory_neurons" = "#e6230b",
                               "Oligodendrocytes" = "#916e00",
                               "B_cells" = "#01b7a6",
                               "T_cells" = "#ff9d47",
                               "White_blood_cells" = "#7ca0ff")

pd = read.csv(paste0(web_path, "/Fig4d_umap_windows.csv"))

p1 = ggplot() +
    geom_point(data = pd, aes(UMAP_1, UMAP_2, color = promoter_enrich), size = 0.03) +
    labs(x="UMAP 1", y="UMAP 2", title="Promoter enrichment (547,317 windows)") +
    theme_classic(base_size = 10) +
    theme(plot.title = element_text(hjust = 0.5)) +
    theme(legend.position="none") +
    theme(axis.text.x = element_text(color="black"), axis.text.y = element_text(color="black")) +
    scale_color_viridis_c(option = "viridis")

p2 = ggplot() +
    geom_point(data = pd, aes(UMAP_1, UMAP_2, color = cell_class) , size = 0.03) +
    labs(x="UMAP 1", y="UMAP 2", title="Cell classes (547,317 windows)") +
    theme_classic(base_size = 10) +
    theme(plot.title = element_text(hjust = 0.5)) +
    theme(legend.position="none") +
    theme(axis.text.x = element_text(color="black"), axis.text.y = element_text(color="black")) +
    scale_color_manual(values=window_cluster_color_plate)

ggsave("Fig4d_umap_windows.png", p1 + p2, height = 5, width = 10, dpi = 300)


#####################################################################################################################################################################################################################################################################################################################################
### Fig. 4e: We randomly sampled 500 windows from each of the 32 non-promoter clusters and computed Pearson's r between their 36-value predicted accessibility vectors vs. those of windows from the same cell class cluster (red), different cell class cluster (green), or orthologous regions of other mammalian genomes (purple).

web_path = "https://shendure-web.gs.washington.edu/content/members/cxqiu/public/backup/jax_atac/download/reproducing_figures"

library(ggplot2)
library(dplyr)

category = data.frame(category_ID = c(1,2,3),
                      category = c("Diff_cluster", "Diff_species", "Same_cluster"))

pd = read.csv(paste0(web_path, "/Fig4e_corr_between_windows.csv"))

pd = pd %>% left_join(category, by = "category_ID")

p = ggplot(pd, aes(x = corr, fill = category)) +
    geom_density(alpha = 0.5) +
    scale_fill_manual(values = c("Diff_cluster" = "#cb6751", "Same_cluster" = "#7aa457", "Diff_species" = "#9e6ebd")) +
    labs(x="Correlation coefficient between windows", y="Density", title="") +
    theme_classic(base_size = 10) +
    theme(plot.title = element_text(hjust = 0.5)) +
    theme(legend.position="none") +
    theme(axis.text.x = element_text(color="black"), axis.text.y = element_text(color="black"))

ggsave("Fig4e_corr_between_windows.pdf", p, height = 5, width = 5)


###########################################################################################################################
### Fig. 4f: Scatterplots of observed vs. predicted chromatin accessibility across 8,223 peaks held out during fine-tuning.

web_path = "https://shendure-web.gs.washington.edu/content/members/cxqiu/public/backup/jax_atac/download/reproducing_figures"

library(ggplot2)
library(dplyr)

pd = read.csv(paste0(web_path, "/Fig4f_held_out_aware.csv"))

p = ggplot(data = pd, aes(x=obs, y=pre)) +
    geom_hex(bins = 70) +
    geom_abline(intercept = 2, slope = 1, color = "red") +
    geom_abline(intercept = 3, slope = 1, color = "red") +
    geom_abline(intercept = -3, slope = 1, color = "red") +
    scale_x_continuous(limits = c(-0.3, 5.25), expand = expansion(mult = 0, add = 0)) +
    scale_y_continuous(limits = c(-0.3, 5.25), expand = expansion(mult = 0, add = 0)) +
    labs(x="Log (observation + 1)", y="Log (prediction + 1)", title="") +
    scale_fill_continuous(type = "viridis", trans = "log10") +
    theme_classic(base_size = 10) +
    theme(plot.title = element_text(hjust = 0.5)) +
    theme(axis.text.x = element_text(color="black"), axis.text.y = element_text(color="black")) 

ggsave("Fig4f_held_out_aware.pdf", p, height = 5, width = 5)


######################################################################################################################################
### Fig. 4g: Scatterplots of observed vs. predicted chromatin accessibility across 1.25M x 100 bp windows spanning mouse chromosome 9.

web_path = "https://shendure-web.gs.washington.edu/content/members/cxqiu/public/backup/jax_atac/download/reproducing_figures"

library(ggplot2)
library(dplyr)

pd = read.csv(paste0(web_path, "/Fig4g_pred_chr9_aware.csv"))

p = ggplot(data = pd, aes(x=log_obs, y=log_pre)) +
    geom_hex(bins = 70) +
    geom_abline(intercept = 2, slope = 1, color = "red") +
    geom_abline(intercept = 3, slope = 1, color = "red") +
    geom_abline(intercept = -3, slope = 1, color = "red") +
    scale_x_continuous(limits = c(-0.3, 5.25), expand = expansion(mult = 0, add = 0)) +
    scale_y_continuous(limits = c(-0.3, 5.25), expand = expansion(mult = 0, add = 0)) +
    labs(x="Log (observation + 1)", y="Log (prediction + 1)", title="") +
    scale_fill_continuous(type = "viridis", trans = "log10") +
    theme_classic(base_size = 10) +
    theme(plot.title = element_text(hjust = 0.5)) +
    theme(axis.text.x = element_text(color="black"), axis.text.y = element_text(color="black")) 

ggsave("Fig4g_pred_chr9_aware.pdf", p, height = 5, width = 5)


########################################################################################################################################################################################################################################################################################
### Fig. 4h: A plot similar to panel g, but instead of coloring each bin by the counts of 100 bp windows, it is colored by the difference between the percentage of windows overlapping tandem repeats and the percentage overlapping promoters (-1,500 bp to +500 bp relative to TSSs).

web_path = "https://shendure-web.gs.washington.edu/content/members/cxqiu/public/backup/jax_atac/download/reproducing_figures"

library(gplots)
library(reshape2)
library(viridis)
library(RColorBrewer)
Colors=rev(brewer.pal(11,"Spectral"))
Colors=colorRampPalette(Colors)(120)

dat = read.csv(paste0(web_path, "/Fig4h_pred_chr9_overlap_promoter_TRF_aware.csv"))

p = ggplot(dat, aes(x = hx_mid, y = hy_mid, color = score)) +
    geom_point(size = 1) +
    scale_color_gradientn(
        colors = Colors,
        limits = c(-1, 1),
        name = "TRF vs. Promoters"
    ) +
    geom_abline(intercept = 2, slope = 1, linetype = "dashed", color = "red") +
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

ggsave("Fig4h_pred_chr9_overlap_promoter_TRF_aware.pdf", p, height = 5, width = 5)












