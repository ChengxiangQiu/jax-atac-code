---
title: Figure 3
kernelspec:
  name: ir
  display_name: R
---

Evolution-naive modeling of chromatin accessibility across 36 Level-2 cell classes.


+++

## Fig. 3b

Scatterplots of observed (normalized Tn5 cut-site counts; natural log-scaled; x-axis) vs. predicted (evolution-naive CREsted model; natural log-scaled; y-axis) chromatin accessibility across across 9,583 held-out peaks for all 36 Level-2 cell classes.


```{code-cell #Fig3b_held_out} r
web_path = "https://shendure-web.gs.washington.edu/content/members/cxqiu/public/backup/jax_atac/download/reproducing_figures"

library(ggridges)
library(ggplot2)
library(dplyr)

df = read.csv(paste0(web_path, "/Fig3b_held_out.csv"))

p = ggplot(data = df, aes(x=obs, y=pre)) +
    geom_hex(bins = 70) +
    geom_abline(intercept = 2, slope = 1, linetype = "dashed", color = "red") +
    geom_abline(intercept = 3, slope = 1, color = "red") +
    geom_abline(intercept = -3, slope = 1, color = "red") +
    scale_x_continuous(limits = c(-0.3, 5.25), expand = expansion(mult = 0, add = 0)) +
    scale_y_continuous(limits = c(-0.3, 5.25), expand = expansion(mult = 0, add = 0)) +
    labs(x="Log (observation + 1)", y="Log (prediction + 1)", title=paste0("Corr = ", round(cor(df$obs, df$pre), 2))) +
    scale_fill_continuous(type = "viridis", trans = "log10") +
    theme_classic(base_size = 10) +
    theme(legend.position="none") +
    theme(plot.title = element_text(hjust = 0.5)) +
    theme(axis.text.x = element_text(color="black"), axis.text.y = element_text(color="black"))

ggsave("Fig3b_held_out.pdf", p, height = 5, width = 5)
p
```

+++

## Fig. 3c

Pairwise Pearson correlations were computed across 36 Level-2 cell classes between normalized Tn5 cut-site counts from scATAC-seq data (rows) and CREsted model predictions (columns) at 9,583 held-out peaks.


```{code-cell #Fig3c_pairwise_pearson_corr} r
web_path = "https://shendure-web.gs.washington.edu/content/members/cxqiu/public/backup/jax_atac/download/reproducing_figures"

library(gplots)
library(reshape2)
library(viridis)

library(RColorBrewer)
Colors=rev(brewer.pal(11,"Spectral"))
Colors=colorRampPalette(Colors)(120)

dat = read.csv(paste0(web_path, "/Fig3c_pairwise_pearson_corr.csv"), row.names=1)

# pdf("Fig3c_pairwise_pearson_corr.pdf")
heatmap.2(as.matrix(dat),
          col=Colors,
          scale="none",
          Rowv = F,
          Colv = F,
          key=T,
          density.info="none",
          trace="none",
          cexRow = 0.5,
          cexCol = 0.5,
          margins = c(5,5))
# dev.off()
```

+++

## Fig. 3d

Observed (top) and predicted (bottom) chromatin accessibility is plotted for each of the 36 Level-2 cell classes, for three representative regions from the test set.


```{code-cell} r
web_path = "https://shendure-web.gs.washington.edu/content/members/cxqiu/public/backup/jax_atac/download/reproducing_figures"

library(ggridges)
library(ggplot2)
library(dplyr)
library(patchwork)

celltype_L2_color_plate = c("Adipocyte_cells_Cyp2e1_L2" = "#7f3e39",
                            "Adipocyte_cells_L2" = "#5e7fbf",
                            "Brain_capillary_endothelial_cells_L2" = "#b25525",
                            "Endocardial_cells_L2" = "#c6a851",
                            "Endothelium_L2" = "#00a34e",
                            "Glomerular_endothelial_cells_L2" = "#7a66d6",
                            "Liver_sinusoidal_endothelial_cells_L2" = "#748132",
                            "Lymphatic_vessel_endothelial_cells_L2" = "#ca47a7",
                            "Epithelial_cells_L2" = "#af9fb6",
                            "Gut_epithelial_cells_L2" = "#ff007a",
                            "Kidney_L2" = "#65c17d",
                            "Lung_and_airway_L2" = "#02b0d1",
                            "Definitive_erythroid_L2" = "#dc453e",
                            "Erythroid_progenitors_L2" = "#cc3f74",
                            "Primitive_erythroid_L2" = "#ffa9a1",
                            "Hepatocytes_L2" = "#185700",
                            "Lateral_plate_and_intermediate_mesoderm_L2" = "#d78860",
                            "Mesoderm_L2" = "#bb46c5",
                            "Cardiomyocytes_L2" = "#643e8c",
                            "Skeletal_muscle_cells_L2" = "#ffa1f5",
                            "NMPs_and_spinal_cord_progenitors_L2" = "#cf6a79",
                            "Melanocyte_cells_L2" = "#88b439",
                            "Neural_crest_PNS_glia_L2" = "#fff167",
                            "Olfactory_ensheathing_cells_L2" = "#b46aa2",
                            "CNS_neurons_L2" = "#e5c000",
                            "Eye_L2" = "#00d450",
                            "Intermediate_neuronal_progenitors_L2" = "#2e0ab7",
                            "Neural_crest_PNS_neurons_L2" = "#b5ce92",
                            "Neuroectoderm_and_glia_L2" = "#f96100",
                            "Corticofugal_neurons_L2" = "#8d6a32",
                            "Olfactory_sensory_neurons_L2" = "#e6230b",
                            "Committed_oligodendrocyte_precursors_L2" = "#408352",
                            "Oligodendrocyte_progenitor_cells_L2" = "#916e00",
                            "B_cells_L2" = "#01b7a6",
                            "T_cells_L2" = "#ff9d47",
                            "White_blood_cells_L2" = "#7ca0ff")

dat = read.csv(paste0(web_path, "/Fig3d_bar_plot_selected_sites.csv"))

dat$celltype_name = factor(dat$celltype_name, levels = unique(dat$celltype_name))

plot_tmp = theme_classic(base_size = 12) +
    theme(legend.position="none") +
    theme(plot.title = element_text(hjust = 0.5)) +
    theme(axis.text.x = element_text(color="black"), axis.text.y = element_text(color="black"))

target_site_list = c("chr9:15561725-15563839", "chr9:92825930-92828044", "chr9:70276731-70278845")
for(target_site in target_site_list){
    print(target_site)
    df = dat[dat$target_site == target_site,]
    print(head(df))
    p1 = ggplot(data = df, aes(x = celltype_name, y = pre, fill = celltype_name)) +
        geom_bar(stat="identity") + plot_tmp + scale_fill_manual(values=celltype_L2_color_plate)
    p2 = ggplot(data = df, aes(x = celltype_name, y = obs, fill = celltype_name)) +
        geom_bar(stat="identity") + plot_tmp + scale_fill_manual(values=celltype_L2_color_plate)
    ggsave(paste0("Fig3d_bar_plot_selected_sites.", target_site, ".pdf"), p1 / p2, width = 6, height = 3)
}
```

+++

## Fig. 3e

Pairwise overlaps (bp) between 100 bp windows with high predicted accessibility and specificity (both a Z-transformed predicted score ≥ 3 and a Z-transformed specificity ≥ 3; rows) and the top 10,000 specific peaks used for training (columns) across all 36 Level-2 cell classes, normalized by the total size of each class’s training peaks.


```{code-cell #Fig3e_overlap_pre_train} r
web_path = "https://shendure-web.gs.washington.edu/content/members/cxqiu/public/backup/jax_atac/download/reproducing_figures"

library(gplots)
library(reshape2)
library(viridis)

library(RColorBrewer)
Colors=rev(brewer.pal(11,"Spectral"))
Colors=colorRampPalette(Colors)(120)

dat = read.csv(paste0(web_path, "/Fig3e_overlap_pre_train.csv"), row.names=1)

# pdf("Fig3e_overlap_pre_train.pdf")
heatmap.2(as.matrix(dat),
          col=viridis(100),
          scale="none",
          Rowv = F,
          Colv = F,
          key=T,
          density.info="none",
          trace="none",
          cexRow = 0.5,
          cexCol = 0.5,
          margins = c(5,5))
# dev.off()
```

+++

## Fig. 3g

Scatterplots of observed (normalized Tn5 cut-site counts; natural log-scaled; x-axis) vs. predicted (evolution-naive CREsted model; natural log-scaled; y-axis) chromatin accessibility across 1.25 million x 100 bp windows spanning mouse chromosome 9 for all 36 Level-2 cell classes.


```{code-cell #Fig3g_pred_chr9} r
web_path = "https://shendure-web.gs.washington.edu/content/members/cxqiu/public/backup/jax_atac/download/reproducing_figures"

library(gplots)
library(reshape2)
library(viridis)

dat = read.csv(paste0(web_path, "/Fig3g_pred_chr9.csv"))
dat$chr = "chr9"
dat$start = seq(1458, 124593458, 100)
dat$end = seq(1558, 124593558, 100)

p = ggplot(data = dat, aes(x=log_obs, y=log_pre)) +
    geom_hex(bins = 70) +
    geom_abline(intercept = 2, slope = 1, linetype = "dashed", color = "red") +
    geom_abline(intercept = 3, slope = 1, color = "red") +
    geom_abline(intercept = -3, slope = 1, color = "red") +
    scale_x_continuous(limits = c(-0.3, 5.25), expand = expansion(mult = 0, add = 0)) +
    scale_y_continuous(limits = c(-0.3, 5.25), expand = expansion(mult = 0, add = 0)) +
    labs(x="Log (observation + 1)", y="Log (prediction + 1)", title=paste0("Corr = ", round(cor(dat$log_obs, dat$log_pre), 2))) +
    scale_fill_continuous(type = "viridis", trans = "log10") +
    theme_classic(base_size = 10) +
    theme(legend.position="none") +
    theme(plot.title = element_text(hjust = 0.5)) +
    theme(axis.text.x = element_text(color="black"), axis.text.y = element_text(color="black"))

ggsave("Fig3g_pred_chr9.pdf", p, height = 5, width = 5)
p
```

+++

## Fig. 3h

A plot similar to panel g, but instead of coloring each bin by the counts of 100 bp windows, it is colored by the difference between the percentage of windows overlapping tandem repeats and the percentage overlapping promoters (−1,500 bp to +500 bp relative to TSSs).


```{code-cell #Fig3h_pred_chr9_overlap_promoter_TRF} r
web_path = "https://shendure-web.gs.washington.edu/content/members/cxqiu/public/backup/jax_atac/download/reproducing_figures"

library(gplots)
library(reshape2)
library(viridis)
library(RColorBrewer)
Colors=rev(brewer.pal(11,"Spectral"))
Colors=colorRampPalette(Colors)(120)

dat = read.csv(paste0(web_path, "/Fig3h_pred_chr9_overlap_promoter_TRF.csv"))

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

ggsave("Fig3h_pred_chr9_overlap_promoter_TRF.pdf", p, height = 5, width = 5)
p
```
