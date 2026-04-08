---
title: Supplementary Figure 7
kernelspec:
  name: ir
  display_name: R
---

Correlation between observed and predicted chromatin accessibility for each cell class.

```{code-cell #Supplementary_Fig7_obs_vs_pre_split_cell_class} r
web_path = "https://shendure-web.gs.washington.edu/content/members/cxqiu/public/backup/jax_atac/download/reproducing_figures"

library(ggplot2)
library(magrittr)

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

celltype_list = c("Mesoderm_L2",
                  "NMPs_and_spinal_cord_progenitors_L2",
                  "Neuroectoderm_and_glia_L2",
                  "Cardiomyocytes_L2",
                  "CNS_neurons_L2",
                  "White_blood_cells_L2",
                  "Neural_crest_PNS_neurons_L2",
                  "Lateral_plate_and_intermediate_mesoderm_L2",
                  "Skeletal_muscle_cells_L2",
                  "Epithelial_cells_L2",
                  "Endothelium_L2",
                  "Neural_crest_PNS_glia_L2",
                  "Definitive_erythroid_L2",
                  "Eye_L2",
                  "Primitive_erythroid_L2",
                  "Hepatocytes_L2",
                  "Liver_sinusoidal_endothelial_cells_L2",
                  "Kidney_L2",
                  "Endocardial_cells_L2",
                  "Gut_epithelial_cells_L2",
                  "Olfactory_sensory_neurons_L2",
                  "Brain_capillary_endothelial_cells_L2",
                  "T_cells_L2",
                  "Lymphatic_vessel_endothelial_cells_L2",
                  "Lung_and_airway_L2",
                  "Erythroid_progenitors_L2",
                  "Intermediate_neuronal_progenitors_L2",
                  "Melanocyte_cells_L2",
                  "Glomerular_endothelial_cells_L2",
                  "Adipocyte_cells_L2",
                  "Committed_oligodendrocyte_precursors_L2",
                  "Oligodendrocyte_progenitor_cells_L2",
                  "Olfactory_ensheathing_cells_L2",
                  "Corticofugal_neurons_L2",
                  "B_cells_L2",
                  "Adipocyte_cells_Cyp2e1_L2")

df = read.csv(paste0(web_path, "/Supplementary_Fig7_obs_vs_pre_split_cell_class.csv"))

df$celltype_name = factor(df$celltype_name, levels = as.vector(celltype_list))

p = df %>%
    ggplot(aes(obs, pre, color = celltype_name)) +
    geom_point(size = 1.5, alpha = 0.6) +
    labs(x="Observation", y="Prediction") +
    theme_classic(base_size = 12) +
    scale_color_manual(values=celltype_L2_color_plate) +
    theme(
        legend.position = "none",
        plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(color = "black", size = 10, family = "Arial"),
        axis.text.y = element_text(color = "black", size = 10, family = "Arial"),
        axis.title.x = element_text(size = 12, family = "Arial"),
        axis.title.y = element_text(size = 12, family = "Arial")
    ) +
    facet_wrap(~celltype_name, nrow = 6, ncol = 6, scales = "free")
ggsave("Supplementary_Fig7_obs_vs_pre_split_cell_class.png", p, dpi = 300, width = 20, height = 20)
p
```
