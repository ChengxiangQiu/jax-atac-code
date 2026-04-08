---
title: Figure 6
kernelspec:
  name: ir
  display_name: R
---

Evolutionary augmentation generalizes genomewide regulatory grammar inference across Mammalia.

+++

## Fig. 6b

Models were trained on mouse windows and orthologs across 1 to 241 species, scaling data from 0.3M to 58.5M sequences, with performance evaluated on held-out mouse data.

```{code-cell #Fig6b_evaluate_multispecies_model} r
web_path = "https://shendure-web.gs.washington.edu/content/members/cxqiu/public/backup/jax_atac/download/reproducing_figures"

library(ggplot2)
library(dplyr)
library(viridis)

dat = read.csv(paste0(web_path, "/Fig6b_evaluate_multispecies_model.csv"))
dat$num_species = factor(dat$num_species)

p = ggplot(dat, aes(x = epoch, y = pearson_r, color = num_species, group = num_species)) +
    geom_point() +
    geom_smooth(method = 'loess', formula = 'y ~ x', se = FALSE) +
    theme_classic(base_size = 10)

ggsave("Fig6b_evaluate_multispecies_model.pdf", p, width = 6, height = 5)
p
```

+++

## Fig. 6d

We examined predicted chromatin accessibility using the STEAM-v1 model across a 200 kb region centered on the Afp TSS.

```{code-cell #Fig6d_Afp_Hepatocytes_prediction} r
web_path = "https://shendure-web.gs.washington.edu/content/members/cxqiu/public/backup/jax_atac/download/reproducing_figures"

library(ggplot2)
library(dplyr)
library(viridis)

cluster_color_plate = c(
    "cluster_0" =  "#7f7f7f",
    "cluster_1" =  "#ff7f0e",
    "cluster_2" =  "#2ca02c",
    "cluster_3" =  "#d62728",
    "cluster_4" =  "#9467bd",
    "cluster_5" =  "#8c564b",
    "cluster_6" =  "#e377c2",
    "cluster_7" =  "#ff9896",
    "cluster_8" =  "#bcbd22",
    "cluster_9" =  "#17becf",
    "cluster_10" = "#1f77b4",
    "cluster_11" = "#ffbb78")

dat = read.csv(paste0(web_path, "/Fig6d_Afp_Hepatocytes_enhancers.csv"))

p = ggplot(dat, aes(x_axis, y_axis, fill= cluster)) +
    geom_tile() +
    scale_fill_manual(values=cluster_color_plate) +
    scale_x_discrete(expand = c(0, 0)) +
    scale_y_discrete(expand = c(0, 0)) +
    theme_void() +
    theme(legend.position="none")

ggsave("Fig6d_Afp_Hepatocytes_enhancers.png", p, width = 5, height = 5, dpi = 300)

dat = read.csv(paste0(web_path, "/Fig6d_Afp_Hepatocytes_prediction.csv"))

p = ggplot(dat, aes(x_axis, y_axis, fill= score)) +
    geom_tile() +
    scale_fill_viridis(discrete=FALSE) +
    scale_x_discrete(expand = c(0, 0)) +
    scale_y_discrete(expand = c(0, 0)) +
    theme_void() +
    theme(legend.position="none",
          panel.background = element_rect(fill = "#440154", color = NA))

ggsave("Fig6d_Afp_Hepatocytes_prediction.png", p, width = 5, height = 5, dpi = 300)
p
```

+++

## Fig. 6e

Across nine clusters, lifted regions in species lacking enhancers defined 570 non-enhancers (vs. 591 enhancers), compared by STEAM-v1 scores and Afp TSS distance.

```{code-cell #Fig6e_compare_enhancers_and_non_enhancers} r
web_path = "https://shendure-web.gs.washington.edu/content/members/cxqiu/public/backup/jax_atac/download/reproducing_figures"

library(ggplot2)
library(dplyr)
library(patchwork)

cluster_color_plate = c(
    "cluster_1"  = "#ff7f0e",
    "cluster_2"  = "#2ca02c",
    "cluster_3"  = "#d62728",
    "cluster_4"  = "#9467bd",
    "cluster_5"  = "#8c564b",
    "cluster_6"  = "#e377c2",
    "cluster_7"  = "#ff9896",
    "cluster_8"  = "#bcbd22",
    "cluster_9"  = "#17becf",
    "cluster_10" = "#1f77b4",
    "cluster_11" = "#ffbb78")

df = read.csv(paste0(web_path, "/Fig6e_compare_enhancers_and_non_enhancers.csv"), sep = "\t")

# Numbered y-axis labels: "1" through "11"
cluster_levels = rev(names(cluster_color_plate))
cluster_labels = as.character(rev(seq_along(cluster_color_plate)))
df$cluster = factor(df$cluster, levels = cluster_levels)
df$detected = df$if_enhancer_detected == "yes"

# Count total orthologs per cluster for right-side labels
counts = df %>% count(cluster) %>% arrange(cluster)
count_labels = setNames(as.character(counts$n), counts$cluster)

base_plot = function(df, x_col) {
    ggplot(df, aes(x = .data[[x_col]], y = cluster, color = cluster)) +
        geom_jitter(data = df %>% filter(!detected),
                    size = 1.2, alpha = 0.5, height = 0.35, shape = 16) +
        geom_jitter(data = df %>% filter(detected),
                    aes(fill = cluster), size = 1.4, alpha = 0.85, height = 0.35,
                    shape = 21, color = "black", stroke = 0.3) +
        scale_color_manual(values = cluster_color_plate) +
        scale_fill_manual(values = cluster_color_plate) +
        scale_y_discrete(labels = setNames(cluster_labels, cluster_levels),
                         sec.axis = dup_axis(labels = count_labels,
                                             name = "Total orthologs")) +
        theme_classic(base_size = 10) +
        theme(legend.position = "none",
              panel.grid.major.x = element_line(color = "grey90")) +
        labs(y = "Synteny group")
}

p1 = base_plot(df, "distance") +
    scale_x_continuous(labels = function(x) x / 1000,
                       breaks = seq(-100000, 100000, 25000)) +
    labs(x = "Distance to Afp TSS (kb)")

p2 = base_plot(df, "log_score") +
    labs(x = expression(Log[2]~(Predicted~score + 1)))

p3 = base_plot(df, "log_phred_score") +
    labs(x = "log_phred_score")

p = p1 / p2 / p3

ggsave("Fig6e_compare_enhancers_and_non_enhancers.pdf", p, width = 11, height = 9)
p
```

+++

## Fig. 6f

Zoom-in -10k to +2k of Afp TSS region.

```{code-cell #Fig6f_Afp_Hepatocytes_prediction_subregion} r
web_path = "https://shendure-web.gs.washington.edu/content/members/cxqiu/public/backup/jax_atac/download/reproducing_figures"

library(ggplot2)
library(dplyr)

cluster_color_plate = c(
    "cluster_0" =  "#7f7f7f",
    "cluster_1" =  "#ff7f0e",
    "cluster_2" =  "#2ca02c",
    "cluster_3" =  "#d62728",
    "cluster_4" =  "#9467bd",
    "cluster_5" =  "#8c564b",
    "cluster_6" =  "#e377c2",
    "cluster_7" =  "#ff9896",
    "cluster_8" =  "#bcbd22",
    "cluster_9" =  "#17becf",
    "cluster_10" = "#1f77b4",
    "cluster_11" = "#ffbb78")

dat = read.csv(paste0(web_path, "/Fig6d_Afp_Hepatocytes_enhancers.csv"))

p = ggplot(dat[dat$x_axis > (-10000) & dat$x_axis < 2000,], aes(x_axis, y_axis, fill= cluster)) +
    geom_tile() +
    scale_fill_manual(values=cluster_color_plate) +
    scale_x_continuous(limits = c(-10000, 2000), expand = c(0, 0)) +
    scale_y_discrete(expand = c(0, 0)) +
    theme_void() +
    theme(legend.position="none")

ggsave("Fig6d_Afp_Hepatocytes_enhancers.png", p, width = 5, height = 5, dpi = 300)

dat = read.csv(paste0(web_path, "/Fig6d_Afp_Hepatocytes_prediction.csv"))

p = ggplot(dat[dat$x_axis > (-100) & dat$x_axis < 20,], aes(x_axis, y_axis, fill= score)) +
    geom_tile() +
    scale_fill_viridis(discrete=FALSE) +
    scale_x_continuous(limits = c(-100, 20), expand = c(0, 0)) +
    scale_y_discrete(expand = c(0, 0)) +
    theme_void() +
    theme(legend.position="none",
          panel.background = element_rect(fill = "#440154", color = NA))

ggsave("Fig6f_Afp_Hepatocytes_prediction_subregion.png", p, width = 5, height = 5, dpi = 300)
p
```
