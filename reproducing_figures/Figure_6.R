
############################################################################################################
### Figure 6. Evolutionary augmentation generalizes genomewide regulatory grammar inference across Mammalia.

######################################################################################################################################################################################
### Fig. 6b: Models were trained on mouse windows and orthologs across 1 to 241 species, scaling data from 0.3M to 58.5M sequences, with performance evaluated on held-out mouse data.

web_path = "https://shendure-web.gs.washington.edu/content/members/cxqiu/public/backup/jax_atac/download/reproducing_figures"

library(ggplot2)
library(dplyr)

dat = read.csv(paste0(web_path, "/Fig6b_evaluate_multispecies_model.csv"))
dat$num_species = factor(dat$num_species)

p = ggplot(dat, aes(x = epoch, y = pearson_r, color = num_species, group = num_species)) +
    geom_point() +
    geom_smooth(method = 'loess', formula = 'y ~ x', se = FALSE) +
    theme_classic(base_size = 10)

ggsave("Fig6b_evaluate_multispecies_model.pdf", p, width = 6, height = 5)


###################################################################################################################################
### Fig. 6d: We examined predicted chromatin accessibility using the STEAM-v1 model across a 200 kb region centered on the Afp TSS.

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


###############################################################################################################################################################################
### Fig. 6e: Across nine clusters, lifted regions in species lacking enhancers defined 570 non-enhancers (vs. 591 enhancers), compared by STEAM-v1 scores and Afp TSS distance.

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

df = pd.read_csv(f"{web_path}/Fig6e_compare_enhancers_and_non_enhancers.csv", sep="\t")

df["cluster_n"] = df["cluster"].str.replace("cluster_", "").astype(int)
clusters = sorted(df["cluster_n"].unique())
cluster_to_y = {c: i for i, c in enumerate(clusters)}

cluster_color_plate = {
"cluster_1": "#ff7f0e",
"cluster_2": "#2ca02c",
"cluster_3": "#d62728",
"cluster_4": "#9467bd",
"cluster_5": "#8c564b",
"cluster_6": "#e377c2",
"cluster_7": "#ff9896",
"cluster_8": "#bcbd22",
"cluster_9": "#17becf",
"cluster_10": "#1f77b4",
"cluster_11": "#ffbb78"}

metrics = ["distance", "log_score", "log_phred_score"]
fig, axes = plt.subplots(3, 1, figsize=(11, 9), sharey=True)
rng = np.random.default_rng(0)

for ax, metric in zip(axes, metrics):
    for c in clusters:
        sub = df[df["cluster_n"] == c]
        x = sub[metric].values
        y = np.full(len(x), cluster_to_y[c]) + rng.uniform(-0.35, 0.35, len(x))
        det = (sub["if_enhancer_detected"] == "yes").values
        cluster_name = f"cluster_{c}"
        col = cluster_color_plate.get(cluster_name, "grey")
        ax.scatter(x[~det], y[~det], s=14, alpha=0.5, color=col, edgecolors="none")
        ax.scatter(x[det], y[det], s=16, alpha=0.85, color=col, edgecolors="black", linewidths=0.5)
    ax.set_yticks(range(len(clusters)))
    ax.set_yticklabels([f"cluster_{c}" for c in clusters])
    ax.set_xlabel(metric)
    ax.invert_yaxis()
    ax.grid(axis="x", alpha=0.3)
    ax.set_axisbelow(True)

axes[0].set_title("Per-species ortholog scores by cluster  (n per cluster shown at right)")

for c in clusters:
    n = (df["cluster_n"] == c).sum()
    col = cluster_color_plate.get(f"cluster_{c}", "grey")
    axes[0].annotate(f"n={n}", xy=(1.01, cluster_to_y[c]),
                 xycoords=("axes fraction", "data"),
                 va="center", fontsize=8, color=col)

plt.tight_layout()
out = "Fig6e_compare_enhancers_and_non_enhancers.pdf"
plt.savefig(out, bbox_inches="tight")
print(f"wrote {out}")


###################################################
### Fig. 6f: Zoom-in -10k to +2k of Afp TSS region.

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




