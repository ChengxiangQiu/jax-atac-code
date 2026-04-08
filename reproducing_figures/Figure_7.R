
################################################################################
### Figure 7. STEAM-v1 predicts human and mouse genome-wide enhancer landscapes. 

#############################################################################################################################################################################################################
### Fig. 7a: Genome-wide STEAM-v1 enhancer predictions were merged across cell classes. Histogram of merged regions by the number of cell classes for which they are predicted enhancers, in mouse and human.

web_path = "https://shendure-web.gs.washington.edu/content/members/cxqiu/public/backup/jax_atac/download/reproducing_figures"

library(ggplot2)
library(dplyr)

df = read.csv(paste0(web_path, "/Fig7a_hist_num_enhancer_cell_class.csv"))

df$num_celltypes = factor(df$num_celltypes, levels =1:10)
df$species = factor(df$species, levels = c("mouse", "human"))

p = ggplot(data=df, aes(x=num_celltypes, y=num_enhancers, fill = species)) +
    geom_bar(stat="identity", position=position_dodge()) +
    theme_classic(base_size = 10) +
    scale_fill_manual(values=c("human" = "#658bca", "mouse" = "#d79c36"))

ggsave("Fig7a_hist_num_enhancer_cell_class.pdf", p, width = 5, height = 5)


###################################################################################################################################################################
### Fig. 7b: For each category in panel a, the percentage of merged regions overlapping a promoter (+/- 2.5 kb around an annotated TSS) is shown for mouse and human.

web_path = "https://shendure-web.gs.washington.edu/content/members/cxqiu/public/backup/jax_atac/download/reproducing_figures"

library(ggplot2)
library(dplyr)

df = read.csv(paste0(web_path, "/Fig7b_hist_enhancer_enrich_promoter.csv"))

df$num_celltype = factor(df$num_celltype, levels =1:10)
df$species = factor(df$species, levels = c("mouse", "human"))

p = ggplot(data=df, aes(x=num_celltype, y=frac, fill = species)) +
    geom_bar(stat="identity", position=position_dodge()) +
    geom_hline(yintercept = 9) +
    geom_hline(yintercept = 10) +
    theme_classic(base_size = 10) +
    scale_fill_manual(values=c("human" = "#658bca", "mouse" = "#d79c36"))

ggsave("Fig7b_hist_enhancer_enrich_promoter.pdf", p, width = 5, height = 5)


###################################################################################################################################################################
### Fig. 7c: For each category in panel a, the percentage of merged regions overlapping a promoter (+/- 2.5 kb around an annotated TSS) is shown for mouse and human.

web_path = "https://shendure-web.gs.washington.edu/content/members/cxqiu/public/backup/jax_atac/download/reproducing_figures"

library(gplots)
library(viridis)

df = read.csv(paste0(web_path, "/Fig7c_heatmap_mouse.csv"))

pdf("Fig7c_heatmap_mouse.pdf", 8, 8)
p = heatmap.2(as.matrix(df), 
              col=viridis(100), 
              breaks = seq(0, 0.66, length.out = 101),
              scale="none", 
              Rowv = F, 
              Colv = F, 
              key=T, 
              density.info="none", 
              trace="none", 
              cexRow = 0.5, 
              cexCol = 0.5,
              margins = c(5,5))
dev.off()

df = read.csv(paste0(web_path, "/Fig7c_heatmap_human.csv"))

pdf("Fig7c_heatmap_human.pdf", 8, 8)
p = heatmap.2(as.matrix(df), 
              col=viridis(100), 
              breaks = seq(0, 0.66, length.out = 101),
              scale="none", 
              Rowv = F, 
              Colv = F, 
              key=T, 
              density.info="none", 
              trace="none", 
              cexRow = 0.5, 
              cexCol = 0.5,
              margins = c(5,5))
dev.off()


###################################################################################################################################################################
### Fig. 7d: Human and mouse predicted hepatocyte enhancers were reciprocally lifted over: mouse-only (n = 4,948), human-mouse shared (n = 2,597), and human-only (n = 5,790) predictions. Observed chromatin accessibility from human fetal liver hepatoblasts and mouse hepatocytes (this study) was extracted at averaged base resolution, and the distribution of log2-scaled human/mouse fold-changes is plotted for the three enhancer sets 

web_path = "https://shendure-web.gs.washington.edu/content/members/cxqiu/public/backup/jax_atac/download/reproducing_figures"

library(ggplot2)
library(dplyr)

df = read.csv(paste0(web_path, "/Fig7d_human_mouse_enhancer_overlap.csv"))

df$category = factor(df$category, levels = c("mouse_only", "human_mouse_share", "human_only"))

p = ggplot(df[df$celltype == "Hepatocytes",], aes(x = category, y = log2_human_mouse_obs_ratio, fill = celltype)) +
    geom_boxplot() + 
    labs(x="", y="Log2 (Hum Obs / Mus Obs)", title="") +
    theme_classic(base_size = 10) +
    theme(plot.title = element_text(hjust = 0.5)) +
    theme(legend.position="none") +
    theme(axis.text.x = element_text(color="black"), axis.text.y = element_text(color="black")) +
    scale_fill_manual(values=c("Hepatocytes" = "#7aa457", "Erythroid_cells" = "#cb6751"))

ggsave("Fig7d_human_mouse_enhancer_overlap_hepatocyte.pdf", p, width = 5, height = 5)


###################################################################################################################################################################
### Fig. 7e: Same as panel d for erythroid enhancers: mouse-only (n = 3,492), human-mouse shared (n = 1,045), and human-only (n = 4,353) predicted enhancers with successful cross-species liftover. Observed chromatin accessibility from human fetal liver erythroblasts and mouse definitive erythroid cells (this study) was used.

web_path = "https://shendure-web.gs.washington.edu/content/members/cxqiu/public/backup/jax_atac/download/reproducing_figures"

library(ggplot2)
library(dplyr)

df = read.csv(paste0(web_path, "/Fig7d_human_mouse_enhancer_overlap.csv"))

df$category = factor(df$category, levels = c("mouse_only", "human_mouse_share", "human_only"))

p = ggplot(df[df$celltype == "Erythroid_cells",], aes(x = category, y = log2_human_mouse_obs_ratio, fill = celltype)) +
    geom_boxplot() + 
    labs(x="", y="Log2 (Hum Obs / Mus Obs)", title="") +
    theme_classic(base_size = 10) +
    theme(plot.title = element_text(hjust = 0.5)) +
    theme(legend.position="none") +
    theme(axis.text.x = element_text(color="black"), axis.text.y = element_text(color="black")) +
    scale_fill_manual(values=c("Hepatocytes" = "#7aa457", "Erythroid_cells" = "#cb6751"))

ggsave("Fig7e_human_mouse_enhancer_overlap_erythroid.pdf", p, width = 5, height = 5)



###################################################################################################################################################################
### Fig. 7f: Human-specific hepatocyte (n = 3,239; green) and erythroid (n = 3,842; red) predicted enhancers failing liftover to the mouse genome are plotted by observed chromatin accessibility in human fetal liver erythroblasts (x-axis) vs. hepatoblasts (y-axis), log2-scaled.

web_path = "https://shendure-web.gs.washington.edu/content/members/cxqiu/public/backup/jax_atac/download/reproducing_figures"

library(ggplot2)
library(dplyr)

df = read.csv(paste0(web_path, "/Fig7f_human_only.csv"))

p = ggplot(df[sample(1:nrow(df)),], aes(x = log2_Erythroid_cells, y = log2_Hepatocytes, color = celltype)) +
    geom_point(alpha = 0.6) +
    labs(x="Log2(Erythroid Obs)", y="Log2(Hepatocyte Obs)", title="") +
    theme_classic(base_size = 10) +
    scale_color_manual(values=c("Hepatocytes" = "#7aa457", "Erythroid_cells" = "#cb6751")) +
    coord_cartesian(xlim = c(-15, -2), ylim = c(-15, -2)) +
    theme(plot.title = element_text(hjust = 0.5)) +
    theme(legend.position="none") +
    theme(axis.text.x = element_text(color="black"), axis.text.y = element_text(color="black")) 

ggsave("Fig7f_human_only_scatter.pdf", p, width = 5, height = 5)



###################################################################################################################################################################
### Fig. 7g: For the same sets of human-specific enhancer predictions failing liftover to the mouse genome shown in panel f, the log2-scaled ratio of observed human hepatoblast to erythroblast accessibility is compared between hepatocyte-predicted (left) and erythroid-predicted (right) enhancers.

web_path = "https://shendure-web.gs.washington.edu/content/members/cxqiu/public/backup/jax_atac/download/reproducing_figures"

library(ggplot2)
library(dplyr)

df = read.csv(paste0(web_path, "/Fig7f_human_only.csv"))

df$celltype = factor(df$celltype, levels = c("Hepatocytes", "Erythroid_cells"))

p = ggplot(df, aes(x = celltype, y = log2_fc, fill = celltype)) +
    geom_boxplot() +
    labs(x="", y="Log2(Hepatocytes/Erythroid_cells obs)", title="") +
    scale_fill_manual(values=c("Hepatocytes" = "#7aa457", "Erythroid_cells" = "#cb6751")) +
    theme_classic(base_size = 10) +
    theme(plot.title = element_text(hjust = 0.5)) +
    theme(legend.position="none") +
    theme(axis.text.x = element_text(color="black"), axis.text.y = element_text(color="black")) 

ggsave("Fig7f_human_only_boxplot.pdf", p, width = 5, height = 5)














