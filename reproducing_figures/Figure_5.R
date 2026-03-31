
###################################################################################################################
### Figure 5. Evolution-aware modeling identifies distal enhancers that predict cell-type-specific gene expression.

##########################################################################################################################################
### Fig. 5b: Distribution of mean phyloP scores (60 placental mammals), computed by averaging nucleotide-level scores within each element.

web_path = "https://shendure-web.gs.washington.edu/content/members/cxqiu/public/backup/jax_atac/download/reproducing_figures"

library(ggplot2)
library(dplyr)

df = read.csv(paste0(web_path, "/Fig5b_phyloP_distribution.csv"))

p = ggplot(df, aes(x = phyloP, fill = category)) +
    geom_density(alpha = 0.5) +
    scale_fill_manual(values = c("protein_coding" = "#cb6751", "enhancer" = "#7aa457", "background" = "#9e6ebd")) +
    labs(x="PhyloP scores", y="Density", title="") +
    theme_classic(base_size = 10) +
    theme(plot.title = element_text(hjust = 0.5)) +
    theme(legend.position="none") +
    theme(axis.text.x = element_text(color="black"), axis.text.y = element_text(color="black")) +
    coord_cartesian(xlim = c(-3, 3)) +
    facet_wrap(~category, scales = "free_y", ncol = 1)

ggsave("Fig5b_phyloP_distribution.pdf", p, width = 6, height = 4.5)


####################################################################################################################################
### Fig. 5d: Relationship between cell-class-specific enhancer prediction scores vs. cell-class-specific expression of nearby genes.

web_path = "https://shendure-web.gs.washington.edu/content/members/cxqiu/public/backup/jax_atac/download/reproducing_figures"

library(ggplot2)
library(dplyr)

dat = read.table(paste0(web_path, "/Fig5d_predict_gene_exp.csv"), sep='\t', header=T)

bin_label_color = c("#440154", "#46327e", "#365c8d", "#277f8e", "#1fa187", "#4ac16d", "#a0da39", "#e65d2f")
names(bin_label_color) = names(table(dat$bin_label))

p = ggplot() +
    geom_point(data = dat, aes(x = distance_index, y = delta_exp, color = bin_label)) +
    geom_smooth(data = dat, aes(x = distance_index, y = delta_exp, group = bin_label, color = bin_label), method = 'loess', formula = 'y ~ x', se = FALSE) +
    labs(
        x = "Distance from gene TSS to peak (kb)",
        y = "Cell class-specific gene expression metric",
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

ggsave("Fig5d_predict_gene_exp.pdf", p, height = 5, width = 6)












