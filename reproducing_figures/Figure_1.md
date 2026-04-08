---
title: Figure 1
kernelspec:
  name: ir
  display_name: R
---

A single-cell chromatin-accessibility timelapse of mouse development from E10 to P0.

+++

## Fig. 1b

The log2 scaled number of cells profiled per embryo, sorted by staging bin.

```{code-cell #Fig1b_cell_num_per_timepoint} r
web_path = "https://shendure-web.gs.washington.edu/content/members/cxqiu/public/backup/jax_atac/download/reproducing_figures"

library(ggplot2)
library(scales)
library(dplyr)

day_color_plate = c("#5E4FA2", "#515FA9", "#446FB1", "#387FB9", "#3990B9", "#48A0B2", "#57B1AB",
                    "#65C2A5", "#79C9A4", "#8DD1A4", "#A1D9A4", "#B3E0A2", "#C4E79E", "#D5EE9B",
                    "#E5F498", "#EDF7A3", "#F4FAAE", "#FBFDB9", "#FEFAB7", "#FEF1A8", "#FEE899",
                    "#FEE08B", "#FDD17F", "#FDC373", "#FDB567", "#FBA45C", "#F99254", "#F67F4B",
                    "#F46D43", "#EB5F46", "#E25249", "#D9444D", "#CD354D", "#BD2349", "#AD1245", "#9E0142")
day_list = c("E10.0", "E10.25", "E10.5", "E10.75", "E11.0", "E11.25", "E11.5",
             "E11.75", "E12.0", "E12.25", "E12.5", "E12.75", "E13.0", "E13.25",
             "E13.5", "E13.75", "E14.0", "E14.25", "E14.375", "E14.75", "E15.0",
             "E15.25", "E15.5", "E15.75", "E16.0", "E16.25", "E16.5", "E16.75",
             "E17.0", "E17.25", "E17.5", "E18.0", "E18.25", "E18.5", "E18.75", "P0")
names(day_color_plate) = day_list

dat = read.csv(paste0(web_path, "/Fig1b_cell_num_per_timepoint.csv"))
dat$day = factor(dat$day, levels = rev(names(day_color_plate)))

p = ggplot(data = dat, aes(day, cell_num, fill = day)) +
    geom_bar(stat="identity") +
    coord_flip() +
    scale_fill_manual(values = day_color_plate) +
    scale_y_continuous(trans = log2_trans(), breaks = trans_breaks("log2", function(x) 2^x), labels = trans_format("log2", math_format(2^.x))) +
    geom_text(aes(label = scales::comma(cell_num)),
              hjust = -0.1,
              position = position_dodge(width = 1),
              inherit.aes = TRUE,
              size = 3) +
    labs(x = "Embryonic Day", y = "Cell number (log2)") +
    theme_classic(base_size = 15) +
    theme(legend.position="none") +
    theme(axis.text.x = element_text(color="black"), axis.text.y = element_text(color="black"))

ggsave("Fig1b_cell_num_per_timepoint.pdf", p, height = 5, width = 5)
p
```

+++

## Fig. 1c

Embeddings of pseudo-bulk ATAC-seq profiles of 36 embryos in PCA space, plotting PC1 vs. PC2.

```{code-cell #Fig1c_pseudobulk_PCA} r
web_path = "https://shendure-web.gs.washington.edu/content/members/cxqiu/public/backup/jax_atac/download/reproducing_figures"

library(ggplot2)
library(scales)
library(dplyr)

day_color_plate = c("#5E4FA2", "#515FA9", "#446FB1", "#387FB9", "#3990B9", "#48A0B2", "#57B1AB",
                    "#65C2A5", "#79C9A4", "#8DD1A4", "#A1D9A4", "#B3E0A2", "#C4E79E", "#D5EE9B",
                    "#E5F498", "#EDF7A3", "#F4FAAE", "#FBFDB9", "#FEFAB7", "#FEF1A8", "#FEE899",
                    "#FEE08B", "#FDD17F", "#FDC373", "#FDB567", "#FBA45C", "#F99254", "#F67F4B",
                    "#F46D43", "#EB5F46", "#E25249", "#D9444D", "#CD354D", "#BD2349", "#AD1245", "#9E0142")
day_list = c("E10.0", "E10.25", "E10.5", "E10.75", "E11.0", "E11.25", "E11.5",
             "E11.75", "E12.0", "E12.25", "E12.5", "E12.75", "E13.0", "E13.25",
             "E13.5", "E13.75", "E14.0", "E14.25", "E14.375", "E14.75", "E15.0",
             "E15.25", "E15.5", "E15.75", "E16.0", "E16.25", "E16.5", "E16.75",
             "E17.0", "E17.25", "E17.5", "E18.0", "E18.25", "E18.5", "E18.75", "P0")
names(day_color_plate) = day_list

dat = read.csv(paste0(web_path, "/Fig1c_pseudobulk_PCA.csv"))

p = ggplot(data = dat) +
    geom_point(aes(x = PC_2, y = PC_1), color = "grey80", size=5) +
    geom_point(aes(x = PC_2, y = PC_1, color = day), size=4.5) +
    theme_classic(base_size = 12) +
    scale_color_manual(values = day_color_plate) +
    labs(x = "PC_2 (12.3%)", y = "PC_1 (76.8%)") +
    theme(legend.position="none") +
    theme(axis.text.x = element_text(color="black"), axis.text.y = element_text(color="black"))

ggsave("Fig1c_pseudobulk_PCA.pdf", p, height = 5, width = 5)
p
```

+++

## Fig. 1d

3D UMAP visualization of co-embedded cells from published scRNA-seq and scATAC-seq (this study) data.

```{code-cell #Fig1d_coembedding_RNA_ATAC} r
web_path = "https://shendure-web.gs.washington.edu/content/members/cxqiu/public/backup/jax_atac/download/reproducing_figures"

library(plotly)
library(htmlwidgets)

t1 = list(family = 'Helvetica',
          size = 25,
          color = "black")
t2 = list(family = 'Helvetica',
          size = 15,
          color = "grey")


df = read.csv(paste0(web_path, "/Fig1d_coembedding_RNA_ATAC.csv"))

fig = plot_ly(df[sample(1:nrow(df), 200000),], x=~UMAP_1, y=~UMAP_2, z=~UMAP_3, size = I(30), color = ~group) %>%
    layout(scene = list(xaxis=list(title = list(text ='UMAP_1', font = t1), tickfont = t2),
                        yaxis=list(title = list(text ='UMAP_2', font = t1), tickfont = t2),
                        zaxis=list(title = list(text ='UMAP_3', font = t1), tickfont = t2)))

saveWidget(fig, "Fig1d_coembedding_RNA_ATAC.html")
fig
```

+++

## Fig. 1e

For scATAC-seq cells from each staging bin (x-axis), the density of staging-bins-of-origin of their nearest scRNA-seq derived neighbors (n = 10) in the scGLUE-derived co-embedding is plotted (y-axis).

```{code-cell #Fig1e_RNA_ATAC_timepoint_match} r
web_path = "https://shendure-web.gs.washington.edu/content/members/cxqiu/public/backup/jax_atac/download/reproducing_figures"

library(ggridges)
library(ggplot2)
library(dplyr)
library(viridis)

dat = read.csv(paste0(web_path, "/Fig1e_RNA_ATAC_timepoint_match.csv"))

rna_day = as.vector(dat$rna_day)
rna_day[dat$rna_day == "P0"] = "E19"
dat$rna_day = as.numeric(as.vector(gsub("E","",rna_day)))

p = ggplot(dat, aes(x = rna_day, y = atac_day, fill = atac_day)) +
    geom_density_ridges(bandwidth = 1) +
    theme_classic(base_size = 10) +
    scale_fill_viridis(discrete=TRUE) +
    theme(legend.position = "none") +
    coord_flip() +
    labs(x = "Staging bins of scRNA-seq nearest neighbors", y = "scATAC-seq staging bin") +
    theme(axis.text.x = element_text(color="black", angle = 90, vjust = 1, hjust=1), axis.text.y = element_text(color="black"))

ggsave("Fig1e_RNA_ATAC_timepoint_match.pdf", p, height = 5, width = 5)
p
```
