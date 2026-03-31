
######################################################################################################
### Figure 2. Iterative annotation of a single cell mouse developmental chromatin accessibility atlas.

##################################
### Fig. 2a: 3D UMAP visualization of the entire dataset (left), epithelial cells (middle), and kidney cells (right).

web_path = "https://shendure-web.gs.washington.edu/content/members/cxqiu/public/backup/jax_atac/download/reproducing_figures"

library(dplyr)
library(plotly)
library(htmlwidgets)
library(gplots)

t1 = list(family = 'Helvetica',
          size = 25,
          color = "black")
t2 = list(family = 'Helvetica',
          size = 15,
          color = "grey")

celltype_L1_color_plate = c("Neuroectoderm" = "#f96100",
                            "Mesoderm" = "#bb46c5",
                            "Epithelial_cells" = "#af9fb6",
                            "Erythroid_cells" = "#dc453e",
                            "Muscle_cells" = "#ffa1f5",
                            "White_blood_cells" = "#7ca0ff",
                            "Hepatocytes" = "#185700",
                            "Endothelium" = "#00a34e",
                            "NMPs_and_spinal_cord_progenitors" = "#cf6a79",
                            "Adipocytes" = "#5e7fbf",
                            "Neural_crest_PNS_glia" = "#fff167",
                            "Olfactory_sensory_neurons" = "#e6230b",
                            "Oligodendrocytes" = "#916e00")

celltype_L2_color_plate = c("Epithelial_cells_L2" = "#af9fb6",
                            "Gut_epithelial_cells_L2" = "#ff007a",
                            "Kidney_L2" = "#65c17d",
                            "Lung_and_airway_L2" = "#02b0d1")

### Left:
dat = read.csv(paste0(web_path, "/Fig2a_3D_UMAP_whole_dataset.csv"))
fig = plot_ly(dat, x = ~UMAP_1, y = ~UMAP_2, z = ~UMAP_3, size = I(30), color = ~celltype_L1, colors = celltype_L1_color_plate) %>% 
    layout(scene = list(xaxis=list(title = list(text ='UMAP_1', font = t1), tickfont = t2),
                        yaxis=list(title = list(text ='UMAP_2', font = t1), tickfont = t2),
                        zaxis=list(title = list(text ='UMAP_3', font = t1), tickfont = t2),
                        camera = list(eye = list(x = -0.8, y = 2, z = 1.5))))

saveWidget(fig, "Fig2a_3D_UMAP_whole_dataset_Left.html")

### Middle:
dat = read.csv(paste0(web_path, "/Fig2a_3D_UMAP_epithelial_cells.csv"))
fig = plot_ly(dat, x = ~UMAP_1, y = ~UMAP_2, z = ~UMAP_3, size = I(30), color = ~celltype_L2, colors = celltype_L2_color_plate) %>% 
    layout(scene = list(xaxis=list(title = list(text ='UMAP_1', font = t1), tickfont = t2),
                        yaxis=list(title = list(text ='UMAP_2', font = t1), tickfont = t2),
                        zaxis=list(title = list(text ='UMAP_3', font = t1), tickfont = t2),
                        camera = list(eye = list(x = -0.8, y = 2, z = 1.5))))

saveWidget(fig, "Fig2a_3D_UMAP_whole_dataset_Middle.html")

### Right:
dat = read.csv(paste0(web_path, "/Fig2a_3D_UMAP_kidney_cells.csv"))
fig = plot_ly(dat, x = ~UMAP_1, y = ~UMAP_2, z = ~UMAP_3, size = I(30), color = ~celltype_L3) %>% 
    layout(scene = list(xaxis=list(title = list(text ='UMAP_1', font = t1), tickfont = t2),
                        yaxis=list(title = list(text ='UMAP_2', font = t1), tickfont = t2),
                        zaxis=list(title = list(text ='UMAP_3', font = t1), tickfont = t2),
                        camera = list(eye = list(x = -0.8, y = 2, z = 1.5))))

saveWidget(fig, "Fig2a_3D_UMAP_whole_dataset_Right.html")


################################################################################################################################################################
### Fig. 2c: Heatmap of combined regression coefficients (NNLS; row-scaled) for 13 major scRNA-seq cell clusters (rows) vs. 13 Level-1 scATAC-seq cell lineages.

web_path = "https://shendure-web.gs.washington.edu/content/members/cxqiu/public/backup/jax_atac/download/reproducing_figures"

library(gplots)
library(reshape2)
library(viridis)

dat = read.csv(paste0(web_path, "/Fig2c_NNLS_RNA_ATAC.csv"), row.names=1)

pdf("Fig2c_NNLS_RNA_ATAC.pdf")
heatmap.2(as.matrix(dat), 
          col=viridis, 
          scale="row", 
          Rowv = FALSE, 
          Colv = FALSE, 
          key=T, 
          density.info="none", 
          trace="none", 
          cexRow = 0.5, 
          cexCol = 0.5,
          margins = c(5,5))
dev.off()


####################################################################################################################################################################
### Fig. 2d: Overlap of 36 peak sets from Level-2 scATAC-seq cell classes (rows) vs. 12 peak sets from bulk ATAC-seq of dissected mouse embryonic tissues by ENCODE.

web_path = "https://shendure-web.gs.washington.edu/content/members/cxqiu/public/backup/jax_atac/download/reproducing_figures"

library(gplots)
library(reshape2)
library(viridis)

library(RColorBrewer)
Colors=rev(brewer.pal(11,"Spectral"))
Colors=colorRampPalette(Colors)(120)

dat = read.csv(paste0(web_path, "/Fig2d_overlap_ENCODE.csv"), row.names=1)

pdf("Fig2d_overlap_ENCODE.pdf")
heatmap.2(as.matrix(dat), 
          col=Colors, 
          scale="row", 
          Rowv = F, 
          Colv = F, 
          key=T, 
          density.info="none", 
          trace="none", 
          cexRow = 0.5, 
          cexCol = 0.5,
          margins = c(5,5))
dev.off()









