
################################################
### Evaluating evolution-naive model performance

### Support data can be downloaded from:
### https://shendure-web.gs.washington.edu/content/members/cxqiu/public/backup/jax_atac/download/

### Please contact Chengxiang (CX) Qiu for any questions!
### cxqiu@uw.edu or chengxiang.qiu@dartmouth.edu


########################################
### Step-1: output the prediction matrix
### <Python script>

import sys, os
import anndata as ad
import crested
import numpy as np
import matplotlib
import pandas as pd
import keras

work_path = ""
adata = ad.read_h5ad(os.path.join(work_path, "celltype_L2_cut_norm", "data_celltype_L2_top3K.h5ad"))
print(adata.var["split"].value_counts())

genome = crested.Genome(
    "mm10.fa",
    "chromosome_sizes.txt"
)
crested.register_genome(
    genome
)
print(genome.fetch("chr1", 10000000, 10000010))

datamodule = crested.tl.data.AnnDataModule(
    adata,
    batch_size=256,
)

# load an existing model
evaluator = crested.tl.Crested(data=datamodule)
model_path = os.path.join(work_path, "celltype_L2_cut_norm", "celltype_L2", "finetuned_model", "checkpoints", "07.keras")

evaluator.load_model(
    model_path,
    compile=True,
)

evaluator.test()

model = keras.models.load_model(model_path, compile=False)
predictions = crested.tl.predict(adata, model)
adata.layers["celltype_L2"] = predictions.T

adata.write_h5ad(os.path.join(work_path, "celltype_L2_cut_norm", "data_celltype_L2_top3K_prediction.h5ad"))

print(adata.X.shape)
print(adata.layers["celltype_L2"].shape)
print(adata.var.shape)

output_dir = work_path + "/celltype_L2_cut_norm/celltype_L2"

X_array = adata.X.toarray() if hasattr(adata.X, "toarray") else adata.X
X_df = pd.DataFrame(np.round(X_array.T, 3), index=adata.var_names, columns=adata.obs_names)
X_df.to_csv(os.path.join(output_dir, "celltype_L2_3K.obs.txt"), sep="\t")

layer_array = adata.layers["celltype_L2"].toarray() if hasattr(adata.layers["celltype_L2"], "toarray") else adata.layers["celltype_L2"]
layer_df = pd.DataFrame(np.round(layer_array.T, 3), index=adata.var_names, columns=adata.obs_names)
layer_df.to_csv(os.path.join(output_dir, "celltype_L2_3K.pre.txt"), sep="\t")

adata.var.to_csv(os.path.join(output_dir, "celltype_L2_3K.var.txt"), sep="\t")


###############################################################################
### Step-2: making barplot of predictions across cell classes for specific site

work_path = ""
web_path = "https://shendure-web.gs.washington.edu/content/members/cxqiu/public/backup/jax_atac/download"
source("help_code/utils.R")

celltype_list = read.table(paste0(web_path, "/cell_class_id.txt"), sep="\t")
colnames(celltype_list) = c("celltype_name", "celltype_id")
celltype_list$celltype_id = paste0("celltype_", celltype_list$celltype_id)

obs = read.table("celltype_L2_3K.obs.txt"), header=T, row.names=1, as.is=T)
obs = obs[,as.vector(celltype_list$celltype_id)]
obs = as.matrix(obs)

pre = read.table("celltype_L2_3K.pre.txt", header=T, row.names=1, as.is=T)
pre = pre[,as.vector(celltype_list$celltype_id)]
pre = as.matrix(pre)

dat = read.table("celltype_L2_3K.var.txt", header=T, row.names=1, as.is=T)

keep = dat$split == "test"
obs = obs[keep,]
pre = pre[keep,]
dat = dat[keep,]

x = cor(t(obs), t(pre))
x = diag(x)
x = x[order(x, decreasing = T)]

target_site_list = c("chr9:15561725-15563839", "chr9:92825930-92828044", "chr9:70276731-70278845")
for(target_site in target_site_list){
    print(target_site)
    target_site_out = gsub('[:|-]', '_', target_site)
    
    df = data.frame(celltype_id = colnames(obs),
                    obs = as.vector(obs[target_site,]),
                    pre = as.vector(pre[target_site,]))
    df = df %>% left_join(celltype_list, by = "celltype_id")
    df$celltype_name = factor(df$celltype_name, levels = as.vector(celltype_list$celltype_name))
    df$target_siet = target_site
    
    p1 = ggplot(data = df, aes(x = celltype_name, y = pre, fill = celltype_name)) +
        geom_bar(stat="identity") + plot_tmp + scale_fill_manual(values=celltype_L2_color_plate)
    p2 = ggplot(data = df, aes(x = celltype_name, y = obs, fill = celltype_name)) +
        geom_bar(stat="identity") + plot_tmp + scale_fill_manual(values=celltype_L2_color_plate)
    ggsave(paste0(target_site_out, ".pdf"), p1 / p2, width = 6, height = 3)
}


#############################################
### Step-3: making heatmap cross cell classes

work_path = ""
web_path = "https://shendure-web.gs.washington.edu/content/members/cxqiu/public/backup/jax_atac/download"
source("help_code/utils.R")

celltype_list = read.table(paste0(web_path, "/cell_class_id.txt"), sep="\t")
colnames(celltype_list) = c("celltype_name", "celltype_id")
celltype_list$celltype_id = paste0("celltype_", celltype_list$celltype_id)

obs = read.table("celltype_L2_3K.obs.txt"), header=T, row.names=1, as.is=T)
obs = obs[,as.vector(celltype_list$celltype_id)]
obs = as.matrix(obs)

pre = read.table("celltype_L2_3K.pre.txt", header=T, row.names=1, as.is=T)
pre = pre[,as.vector(celltype_list$celltype_id)]
pre = as.matrix(pre)

dat = read.table("celltype_L2_3K.var.txt", header=T, row.names=1, as.is=T)

keep = dat$split == "test"
obs = obs[keep,]
pre = pre[keep,]
dat = dat[keep,]

res = cor(obs, pre)

library("gplots")
library(viridis)
library(RColorBrewer)
Colors=rev(brewer.pal(11,"Spectral"))
Colors=colorRampPalette(Colors)(120)
pdf("correlation_cross_celltypes.pdf", 10, 10)
p = heatmap.2(as.matrix(res), 
              col=Colors, 
              scale="none", 
              Rowv = T, 
              Colv = T, 
              key=T, 
              density.info="none", 
              trace="none", 
              cexRow = 0.5, 
              cexCol = 0.5,
              margins = c(5,5))
dev.off()

x = colnames(res)[p$rowInd]
x = c(x[1:2], 'celltype_26',x[3:7],x[9:36])

pdf("correlation_cross_celltypes_reorder.pdf", 10, 10)
heatmap.2(as.matrix(res[x,x]), 
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
dev.off()



###########################################################
### Step-4: making scatter plot for individual cell classes

work_path = ""
web_path = "https://shendure-web.gs.washington.edu/content/members/cxqiu/public/backup/jax_atac/download"
source("help_code/utils.R")

celltype_list = read.table(paste0(web_path, "/cell_class_id.txt"), sep="\t")
colnames(celltype_list) = c("celltype_name", "celltype_id")
celltype_list$celltype_id = paste0("celltype_", celltype_list$celltype_id)

obs = read.table("celltype_L2_3K.obs.txt"), header=T, row.names=1, as.is=T)
obs = obs[,as.vector(celltype_list$celltype_id)]
obs = as.matrix(obs)

pre = read.table("celltype_L2_3K.pre.txt", header=T, row.names=1, as.is=T)
pre = pre[,as.vector(celltype_list$celltype_id)]
pre = as.matrix(pre)

dat = read.table("celltype_L2_3K.var.txt", header=T, row.names=1, as.is=T)

keep = dat$split == "test"
obs = obs[keep,]
pre = pre[keep,]
dat = dat[keep,]

obs = log1p(obs)
pre = log1p(pre)

print(round(diag(cor(obs, pre)), 2))

df = data.frame(peak_id = rep(rownames(dat), ncol(obs)),
                obs = c(obs), pre = c(pre), 
                celltype_id = rep(celltype_list$celltype_id, each = nrow(obs)),
                celltype_name = rep(celltype_list$celltype_name, each = nrow(obs)))
df$celltype_name = factor(df$celltype_name, levels = as.vector(celltype_list$celltype_name))

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
ggsave("corr_each_celltype.png", dpi = 300, width = 20, height = 20)


##############################################################################
### Step-5: merging all cell types to make a scatter plot (presenting density)

p = ggplot(data = df, aes(x=obs, y=pre) ) +
    geom_hex(bins = 70) +
    geom_abline(intercept = 2, slope = 1, linetype = "dashed", color = "red") +
    geom_abline(intercept = -3, slope = 1, linetype = "dashed", color = "red") +
    scale_x_continuous(limits = c(-0.3, 5.25), expand = expansion(mult = 0, add = 0)) +
    scale_y_continuous(limits = c(-0.3, 5.25), expand = expansion(mult = 0, add = 0)) +
    labs(x="Log (observation + 1)", y="Log (prediction + 1)", title="") +
    scale_fill_continuous(type = "viridis", trans = "log10") +
    theme_classic(base_size = 10) +
    theme(legend.position="none") +
    theme(plot.title = element_text(hjust = 0.5)) +
    theme(axis.text.x = element_text(color="black"), axis.text.y = element_text(color="black")) 
ggsave("held_out_set.pdf", p, height = 5, width = 5)

print(cor(df$obs, df$pre))


