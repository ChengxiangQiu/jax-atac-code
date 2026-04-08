
###################################
### Train the evolution-aware model

### Support data can be downloaded from:
### https://shendure-web.gs.washington.edu/content/members/cxqiu/public/backup/jax_atac/download/

### Please contact Chengxiang (CX) Qiu for any questions!
### cxqiu@uw.edu or chengxiang.qiu@dartmouth.edu


############################################################
### Step-1: Preparing the candidate window list for training

work_path = ""
web_path = "https://shendure-web.gs.washington.edu/content/members/cxqiu/public/backup/jax_atac/download"
source("help_code/utils.R")

window_list_uniq = readRDS(paste0(work_path, "/corr_Mus_musculus/window_list_uniq.rds"))
pd = readRDS(paste0(work_path, "/umap_Mus_musculus/pd.rds"))
print(sum(window_list_uniq$window_ID == pd$window_ID))
window_list_uniq$celltype = as.vector(pd$anno_L2)

window_list_uniq_no_promoter = window_list_uniq[window_list_uniq$celltype != "Promoters",]

options(scipen = 999)
write.table(window_list_uniq_no_promoter[,c("chr", "start", "end", "celltype")],
            paste0(work_path, "/candidate_region.bed"), row.names=F, col.names=F, sep="\t", quote=F)
### n = 357,282 windows

### excluding windows which are overlapped with ENCODE blacklist
### n = 354,450 windows

### https://shendure-web.gs.washington.edu/content/members/cxqiu/public/backup/jax_atac/download/candidate_windows_354k.bed

dat = read.table(paste0(web_path, "/candidate_windows_354k.bed"))
celltype_list = names(table(dat$V4))
for(celltype in celltype_list){
    print(celltype)
    dat_sub = dat[dat$V4 == celltype,]
    write.table(dat_sub, paste0(work_path, "/", celltype, ".bed"), row.names=F, col.names=F, sep="\t", quote=F)
}

celltype_list=(Adipocyte_cells Adipocyte_cells_Cyp2e1 B_cells Brain_capillary_endothelial_cells Cardiomyocytes CNS_neurons Corticofugal_neurons Endocardial_cells Endothelium Epithelial_cells Erythroid_cells Eye Glia Glomerular_endothelial_cells Gut_epithelial_cells Hepatocytes Intermediate_neuronal_progenitors Kidney Lateral_plate_and_intermediate_mesoderm Liver_sinusoidal_endothelial_cells Lung_and_airway Lymphatic_vessel_endothelial_cells Melanocyte_cells Mesoderm Neural_crest_PNS_neurons Neuroectoderm_and_glia Olfactory_ensheathing_cells Olfactory_neurons Oligodendrocytes Skeletal_muscle_cells T_cells White_blood_cells)
for celltype in "${celltype_list[@]}"; do
    echo "${celltype}"
    bedtools sort -i "${celltype}.bed" | bedtools merge -i - > "${celltype}.merge.bed"
done

dat = NULL
for(celltype in celltype_list){
    print(celltype)
    dat_sub = read.table(paste0(work_path, "/", celltype, ".merge.bed"))
    dat_sub$celltype = celltype
    dat = rbind(dat, dat_sub)
}
### 74,621 regions after merging windows within each cell class; 35,445,000



##############################################
### Step-2: Training the evolution-aware model

import anndata as ad
import crested
import numpy as np
import pandas as pd
import os, sys
import matplotlib

work_path = ""

bigwigs_folder = "BigWig_cut_site_norm"
regions_file = "candidate_windows_354k.bed"

# Set the genome (this only includes regular chrs)
genome = crested.Genome(
    "mm10.fa",
    "chromosome_sizes.txt"
)
crested.register_genome(
    genome
)

print(genome.fetch("chr1", 10000000, 10000010))

adata = crested.import_bigwigs(
    bigwigs_folder=bigwigs_folder,
    regions_file=regions_file,
    target_region_width=1000,
    target="count",
)
adata

celltype_list = ["Adipocyte_cells","Adipocyte_cells_Cyp2e1","B_cells","Brain_capillary_endothelial_cells","CNS_neurons","Cardiomyocytes","Corticofugal_neurons","Endocardial_cells","Endothelium","Epithelial_cells","Erythroid_cells","Eye","Glia","Glomerular_endothelial_cells","Gut_epithelial_cells","Hepatocytes","Intermediate_neuronal_progenitors","Kidney","Lateral_plate_and_intermediate_mesoderm","Liver_sinusoidal_endothelial_cells","Lung_and_airway","Lymphatic_vessel_endothelial_cells","Melanocyte_cells","Mesoderm","Neural_crest_PNS_neurons","Neuroectoderm_and_glia","Olfactory_ensheathing_cells","Olfactory_neurons","Oligodendrocytes","Skeletal_muscle_cells","T_cells","White_blood_cells"]
celltype_dict = {}
with open(f"{web_path}/celltype_list_convert.txt") as f:
    for line in f:
        l = line.rstrip().split('\t')
        celltype_dict[l[0]] = celltype_dict.get(l[0], [])
        celltype_dict[l[0]].append(l[2])

X_new = []
for celltype in celltype_list:
    print(celltype)
    mask = adata.obs.index.isin(celltype_dict[celltype])
    X_sub = adata.X[mask]
    if X_sub.shape[0] > 0:
        X_new.append(np.asarray(X_sub).max(axis=0))

X_new = np.array(X_new)

adata_new = ad.AnnData(X=X_new, obs=pd.DataFrame(index=celltype_list, data={"file_path": celltype_list}), var=adata.var)

adata = adata_new

crested.pp.train_val_test_split(
    adata, strategy="chr", val_chroms=["chr8", "chr10"], test_chroms=["chr9", "chr18"]
)

print(adata.var["split"].value_counts())
adata.var.head(3)

crested.pp.change_regions_width(
    adata, 2114
)

crested.pp.normalize_peaks(
    adata, top_k_percent=0.03
)

adata.write_h5ad(os.path.join(work_path, "data_window_cluster.h5ad"))

### subset 10K top peaks in each cell type
def compute_softmax_stable(x):
    x = np.array(x, dtype=np.float64)
    exp_x = np.exp(x - np.max(x))
    return exp_x / np.sum(exp_x)

X_softmax = np.apply_along_axis(compute_softmax_stable, axis=0, arr=adata.X)

topk = 10000
selected_columns = set()
for row in X_softmax:
    top_idx = np.argsort(row)[-topk:]
    selected_columns.update(top_idx)

selected_columns = sorted(list(selected_columns))
adata_sub = adata[:, selected_columns].copy()
adata_sub.write_h5ad(os.path.join(work_path, "data_window_cluster_top10K.h5ad"))

adata = ad.read_h5ad(os.path.join(work_path, "data_window_cluster_top10K.h5ad"))
adata

print(adata.var["split"].value_counts())
adata.var.head(3)

datamodule = crested.tl.data.AnnDataModule(
    adata,
    batch_size=256,
    max_stochastic_shift=3,
    always_reverse_complement=True,
)

model_architecture = crested.tl.zoo.dilated_cnn(
    seq_len=2114, num_classes=len(list(adata.obs_names))
)

config = crested.tl.default_configs(
    "peak_regression"
)
print(config)

import keras

optimizer = keras.optimizers.Adam(learning_rate=1e-3)
loss = crested.tl.losses.CosineMSELogLoss(max_weight=100, multiplier=1)
metrics = [
    keras.metrics.MeanAbsoluteError(),
    keras.metrics.MeanSquaredError(),
    keras.metrics.CosineSimilarity(axis=1),
    crested.tl.metrics.PearsonCorrelation(),
    crested.tl.metrics.ConcordanceCorrelationCoefficient(),
    crested.tl.metrics.PearsonCorrelationLog(),
    crested.tl.metrics.ZeroPenaltyMetric(),
]

alternative_config = crested.tl.TaskConfig(optimizer, loss, metrics)
print(alternative_config)

trainer = crested.tl.Crested(
    data=datamodule,
    model=model_architecture,
    config=alternative_config,
    project_name="window_cluster",
    run_name="basemodel",
    logger="wandb",
    seed=7,
)

trainer.fit(
    epochs=60,
    learning_rate_reduce_patience=3,
    early_stopping_patience=6,
    model_checkpointing_best_only=False
)


####################################################
### Step-3: Finetuning on cell type-specific regions

X_softmax = np.apply_along_axis(compute_softmax_stable, axis=0, arr=adata.X)

topk = 3000
selected_columns = set()
for row in X_softmax:
    top_idx = np.argsort(row)[-topk:]
    selected_columns.update(top_idx)

selected_columns = sorted(list(selected_columns))
adata_sub = adata[:, selected_columns].copy()
adata_sub.write_h5ad(os.path.join(work_path, "data_window_cluster_top3K.h5ad"))

adata = ad.read_h5ad(os.path.join(work_path, "data_window_cluster_top3K.h5ad"))
adata

datamodule = crested.tl.data.AnnDataModule(
    adata,
    batch_size=64,
    max_stochastic_shift=3,
    always_reverse_complement=True,
)

model_architecture = keras.models.load_model(
    f"{work_path}/mouse_fake_track_14/window_cluster/basemodel/checkpoints/15.keras",
    compile=False,
)

optimizer = keras.optimizers.Adam(learning_rate=1e-4)  # Lower LR!
loss = crested.tl.losses.CosineMSELogLoss(max_weight=100, multiplier=1)
metrics = [
    keras.metrics.MeanAbsoluteError(),
    keras.metrics.MeanSquaredError(),
    keras.metrics.CosineSimilarity(axis=1),
    crested.tl.metrics.PearsonCorrelation(),
    crested.tl.metrics.ConcordanceCorrelationCoefficient(),
    crested.tl.metrics.PearsonCorrelationLog(),
    crested.tl.metrics.ZeroPenaltyMetric(),
]

alternative_config = crested.tl.TaskConfig(optimizer, loss, metrics)
print(alternative_config)

trainer = crested.tl.Crested(
    data=datamodule,
    model=model_architecture,
    config=alternative_config,
    project_name="window_cluster",
    run_name="finetuned_model",
    logger="wandb",
)

trainer.fit(
    epochs=60,
    learning_rate_reduce_patience=3,
    early_stopping_patience=6,
    model_checkpointing_best_only=False
)


##############################
### Step-4: Evaluate the model

adata = ad.read_h5ad(os.path.join(work_path, "data_window_cluster_top3K.h5ad"))

datamodule = crested.tl.data.AnnDataModule(
    adata,
    batch_size=256,
)

# load an existing model
evaluator = crested.tl.Crested(data=datamodule)
model_path = f"{web_path}/CREsted_model/evolution_aware_model.keras"

evaluator.load_model(
    model_path,
    compile=True,
)

# evaluate the model on the test set
evaluator.test()




